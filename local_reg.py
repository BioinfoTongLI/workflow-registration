#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
from pathlib import Path
from ome_zarr.reader import Reader
from ome_zarr.io import parse_url
import xml.etree.ElementTree as ET
import numpy as np
# import dask_stitch as ds
from bigstream import affine
import tifffile as tf

ns = {"ome":"http://www.openmicroscopy.org/Schemas/OME/2016-06"}

import zarr

from bigstream import transform

def get_ch_names(root):
    image_md = root.find('ome:Image', ns)
    pixels = image_md.find('ome:Pixels', ns)
    return [x.attrib["Name"] for x in pixels if "Name" in x.attrib]

def main(zarr_root):
    zr = zarr.open(zarr_root + "/0/0", mode='r')
    print(zr)
    # img_zarr = Path(zarr_root + "/0")
    # reader = Reader(parse_url(img_zarr))
    # raw_data = list(reader())[0].data[0].squeeze()
    root = ET.parse(zarr_root + "/OME/METADATA.ome.xml").getroot()
    ch_names = np.array(get_ch_names(root))
    dapi_indexes = [i for i, ch_n in enumerate(ch_names) if "DAPI" in ch_n]
    fix_img = zr[0, dapi_indexes[0], :]
    mov_img = zr[0, dapi_indexes[1], :]

    spacing = np.array([1, 0.15, 0.15])
    # Note use of mov_lowres_aligned as moving image rather than mov_lowres_data
    # Note also that fix_lowres_spacing is used as the "moving" voxel spacing here
    local_affines = affine.prepare_piecewise_ransac_affine(
        fix_img, mov_img,
        spacing, spacing,
        min_radius=6, max_radius=20, match_threshold=0.75,
        blocksize=[1, 128, 128],
    )
    print(local_affines)

    # function for executing prepared computations
    from bigstream.cluster import execute

    # execute the prepared computation, specifying resource limitations
    local_affines = execute(
        local_affines,
        n_workers=None,
        threads_per_worker=None,
        memory_limit=None,
        config={},
    )

    # now we have a numpy array
    print(type(local_affines))
    print(local_affines.shape)

    # apply the local affines to the moving image
    #   Note we're using mov_lowres_data again - it's better
    #   to provide the global and local affines together. They
    #   are composed into a single transform - that way the moving
    #   image is only resampled one time.
    mov_aligned = transform.prepare_apply_local_affines(
        fix_img, mov_img,
        spacing, spacing,
        local_affines,
        blocksize=[1, 128, 128],
        # global_affine=global_affine,
    )

    # prepared computation, so not a numpy array yet
    print(type(mov_aligned))
    print(mov_aligned.shape)

    # execute using defaults
    # CONSIDER if you need to change resource parameters as discussed in "Executing a prepared compuation" section above
    mov_aligned = execute(mov_aligned)
    print(mov_aligned.shape + " final")
    tf.imwrite(f"{zarr_root}_local_registered.tif", np.squeeze(mov_aligned))


if __name__ == "__main__":
    fire.Fire(main)
