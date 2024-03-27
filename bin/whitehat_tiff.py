#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
import numpy as np
import cupy as cp
from dask_image.imread import imread
import tifffile as tf
from cucim.skimage.morphology import white_tophat, disk
from ome_types import from_tiff, to_xml


def white_tophat_cp(chunk, **kwargs):
    return white_tophat(cp.array(chunk), **kwargs).get()


def main(stem, tiff, spot_radius, method="max", enhance=True):
    img = imread(tiff).squeeze()
    dapi_ch = img[[0]]
    coding_chs = img[1:]
    # rechunked_img = img.rechunk({0:1, 1:"auto", 2:"auto"})
    rechunked_coding_chs = coding_chs.rechunk({0: 1, 1: 4000, 2: 4000})
    print(rechunked_coding_chs)

    if enhance:
        footprint = np.expand_dims(disk(spot_radius), 0)
        print(footprint.shape)
        hat_enhenced = rechunked_coding_chs.map_overlap(
            white_tophat_cp,
            depth=(0, spot_radius * 3, spot_radius * 3),
            footprint=footprint,
            dtype=np.float16,
        )
    else:
        hat_enhenced = rechunked_coding_chs
        print("not running enhancement")

    if method == "std":
        fake_anchor = np.std(hat_enhenced, 0)
    else:
        fake_anchor = np.max(hat_enhenced, 0)
    print(fake_anchor)
    all_stakced = np.vstack([dapi_ch, np.expand_dims(fake_anchor, 0), hat_enhenced])
    print(all_stakced)

    # ome_md = from_tiff(tiff)
    # print(ome_md)
    # ome_md.images[0].pixels.size_c = all_stakced.shape[0]
    # print(ome_md.images[0].pixels.size_c)
    # omexml = to_xml(ome_md)
    # print(omexml)
    tf.imwrite(
        f"{stem}_whitehat_enhanced_radius_{spot_radius}.ome.tif",
        all_stakced.compute(),
        # description=omexml.encode("utf-8"),
    )


if __name__ == "__main__":
    fire.Fire(main)
