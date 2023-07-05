#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2023 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
import numpy as np
from aicsimageio import AICSImage
from aicsimageio.writers.ome_tiff_writer import OmeTiffWriter
from aicsimageio.types import PhysicalPixelSizes


def main(ref, ref_index, out, *movings):
    all_tifs = [p for p in movings]
    all_tifs.insert(ref_index, ref)
    stack = []
    channel_names = []
    for p in all_tifs:
        img = AICSImage(p)
        stack.append(img.get_image_dask_data())
        dims = img.dims
        dim_order = img.dims.order
        channel_names += img.channel_names
        pixel_sizes = img.physical_pixel_sizes
    stack = np.concatenate(stack, axis=1)
    writer = OmeTiffWriter()
    print(stack.shape)
    print(dim_order)
    pixel_sizes = PhysicalPixelSizes(1, pixel_sizes.X, pixel_sizes.Y)
    writer.save(
        stack, out,
        dim_order=[dim_order],
        channel_names=[channel_names],
        image_names=["registered_and_stacked"],
        physical_pixel_sizes=[pixel_sizes],
    )


if __name__ == "__main__":
    fire.Fire(main)
