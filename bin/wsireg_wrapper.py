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
from wsireg.wsireg2d import WsiReg2D


def main(stem, *images):
    # initialize registration graph
    reg_graph = WsiReg2D(stem, stem)

    print(images, type(images))

    cycle_names = []
    for i, img in enumerate(images):
        cycle_name = f"cycle{i + 1}"
        cycle_names.append(cycle_name)
        # add registration images (modalities)
        reg_graph.add_modality(
            cycle_name,
            f"{img}",
            image_res=0.65,
            channel_names=["DAPI", "A", "G", "C", "T"],
            channel_colors=["blue", "green", "red", "yellow", "cyan"],
            preprocessing={
                "image_type": "FL",
                "ch_indices": [0],
                "as_uint8": True,
                "contrast_enhance": True,
            },
        )

    # reg_graph.add_modality(
        # "modality_brightfield",
        # "./data/im2.tiff",
        # image_res=1,
        # preprocessing={"image_type": "BF", "as_uint8": True, "invert_intensity": True},
    # )


    for cycleName in cycle_names[1:]:
        print(f"{cycleName} to ref added")
        reg_graph.add_reg_path(
            cycle_names[0],
            cycleName,
            thru_modality=None,
            reg_params=["rigid", "affine"],
        )

    reg_graph.register_images()
    reg_graph.save_transformations()
    reg_graph.transform_images(file_writer="ome.tiff")


if __name__ == "__main__":
    fire.Fire(main)
