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
from wsireg.parameter_maps.reg_params import DEFAULT_REG_PARAM_MAPS
from ome_types import from_tiff
from dask_image.imread import imread


def main(*images, stem="WSIReg", target_ch_index=0):
    ref_img = sorted(images)[0]
    ome_md = from_tiff(ref_img, parser="lxml")
    pixelsize = ome_md.images[0].pixels.physical_size_x
    print(pixelsize)

    # initialize registration graph
    reg_graph = WsiReg2D(stem, stem)
    cycle_names = []
    for i, img in enumerate(images):
        cycle_name = f"cycle{i + 1}"
        cycle_names.append(cycle_name)
        # add registration images (modalities)
        reg_graph.add_modality(
            cycle_name,
            f"{img}",
            image_res=pixelsize,
            channel_names=["DAPI", "anchor", "A", "G", "C", "T"],
            # channel_colors=["blue", "green", "red", "yellow", "cyan"],
            preprocessing={
                "image_type": "FL",
                "ch_indices": [target_ch_index],
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
        print(f"{cycleName} to {cycle_names[0]} added")
        reg_graph.add_reg_path(
            cycleName,
            cycle_names[0],
            thru_modality=None,
            # reg_params=["rigid", "affine"],
            # reg_params=[RegModel.rigid, RegModel.affine],
            reg_params=[
                DEFAULT_REG_PARAM_MAPS["rigid"],
                DEFAULT_REG_PARAM_MAPS["affine"],
                DEFAULT_REG_PARAM_MAPS["nl"],
            ],
        )

    reg_graph.register_images()
    # reg_graph.save_transformations()
    reg_graph.transform_images(file_writer="ome.tiff")


if __name__ == "__main__":
    fire.Fire(main)
