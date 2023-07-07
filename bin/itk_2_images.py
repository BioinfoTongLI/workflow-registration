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
import itk

# from ome_types import from_tiff
# from dask_image.imread import imread


def main(ref, moving, target_ch_index=0):
    # ome_md = from_tiff(ref_img, parser="lxml")
    # pixelsize = ome_md.images[0].pixels.physical_size_x
    # print(pixelsize)

    fixed_image = itk.imread(ref)
    print(fixed_image)
    moving_image = itk.imread(moving)
    print(moving_image)

    parameter_object = itk.ParameterObject.New()
    default_rigid_params = parameter_object.GetDefaultParameterMap("rigid")
    parameter_object.AddParameterMap(default_rigid_params)

    registered_image, params = itk.elastix_registration_method(
        fixed_image,
        moving_image,
        parameter_object=parameter_object,
        log_to_console=False,
    )

    print(params)
    itk.imwrite(registered_image, "itk_reg_stack.tif")
    # print(tmats.tsv)


if __name__ == "__main__":
    fire.Fire(main)
