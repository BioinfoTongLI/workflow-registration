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
from microaligner import FeatureRegistrator, transform_img_with_tmat
from dask_image.imread import imread


def main(*imgs):
    freg = FeatureRegistrator()
    freg.ref_img = imread(imgs[0])[0]
    freg.mov_img = imread(imgs[1])[0]
    transformation_matrix = freg.register()
    print(transformation_matrix)


if __name__ == "__main__":
    fire.Fire(main)
