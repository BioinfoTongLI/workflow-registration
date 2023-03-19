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
import pandas as pd
from dask_image.imread import imread


def main(anchor_peaks, frame, tiff):
    anchors = pd.read_csv(anchor_peaks, sep="\t", index_col=0)
    img = imread(tiff)
    coding_chs = img[2:]
    assert coding_chs.shape[0] == 4
    peak_profile = coding_chs.compute()[:, anchors.y.astype(int), anchors.x.astype(int)]
    profile_df = pd.DataFrame(peak_profile.astype(int).T, columns=["A", "G", "C", "T"])
    df_to_save = pd.merge(anchors, profile_df, left_index=True, right_index=True)
    df_to_save.to_csv(f"anchor_profile_of_frame_{frame}.tsv", sep="\t")


if __name__ == "__main__":
    fire.Fire(main)
