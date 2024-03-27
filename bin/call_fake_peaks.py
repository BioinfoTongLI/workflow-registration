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
import pandas as pd
from dask_image.imread import imread
import trackpy as tp


def main(stem, image, frame, peak_ch_index, spot_radius=3, percentile=90):
    img = imread(image)
    anchor = img[peak_ch_index]
    coding_chs = img[2:]
    print(coding_chs)
    peaks = tp.locate(
        anchor.compute(), spot_radius * 2 - 1, percentile=percentile, separation=2
    )
    peaks["frame"] = int(frame)
    df_to_save = peaks
    if frame != 0:
        peak_profile = coding_chs.compute()[:, peaks.y.astype(int), peaks.x.astype(int)]
        profile_df = pd.DataFrame(
            peak_profile.astype(int).T, columns=["A", "G", "C", "T"]
        )
        df_to_save = pd.merge(peaks, profile_df, left_index=True, right_index=True)
    df_to_save["y_int"] = df_to_save.y.astype(int)
    df_to_save["x_int"] = df_to_save.x.astype(int)
    df_to_save.to_csv(
        f"{stem}_peaks_radius_{spot_radius}_percentile_{percentile}_with_profile.tsv",
        sep="\t",
    )


if __name__ == "__main__":
    fire.Fire(main)
