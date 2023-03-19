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
import numpy as np
import re


def main(*peak_profiles_ps):
    df_list = []
    frames = []
    for p in peak_profiles_ps:
        df = pd.read_csv(p, sep="\t", index_col=0)
        frame = re.search(".*frame_(\d).tsv", p).group(1)
        df["frame"] = frame
        frames.append(frame)
        df.rename(
            {
                "A": f"{frame}_A",
                "G": f"{frame}_G",
                "C": f"{frame}_C",
                "T": f"{frame}_T",
            },
            axis=1,
            inplace=True,
        )
        df_list.append(df)
    frames = sorted(np.unique(frames))
    n_frame = len(frames)
    full_df = pd.concat(df_list, axis=1)
    # print(full_df)
    target_channels = []
    for x in ["A", "G", "C", "T"]:
        for f in frames:
            target_channels.append(f"{f}_{x}")
    profile = full_df[target_channels]
    print(profile)
    # print(profile.values.shape)
    reshaped_profile = profile.values.reshape([profile.shape[0], n_frame, 4])
    swapped = np.swapaxes(reshaped_profile, 1, 2)
    print(swapped)
    np.save("anchor_profile_for_decoding.npy", swapped, allow_pickle=True)


if __name__ == "__main__":
    fire.Fire(main)
