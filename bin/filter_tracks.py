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
import trackpy as tp


def main(stem, tracks):
    tracks_df = pd.read_csv(tracks, sep="\t", index_col=0)
    n_cycle = tracks_df.frame.max()
    print(tracks_df.shape, n_cycle)
    fully_connected = tp.filter_stubs(tracks_df, n_cycle)
    print(fully_connected)
    # tracks_df.to_csv(f"{stem}_filtered.tsv", sep="\t")


if __name__ == "__main__":
    fire.Fire(main)
