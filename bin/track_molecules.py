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


def main(*tsvs, search_range=6):
    tp.linking.Linker.MAX_SUB_NET_SIZE = 200
    full_df = pd.concat([pd.read_csv(tsv, sep="\t", index_col=0) for tsv in tsvs])
    tracked_df = tp.link(full_df, search_range=search_range)
    tracked_df["search_range"] = search_range
    tracked_df.to_csv(f"tracked_peaks_search_range_{search_range}.tsv", sep="\t")


if __name__ == "__main__":
    fire.Fire(main)
