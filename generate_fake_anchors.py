#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Take ome.tif as input and generate fake anchors channels for decoding cycles
"""
import argparse
import tifffile as tf
import numpy as np
import re
from pathlib import Path
from apeer_ometiff_library import omexmlClass, io
import copy
import pysnooper
import xml.etree.ElementTree as ET


def sort_channels(node, tag="ID"):
    return sorted(node, key=lambda child: int(child.get(tag).split(":")[-1]))


def normalize(stack, quantile):
    img_min = np.nanmin(stack, axis=(0, 1), keepdims=True)
    img_max = np.nanpercentile(stack, q=quantile, axis=(0, 1), keepdims=True)
    normalized = ((stack - img_min) / (img_max - img_min) * 2 * 16).astype(np.uint16)
    return normalized


def update_ome_xml_for_fake_channels(ome_string):
    md = omexmlClass.OMEXML(ome_string)
    pixels = md.image().Pixels

    last_cycle_ind = None
    n_ch_added = 0
    origin_n_ch = pixels.get_channel_count()
    for i in range(origin_n_ch):
        # pixels.remove_channel(i)
        current_name = pixels.Channel(i).get_Name()
        current_cycle_ind = re.search("c0(\d).*", current_name).group(1)
        if not last_cycle_ind:
            last_cycle_ind = current_cycle_ind

        if last_cycle_ind != current_cycle_ind:
            pixels.append_channel(i + n_ch_added, "c0%s" % last_cycle_ind + " anchor")
            n_ch_added += 1

        pixels.Channel(i).set_ID("Channel:0:" + str(i + n_ch_added))

        # for k in pixels.Channel(i).node.attrib:
        # print(k, pixels.Channel(i).node.get(k))
        # print(i, pixels.Channel(i).node.attrib)
        last_cycle_ind = current_cycle_ind
    pixels.append_channel(origin_n_ch + n_ch_added, "c0%s" % last_cycle_ind + " anchor")
    pixels.populate_TiffData()
    return md


def xml_cleanup(ome_xml_str):
    ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    root = ET.ElementTree(ET.fromstring(ome_xml_str)).getroot()
    pixels_et = root[2][3]  # Image -> Pixels
    all_channels = pixels_et.findall("ome:Channel", ns)
    all_tiffdata = pixels_et.findall("ome:TiffData", ns)
    pixel_attrs = pixels_et.attrib
    pixels_et.clear()  # clear the Pixel
    pixels_et.attrib = pixel_attrs
    for ch in sort_channels(all_channels):
        pixels_et.append(ch)
    for td in all_tiffdata:
        pixels_et.append(td)
    root[3].clear()  # remove the pyramid metadata
    return root


def map_old_new_chs(ome_string, new_md):
    ori_ch_names = omexmlClass.OMEXML(ome_string).image().Pixels.get_channel_names()
    ori_ch_names_i = {ch: i for i, ch in enumerate(ori_ch_names)}

    pixels = new_md.image().Pixels
    new_i_ch_name = {
        int(pixels.Channel(i).node.get("ID").split(":")[-1]): pixels.Channel(i).Name
        for i in range(pixels.get_channel_count())
    }

    maps = {}
    for i in sorted(new_i_ch_name.keys()):
        ch_name = new_i_ch_name[i]
        if "anchor" not in ch_name:
            maps[i] = (ch_name, ori_ch_names_i[new_i_ch_name[i]])
        else:
            maps[i] = (ch_name, None)
    return maps


@pysnooper.snoop()
def main(args):
    ignore_dapi = True
    with tf.TiffFile(args.ome_tif, "r", is_ome=True) as fh:
        print(fh)
        ome_string = fh.ome_metadata
        shape = fh.pages[0].shape

    updated_md = update_ome_xml_for_fake_channels(ome_string)

    clean_xml = xml_cleanup(updated_md.to_xml())
    meta = ET.tostring(clean_xml, encoding="ascii")
    ch_map = map_old_new_chs(ome_string, omexmlClass.OMEXML(meta))
    # print(ch_map)

    tmp_anchor = np.zeros(shape, dtype=np.uint16)
    new_name = Path(args.ome_tif).stem + "_with_anchors.ome.tif"

    if args.known_anchor:
        known_anchor_cyc = args.known_anchor.split(" ")[0]

    writer = tf.TiffWriter(new_name, bigtiff=True)
    with tf.TiffFile(args.ome_tif, "r", is_ome=True) as fh:
        for i in sorted(ch_map.keys()):
            print(i, ch_map[i])
            old_ind = ch_map[i][1]
            cur_ch_name = ch_map[i][0]
            if old_ind is not None:
                tmp_array = fh.pages[old_ind].asarray().squeeze()
                writer.save(tmp_array, photometric="minisblack", description=meta)

                if (
                    args.known_anchor != ""
                    and cur_ch_name != args.known_anchor
                    and cur_ch_name.startswith(known_anchor_cyc)
                ):
                    continue

                if ignore_dapi and "DAPI" in cur_ch_name:
                    continue
                print(cur_ch_name + " max projected")
                tmp_anchor = np.max(np.array([normalize(tmp_array, 99.5), tmp_anchor]), axis=0)
                # tmp_anchor = np.max(np.array([normalize(tmp_array, 98), tmp_anchor]), axis=0)
            else:
                print("save anchor, and re-intialize anchor image")
                writer.save(tmp_anchor, photometric="minisblack", description=meta)

                tmp_anchor = np.zeros(shape, dtype=np.uint16)

    writer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-ome_tif", type=str, required=True)
    parser.add_argument("-known_anchor", type=str, default=None)

    args = parser.parse_args()

    main(args)
