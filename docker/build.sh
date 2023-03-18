#! /bin/sh
#
# build.sh
# Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
#


docker build -t microaligner -f Dockerfile.microaligner .
singularity build /lustre/scratch117/cellgen/team283/imaging_sifs/microaligner.sif docker-daemon://microaligner:latest
