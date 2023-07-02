#! /bin/sh
#
# build.sh
# Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
#


docker build -t registration:itk -f Dockerfile.itk .
docker build -t registration:microaligner -f Dockerfile.microaligner .
#singularity build /lustre/scratch126/cellgen/team283/imaging_sifs/microaligner.sif docker-daemon://microaligner:latest
