#! /bin/sh
#
# build.sh
# Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
#


docker build -t bioinfotongli/registration:itk -f Dockerfile.itk .
docker build -t bioinfotongli/registration:microaligner -f Dockerfile.microaligner .
docker build -t bioinfotongli/registration:wsireg -f Dockerfile.wsireg .
#singularity build /lustre/scratch126/cellgen/team283/imaging_sifs/microaligner.sif docker-daemon://microaligner:latest
