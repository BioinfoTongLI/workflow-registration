#! /bin/sh
#
# build.sh
# Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
#


docker build -t bioinfotongli/registration:microaligner -f Dockerfile.microaligner .
docker build -t bioinfotongli/registration:wsireg -f Dockerfile.wsireg .
docker build -t bioinfotongli/registration:qc -f Dockerfile.QC .
#singularity build /lustre/scratch126/cellgen/team283/imaging_sifs/microaligner.sif docker-daemon://microaligner:latest
