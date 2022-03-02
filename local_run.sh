#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

DATE_WITH_TIME=`date "+%Y%m%d%H%M"`
TMP_NF_WORK="/tmp/${USER}_${DATE_WITH_TIME}_registration_work"

NXF_OPTS='-Dleveldb.mmap=false' NXF_WORK=$TMP_NF_WORK nextflow run /lustre/scratch117/cellgen/team283/tl10/workflow-registration/main.nf \
	-params-file $1 \
	-profile standard,singularity \
	-resume
