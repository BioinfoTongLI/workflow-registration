#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

MOUNT_POINT='/tmp/work/'

DATE_WITH_TIME=`date "+%Y%m%d%H%M"`
TMP_NF_WORK="$MOUNT_POINT/${USER}_${DATE_WITH_TIME}_registration_work"

NXF_OPTS='-Dleveldb.mmap=false' NXF_VER=22.04.5 NXF_WORK=$TMP_NF_WORK nextflow -trace nextflow.executor run /lustre/scratch126/cellgen/team283/tl10/workflow-registration/main.nf \
	-params-file $1 \
	-profile local \
	-with-report
	#-resume
