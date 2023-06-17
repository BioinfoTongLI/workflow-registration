#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

TMP_NF_WORK="/tmp/work/${USER}_registration_work"

NXF_OPTS='-Dleveldb.mmap=false' NXF_WORK=$TMP_NF_WORK nextflow -trace nextflow.executor run /scratch/iss_decoding/nf/workflow-registration/main.nf \
	-params-file $1 \
	-profile local \
	-resume
	#-with-report
