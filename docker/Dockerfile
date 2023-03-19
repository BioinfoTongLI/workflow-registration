FROM nvidia/cuda:11.2.0-cudnn8-runtime-ubuntu20.04

RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
	python3-dev \
   	python3-pip \
    	python3-wheel \
    	python3-setuptools \
    	python-is-python3 \
    	procps \
    	libgomp1 ffmpeg libsm6 libxext6 -y && \
    	apt-get clean && \
    	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN python3 -m pip install --upgrade pip \
	itk==v5.3rc04 \
	fire==0.4.0 \
      	dask_image==2022.9.0 \
      	trackpy==0.5.0 \
      	cucim==22.10.0 \
      	wsireg==0.3.7 \
	ome-types==0.3.1 \
	cupy-cuda11x \
	xmlschema==1.4.1 \
	--no-cache-dir
