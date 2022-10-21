FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps ffmpeg libsm6 libxext6 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#RUN python3 -m pip install --upgrade pip itk-elastix==0.12.0 wsireg fire dask_image --no-cache-dir
RUN python3 -m pip install --upgrade pip itk==v5.3rc04 wsireg fire dask_image --no-cache-dir
