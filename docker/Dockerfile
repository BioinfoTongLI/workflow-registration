FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps ffmpeg libsm6 libxext6 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN python3 -m pip install --upgrade pip wsireg fire dask_image --no-cache-dir
