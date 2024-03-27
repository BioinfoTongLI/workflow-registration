FROM python:3.9-slim

RUN apt-get --allow-releaseinfo-change update && apt-get install --fix-missing  libgl1-mesa-glx procps git -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY requirements.txt /tmp/requirements.txt
RUN python3 -m pip install --no-cache-dir -r /tmp/requirements.txt
