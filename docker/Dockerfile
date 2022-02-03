FROM continuumio/miniconda3:4.10.3
RUN apt-get --allow-releaseinfo-change update && apt-get install --fix-missing  libgl1-mesa-glx procps -y
RUN git clone https://github.com/BayraktarLab/opt_flow_reg.git
RUN git clone https://github.com/BayraktarLab/feature_reg.git
COPY ./environment.yml /environment.yml
RUN conda update conda
RUN conda env update -n base --file /environment.yml
