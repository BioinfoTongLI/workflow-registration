FROM continuumio/miniconda3:4.9.2
RUN conda install numpy pandas tifffile scikit-image dask
RUN pip install opencv-contrib-python
RUN git clone https://github.com/BayraktarLab/opt_flow_reg.git && cd opt_flow_reg && git pull && cd ..
#RUN git clone https://github.com/VasylVaskivskyi/opt_flow_reg.git
#RUN cd opt_flow_reg && git checkout development && git pull && cd ..
RUN git clone https://github.com/BayraktarLab/feature_reg.git && cd feature_reg && git pull && cd ..
#RUN git clone https://github.com/VasylVaskivskyi/image_registrator.git
#RUN cd image_registrator && git checkout development && git pull && cd ..
RUN apt-get --allow-releaseinfo-change update && apt-get install --fix-missing  libgl1-mesa-glx procps -y
