FROM continuumio/miniconda3:4.9.2
RUN apt update & apt install libgl1-mesa-glx -y
RUN conda install numpy pandas tifffile scikit-image dask
RUN pip install opencv-contrib-python
RUN git clone https://github.com/VasylVaskivskyi/opt_flow_reg.git
RUN git clone https://github.com/VasylVaskivskyi/image_registrator.git
