FROM brainlife/fsl:6.0.4-patched2

## MAINTAINER Marco Pereañez <marco.pereanez@mssm.edu> 

RUN apt update

# install Python3 + NumPy
RUN apt -y install python3.8 python3-pip

RUN pip3 install numpy

# install DICOM to NIFTI
RUN apt -y install dcm2niix

# install ants / fsl / other requirements
RUN apt -y install ants 

## add N4BiasFieldCorrection to path
ENV ANTSPATH=/usr/lib/ants
ENV PATH=$PATH:/usr/lib/ants

## run distributed script to set up fsl
RUN . /usr/local/fsl/etc/fslconf/fsl.sh


RUN apt -y install fsl-first-data fsl-atlases 

# ## add better eddy functions
# RUN wget https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.11/centos6/eddy_cuda8.0
# RUN mv eddy_cuda8.0 /usr/local/bin/eddy_cuda
# RUN wget https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.11/centos6/eddy_openmp -P /usr/local/bin
# RUN chmod +x /usr/local/bin/eddy_cuda
# RUN chmod +x /usr/local/bin/eddy_openmp

## install mrtrix3 requirements
RUN apt update && apt -y install git g++ python libeigen3-dev zlib1g-dev libqt5opengl5-dev libqt5svg5-dev libgl1-mesa-dev libfftw3-dev libtiff5-dev libpng-dev

## install and compile mrtrix3
RUN git clone https://github.com/MRtrix3/mrtrix3.git
RUN cd mrtrix3 && ./configure -nogui && ./build

## manually add to path
ENV PATH=$PATH:/mrtrix3/bin

# install Tree
RUN apt -y install tree

## make it work under singularity 
# RUN ldconfig && mkdir -p /N/u /N/home /N/dc2 /N/soft /mnt/scratch

## https://wiki.ubuntu.com/DashAsBinSh 
#/RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Configure entrypoint
# ENTRYPOINT ["/usr/bin/python3 ${FLYWHEEL}/run.py"]
ENTRYPOINT ["/bin/bash"]

## Make directory for flywheel spec (v0)
ENV FLYWHEEL /flywheel/v0
RUN mkdir -p ${FLYWHEEL}
VOLUME ${FLYWHEEL}
WORKDIR ${FLYWHEEL}
COPY manifest.json ${FLYWHEEL}/manifest.json
COPY run.py ${FLYWHEEL}/run.py

