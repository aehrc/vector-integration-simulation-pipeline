# Dockerfile for viral integration simulation pipeline

FROM ubuntu:20.04

#  $ docker build . -t szsctt/intvi_sim:latest -t szsctt/intvi_sim:1
#  $ docker run --rm -it szsctt/intvi_sim:latest /bin/bash
#  $ docker push szsctt/intvi_sim:latest
#  $ docker push szsctt/intvi_sim:1


ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git 

# https://hub.docker.com/r/continuumio/miniconda/dockerfile
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc &&\
    /opt/conda/bin/conda update conda python>3 -y &&\
    /opt/conda/bin/conda clean --all -y 
ENV PATH /opt/conda/bin:$PATH

# install conda stuff
ADD scripts/consolidate_envs.py /opt/intvi_simulation/scripts/
ADD envs /opt/intvi_simulation/envs/
RUN /opt/conda/bin/conda install -n base -c anaconda pip pyyaml=5.3 -y &&\
	python3 /opt/intvi_simulation/scripts/consolidate_envs.py /opt/intvi_simulation/envs/*yml /opt/intvi_simulation/envs/sim.yml &&\
	/opt/conda/bin/conda env update -n base -f /opt/intvi_simulation/envs/sim.yml &&\
	/opt/conda/bin/conda clean --all -y 	

# include intvi_simulation scripts, etc
ADD scripts /opt/intvi_simulation/scripts/
ADD snakemake_rules /opt/intvi_simulation/snakemake_rules
ADD Snakefile /opt/intvi_simulation/Snakefile

# add test files
ADD test /opt/intvi_simulation/test

WORKDIR /opt/intvi_simulation