# Dockerfile for viral integration simulation pipeline

FROM mambaorg/micromamba:0.19.1

USER root

#  $ docker build . -t szsctt/simvi:latest -t szsctt/simvi:1
#  $ docker run --rm -it szsctt/simvi:latest /bin/bash
#  $ docker push szsctt/simvi:latest
#  $ docker push szsctt/simvi:1


ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git 

# install conda stuff
ADD scripts/consolidate_envs.py /opt/simvi/scripts/
ADD envs /opt/simvi/envs/
RUN micromamba install -n base -c anaconda pip pyyaml=5.3 python=3 -y &&\
	python3 /opt/simvi/scripts/consolidate_envs.py /opt/simvi/envs/*yml /opt/simvi/envs/simvi.yml
	
RUN micromamba env update -n base -f /opt/simvi/envs/simvi.yml &&\
	micromamba clean --all -y 	

# include simvi scripts, etc
ADD scripts /opt/simvi/scripts/
ADD snakemake_rules /opt/simvi/snakemake_rules
ADD Snakefile /opt/simvi/Snakefile

# add test files
ADD test /opt/simvi/test

WORKDIR /opt/simvi

CMD test/runMe.sh
