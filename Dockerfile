#############################################################
# Dockerfile tools for pyBsHap
#############################################################
FROM continuumio/miniconda3:latest
MAINTAINER Rahul Pisupati <rahul.pisupati@gmi.oeaw.ac.at>
LABEL authors="rahul.pisupati@gmi.oeaw.ac.at" \
    description="Docker image containing all requirements for the nf-core/methylpy pipeline"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN /bin/bash -c "source activate /opt/conda/envs/pybshap" && pip install scikit-allel==0.20.3 && pip install git+https://github.com/Gregor-Mendel-Institute/pyBsHap.git
ENV PATH /opt/conda/envs/pybshap/bin:$PATH
