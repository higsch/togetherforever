FROM nfcore/base
MAINTAINER Matthias Stahl <matthias.stahl@ki.se>

COPY ./piDeepNet /piDeepNet
COPY ./environment.yml /
COPY ./environment-percolator.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env create -f /environment-percolator.yml && conda clean -a
ENV PATH /opt/conda/envs/togetherforever/bin:$PATH
ENV PATH /opt/conda/envs/togetherforever-percolator/bin:$PATH
