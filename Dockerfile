FROM nfcore/base
MAINTAINER Matthias Stahl <matthias.stahl@ki.se>

COPY ./piDeepNet /piDeepNet
COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/togetherforever/bin:$PATH
