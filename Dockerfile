FROM nfcore/base
LABEL authors="Matthias Stahl <matthias.stahl@ki.se>"
LABEL description="Docker image containing all requirements for the togetherforever pipeline."

COPY ./bin/*.py /bin/
COPY ./piDeepNet /piDeepNet
COPY ./environment.yml /
COPY ./environment-percolator.yml /

RUN apt-get update && apt-get -y install gcc

RUN conda update conda
RUN conda env create -f /environment-percolator.yml && conda clean -a
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/togetherforever/bin:$PATH
ENV PATH /opt/conda/envs/togetherforever-percolator/bin:$PATH

RUN Rscript -e "devtools::install_git('git://github.com/ypriverol/pIR.git')"
RUN Rscript -e "devtools::install_git('git://github.com/dosorio/Peptides.git')"

