FROM nfcore/base:1.7
LABEL authors="Alex Lemenze" \
      description="Docker image containing all requirements for alemenze/magicuniqueamplicons pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN apt-get install -y libltdl7 
ENV PATH /opt/conda/envs/magicuniqueamplicons/bin:$PATH
