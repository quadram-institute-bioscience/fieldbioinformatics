FROM debian:sid-slim
LABEL authors="Thanh Le Viet thanh.le-viet@quadram.ac.uk" \
      description="Docker image containing all requirements for the artic pipeline with coronahit protocol"

RUN apt-get update &&\
    apt-get install -y wget g++ git make procps samtools=1.11-1 &&\
    apt-get clean -y &&\
    rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    bash ./Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b  &&\
    rm Miniconda3-latest-Linux-x86_64.sh

RUN /miniconda/bin/conda install -c conda-forge mamba
COPY . /fieldbioinformatics
RUN /miniconda/bin/mamba env create -f /fieldbioinformatics/environment.yml &&\
    cd /fieldbioinformatics && /miniconda/envs/artic/bin/python -m pip install . &&\
    /miniconda/bin/mamba clean -y --all &&\
    rm /miniconda/envs/artic/bin/samtools &&\
    ln -s /usr/bin/samtools /miniconda/envs/artic/bin/samtools

ENV PATH=/miniconda/envs/artic/bin:$PATH
