FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# basic OS packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    curl \
    git \
    vim \
    bzip2 \
    xz-utils \
    && rm -rf /var/lib/apt/lists/*

# install miniconda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm -f /tmp/miniconda.sh

ENV PATH="/opt/conda/bin:${PATH}"

# copy only env file first (better layer caching)
WORKDIR /app
COPY env.yaml /app/env.yaml

# create env (IMPORTANT: create by NAME so location is standard in container)
RUN conda env create -n virusscan -f /app/env.yaml && \
    conda clean -afy

# make env available without "conda activate"
ENV PATH="/opt/conda/envs/virusscan/bin:${PATH}"

# copy your pipeline code
COPY . /app

CMD ["/bin/bash"]