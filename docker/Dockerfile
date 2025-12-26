# we start from the ubuntu:18.04 base image, which will work for most applications
FROM ubuntu:18.04

# installing some extra tools that we don't get by default
RUN apt-get update && apt-get install -y vim wget curl git

# installing miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:$PATH"

# copy the env.yaml file into the docker image
# we do this before copying the project directory so docker doesn't try and rebuild everything every time
COPY ./env.yaml /app/env.yaml

# Install additional Python packages using pip before creating the conda environment
#RUN pip install pysam==0.16.0.1

# create the conda environment within the container
RUN conda env create --file /app/env.yaml

# copy the project directory into /app
COPY . /app

# set the PATH variable so it can find the conda environment
# MAKE SURE TO REPLACE <ENV_NAME> WITH THE NAME OF YOUR CONDA ENVIRONMENT
ENV PATH="$PATH:/miniconda/envs/virusscan/bin:/miniconda/bin"

# Install additional Python packages using pip within the conda environment: added
#RUN /bin/bash -c "source activate virusscan && conda install -c bioconda pysam==0.16.0.1"
# Activate the conda environment and install additional Python packages using pip
# RUN echo "source activate virusscan" >> ~/.bashrc
# RUN /bin/bash -c "source activate virusscan && pip install about-time==3.1.1 numpy==1.21.6 numpy-groupies==0.9.14 pandas==1.3.5 pybedtools==0.8.1 pysam==0.16.0.1 rpy2==2.9.4 scipy==1.7.3"

#Run pip install about-time==3.1.1 numpy==1.21.6 numpy-groupies==0.9.14 pandas==1.3.5 pybedtools==0.8.1 pysam==0.16.0.1 rpy2==2.9.4 scipy==1.7.3

# command that is automatically ran when running a docker image
CMD /bin/bash