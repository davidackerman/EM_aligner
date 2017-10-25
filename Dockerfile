FROM debian:jessie as mcr

RUN apt-get update && apt-get install -y \
    wget \
    tar \
    unzip \
    xorg \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /mcr-install && \
    cd /mcr-install && \
    wget -nv http://ssd.mathworks.com/supportfiles/downloads/R2016b/deployment_files/R2016b/installers/glnxa64/MCR_R2016b_glnxa64_installer.zip && \
    unzip MCR_R2015b_glnxa64_installer.zip && \
    mkdir /opt/mcr && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf /mcr-install

ENV LD_LIBRARY_PATH /opt/mcr/v90/runtime/glnxa64:/opt/mcr/v90/bin/glnxa64:/opt/mcr/v90/sys/os/glnxa64
COPY template_production_scripts/pipeline_scripts/bin /usr/bin/EM_Aligner

