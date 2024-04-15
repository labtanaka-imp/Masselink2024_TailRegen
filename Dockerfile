################## BASE IMAGE ########################
FROM ubuntu:20.04

################## INSTALLATION ######################
ARG DEBIAN_FRONTEND=noninteractive
ARG dependencies="curl wget cmake unzip build-essential git default-jre python3 python3-pip libbz2-dev liblzma-dev zlib1g-dev libncurses5-dev libncursesw5-dev"

###Versions
ENV Seqkit_Version=2.3.1
ENV Samtools_Version=1.16

###Install dependencies and upgrade
RUN mkdir -p /Apps/ \
    && apt-get update \
    && apt-get install -y $dependencies  \
    && apt-get upgrade -y

##Install Seqkit
RUN cd /Apps/ \
    && wget -O seqkit.tar.gz https://github.com/shenwei356/seqkit/releases/download/v${Seqkit_Version}/seqkit_darwin_amd64.tar.gz \
    && tar xf seqkit.tar.gz \
    && mv seqkit /usr/local/bin/ \
    && rm seqkit.tar.gz

##Install Samtools
###Install htslib
RUN cd /Apps/ \
    && wget https://github.com/samtools/htslib/releases/download/${Samtools_Version}/htslib-${Samtools_Version}.tar.bz2 \
    && tar xvjf htslib-${Samtools_Version}.tar.bz2 \
    && cd htslib-${Samtools_Version} \
    && ./configure \
    && make \
    && make install \
    && rm /Apps/htslib-${Samtools_Version}.tar.bz2

###Install samtools
RUN cd /Apps/ \
    && wget https://github.com/samtools/samtools/releases/download/${Samtools_Version}/samtools-${Samtools_Version}.tar.bz2 \
    && tar -xvjf samtools-${Samtools_Version}.tar.bz2 \
    && cd samtools-${Samtools_Version} \
    && ./configure \
    && make \
    && make install \
    && rm /Apps/samtools-${Samtools_Version}.tar.bz2