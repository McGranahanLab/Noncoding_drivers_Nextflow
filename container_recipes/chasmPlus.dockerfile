FROM ubuntu:22.04
MAINTAINER Maria Litovchenko m.litovchenko@ucl.ac.uk

RUN apt-get -y update \
    && apt-get -y upgrade \
    && apt-get install -y git libc6 python3-pip \
    && pip3 install numpy==1.21.6 pandas==1.2.5 scipy==1.11.3 open-cravat==2.4.2 \
    && oc module install-base \
    && oc module install -y chasmplus chasmplus_LUAD chasmplus_LUSC