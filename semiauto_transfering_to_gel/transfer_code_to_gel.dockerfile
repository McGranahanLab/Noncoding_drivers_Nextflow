From ubuntu:22.04
MAINTAINER Maria Litovchenko m.litovchenko@ucl.ac.uk

ARG DEBIAN_FRONTEND=noninteractive 
ARG DEBCONF_NONINTERACTIVE_SEEN=true

COPY tmp_dir/* .