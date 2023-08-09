From ubuntu:22.04
MAINTAINER Maria Litovchenko m.litovchenko@ucl.ac.uk

ARG DEBIAN_FRONTEND=noninteractive 
ARG DEBCONF_NONINTERACTIVE_SEEN=true

# install R v 4.3.0
# BioConductor packages, version is controlled though BiocManager version (3.17)
# biomaRt v2.56.1
# GenomicFeatures v1.52.2, GenomicRanges v1.52.0 
# rtracklayer v1.60.1
# TxDb.Hsapiens.UCSC.hg19.knownGene v3.2.2, TxDb.Hsapiens.UCSC.hg38.knownGene v3.17.0
# installation of modified dNdScv package from github
RUN apt-get update && apt-get install -y git autoconf gcc git make ssh \
                wget vim build-essential software-properties-common \
                ca-certificates libssl-dev libcurl4-gnutls-dev \
                libatlas-base-dev libffi-dev libxml2-dev \
                libncurses5-dev libsm6 libxrender1 libfontconfig1 \
                libxt6 libtcl8.6 libtk8.6 glibc-source libharfbuzz-dev \
                libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev \
                libjpeg-dev fontconfig cmake libfontconfig1-dev xclip libtiff5 \
    && cp /etc/apt/sources.list /etc/apt/sources.list~ \
    && sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list \
    && apt-get update \
    && apt-get -y build-dep r-base-dev \
    && wget -c https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz \
    && tar -xf R-4.3.0.tar.gz \
    && cd R-4.3.0 \
    && ./configure \
    && make -j9 \
    && make install \
    && cd ../ \
    && Rscript -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('devtools', version = '2.4.5', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('argparse', version = '2.2.2', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('BiocManager', version = '1.30.18', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('box', version = '1.1.3', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('data.table', version = '1.14.0', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('plyr', version = '1.8.8', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('poilog', version = '0.4.2', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('seqinr', version = '4.2-30', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "BiocManager::install(c('biomaRt'), version = '3.17')" \
    && Rscript -e "BiocManager::install(c('GenomicFeatures', 'GenomicRanges'), version = '3.17')" \
    && Rscript -e "BiocManager::install(c('rtracklayer'), version = '3.17')" \
    && Rscript -e "BiocManager::install(c('TxDb.Hsapiens.UCSC.hg19.knownGene', 'TxDb.Hsapiens.UCSC.hg38.knownGene'), version = '3.17')" \
    && git clone https://github.com/marialitovchenko/dNdScv_0.1.0_indel.git \
    && Rscript -e "install.packages('dNdScv_0.1.0_indel', repos = NULL, type = 'source')"