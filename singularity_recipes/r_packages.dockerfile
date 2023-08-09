FROM ubuntu:22.04
MAINTAINER Maria Litovchenko m.litovchenko@ucl.ac.uk

ARG DEBIAN_FRONTEND=noninteractive 
ARG DEBCONF_NONINTERACTIVE_SEEN=true

# install R v 4.3.0, CRAN packages & BioConductor packages
# For BioConductor packages, version is controlled though BiocManager 
# version (3.17): GenomicFeatures v1.52.2, GenomicRanges v1.52.0 
# org.Hs.eg.db v3.17.0, BSgenome.Hsapiens.UCSC.hg19 v1.4.3
# TxDb.Hsapiens.UCSC.hg19.knownGene v3.2.2
# VariantAnnotation v1.46.0, EmpiricalBrownsMethod v1.28.0
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
    && Rscript -e "devtools::install_version('dplyr', version = '1.1.3', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('plyr', version = '1.8.8', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('poolr', version = '1.1-1', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('reshape2', version = '1.4.4', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('R.utils', version = '2.12.2', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('strex', version = '1.6.0', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "BiocManager::install(c('GenomicFeatures', 'GenomicRanges'), version = '3.17')" \
    && Rscript -e "BiocManager::install(c('maftools', 'plyranges', 'rtracklayer'), version = '3.17')" \
    && Rscript -e "BiocManager::install(c('org.Hs.eg.db', 'BSgenome.Hsapiens.UCSC.hg19', 'TxDb.Hsapiens.UCSC.hg19.knownGene'), version = '3.17')" \
    && Rscript -e "BiocManager::install(c('VariantAnnotation', 'EmpiricalBrownsMethod'), version = '3.17')"  \
    && Rscript -e "Sys.setenv(TAR = '/bin/tar');devtools::install_github('im3sanger/dndscv@1b39f8267d320a9fbbd7547c7e71e2d5e133ba3e')" \
    && Rscript -e "devtools::install_github('NKI-CCB/DISCOVER/R@r_v0.9.4')"