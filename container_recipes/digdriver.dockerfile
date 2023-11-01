From ubuntu:22.04
MAINTAINER Maria Litovchenko m.litovchenko@ucl.ac.uk

ARG DEBIAN_FRONTEND=noninteractive 
ARG DEBCONF_NONINTERACTIVE_SEEN=true

# python 3.10 is default for ubuntu:22.04
# digdriver needs bedtools
# install R v 4.3.0
# BioConductor packages, version is controlled though BiocManager version (3.17)
# BiocFileCache v2.8.0, biomaRt v2.56.1, Biostrings v2.68.1
# GenomeInfoDb v1.16.3, GenomicFeatures v1.52.2, GenomicRanges v1.52.0 
# rtracklayer v1.60.1
# installation of dNdScv package from github
# The DIGDriver code us changed during installation to be able to give our own
# rda object
RUN apt-get update && apt-get install -y git autoconf gcc git make ssh \
                wget vim build-essential software-properties-common \
                ca-certificates libssl-dev libcurl4-gnutls-dev \
                libatlas-base-dev libffi-dev libxml2-dev \
                libncurses5-dev libsm6 libxrender1 libfontconfig1 \
                libxt6 libtcl8.6 libtk8.6 glibc-source libharfbuzz-dev \
                libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev \
                libjpeg-dev fontconfig cmake libfontconfig1-dev xclip libtiff5 \
    && apt-get install -y python3-pip \
    && apt-get install -y bedtools=2.30.0+dfsg-2ubuntu0.1 \
    && pip3 install h5py==3.9.0 numpy==1.24.4 pandas==2.0.3 pybbi==0.3.5 \
                 pybedtools==0.9.1 pysam==0.21.0 scipy==1.10.1 \
                 statsmodels==0.14.0 tables==3.8.0 \
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
    && apt-get install -y r-cran-devtools \
    && Rscript -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('devtools', version = '2.4.5', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('argparse', version = '2.2.2', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('BiocManager', version = '1.30.22', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('data.table', version = '1.14.0', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('httr', version = '1.4.7', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('plyr', version = '1.8.8', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('poilog', version = '0.4.2', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('RCurl', version = '1.98-1.12', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('seqinr', version = '4.2-30', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('openssl', version = '2.1.1', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('XML', version = '3.99-0.14', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_version('xml2', version = '1.3.5', repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "BiocManager::install(version = '3.17', ask = F)" \
    && Rscript -e "BiocManager::install(c('BiocFileCache', 'biomaRt', 'Biostrings'), version = '3.17')" \
    && Rscript -e "BiocManager::install(c('GenomeInfoDb', 'GenomicFeatures', 'GenomicRanges'), version = '3.17')" \
    && Rscript -e "BiocManager::install(c('rtracklayer'), version = '3.17')" \
    && Rscript -e "Sys.setenv(TAR = '/bin/tar');devtools::install_github('im3sanger/dndscv@1b39f8267d320a9fbbd7547c7e71e2d5e133ba3e')" \
    && git clone https://github.com/maxwellsh/DIGDriver \
    && cd DIGDriver/ \
    && git reset --hard 5bb565a \
    && cd scripts/ \
    && sed -i "s/# parse_d.add_argument('refdb'/parse_d.add_argument('refdb'/g" DigPreprocess.py \
    && sed -i "s/args.fmut, refdb,/args.fmut, args.refdb,/g" DigPreprocess.py \
    && cd ../ \
    && python3 setup.py install

ENTRYPOINT ["python3"]
