FROM ubuntu:18.04
MAINTAINER Maria Litovchenko m.litovchenko@ucl.ac.uk
 
# python 3.6 is default for ubuntu:18.04
RUN apt-get -y update \
    && apt-get -y upgrade \
    && apt-get install -y procps libcurl4-openssl-dev libssl-dev gcc libz-dev\
    && apt-get install -y python3-pip \
    && pip3 install ago==0.0.9 appdirs==1.4.4 bgcache==0.1 bgconfig==0.8 \
                 bgdata==2.0.2 bglogs==0.6 bgparsers==0.9 bgreference==0.5 \
                 bgsignature==0.2 bokeh==0.12.4 brotlipy==0.7.0 \
                 certifi==2020.6.20 cffi==1.14.0 chardet==3.0.4 click==6.7 \
                 conda==4.3.16 configobj==5.0.6 \
                 cycler==0.10.0 Cython==0.28.5 daiquiri==2.1.1. dill==0.3.2 \
                 homura==0.1.5 humanize==2.6.0 intervaltree==2.1.0 \
                 Jinja2==2.11.2 joblib==0.16.0 kiwisolver==1.2.0 \
                 MarkupSafe==1.1.1 matplotlib==2.2.3 numpy==1.15.2 \
                 oncodrivefml==2.3.0 oncodriveclustl==1.1.1 \
                 pandas==0.23.4 patsy==0.5.1 pycosat==0.6.3 pycurl==7.43.0.6 \
                 pyparsing==2.4.7 PySocks==1.7.1 pytabix==0.0.2 \
                 python-dateutil==2.8.1 python-json-logger==2.0.0 \
                 pytz==2020.1 PyYAML==5.3.1 ruamel-yaml==0.15.87 \
                 scikit-learn==0.23.2 scipy==1.1.0 six==1.15.0 \
                 sortedcontainers==2.2.2 statsmodels==0.9.0 \
                 threadpoolctl==2.1.0 tornado==6.0.4 tqdm==4.42.1 \
    && echo "export LC_ALL=C.UTF-8" >> ~/.bashrc \
    && echo "export LANG=C.UTF-8" >> ~/.bashrc

CMD source /root/.bashrc