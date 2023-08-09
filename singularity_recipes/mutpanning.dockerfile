From ubuntu:18.04
MAINTAINER Maria Litovchenko m.litovchenko@ucl.ac.uk

RUN apt-get update && apt-get upgrade -y && \
    apt-get install unzip && \
    apt-get install -y wget openjdk-8-jdk python  && \
    cd bin/  && \
    wget https://github.com/vanallenlab/MutPanningV2/archive/master.zip  && \
    unzip master.zip && \
    rm master.zip && \
    mv MutPanningV2-master MutPanningV2 && \
    cd MutPanningV2/ && \
    wget https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.zip && \
    unzip Hg19.zip && \
    rm Hg19.zip && \
    cd ../ && \
    ff="MutPanningV2/"  && \
    javac -classpath $ff/commons-math3-3.6.1.jar:$ff/jdistlib-0.4.5-bin.jar $ff/*.java && \
    chmod -R 775 MutPanningV2/Hg19 