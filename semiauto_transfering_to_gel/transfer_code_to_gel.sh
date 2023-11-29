#!/bin/bash

CODE_DIR='/Users/maria/Desktop/BitBucket/Noncoding_drivers_Nextflow/'
TMP_DIR=`pwd`'/tmp_dir'
dockerfile='transfer_code_to_gel.dockerfile'
TAG='transfer_code'

mkdir -p $TMP_DIR

chmod +x zip_the_code.sh
./zip_the_code.sh $CODE_DIR $TMP_DIR


docker build . -t noncoding_driver_pipeline:$TAG \
               -f $dockerfile --no-cache --progress=plain
docker tag noncoding_driver_pipeline:$TAG \
           marialitovchenko/noncoding_driver_pipeline:$TAG
docker push marialitovchenko/noncoding_driver_pipeline:$TAG

rm -r $TMP_DIR