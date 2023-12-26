#!/bin/bash

CODE_DIR='/Users/Masha/Desktop/BitBucket/Noncoding_drivers_Nextflow/'
TMP_DIR=`pwd`'/tmp_dir'
MD5_FILE=$TMP_DIR'/md5sums_originals.txt'
dockerfile='transfer_code_to_gel.dockerfile'
TAG='transfer_code'

mkdir -p $TMP_DIR

chmod +x zip_the_code.sh
./zip_the_code.sh $CODE_DIR $TMP_DIR
chmod +x md5sum_code.sh
./md5sum_code.sh $CODE_DIR $MD5_FILE


docker build . -t noncoding_driver_pipeline:$TAG \
               -f $dockerfile --no-cache --progress=plain
docker tag noncoding_driver_pipeline:$TAG \
           marialitovchenko/noncoding_driver_pipeline:$TAG
docker push marialitovchenko/noncoding_driver_pipeline:$TAG

mv $TMP_DIR'/md5sums_originals.txt' .

rm -r $TMP_DIR

# to download on GEL
# module load tools/singularity/3.8.3
# singularity pull docker://marialitovchenko/noncoding_driver_pipeline:transfer_code
# BASE_DIR=`pwd`'/'
# singularity shell --bind $BASE_DIR:$BASE_DIR noncoding_driver_pipeline_transfer_code.sif 
# cp /*.zip .
# exit
# unzip *.zip