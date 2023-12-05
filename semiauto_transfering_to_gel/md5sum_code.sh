#!/bin/bash

BASE_DIR=$1
OUT_FILE=$2

rm $OUT_FILE

md5sum-lite $BASE_DIR'/'bin/*.R \
            $BASE_DIR'/'conf/*.config \
            $BASE_DIR'/'container_recipes/*.d* \
            $BASE_DIR'/'DOCUMENTATION.Rmd \
            $BASE_DIR'/'GEL_modification.Rmd \
            $BASE_DIR'/'launch_pipeline_helix.sh \
            $BASE_DIR'/'main.nf \
            $BASE_DIR'/'nextflow.config \
            $BASE_DIR'/'modules/*.nf \
            $BASE_DIR'/'preprocessing_scripts/*.R \
            $BASE_DIR'/'subworkflows/*.nf \
            $BASE_DIR'/'preprocessing_scripts_gel_specific/*.R > $OUT_FILE