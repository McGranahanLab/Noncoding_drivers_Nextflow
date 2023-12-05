#!/bin/bash

BASE_DIR=$1
OUT_DIR=$2

zip $OUT_DIR'/bin.zip'                                $BASE_DIR'/'bin/*.R
zip $OUT_DIR'/conf.zip'                               $BASE_DIR'/'conf/*.config
zip $OUT_DIR'/container_recipes.zip'                  $BASE_DIR'/'container_recipes/*.d*
zip $OUT_DIR'/launch_and_doc.zip'                     $BASE_DIR'/'DOCUMENTATION.Rmd \
                                                      $BASE_DIR'/'GEL_modification.Rmd \
                                                      $BASE_DIR'/'launch_pipeline_helix.sh \
                                                      $BASE_DIR'/'main.nf \
                                                      $BASE_DIR'/'nextflow.config
zip $OUT_DIR'/modules.zip'                            $BASE_DIR'/'modules/*.nf
zip $OUT_DIR'/preprocessing_scripts.zip'              $BASE_DIR'/'preprocessing_scripts/*.R
zip $OUT_DIR'/subworkflows.zip'                       $BASE_DIR'/'subworkflows/*.nf
zip $OUT_DIR'/preprocessing_scripts_gel_specific.zip' $BASE_DIR'/'preprocessing_scripts_gel_specific/*.R