#!/bin/bash

#BSUB -L /bin/bash
#BSUB -q inter
#BSUB -P re_gecip_cancer_lung
#BSUB -e noncoding_driver_pipeline_YYYY-MM-DD.err
#BSUB -o noncoding_driver_pipeline_YYYY-MM-DD.out
#BSUB -J nextflow_master
#BSUB -R "rusage[mem=32000] span[hosts=1]"
#BSUB -M 32000
#BSUB -n 1
#BSUB -cwd "."

module load lang/Java/19.0
module load tools/singularity/3.8.3

export TMPDIR=/re_scratch/mlitovchenko/
mkdir -p $TMPDIR

export NFX_OPTS="-Xms=8g -Xmx=24g"

./nextflow-23.04.2-all run main.nf -profile genomicsEngland \
                                   -entry CALL_DE_NOVO_CANCER_DRIVERS -resume
./nextflow-23.04.2-all run main.nf -profile genomicsEngland \
                                   -entry POSTPROCESSING -resume

#./nextflow-23.04.2-all run main.nf -entry POSTPROCESSING -profile local -resume