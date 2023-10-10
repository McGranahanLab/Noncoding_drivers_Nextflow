#!/bin/bash

#BSUB -L /bin/bash
#BSUB -q inter
#BSUB -P re_gecip_cancer_lung
#BSUB -e nextflow_YYYY-MM-DD.err
#BSUB -o nextflow_YYYY-MM-DD.out
#BSUB -J nextflow
#BSUB -R "rusage[mem=8000] span[hosts=1]"
#BSUB -n 1
#BSUB -cwd "."

module load lang/Java/19.0
module load tools/singularity/3.8.3

./nextflow-23.04.2-all run main.nf -profile genomicsEngland