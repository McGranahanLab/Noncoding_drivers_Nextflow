#!/bin/bash

# To produce gdc_sample_sheet.2023-11-15.tsv and 
# gdc_manifest_20231115_125954.txt, go to 
# https://portal.gdc.cancer.gov/repository and select data based on following
# filters: 
# files.access in ["open"] and 
# files.analysis.workflow_type in ["ASCAT3","Aliquot Ensemble Somatic Variant Merging and Masking","STAR - Counts"] and
# files.cases.primary_site in ["bronchus and lung"] and 
# files.data_type in ["Allele-specific Copy Number Segment","Gene Expression Quantification","Masked Somatic Mutation"]

Rscript --vanilla select_tumors_from_TCGA.R \
        --sample_sheet gdc_sample_sheet.2023-11-15.tsv \
        --manifest gdc_manifest_20231115_125954.txt \
        --sample_type "Primary Tumor" \
        --data_type "Allele-specific Copy Number Segment" "Masked Somatic Mutation" "Gene Expression Quantification" \
        --n_samples 100 \
        --output_sample_sheet gdc_sample_sheet_test_cancer_drivers_pipeline.tsv \
        --output_manifest gdc_manifest_test_cancer_drivers_pipeline.txt

# download executable of GDC Data Transfer Tool
wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
unzip gdc-client_v1.6.1_Ubuntu_x64.zip
rm gdc-client_v1.6.1_Ubuntu_x64.zip

# download the data for selected samples
mkdir test_dataset/
./gdc-client download -m gdc_manifest_test_cancer_drivers_pipeline.txt -d test_dataset/
