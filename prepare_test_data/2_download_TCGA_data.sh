#!/bin/bash

# To produce gdc_sample_sheet.2023-11-15.tsv and 
# gdc_manifest_20231115_125954.txt, go to 
# https://portal.gdc.cancer.gov/repository and select data based on following
# filters: 
# files.access in ["open"] and 
# files.analysis.workflow_type in ["ASCAT3","Aliquot Ensemble Somatic Variant Merging and Masking","STAR - Counts"] and
# files.cases.primary_site in ["bronchus and lung"] and 
# files.data_type in ["Allele-specific Copy Number Segment","Gene Expression Quantification","Masked Somatic Mutation"]

SAMPLE_SHEET='gdc_sample_sheet_test_cancer_drivers_pipeline.tsv'
MANIFEST='gdc_manifest_test_cancer_drivers_pipeline.txt'
TEST_DATA_FOLDER='test_dataset/'
FPKM_TABLE_PER_SAMPLE='../data/assets_raw/TCGA_FPKM_per_sample.csv'
FPKM_TABLE_PER_TUMOR_SUBTYPE='../data/assets_raw/TCGA_FPKM_per_tumor_subtype.csv'

Rscript --vanilla 0_select_tumors_from_TCGA.R \
        --sample_sheet gdc_sample_sheet.2023-11-15.tsv \
        --manifest gdc_manifest_20231115_125954.txt \
        --sample_type "Primary Tumor" \
        --data_type "Allele-specific Copy Number Segment" "Masked Somatic Mutation" "Gene Expression Quantification" \
        --n_samples 100 \
        --output_sample_sheet $SAMPLE_SHEET \
        --output_manifest $MANIFEST

# download executable of GDC Data Transfer Tool
wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
unzip gdc-client_v1.6.1_Ubuntu_x64.zip
rm gdc-client_v1.6.1_Ubuntu_x64.zip

# download the data for selected samples
mkdir $TEST_DATA_FOLDER
./gdc-client download -m $MANIFEST -d $TEST_DATA_FOLDER

# produce FPKM tables: one containing expression data for every sample and one 
# containing expression data per tumor subtype (median)
Rscript --vanilla 1_create_FPKM_table_TCGA.R --sample_sheet $SAMPLE_SHEET \
        --folder $TEST_DATA_FOLDER \
        --data_type "Gene Expression Quantification" \
        --tumor_subtype_column 'Project ID' --do_aggregate T \
        --aggregate_name PANLUNG  --output_per_sample $FPKM_TABLE_PER_SAMPLE \
        --output_per_tumor_subtype $FPKM_TABLE_PER_TUMOR_SUBTYPE