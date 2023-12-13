#!/usr/bin/env Rscript
# FILE: prepare_CN_on_GEL.R ---------------------------------------------------
#
# DESCRIPTION: creates tables with copy number segments from ASCAT, one per
# participant id.
#
# USAGE: run in R studio
# OPTIONS:
# EXAMPLE: 
# REQUIREMENTS: R 4.0.2, data.table
#               
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  05.12.2023
# REVISION: 05.12.2023

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

# Inputs ----------------------------------------------------------------------
ascatPath <- '/re_gecip/cancer_lung/analysisResults/GEL/release_v8/3.chromosomeAberrations/K.ASCAT/hg38_withSVs/combined_ascatTable_hg19_liftedOver_by_Maria.rds'
outputDir <- '/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/data/ascat/'
dir.create(outputDir, recursive = T)

# Read ASCAT table ------------------------------------------------------------
ascat <- as.data.table(readRDS(ascatPath))
cols_to_get <- c('patient', 'chr', 'startpos', 'endpos', 'nMajor', 'nMinor', 
                 'Ploidy')
ascat <- ascat[, intersect(colnames(ascat), cols_to_get), with = F]
setnames(ascat, c('patient', 'startpos', 'endpos'),
         c('participant_id', 'start', 'end'))
ascat <- split(ascat, by = 'participant_id')

# Output to files -------------------------------------------------------------
lapply(names(ascat), 
       function(x) write.table(ascat[[x]], 
                               paste0(outputDir, '/', x, '.ascat.csv'), 
                               append = F, quote = F, sep = '\t', 
                               row.names = F, col.names = T))