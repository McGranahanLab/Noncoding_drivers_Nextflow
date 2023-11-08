#!/usr/bin/env Rscript
# FILE: preprocess_cgc.R ------------------------------------------------------
#
# DESCRIPTION: A script to process table with known cancer associated genes
#              downloaded from CGC (https://cancer.sanger.ac.uk/cosmic)
#
# USAGE: Rscript --vanilla preprocess_cgc.R 
#
# OPTIONS: Please run to see options and help:
#          Rscript --vanilla preprocess_cgc.R -h
#
# REQUIREMENTS: argparse, data.table
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  08.11.2023
# REVISION: 08.11.2023

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

options(scipen = 999)

# Test arguments --------------------------------------------------------------
args <- list(input = 'data/assets_raw/COSMIC_GeneCensus_all_22102020.csv',
             output = 'data/assets/cgc_knownCancerGenes.csv')

# Read CGC & extract ensembl gene id ------------------------------------------
CGC <- fread(args$input, header = T, stringsAsFactors = F, 
             select = c('Gene Symbol', 'Hallmark', 'Role in Cancer',
                        'Synonyms', 'Tumour Types(Somatic)'), 
             col.names = c('gene_name', 'is_hallmark', 'known_cancer_biotype',
                           'gene_id', 'known_in_tumor_subtype'))

# extract ensembl gene id
CGC[, gene_id := gsub('.*ENS', 'ENS', gene_id)]
CGC[, gene_id := gsub('[.].*', '', gene_id)]
CGC <- CGC[!duplicated(CGC)]

# in RefSeq GTF we unfortunately do not have ensembl gene ids
# CGC[, gene_id := gene_name]

CGCnoID <- CGC[is.na(gene_id) | gene_id == '']
message('[', Sys.time(), '] Could not find gene_id for ', nrow(CGCnoID), ' ',
        'CGC genes: ', paste0(CGCnoID$gene_name, collapse = ', '), '. Please ',
        'consider manual annotation.')

CGC[, is_known_cancer := T]
if (!is.null(args$known_db_to_use)) {
  genes_in_db <- sapply(known_cancer$db_name, 
                        function(x) any(unlist(strsplit(x, ', ')) %in% 
                                          args$known_db_to_use))
  known_cancer <- known_cancer[genes_in_db]
} 

CGC[, is_hallmark := ifelse(is_hallmark == 'Yes', T, F)]

# Write to file ---------------------------------------------------------------
CGC[, db_name := 'CGC']
setcolorder(CGC, c('gene_id', 'gene_name', 'is_known_cancer', 'is_hallmark',
                   'known_in_tumor_subtype', 'known_cancer_biotype',
                   'db_name'))
write.table(CGC, args$output, append = F, quote = F, sep = '\t', row.names = F,
            col.names = T)