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
# REVISION: 13.11.2023

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

options(scipen = 999)

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'preprocess_cgc.R')

inputHelp <- paste('Path to a file with known cancer genes from CGC.',
                   'Navigate to https://cancer.sanger.ac.uk/cosmic/download,',
                   'login and download file from "Cancer Gene Census"',
                   'section. Unzip it and find Cosmic_CancerGeneCensus_vXX_GRChXX.tsv',
                   'file. Genome version does not matter. Tested on version',
                   '98. Essential columns: GENE_SYMBOL, SOMATIC,',
                   'ROLE_IN_CANCER, SYNONYMS, TUMOUR_TYPES_SOMATIC')
parser$add_argument("-i", "--input", required = T, type = 'character', 
                    help = inputHelp)

ensemblHelp <- paste('Whether or not ensembl gene IDs should be used in',
                     'gene_id column. If set to FALSE, then gene_id will be',
                     'equal to gene_name')
parser$add_argument("-e", "--use_ensembl", required = F, type = 'character', 
                    help = ensemblHelp, default = T, choices = c('T', 'F'))

outputHelp <- 'Path to output file'
parser$add_argument("-o", "--output", required = T, type = 'character', 
                    help = inputHelp)

args <- parser$parse_args()
args$use_ensembl <- as.logical(args$use_ensembl)
timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')

if (!dir.exists(dirname(args$ouput))) {
  message('[', Sys.time(), '] Directory does not exist: ', dirname(args$ouput),
          '. Will create it.')
  dir.create(dirname(args$ouput), recursive = T)
}

# Test arguments --------------------------------------------------------------
# args <- list(input = '../data/assets_raw/Cosmic_CancerGeneCensus_v98_GRCh37.tsv',
#              output = '../data/assets/cgc_knownCancerGenes.csv',
#              use_ensembl = F)

# Read CGC & extract ensembl gene id ------------------------------------------
CGC <- fread(args$input, header = T, stringsAsFactors = F, 
             select = c('GENE_SYMBOL', 'SOMATIC', 'ROLE_IN_CANCER', 'SYNONYMS',
                        'TUMOUR_TYPES_SOMATIC'), 
             col.names = c('gene_name', 'somatic', 'known_cancer_biotype',
                           'gene_id', 'known_in_tumor_subtype'))

# select only somatic genes
CGC <- CGC[somatic == 'y']
CGC[, somatic := NULL]

# extract ensembl gene id
if (args$use_ensembl) {
  CGC[, gene_id := gsub('.*ENS', 'ENS', gene_id)]
  CGC[, gene_id := gsub('[.].*', '', gene_id)]
} else {
  CGC[, gene_id := gene_name]
}
CGC <- CGC[!duplicated(CGC)]
CGCnoID <- CGC[is.na(gene_id) | gene_id == '']
if (nrow(CGCnoID) != 0) {
  message('[', Sys.time(), '] Could not find gene_id for ', nrow(CGCnoID), ' ',
          'CGC genes: ', paste0(CGCnoID$gene_name, collapse = ', '), '. ',
          'Please consider manual annotation.')
}
CGC[, is_known_cancer := T]

CGC[, known_in_tumor_subtype := gsub(', ', ',', known_in_tumor_subtype)]
CGC[, known_cancer_biotype := gsub(', ', ',', known_cancer_biotype)]

# Write to file ---------------------------------------------------------------
CGC[, db_name := 'CGC']
setcolorder(CGC, c('db_name', 'gene_id', 'gene_name', 'is_known_cancer',
                   'known_in_tumor_subtype', 'known_cancer_biotype'))
write.table(CGC, args$output, append = F, quote = F, sep = '\t', row.names = F,
            col.names = T)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')