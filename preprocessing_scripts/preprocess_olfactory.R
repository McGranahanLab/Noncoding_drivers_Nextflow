#!/usr/bin/env Rscript
# FILE: preprocess_olfactory.R ------------------------------------------------
#
# DESCRIPTION: A script to process table with olfactory genes downloaded from
#              Barnes et al 2020 
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7055050/bin/12864_2020_6583_MOESM2_ESM.xlsx)
#
# USAGE: Rscript --vanilla preprocess_olfactory.R 
#
# OPTIONS: Please run to see options and help:
#          Rscript --vanilla preprocess_olfactory.R -h
#
# REQUIREMENTS: argparse, data.table, readxl
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  14.11.2023
# REVISION: 14.11.2023

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))

options(scipen = 999)

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'preprocess_olfactory.R')

inputHelp <- paste('Path to a xls table with olfactory genes published by',
                   'Barnes et al 2020. Navigate to',
                   'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7055050 and',
                   'download Additional file 2 Supplementary File 1')
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

# Test inputs -----------------------------------------------------------------
args <- list(input = '../data/assets_raw/12864_2020_6583_MOESM2_ESM.xlsx',
             use_ensembl = T,
             output = '../data/assets/olfactory_barnes_2020.csv')

# Read in table ---------------------------------------------------------------
result <- as.data.table(read_excel(args$input))
result <- result[, c('Gene symbol', 'Ensembl gene ID'), with = F]
colnames(result) <- c('gene_name', 'gene_id')
if (!args$use_ensembl) {
  result[, gene_id := gene_name]
}
result <- unique(result)

# Write to output -------------------------------------------------------------
write.table(result, args$output, col.names = T, row.names = F, quote = F, 
            sep = '\t')

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')