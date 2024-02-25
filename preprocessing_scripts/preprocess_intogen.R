#!/usr/bin/env Rscript
# FILE: preprocess_intogen.R --------------------------------------------------
#
# DESCRIPTION: A script to process table with genes detected in IntOgen study.
#              Downloaded from 
#              https://www.intogen.org/download?file=IntOGen-Drivers-20230531.zip
#              unzip download\?file\=IntOGen-Drivers-20230531.zip
#              rm -r 2023-05-31_IntOGen-Drivers \
#                    download\?file\=IntOGen-Drivers-20230531.zip
#              mv 2023-05-31_IntOGen-Drivers/Compendium_Cancer_Genes.tsv . 
#
# USAGE: Rscript --vanilla preprocess_intogen.R 
#
# OPTIONS: Please run to see options and help:
#          Rscript --vanilla preprocess_intogen.R -h
#
# REQUIREMENTS: argparse, data.table
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  08.11.2023
# REVISION: 25.02.2023

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

options(scipen = 999)

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'preprocess_intogen.R')

inputHelp <- paste('Path to a file with genes detected as cancer associated in',
                   'IntOgen study. Downloaded from',
                   'https://www.intogen.org/download?file=IntOGen-Drivers-20230531.zip',
                   'and extracted file Compendium_Cancer_Genes.tsv. Essential',
                   'columns: SYMBOL, CANCER_TYPE.')
parser$add_argument("-i", "--input", required = T, type = 'character', 
                    help = inputHelp)

outputHelp <- 'Path to output file'
parser$add_argument("-o", "--output", required = T, type = 'character', 
                    help = inputHelp)

args <- parser$parse_args()
timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')

if (!dir.exists(dirname(args$ouput))) {
  message('[', Sys.time(), '] Directory does not exist: ', dirname(args$ouput),
          '. Will create it.')
  dir.create(dirname(args$ouput), recursive = T)
}

# Test arguments --------------------------------------------------------------
# args <- list(input = '../data/assets_raw/Compendium_Cancer_Genes.tsv',
#               output = '../data/assets/intogene_detectedCancerGenes.csv')

# Read intogene ---------------------------------------------------------------
intogen <- fread(args$input,  header = T, stringsAsFactors = F,
                 select = c("SYMBOL", "CANCER_TYPE"), 
                 col.names = c('gene_name', 'known_in_tumor_subtype'))
intogen <- unique(intogen)

# Write to file ---------------------------------------------------------------
setcolorder(intogen, c('gene_name', 'known_in_tumor_subtype'))
write.table(intogen, args$output, append = F, quote = F, sep = '\t',
            row.names = F, col.names = T)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')