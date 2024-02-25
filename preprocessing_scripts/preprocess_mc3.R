#!/usr/bin/env Rscript
# FILE: preprocess_mc3.R --------------------------------------------------
#
# DESCRIPTION: A script to process table with genes detected in MC3 study.
#              Downloaded from 
#              https://ars.els-cdn.com/content/image/1-s2.0-S009286741830237X-mmc1.xlsx
#              and then S1 is extracted to 
#              1-s2.0-S009286741830237X-mmc1_S1.csv manually.
#
# USAGE: Rscript --vanilla preprocess_mc3.R 
#
# OPTIONS: Please run to see options and help:
#          Rscript --vanilla preprocess_mc3.R -h
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
parser <- ArgumentParser(prog = 'preprocess_mc3.R')

inputHelp <- paste('Path to a csv file with S1 table downloaded from',
                   'https://ars.els-cdn.com/content/image/1-s2.0-S009286741830237X-mmc1.xlsx',
                   'S1 is extracted to csv manually. Essential columns: ',
                   'Gene, Cancer.')
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
# args <- list(input = '../data/assets_raw/1-s2.0-S009286741830237X-mmc1_S1.csv',
#               output = '../data/assets/mc3_detectedCancerGenes.csv')

# Read intogene ---------------------------------------------------------------
mc3 <- fread(args$input,  header = T, stringsAsFactors = F, skip = 2,
                 select = c("Gene", "Cancer"), 
                 col.names = c('gene_name', 'known_in_tumor_subtype'))
mc3 <- unique(mc3)

# Write to file ---------------------------------------------------------------
setcolorder(mc3, c('gene_name', 'known_in_tumor_subtype'))
write.table(mc3, args$output, append = F, quote = F, sep = '\t', row.names = F,
            col.names = T)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')