#!/usr/bin/env Rscript
# FILE: preprocess_gene_name_synonyms.R ---------------------------------------
#
# DESCRIPTION: A script to process table with gene name synonyms downloaded 
#              from EBI: 
# http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/
#
# USAGE: Rscript --vanilla preprocess_gene_name_synonyms.R 
#
# OPTIONS: Please run to see options and help:
#          Rscript --vanilla preprocess_gene_name_synonyms.R -h
#
# REQUIREMENTS: argparse, data.table
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  14.11.2023
# REVISION: 14.11.2023

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

options(scipen = 999)

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'preprocess_gene_name_synonyms.R')

inputHelp <- paste('Path to a file with gene name synonyms downloaded from ',
                   'EBI: http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/',
                   'Essential columns: symbol, alias_symbol, prev_symbol')
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

# Test inputs -----------------------------------------------------------------
# args <- list(input = '../data/assets_raw/hgnc_complete_set_2023-10-01.txt',
#              ouput = '../data/assets/hgnc_complete_set_processed.csv')

# Read in table ---------------------------------------------------------------
hugo <- fread(args$input, header = T, stringsAsFactors = F, 
              select = c('symbol', 'alias_symbol', 'prev_symbol'))
hugo <- apply(hugo, 1,
              function(x) unlist(c(x['symbol'], 
                                   strsplit(x[c('alias_symbol', 
                                                'prev_symbol')], '[|]'))))
hugo <- lapply(1:length(hugo), 
               function(x) data.table(idx = x, gene_name = hugo[[x]]))
hugo <- do.call(rbind, hugo)

# Write to file ---------------------------------------------------------------
write.table(hugo, args$ouput, col.names = T, row.names = F, quote = F, 
            sep = '\t')