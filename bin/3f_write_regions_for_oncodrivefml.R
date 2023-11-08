#!/usr/bin/env Rscript
# FILE: 3f_write_regions_for_oncodrivefml.R -----------------------------------
#
# DESCRIPTION: Formats BED12 file containing genomic regions for one tumor 
#              subtype to OncodriveFML input format.
#
# USAGE: Rscript --vanilla 3f_write_regions_for_oncodrivefml.R \
#                --bed [path to BED12 file with all regions for that subtype] \
#                --cancer_subtype [cancer subtype of interest] \
#                --gr_id [list of genomic regions IDs] \
#                --output [path to folder to write files to] \
#
# OPTIONS: Run 
#          Rscript --vanilla  3f_write_regions_for_oncodrivefml.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: 
# BUGS: --
# NOTES:
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  28.09.2023
# REVISION: 28.09.2023

# Source custom functions -----------------------------------------------------
#' get_script_dir
#' @description Returns parent directory of the currently executing script
#' @author https://stackoverflow.com/questions/47044068/get-the-path-of-current-script
#' @return string, absolute path
#' @note this functions has to be present in all R scripts sourcing any other
#' script. Sourcing R scripts with use of box library is unstable then multiple
#' processes try to execute the source at the same time.
get_script_dir <- function() {
  cArgs <- tibble::enframe(commandArgs(), name = NULL)
  cArgsSplit <- tidyr::separate(cArgs, col = value, into = c("key", "value"),
                                sep = "=", fill = "right")
  cArgsFltr <- dplyr::filter(cArgsSplit, key == "--file")
  
  result <- dplyr::pull(cArgsFltr, value)
  result <- tools::file_path_as_absolute(dirname(result))
  result
}

srcDir <- get_script_dir()
# to spread out multiple processes accessing the same file
Sys.sleep(sample(1:15, 1))
source(paste0(srcDir, '/custom_functions.R'))
source(paste0(srcDir, '/custom_functions_preprocessing.R'))

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(tools)))
options(scipen = 999)

# Input arguments -------------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'write_regions_for_oncodrivefml.R')

bedHelp <- 'A path to BED12 file with all regions for that cancer subtype'
parser$add_argument("-b", "--bed", required = T, type = 'character',
                    default = NULL, help = bedHelp)

subtypeHelp <- paste('A cancer subtype to select from patientsInv table. Only',
                     'mutations from patients with that cancer type will be',
                     'selected. In case an analysis of several cancer types',
                     'needed to be performed please run this script ',
                     'separetedly for each cancer type.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    default = NULL, help = subtypeHelp)

grIdHelp <- paste('IDs of genomic regions to be created. Multiple are ',
                  'accepted')
parser$add_argument("-g", "--gr_id", required = T, type = 'character',
                    nargs = '+', default = NULL, help = grIdHelp)

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output folder")

args <- parser$parse_args()
check_input_arguments_preproc(args, outputType = 'folder')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test arguments --------------------------------------------------------------
# args <- list('bed' = '../TEST/inputs/inputGR-LUAD-hg19.bed',
#              'cancer_subtype' = 'LUAD',
#              'gr_id' = list('5primeUTR', 'CDS', 'lincRNA', 'ss'))

# Software specific parameters ------------------------------------------------
colsToGet <- c("seqnames", "start", "end", "gene_id",
               "strand", "gene_name")
colsOutNames <- c('CHROMOSOME', 'START', 'END', 'ELEMENT', 'STRAND', 'SYMBOL')
printColnames <- T

message('[', Sys.time(), '] Formatting genomic regions for OncodriveFML, ', 
        args$cancer_subtype)
outfileBase <- paste0(args$output, '/inputGR-', args$cancer_subtype, 
                      '-oncodrivefml-')

# READ BED12 file -------------------------------------------------------------
bed <- readBED12(args$bed)
# select only regions of interest
bed <- bed[bed$gr_id %in% args$gr_id]
target_genome_version <- unique(bed$target_genome_version)

# PROCESS bed to software format ----------------------------------------------
mcols(bed) <- mcols(bed)[, c('gr_id', 'gene_id', 'gene_name')]
bed <- as.data.table(bed)
bed <- split(bed, by = 'gr_id')
bed <- lapply(bed, function(x) x[, colsToGet, with = F])
bed <- lapply(bed, setnames, colsToGet, colsOutNames)

# WRITE -----------------------------------------------------------------------
lapply(names(bed),
       function(x) write.table(bed[[x]], col.names = printColnames,
                               row.names = F, quote = F, sep = '\t',
                               file = paste0(args$output, '/', outfileBase, x, 
                                             '-', target_genome_version,
                                             '.csv')))

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')