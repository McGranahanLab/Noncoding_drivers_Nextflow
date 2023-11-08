#!/usr/bin/env Rscript
# FILE: 2c_write_mutations_to_dndscv.R -------------------------------------
#
# DESCRIPTION: Formats MAF file containing mutations for one tumor subtype to
#              dNdScv input format.
#
# USAGE: Rscript --vanilla 2c_write_mutations_to_dndscv.R \
#                --maf [path to MAF file with all mutations for that subtype] \
#                --cancer_subtype [cancer subtype of interest] \
#                --target_genome_version hg19 \
#                --output [path to folder to write files to] \
#
# OPTIONS: Run 
#          Rscript --vanilla  2c_write_mutations_to_dndscv.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: 
# BUGS: --
# NOTES:
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  29.09.2023
# REVISION: 29.09.2023

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
options(scipen = 999)

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'write_mutations_to_dndscv.R')

mafHelp <- 'A path to MAF file with all mutations for that cancer subtype'
parser$add_argument("-m", "--maf", required = T, type = 'character',
                    default = NULL, help = mafHelp)

subtypeHelp <- paste('A cancer subtype to select from patientsInv table. Only',
                     'mutations from patients with that cancer type will be',
                     'selected. In case an analysis of several cancer types',
                     'needed to be performed please run this script ',
                     'separetedly for each cancer type.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    default = NULL, help = subtypeHelp)

targetGenomeHelp <- paste("Genome version, i.e. hg19, to which input variants",
                          "files for software should be brought.",
                          "Default: hg19.")
parser$add_argument("-g", "--target_genome_version", required = F,
                    default = 'hg19', type = 'character',
                    help = targetGenomeHelp)

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output file")

args <- parser$parse_args()
check_input_arguments_preproc(args, outputType = 'file')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test input args -------------------------------------------------------------
#args <- list(maf = 'work/d6/cf644e2aa68b76622493133ccde6f3/inputs/inputMutations-LUSC-hg19.maf',
#             cancer_subtype = 'LUSC', target_genome_version = 'hg19',
#             output = 'test')

# Software specific parameters ------------------------------------------------
colsToGet <- c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 
               'Reference_Allele', 'Tumor_Seq_Allele2')
colsOutNames <- c("sampleID", "chr", "pos", "ref", "mut")
printColnames <- T

# READ MAF file ---------------------------------------------------------------
maf <- fread(args$maf, header = T, stringsAsFactors = F)

message('[', Sys.time(), '] Formatting mutations for dNdScv, ', 
        args$cancer_subtype)

# PROCESS MAF file to suit the software ---------------------------------------
maf <- maf[, intersect(colsToGet, colnames(maf)), with = F]
maf <- maf[order(Chromosome, Start_Position)]

setnames(maf, colsToGet, colsOutNames)
setcolorder(maf, colsOutNames)

write.table(maf, args$output, col.names = printColnames, row.names = F, 
            sep = '\t', quote = F)
message('[', Sys.time(), '] Wrote mutations for dNdScv, ', 
        args$cancer_subtype)

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')