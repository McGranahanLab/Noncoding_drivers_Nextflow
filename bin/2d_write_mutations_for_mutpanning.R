#!/usr/bin/env Rscript
# FILE: 2a_write_mutations_to_mutpanning.R ------------------------------------
#
# DESCRIPTION: Formats MAF file containing mutations for one tumor subtype to
#              MutPanning input format.
#
# USAGE: Rscript --vanilla 2a_write_mutations_to_mutpanning.R \
#                --maf [path to MAF file with all mutations for that subtype] \
#                --cancer_subtype [cancer subtype of interest] \
#                --target_genome_version hg19 \
#                --output [path to folder to write files to] \
#
# OPTIONS: Run 
#          Rscript --vanilla  2a_write_mutations_to_mutpanning.R -h
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

box::use(./custom_functions[...])

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
options(scipen = 999)

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'write_mutations_to_mutpanning.R')

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
                    help = "Path to the output folder")

args <- parser$parse_args()
check_input_arguments(args, outputType = 'folder')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test input args -------------------------------------------------------------
#args <- list(maf = 'work/d6/cf644e2aa68b76622493133ccde6f3/inputs/inputMutations-LUSC-hg19.maf',
#             cancer_subtype = 'LUSC', target_genome_version = 'hg19',
#             output = 'test')

# Software specific parameters ------------------------------------------------
colsToGet <- c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
               'Strand', 'Variant_Classification', 'Variant_Type',
               'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 
               'Tumor_Sample_Barcode')
colsOutNames <- c('Hugo_Symbol', 'Chromosome', 'Start_Position',
                  'End_Position', 'Strand', 'Variant_Classification', 
                  'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 
                  'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode')
printColnames <- T

# READ MAF file ---------------------------------------------------------------
maf <- fread(args$maf, header = T, stringsAsFactors = F)

message('[', Sys.time(), '] Formatting mutations for MutPanning, ', 
        args$cancer_subtype)
outfile <- paste0(args$output, '/mutpanning-inputMutations-',
                  args$cancer_subtype, '-', args$target_genome_version, '.csv')

# PROCESS MAF file to suit the software ---------------------------------------
maf <- maf[, intersect(colsToGet, colnames(maf)), with = F]
maf <- maf[order(Chromosome, Start_Position)]

message('[', Sys.time(), '] MutPanning requires chromosomal names in NCBI ',
        'format. Changed chromosomal names to NCBI format.')
maf[, Chromosome := gsub('chr', '', Chromosome)]

setnames(maf, colsToGet, colsOutNames)
setcolorder(maf, colsOutNames)

write.table(maf, outfile, col.names = printColnames, row.names = F, sep = '\t',
            quote = F)
message('[', Sys.time(), '] Wrote mutations for MutPanning, ', 
        args$cancer_subtype)

# mutpanning also requires txt file with samples, similar to patientsInv
mutpanInvent <- maf[, c('Tumor_Sample_Barcode'), with = F]
mutpanInvent <- unique(mutpanInvent)
colnames(mutpanInvent) <- c('ID')
mutpanInvent[, Sample := ID]
mutpanInvent[, Subtype := args$cancer_subtype]
mutpanInvent[, ConfidenceLevel := 'custom']
mutpanInvent[, Study := 'custom']
mutpanInvent[, Cohort := 'custom']
# output to file
outfile <- paste0('mutpanning-patientsInv-', args$cancer_subtype, '-', 
                  args$target_genome_version, '.csv')
write.table(mutpanInvent, append = F, quote = F, row.names = F, outfile,
            col.names = T, sep = '\t')

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')