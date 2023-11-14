#!/usr/bin/env Rscript
# FILE: preprocess_GTEx.R -----------------------------------------------------
#
# DESCRIPTION: a script to preprocess GTEx data to be used to annotated 
#              potential driver genes with expression.
#
# USAGE: Rscript --vanilla preprocess_GTEx.R
#
# OPTIONS:  
# EXAMPLE:
#
# REQUIREMENTS:  argparse, data.table 
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  03.11.2022
# REVISION: 14.11.2023

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(tools)))

options(scipen = 999)

# Functions -------------------------------------------------------------------
#' parseBED12regName
#' @description Parces name string(s) from BED12 files containing information
#' about genomic regions of interest. Usually those files are used with 
#' driverpower.
#' @author Maria Litovchenko
#' @param nameStr vector of strings
#' @return data table with number of rows = length of nameStr with columns 
#'         target_genome_version, gr_id, gene_id, gene_name, name
#' @export
parseBED12regName <- function(nameStr, sepStr = '--') {
  result <- lapply(nameStr, strsplit, sepStr)
  result <- do.call(rbind, lapply(result, function(x) x[[1]]))
  result <- as.data.table(result)
  result[, name := nameStr]
  colnames(result) <- c('target_genome_version', 'gr_id', 'gene_id',
                        'gene_name', 'name')
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'preprocess_GTEx.R')

gtexInvHelp <- paste('Inventory table matching tumor subtypes used in the',
                     'de-novo driver discovery and tumor subtypes in GTEx',
                     'table. Essential columns: tumor_subtype and',
                     'gtex_tumor_subtype.')
parser$add_argument('-i', '--inventory_gtex', required = T, type = 'character', 
                    help = gtexInvHelp)

gtexHelp <- paste('Path to GTEx expression table downloaded from ',
                  'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz')
parser$add_argument("-g", "--gtex", required = T, type = 'character', 
                    help = gtexHelp)

grHelp <- paste('Path to BED12 files for all tumor subtypes involved in',
                'de-novo driver discovery.')
parser$add_argument("-gr", "--genomic_regions", required = T, nargs = '*',
                    type = 'character', help = gtexHelp)

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

# Test input arguments --------------------------------------------------------
args <- list('inventory_gtex' = '../data/inventory/inventory_gtex.csv',
             'gtex' = '../data/assets_raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',
             'genomic_regions' = list('../TEST/inputs/inputGR-LUAD-hg19.bed',
                                      '../TEST/inputs/inputGR-LUSC-hg19.bed'),
             'use_ensembl' = T, 
             'output' = '../data/assets/GTEx_expression.csv')

# Read GTEx inventory ---------------------------------------------------------
inventory_gtex <- fread(args$inventory_gtex, header = T, stringsAsFactors = F,
                        select = c('tumor_subtype', 'gtex_tumor_subtype'))

# Read GTEx data --------------------------------------------------------------
gtex <- fread(args$gtex, header = T, stringsAsFactors = F, 
              select = c('Name', 'Description', 
                         unique(inventory_gtex$gtex_tumor_subtype)))
setnames(gtex, c("Name", "Description"), c('gene_id', 'gene_name'))
if (identical(colnames(gtex), c('gene_id', 'gene_name'))) {
  gtexTissues <- fread(args$gtex, header = T, stringsAsFactors = F, nrows = 10)
  gtexTissues <- setdiff(c("Name", "Description"), colnames(gtexTissues))
  stop('[', Sys.time(), '] No expression data was read from GTEx table. ',
       'Tissues from GTEx inventory (', 
       paste0(setdiff(unique(inventory_gtex$gtex_tumor_subtype), gtexTissues),
              collapse = ', '), ') are not found in GTEx table: ', 
       paste0(gtexTissues, collapse = ', '))
}
gtex[, gene_id := gsub('[.].*', '', gene_id)]

if (!args$use_ensembl) {
  gtex[, gene_id := gene_name]
}

# Read genomic ranges files & extract gene names to gene ids ------------------
# identify unique genomic regions files
grFilesMd5sum <- sapply(args$genomic_regions, md5sum)
if (length(grFilesMd5sum) != grFilesMd5sum[!duplicated(grFilesMd5sum)]) {
  uniqeGRfiles <- names(grFilesMd5sum[!duplicated(grFilesMd5sum)])
  message('[', Sys.time(), '] Unique genomic regions files: ',
          paste0(uniqeGRfiles, collapse = ', '), '. Other files will not be ',
          'read.')
  args$genomic_regions <- uniqeGRfiles
}

GR <- lapply(c(args$genomic_regions, args$genomic_regions), import)
GR <- unique(unlist(lapply(GR, function(x) unique(x$name))))
GR <- parseBED12regName(GR)
GR <- unique(GR[,.(gene_id, gene_name)])

# Modify gene_name and gene_id of GTEx data to match GRs of interest ----------
# table showing if gene id or gene name from gtex was found in GR
foundInMap <- data.table(found_id = gtex$gene_id %in% GR$gene_id, 
                         found_name = gtex$gene_name %in% GR$gene_name)
# we will only modify entries in gtex if just either gene name or gene id was
# found
changeName <- which(rowSums(foundInMap) == 1 & foundInMap$found_id == T)
setkey(GR, 'gene_id')
updatedNames <- GR[gtex$gene_id[changeName]]
updatedNames <- updatedNames[,.SD[1], by = gene_id]
setkey(updatedNames, 'gene_id')
gtex$gene_name[changeName] <- updatedNames[gtex$gene_id[changeName]]$gene_name

changeID <- which(rowSums(foundInMap) == 1 & foundInMap$found_name == T)
setkey(GR, 'gene_name')
updatedIDs <- GR[gtex$gene_name[changeID]]
updatedIDs <- updatedIDs[,.SD[1], by = gene_name]
setkey(updatedIDs, 'gene_id')
gtex$gene_id[changeID] <- updatedIDs[gtex$gene_id[changeID]]$gene_id

# Write to file ---------------------------------------------------------------
setcolorder(gtex, c('gene_id', 'gene_name'))
write.table(gtex, args$output, append = F, quote = F, sep = '\t',
            row.names = F, col.names = T)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')