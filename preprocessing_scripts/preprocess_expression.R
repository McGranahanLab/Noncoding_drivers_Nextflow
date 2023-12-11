#!/usr/bin/env Rscript
# FILE: preprocess_expression.R -----------------------------------------------
#
# DESCRIPTION: a script to preprocess expression data (GTEx, TCGA) to be used
#              to annotated potential driver genes with expression.
#
# USAGE: Rscript --vanilla preprocess_expression.R
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
# REVISION: 16.11.2023

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(tools)))

options(scipen = 999)

# Functions -------------------------------------------------------------------
#' printArgs
#' @description Prints submitted to the script arguments as messages.
#' @author Maria Litovchenko
#' @param argsList named list representing submitted to the script arguments
#' @return void
#' @export
printArgs <- function(argsList) {
  message('[', Sys.time(), '] Submitted arguments: ')
  for (argName in names(argsList)) {
    oneArg <- argsList[[argName]]
    if (length(oneArg) > 1 & !is.null(names(oneArg))) {
      msg <- paste(apply(data.table(names(oneArg), oneArg), 1, paste, 
                         collapse = ' - '), collapse = ',')
    } else {
      msg <- paste(oneArg, collapse = ', ')
    }
    message('\t', argName, ':', msg)
  }
}

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
parser <- ArgumentParser(prog = 'preprocess_expression.R')

exprInvHelp <- paste('Inventory table matching tumor subtypes used in the',
                     'de-novo driver discovery and tumor subtypes in',
                     'expression table. Essential columns: tumor_subtype, ',
                     'expr_tumor_subtype, expr_db, min_expr')
parser$add_argument('-i', '--inventory_expression', required = T, 
                    type = 'character', help = exprInvHelp)

exprHelp <- paste('Path to expression table. GTEx expression table can be',
                  'used straight away after download from',
                  'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',
                  '. TCGA data first needs to be aggregated via median',
                  'accross the samples per tumor subtype. File must have',
                  'header. First 2 columns of the table must contain gene',
                  'IDs and gene names. Column names of the first 2 columns',
                  'will be ignored. Other columns contain normalized',
                  'expression data, one column per tumor subtype.')
parser$add_argument("-d", "--expression", required = T, type = 'character', 
                    help = exprHelp)

grHelp <- paste('Path to BED12 files containing regions of interest used in',
                'de-novo driver discovery for all tumor subtypes')
parser$add_argument("-gr", "--genomic_regions", required = T, nargs = '*',
                    type = 'character', help = grHelp)

ensemblHelp <- paste('Whether or not ensembl gene IDs should be used in',
                     'gene_id column. If set to FALSE, then gene_id will be',
                     'replaced by values in gene_name column and not used',
                     'to match expression table and genomic regions of',
                     'interest.')
parser$add_argument("-e", "--use_ensembl", required = F, type = 'character', 
                    help = ensemblHelp, default = T, choices = c('T', 'F'))

outputHelp <- 'Path to output file'
parser$add_argument("-o", "--output", required = T, type = 'character', 
                    help = outputHelp)

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
# args <- list('inventory_expression' = '../data/inventory/inventory_expression_gtex.csv',
#              'expression' = '../data/assets_raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',
#              'genomic_regions' = list('../TEST/inputs/inputGR-LUAD-hg19.bed',
#                                       '../TEST/inputs/inputGR-LUSC-hg19.bed'),
#              'use_ensembl' = T, 
#              'output' = '../data/assets/GTEx_expression.csv')
# args <- list('inventory_expression' = '../data/inventory/inventory_expression_tcga.csv',
#              'expression' = '../data/assets_raw/TCGA_FPKM_per_tumor_subtype.csv',
#              'genomic_regions' = list('../TEST/inputs/inputGR-LUAD-hg19.bed',
#                                       '../TEST/inputs/inputGR-LUSC-hg19.bed'),
#              'use_ensembl' = T, 
#              'output' = '../data/assets/TCGA_expression.csv')

# args <- list('inventory_expression' = '../data/inventory/inventory_expression_gtex.csv',
#              'expression' = '../data/assets_raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',
#              'genomic_regions' = list('../completed_runs/06_12_2023/inputs/inputGR-Adenocarcinoma-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Adenocarcinoma_met-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Adenosquamous-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Carcinoid-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Large_cell-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Mesothelioma-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Neuroendocrine_carcinoma-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Panlung-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Panlung_met-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Small_cell-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Squamous_cell-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Squamous_cell_met-hg19.bed'),
#              'use_ensembl' = T, 
#              'output' = '../data/assets/GTEx_expression.csv')
# args <- list('inventory_expression' = '../data/inventory/inventory_expression_tcga.csv',
#              'expression' = '../data/assets_raw/TCGA_FPKM_per_tumor_subtype.csv',
#              'genomic_regions' = list('../completed_runs/06_12_2023/inputs/inputGR-Adenocarcinoma-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Adenocarcinoma_met-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Adenosquamous-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Carcinoid-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Large_cell-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Mesothelioma-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Neuroendocrine_carcinoma-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Panlung-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Panlung_met-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Small_cell-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Squamous_cell-hg19.bed',
#                                       '../completed_runs/06_12_2023/inputs/inputGR-Squamous_cell_met-hg19.bed'),
#              'use_ensembl' = T, 
#              'output' = '../data/assets/TCGA_expression.csv')

# Read & check expression inventory -------------------------------------------
inventory_expr <- fread(args$inventory_expression, header = T,  
                        stringsAsFactors = F,
                        select = c('tumor_subtype', 'expr_tumor_subtype', 
                                   'expr_db', 'min_expr'))
if (ncol(inventory_expr) != 4) {
  stop('[', Sys.time(), '] ', args$inventory_expr, ' does not have all 3 ',
       'columns: tumor_subtype, expr_tumor_subtype, expr_db and min_expr.')
}
if (length(unique(inventory_expr$tumor_subtype)) != nrow(inventory_expr)) {
  stop('[', Sys.time(), '] ', args$inventory_expression, ' has duplicated ',
       'tumor_subtype.')
}
if (any(is.na(inventory_expr$expr_db))) {
  stop('[', Sys.time(), '] ', args$inventory_expression, ' has NA in ',
       'expr_db column')
}
if (length(unique(inventory_expr$expr_db)) > 1) {
  stop('[', Sys.time(), '] ', args$inventory_expression, ' has >1 unique ',
       'value in expr_db column')
}

# Read expression data --------------------------------------------------------
exprDT <- fread(args$expression, header = T, stringsAsFactors = F)
# select only needed columns
exprDTcols <- match(unique(inventory_expr$expr_tumor_subtype),
                    colnames(exprDT))
exprDTcols <- c(1:2, unique(exprDTcols[!is.na(exprDTcols)]))
exprDT <- exprDT[, exprDTcols, with = F]
setnames(exprDT, colnames(exprDT)[1:2], c('gene_id', 'gene_name'))

if (identical(colnames(exprDT), c('gene_id', 'gene_name'))) {
  exprTissues <- fread(args$expression, header = T, stringsAsFactors = F,
                       nrows = 10)
  exprTissues <- colnames(exprTissues)
  stop('[', Sys.time(), '] No expression data was read from expression ',
       'table. Tissues from expression inventory (', 
       paste0(setdiff(unique(inventory_expr$expr_tumor_subtype),
                      exprTissues), collapse = ', '), 
       ') are not found in expression table: ', 
       paste0(exprTissues, collapse = ', '))
}

exprDT[, gene_id := gsub('[.].*', '', gene_id)]
if (args$use_ensembl) {
  exprDT[, gene_id := gene_name]
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

# Modify gene_name and gene_id of expr. data to match GRs of interest ---------
# table showing if gene id or gene name from exprDT was found in GR
foundInGR <- data.table(found_id = exprDT$gene_id %in% GR$gene_id, 
                        found_name = exprDT$gene_name %in% GR$gene_name)
# we will only modify entries in gtex if just either gene name or gene id was
# found
changeName <- which(rowSums(foundInGR) == 1 & foundInGR$found_id == T)
setkey(GR, 'gene_id')
updatedNames <- GR[exprDT$gene_id[changeName]]
updatedNames <- updatedNames[,.SD[1], by = gene_id]
setkey(updatedNames, 'gene_id')
exprDT$gene_name[changeName] <- updatedNames[exprDT$gene_id[changeName]]$gene_name

changeID <- which(rowSums(foundInGR) == 1 & foundInGR$found_name == T)
setkey(GR, 'gene_name')
updatedIDs <- GR[exprDT$gene_name[changeID]]
updatedIDs <- updatedIDs[,.SD[1], by = gene_name]
setkey(updatedIDs, 'gene_id')
exprDT$gene_id[changeID] <- updatedIDs[exprDT$gene_id[changeID]]$gene_id

# Clean up table from duplicates ----------------------------------------------
# previous lines may lead to appearance of duplicates in the expression table
# check it and fix, if needed.
exprDT[, isDupl := .N > 1, by = .(gene_id, gene_name)]
if (any(exprDT$isDupl)) {
  duplRows <- exprDT[isDupl == T]
  message('[', Sys.time(), '] Found ',
          nrow(unique(duplRows[,.(gene_id, gene_name)])), ' duplicated ',
          'entries gene_id-gene_name pairs in the expression table. They ',
          'could be present in the table from the beginning or be created ',
          'by removing "after the dot" part from gene_id, use of ',
          '--use_ensembl flag or merge with gene_id-gene_name pairs from ',
          'scanned genomic regions. Only the entry with the highest mean ',
          'expression will be kept.')
  duplRows[, isDupl := NULL]
  duplRows[, meanExpr := rowMeans(duplRows[, setdiff(colnames(duplRows), 
                                                     c('gene_id', 'gene_name')),
                                           with = F])]
  duplRows <- duplRows[order(gene_id, gene_name, -meanExpr)] 
  duplRows <- duplRows[,.SD[1], by = .(gene_id, gene_name)]
  duplRows[, meanExpr := NULL]
  exprDT <- rbind(exprDT[isDupl == F][, setdiff(colnames(exprDT), 'isDupl'),
                                      with = F], duplRows)
}
# remove rows where expression in all tumor subtypes is NA
allIsNA <- apply(exprDT[, setdiff(colnames(exprDT), c('gene_id', 'gene_name')),
                        with = F], 1, function(x) sum(is.na(x)))
allIsNA <- allIsNA == (ncol(exprDT) - 2)
if (any(allIsNA)) {
  message('[', Sys.time(), '] Found ', sum(allIsNA), ' rows in the ',
          'expression table where all expression values are set to NA. Will ',
          'remove them.')
  exprDT <- exprDT[!allIsNA]
}

# Write to file ---------------------------------------------------------------
exprDT <- exprDT[, c('gene_id', 'gene_name',
                     inventory_expr$expr_tumor_subtype), with = F]
colnames(exprDT) <- c('gene_id', 'gene_name', inventory_expr$tumor_subtype)
write.table(exprDT, args$output, append = F, quote = F, sep = '\t',
            row.names = F, col.names = T)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')