#!/usr/bin/env Rscript

# FILE: run_dndscv.R ----------------------------------------------------------
#
# DESCRIPTION: Runs dNdScv on supplied mutation file
# 
# USAGE: Rscript --vanilla run_dndscv.R 
#                --variants [path to mutation file] \
#                --genomic_regions [path to Rda file] \
#                --with_covariates [whether covariates should be used] \
#                --computeCI [whether confidence interval should be comp.] \
#                --genes_of_interest [path to file with genes of interest] \
#                --genes_of_interest_sep [separator used in file with genes of interest] \
#                --output [path to the output file] \
#                --outputGlobal [path to write global dNdScv rates] \
#                --outputRds [path to write Rds object] 
# 
# REQUIREMENTS: argparse, biomaRt, data.table, dndscv, GenomicRanges, tools,
#               TxDb.Hsapiens.UCSC.hg38.knownGene, rtracklayer
#
# EXAMPLES: 
#
# BUGS: --
# NOTES:
# 
# NOTES:  [ !!! I M P O R T A N T !!! ]
# DO NOT PRE-FILTER VARIANTS BASED ON LOCATION IN CDS! It causes wrong 
# estimation of wspl (splicing sites)
# 
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  15.01.2021
# REVISION: 04.10.2023

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
suppressWarnings(suppressPackageStartupMessages(library(biomaRt)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dndscv)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges)))
suppressWarnings(suppressPackageStartupMessages(library(tools)))
suppressWarnings(suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
options(scipen = 999)

# Functions -------------------------------------------------------------------
#' run_dNdScv_ci95
#' @description Runs dNdScv with computing of confidence intervals. Essentially
#' just a wrapper.
#' @author Maria Litovchenko
#' @param mutsTab data frame with mutations, columns sampleID, chr, pos, ref, 
#' mut
#' @param gois genes of interest, list of genes to restrict the analysis 
#' @param CV covariates (a matrix of covariates -columns- for each gene -rows-)
#' @param refRda object with annotation of a reference genome
#' @param computeCI boolean, indicates, whatever confidence intervals needs to
#' be computed
#' @return list of objects, same as function dndscv
run_dNdScv_ci95 <- function(mutsTab, gois = NULL, CV = NULL, refRda = NULL,
                            computeCI = F) {
  if (is.null(refRda)) {
    message('[', Sys.time(), '] No refRda was submitted, using build in hg19 ',
            'one')
    result <- dndscv(mutsTab, gene_list = gois, outmats = T, cv = CV)
  } else {
    result <- dndscv(mutsTab, gene_list = gois, outmats = T, cv = CV,
                     refdb = refRda)
  }
  result$sel_cv <- as.data.table(result$sel_cv)
  
  # compute confidence interval
  result$ci95 <- NULL
  if (computeCI) {
    result$ci95 <- as.data.table(geneci(result))
    setnames(result$ci95, 'gene', 'gene_name')
    result$sel_cv <- merge(result$sel_cv, result$ci95, by = 'gene_name',  
                           all = T)
  }
  
  result
}

# Parse input arguments -------------------------------------------------------
parser <- ArgumentParser(prog = 'run_dndscv.R')
parser$add_argument('-v', '--variants', nargs = 1, required = T, 
                    type = "character", help = 'Path to variants file')
parser$add_argument('-r', '--genomic_regions', nargs = 1, required = F, 
                    type = "character", default = NULL, 
                    help = 'Path to Rda object')
parser$add_argument('-c', '--with_covariates', nargs = 1, required = F, 
                    type = "logical", default = T, 
                    help = 'Whatever or not covariates should be used')
parser$add_argument('-ci', '--computeCI', nargs = 1, required = F, default = T,
                    type = "logical",
                    help = paste('Whatever confidence intervals for genes ',
                                 'should be computed. Increases run time! ',
                                 'Accepted values: T and F.'))
parser$add_argument('-g', '--genes_of_interest', nargs = 1, default = NULL, 
                    required = F, type = "character",
                    help = paste('Path to file containing  symbols of genes ',
                                 'of interest. They should be contained in',
                                 'column with name "symbol" or ',
                                 '"Gene Symbol".'))
parser$add_argument('-s', '--genes_of_interest_sep', nargs = 1, default = '\t', 
                    required = F, type = "character",
                    help = 'Separator used in genes_of_interest file.')
parser$add_argument('-o', '--output', nargs = 1, required = T,
                    type = "character", help = 'Path to save dNdScv results.')
parser$add_argument('-og', '--outputGlobal', nargs = 1, required = T,
                    type = "character",
                    help = 'Path to save global dNdScv estimations.')
parser$add_argument('-or', '--outputRds', nargs = 1, required = T,
                    type = "character",
                    help = 'Path to save Rds with all dNdScv fields')

args <- parser$parse_args()
args$with_covariates <- as.logical(args$with_covariates)
args$computeCI <- as.logical(args$computeCI)
check_input_arguments_preproc(args, outputType = 'file')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Select genes of interest, if required ---------------------------------------
# load dNdScv data base
if (is.null(args$genomic_regions)) {
  data("refcds_hg19", package = "dndscv")
} else {
  load(args$genomic_regions)
}

# all genes covered in dNdScv object
dndscvAvailGenes <- sapply(RefCDS, function(x) x$gene_name)
rm(RefCDS)
GOIs <- NULL

# select genes
if (!is.null(args$genes_of_interest)) {
  message('[', Sys.time(),'] Gene of interest submitted')
  GOIs <- suppressWarnings(fread(args$genes_of_interest, header = T,
                                 stringsAsFactors = F, fill = T, 
                                 sep = args$genes_of_interest_sep,
                                 select = c('Gene Symbol', 'symbol')))
  
  if (ncol(GOIs) == 0) {
    stop('[', Sys.time(), '] Columns "Gene Symbol" or "symbol" were not found')
  } else {
    GOIs <- unique(unlist(GOIs[, 1]))
    nBefore <- length(GOIs)
    GOIs <- GOIs[GOIs %in% dndscvAvailGenes]
    message('[', Sys.time(),'] Removed ', nBefore - length(GOIs), ' genes ',
            'from genes_of_interest (', 
            round(100 * (nBefore - length(GOIs)) / nBefore), '%) because ',
            'they are not in dNdScv database.')
    
    if (length(GOIs) == 0) {
      stop('[', Sys.time(), '] No genes left to analyse')
    }
  }
}

# Prepare the data ------------------------------------------------------------
# [IMPORTANT NOTE] DO NOT PRE-FILTER VARIANTS BASED ON LOCATION IN CDS! It 
# causes wrong estimation of wspl (splicing sites)
# Read in variants
mutsDF <- fread(args$variants, header = T, sep = '\t')

# After dNdScv will be run, also Rds objects will be saved, because:
# When using a new RefCDS, please exercise caution in interpreting the results.
# For example, low values of theta (e.g. dndsout$nbreg$theta << 1) indicate that
# there is large unexplained variation in the mutation density across genes and 
# may mean that dNdScv is not adequate for this dataset.

# Run dNdScv without covariates -----------------------------------------------
if (!args$with_covariates) {
  results_dNdScvNoCov <- run_dNdScv_ci95(mutsTab = mutsDF, gois = GOIs, 
                                         CV = NULL, refRda = args$genomic_regions,
                                         computeCI = args$computeCI)
  # save Rds
  saveRDS(results_dNdScvNoCov, args$outputRds)
  # write table with global rates
  write.table(results_dNdScvNoCov$globaldnds, args$outputGlobal, append = F, 
              col.names = T, sep = ',', quote = F, row.names = F)
  # table to write to txt file
  write.table(results_dNdScvNoCov$sel_cv, args$output, append = F, 
              col.names = T, sep = ',', quote = F, row.names = F)
}

# Run dNdScv with covariates --------------------------------------------------
if (args$with_covariates) {
  suppressWarnings(rm(RefCDS, gr_genes))
  data("covariates_hg19", package = "dndscv")
  message('[WARNING]: dNdScv covariantes from hg19 will be used!')
  
  # in order to use covariates with custom RefRda, we need to restrict their 
  # rows to genes which are present in RefRda
  if (is.null(args$genomic_regions)) {
    data("refcds_hg19", package = "dndscv")
  } else {
    load(args$genomic_regions)
  }
  
  genesInRefCDS <- sapply(RefCDS, function(x) x$gene_name)
  genesNotInCovs <- genesInRefCDS[genesInRefCDS != 'CDKN2A']
  genesNotInCovs <- genesNotInCovs[!genesNotInCovs %in% rownames(covs)]
  if (length(genesNotInCovs) != 0) {
    stop('[', Sys.time(), '] Genes ', genesNotInCovs, 
         paste(genesNotInCovs, collapse = ', ', ' are listed in the ',
               'submitted --genomic_regions but are not present in dNdScv ',
               'covariates table. Run of dNdScv with covariates is not ',
               'possible.'))
  }
  if (any(genesInRefCDS == 'CDKN2A')) {
    rownames(covs)[rownames(covs) == 'CDKN2A.p14arf'] <- 'CDKN2A'
  }
  covs <- covs[intersect(rownames(covs), genesInRefCDS), ]
  rm(RefCDS, gr_genes)
  
  results_dNdScvCov <- run_dNdScv_ci95(mutsTab = mutsDF, gois = GOIs,
                                       CV = covs, refRda = args$genomic_regions,
                                       computeCI = args$computeCI)
  # save Rds
  saveRDS(results_dNdScvCov, args$outputRds)
  # write table with global rates
  write.table(results_dNdScvCov$globaldnds, args$outputGlobal, append = F, 
              col.names = T, sep = ',', quote = F, row.names = F)
  # table to write to txt file
  write.table(results_dNdScvCov$sel_cv, args$output, append = F, 
              col.names = T, sep = ',', quote = F, row.names = F)
}

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')