#!/usr/bin/env Rscript
# FILE: preprocess_CancerGenomeInterpreter.R ----------------------------------
#
# DESCRIPTION: A script to process table with known cancer driver mutations
#              downloaded from Cancer Genome Interpreter. 
# https://www.cancergenomeinterpreter.org/data/mutations/catalog_of_validated_oncogenic_mutations_20180116.zip
#
# USAGE: Rscript --vanilla preprocess_CancerGenomeInterpreter.R 
#
# OPTIONS: Please run to see options and help:
#          Rscript --vanilla preprocess_CancerGenomeInterpreter.R -h
#
# REQUIREMENTS: argparse, data.table
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  08.11.2023
# REVISION: 25.12.2023

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

options(scipen = 999)

# Functions -------------------------------------------------------------------
#' liftOverGenomicRegion
#' @description Lifts over coordinates from original genome to target genome
#' using chain chainObj
#' @author Maria Litovchenko
#' @param inGR input genomic regions. mcols are not used
#' @param chainObj chain to move genome coordinates from original genome 
#' version to the target one
#' @param min.gapwidt minimum length of a gap between lifted over regions that
#'                    will prevent them from being merged together. 
#'                    Default = Inf (no reduce will be performed).
#' @param min.width minimum width of lifted over genomic regions.
#'                  Default = 0 (no regions will be filtered out). 
#' @return lifted over genomic coordinates
#' @export
liftOverGenomicRegion <- function(inGR, chainObj, min.gapwidt = Inf, 
                                  min.width = 0) {
  message('[', Sys.time(), '] Started liftover')
  if (class(inGR) == 'CompressedGRangesList') {
    inGRunlist <- unlist(inGR)
  } else {
    inGRunlist <- inGR
  }
  # change chromosomal format, if needed
  seqlevelsStyle(inGRunlist) <- seqlevelsStyle(chainObj)[1]
  
  # perform liftover
  result <- liftOver(inGRunlist, chainObj)
  # inform about not lifted regions
  notLifted <- unlist(lapply(width(result), 
                             function(x) identical(x, integer(0))))
  notLiftedPerc <- round(100 * sum(notLifted) / length(inGRunlist), 2)
  message('[', Sys.time(), '] ', sum(notLifted), ' (', notLiftedPerc,
          '%) genomic regions were not lifted over.')
  # inform if a region was split into 2
  splitted <- unlist(lapply(width(result), function(x) length(x) > 1))
  splittedPerc <- round(100 * sum(splitted) / length(splitted), 2)
  message('[', Sys.time(), '] ', sum(splitted), ' (', splittedPerc,
          '%) genomic regions were splitted into 2 and more.')
  
  # merge regions, if min.gapwidt is not infinite
  if (!is.infinite(min.gapwidt)) {
    message('[', Sys.time(), '] Reduction of lifted over genomic regions ',
            'located closer than ', min.gapwidt, 'bp will be performed.')
    result <- reduce(result, min.gapwidt = min.gapwidt)
    # reduce removed all mcols, add them back
    not_empty <- which(!notLifted)
    lo_length <- elementNROWS(result)[!notLifted]
    idx_inOrigin <- rep(not_empty, lo_length)
    result <- unlist(result)
    mcols(result) <- mcols(inGRunlist)[idx_inOrigin, ]
  } else {
    result <- unlist(result)
  }
  
  # remove small regions, min.width is given
  if (min.width > 0) {
    message('[', Sys.time(), '] min.width is given. Summary of region ',
            'lengths before removal of regions with width < ', min.width)
    width_sum <- round(summary(width(result)), 0)
    message(paste(paste0(names(width_sum), collapse = '\t'), '\n', 
                  paste0(width_sum, collapse = '\t')))
    result <- result[width(result) >= min.width]
    message('[', Sys.time(), '] Summary of region lengths after removal of ',
            'regions with width < ', min.width)
    width_sum <- round(summary(width(result)), 0)
    message(paste(paste0(names(width_sum), collapse = '\t'), '\n', 
                  paste0(width_sum, collapse = '\t')))
  }
  
  message('[', Sys.time(), '] Ratio of total region length after and ',
          'before liftover: ', round(sum(width(result)) / sum(width(inGR)), 3))
  message('[', Sys.time(), "] Ratio of regions' numbers after and before ",
          'liftover: ', round(length(result) / length(inGR), 3))
  
  # set style of chromosomal names to the same as input data
  seqlevelsStyle(result) <- seqlevelsStyle(inGR)
  
  message('[', Sys.time(), '] Finished liftover')
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'preprocess_CancerGenomeInterpreter.R')

inputHelp <- paste('Path to a file with known cancer mutations from Cancer',
                   'Genome Interpreter. Navigate to',
                   'https://www.cancergenomeinterpreter.org/data/mutations/',
                   'and download catalog_of_validated_oncogenic_mutations_20180116.zip.',
                   'Unzip it and find catalog_of_validated_oncogenic_mutations.tsv.',
                   'Check file cancer_acronyms.txt to select cancer subtypes',
                   'acronyms corresponding to cancer subtype of your choise.')
parser$add_argument("-i", "--input", required = T, type = 'character', 
                    help = inputHelp)

acronymsHelp <- paste('Acronyms')
parser$add_argument("-a", "--acronyms", required = T, type = 'character', 
                    help = acronymsHelp, nargs = '+')

refGenomeHelp <- paste("Genome version, i.e. hg19, which will be used in", 
                       "de-novo driver discovery analysis. Default: hg19.")
parser$add_argument("-g", "--target_genome_version", required = T,
                    default = 'hg19', type = 'character',
                    help = refGenomeHelp)

refGenomeFormatHelp <- paste("Reference genome version which will be used in",
                             'de-novo driver discovery analysis. One of UCSC',
                             'or NCBI.')
parser$add_argument("-gf", "--reference_genome_format", required = T,
                    type = 'character', choices = c('NCBI', 'UCSC'),
                    help = refGenomeFormatHelp)

mutGenomeHelp <- paste("Genome version, i.e. hg19, of Cancer Genome",
                       "Interpreter files.")
parser$add_argument("-k", "--known_mut_drivers_genome", required = T,
                    default = 'hg19', type = 'character',
                    help = mutGenomeHelp)

outputHelp <- 'Path to output file'
parser$add_argument("-o", "--output", required = T, type = 'character', 
                    help = inputHelp)

args <- parser$parse_args()
args$acronyms <- unlist(args$acronyms)
timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')

if (!dir.exists(dirname(args$ouput))) {
  message('[', Sys.time(), '] Directory does not exist: ', dirname(args$ouput),
          '. Will create it.')
  dir.create(dirname(args$ouput), recursive = T)
}

# Test arguments --------------------------------------------------------------
# args <- list(input = '../data/assets_raw/CancerGenomeInterpreter/catalog_of_validated_oncogenic_mutations.tsv',
#              acronyms = list('L', 'LUAD', 'LUSC', 'NSCLC', 'SCLC',
#                                     'CANCER', 'CANCER-PR'),
#              known_mut_drivers_genome = 'hg19', 
#              reference_genome_format = 'UCSC', target_genome_version = 'hg19',
#              output = '../data/assets/CancerGenomeInterpreter_lung_hg19.tsv')

# Read known driver mutations -------------------------------------------------
message('[', Sys.time(), '] Started reading ', args$input)
knownDriverMuts <- fread(args$input, header = T, stringsAsFactors = F, 
                         select = c('gene', 'gdna', 'context', 
                                    'cancer_acronym', 'source'))
message('[', Sys.time(), '] Finished reading ', args$input)

# select only SNVs as it's problematic to perform proper match to indels
nbefore <- nrow(knownDriverMuts)
knownDriverMuts <- knownDriverMuts[grepl('>', gdna)]
nafter <- nrow(knownDriverMuts)
message('[', Sys.time(), '] Selected ', nafter, ' out of ', nbefore, ' (',
        100 * round(nafter/nbefore, 4), '%) known driver mutations as they ',
        'were SNVs.')

# parse gdna
knownDriverMuts[, chr := gsub(':.*', '', gdna)]
knownDriverMuts[, chr := gsub('chr', '', chr)]
if (args$reference_genome_format == 'UCSC') {
  knownDriverMuts[, chr := paste0('chr', chr)]
}
knownDriverMuts[, start := gsub('.*g.', '', gdna)]
knownDriverMuts[, start := gsub('[a-zA-Z].*', '', start)]
knownDriverMuts[, start := as.integer(start)]
# this will remove occasional indels
knownDriverMuts <- knownDriverMuts[!is.na(start)]
knownDriverMuts[, end := start]
knownDriverMuts[, Reference_Allele := gsub('.*g.', '', gdna)]
knownDriverMuts[, Reference_Allele := gsub('\\d', '', Reference_Allele)]
knownDriverMuts[, Reference_Allele := gsub('>.*', '', Reference_Allele)]
knownDriverMuts[, Tumor_Seq_Allele2 := gsub('.*>', '', gdna)]

knownDriverMuts <- knownDriverMuts[,.(chr, start, end, gene, Reference_Allele,
                                      Tumor_Seq_Allele2, cancer_acronym, 
                                      context, source)]

# Select cancer subtypes of interest ------------------------------------------
knownDriverMuts[, cancer_acronym := lapply(cancer_acronym,
                                           function(x) unlist(strsplit(x, 
                                                                       '_')))]
knownDriverMuts <- knownDriverMuts[sapply(cancer_acronym, 
                                          function(x) any(unlist(x) %in% 
                                                            args$acronyms))]
knownDriverMuts[, cancer_acronym := sapply(cancer_acronym, paste0, 
                                           collapse = '_')]

# Liftover, if needed ---------------------------------------------------------
if (args$target_genome_version != args$known_mut_drivers_genome) {
  knownDriverMuts <- makeGRangesFromDataFrame(knownDriverMuts, 
                                              keep.extra.columns = T)
  message('[', Sys.time(), '] Started lifover of known driver ',
          'mutations from ', args$known_mut_drivers_genome, ' to ', 
          args$target_genome_version)
  chain <- import.chain(args$driverMuts_chain)
  seqlevelsStyle(chain) <- seqlevelsStyle(knownDriverMuts)[1]
  knownDriverMuts <- liftOverGenomicRegion(knownDriverMuts, chain)
  message('[', Sys.time(), '] Finished lifover of known driver ',
          'mutations from ', args$known_mut_drivers_genome, ' to ', 
          args$target_genome_version)
  knownDriverMuts <- as.data.table(knownDriverMuts)
  rm(chain)
  gc()
}

# create keys in the same format as used in mutation table
knownDriverMuts <- as.data.table(knownDriverMuts)
setnames(knownDriverMuts, 'seqnames', 'chr', skip_absent = T)
knownDriverMuts[, key := apply(knownDriverMuts[,.(chr, start, Reference_Allele,
                                                  Tumor_Seq_Allele2)], 
                               1, paste, collapse = ':')]
knownDriverMuts[, key := gsub(' ', '', key)]

# Write to file ---------------------------------------------------------------
setnames(knownDriverMuts, c('gene', 'Reference_Allele', 'Tumor_Seq_Allele2'),
         c('gene_name', 'ref', 'alt'))
setcolorder(knownDriverMuts, c('chr', 'start', 'end', 'ref', 'alt', 'key',
                               'gene_name', 'cancer_acronym'))
write.table(knownDriverMuts, args$output, append = F, quote = F, sep = '\t',
            row.names = F, col.names = T)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')