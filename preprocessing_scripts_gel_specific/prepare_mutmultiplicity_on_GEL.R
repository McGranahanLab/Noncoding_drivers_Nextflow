#!/usr/bin/env Rscript
# FILE: prepare_mutmultiplicity_on_GEL.R --------------------------------------
#
# DESCRIPTION: creates tables with mutation multiplicity values computed by 
# custom script one per participant id.
#
# USAGE: run in R studio
# OPTIONS:
# EXAMPLE: 
# REQUIREMENTS: R 4.0.2, data.table
#               
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  05.12.2023
# REVISION: 05.12.2023

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

# Functions -------------------------------------------------------------------
#' parse_mutation_id
#' @description Parses mutation_id vector to data table of columns 
#'              participant_id, chr, pos, REF
#' @author Maria Litovchenko
#' @param mutIDs vector of mutation ids - either from column mutation_id or 
#'               key, key_orig
#' @param delim character - delimiter used in coding mutation_id
#' @param type one of mutmultiplicity or mutmap
#' @return data table with participant_id, chr, pos, REF
parse_mutation_id <- function(mutIDs, delim = ':', type = 'mutmultiplicity') {
  if (!type %in% c('mutmultiplicity', 'mutrate')) {
    stop('[', Sys.time(), '] Only mutmultiplicity or mutrate is accepted as ',
         'type')
  }
  
  result <- strsplit(unlist(mutIDs), split = delim)
  result <- lapply(result, unlist)
  result <- as.data.table(do.call(rbind, result))
  if (type == 'mutmultiplicity') {
    colnames(result) <- c('participant_id', 'chr', 'pos', 'REF')
    result[, participant_id := gsub('^X', '', participant_id)]
  } else {
    colnames(result) <- c('chr', 'pos', 'REF', 'ALT')
  }
  
  result
}

# Inputs ----------------------------------------------------------------------
mutmultiplicity_dir <- '/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/Timing/mutation_timing/'
# this file is produced during the pipeline run.
mutations_pancan_path <- '/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/completed_runs/06_12_2023/inputs/inputMutations-Panlung-hg19.maf'
outputDir <- '/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/data/mutmultiplicity/'
dir.create(outputDir, recursive = T)

# Read mutation multiplicity data & output them to file -----------------------
timing_RData <- list.files(mutmultiplicity_dir, pattern = '.RData$', 
                           full.names = T)
timing_RData <- unique(timing_RData)

timingDT <- data.table()
for (i in 1:length(timing_RData)) {
  message('[', Sys.time(), '] Loading timing data ...', i, '/',
          length(timing_RData), ' (', round(100 * i / length(timing_RData), 2),
          '%)')
  load(timing_RData[[i]])
  region.earlyLate <- as.data.table(region.earlyLate)
  region.earlyLate <- region.earlyLate[,.(mutation_id, Reference_Base,
                                          Alternate_Base, mut.multi)]
  region.earlyLate <- cbind(parse_mutation_id(region.earlyLate$mutation_id, 
                                              type = 'mutmultiplicity'), 
                            region.earlyLate)
  region.earlyLate <- region.earlyLate[,.(participant_id, chr, pos,
                                          Reference_Base, Alternate_Base,
                                          mut.multi)]
  setnames(region.earlyLate, c('pos', 'Reference_Base', 'Alternate_Base'),
           c('start', 'ref', 'var'))
  region.earlyLate <- as.data.table(apply(region.earlyLate, 2, unlist))
  
  # output to file
  write.table(region.earlyLate, append = F, quote = F, sep = '\t',
              file = paste0(outputDir, '/', 
                            unique(region.earlyLate$participant_id), 
                            '.mutmultiplicity.csv'),
              row.names = F, col.names = T)
  
  rm(region.earlyLate)
}