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
mutations_pancan_path <- '/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/completed_runs/29_11_2023/inputs/inputMutations-Panlung-hg19.maf'
outputDir <- '/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/data/mutmultiplicity/'
dir.create(outputDir, recursive = T)

# Read mutation multiplicity data -------------------------------------------
timing_RData <- list.files(mutmultiplicity_dir, pattern = '.RData$', 
                           full.names = T)
timing_RData <- unique(timing_RData)

timingDT <- data.table()
for (i in 1:length(timing_RData)) {
  message('[', Sys.time(), '] Loading timing data ...', i, '/',
          length(timing_RData), ' (', round(100 * i / length(timing_RData), 2),
          '%)')
  timing_i <- timing_RData[[i]]
  load(timing_i)
  region.earlyLate <- as.data.table(region.earlyLate)
  region.earlyLate <- region.earlyLate[,.(mutation_id, mut.multi)]
  timingDT <- rbind(timingDT, region.earlyLate)
  rm(region.earlyLate, timing_i)
}

timingDT <- as.data.table(apply(timingDT, 2, unlist))
timingDT[, mut.multi := as.numeric(mut.multi)]
timingDT <- cbind(parse_mutation_id(timingDT$mutation_id,  
                                    type = 'mutmultiplicity'), timingDT)
timingDT[, mutation_id := NULL]

# Read mutations to genomic regions maps for pancan ---------------------------
mutations_pancan <- fread(mutations_pancan_path, header = T, 
                          stringsAsFactors = F, 
                          select = c('Tumor_Sample_Barcode', 'key', 
                                     'key_orig'))
setnames(mutations_pancan, 'Tumor_Sample_Barcode', 'participant_id')
mutations_pancan <- unique(mutations_pancan)
mutations_pancan <- cbind(parse_mutation_id(mutations_pancan$key_orig, 
                                            type = 'mutrate'), 
                          mutations_pancan)

# Match mutation multiplicity with key of mutations ---------------------------
timingDT[, chr := gsub('^chr', '', chr)]
if (grepl('^chr', mutations_pancan$chr[1])) {
  timingDT[, chr := paste0('chr', chr)]
}
mutations_pancan[, participant_id := as.character(participant_id)]
timingDT[, participant_id := as.character(participant_id)]

timingDT <- merge(timingDT, mutations_pancan,
                  by = c('participant_id', 'chr', 'pos', 'REF'))
timingDT <- timingDT[,.(participant_id, key, key_orig, mut.multi)]
timingDT <- split(timingDT, by = 'participant_id')

# Output to files -------------------------------------------------------------
lapply(names(timingDT), 
       function(x) write.table(timingDT[[x]], 
                               paste0(outputDir, '/', x, 
                                      '.mutmultiplicity.csv'), 
                               append = F, quote = F, sep = '\t', 
                               row.names = F, col.names = T))