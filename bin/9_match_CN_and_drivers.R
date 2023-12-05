#!/usr/bin/env Rscript
# FILE: match_CN_and_drivers.R -----------------------------------------------
#
# DESCRIPTION: 
#
# USAGE: 
# OPTIONS:
# EXAMPLE: 
# REQUIREMENTS: 
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  17.12.2020
# REVISION: 05.12.2023

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

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges)))
suppressWarnings(suppressPackageStartupMessages(library(plyranges)))

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'match_CN_and_drivers.R')

# Test arguments --------------------------------------------------------------
args <- list(cancer_subtype = 'LUAD',
             inventory_patients = '../data/inventory/inventory_patients_tcga.csv',
             genomic_regions = '../TEST/inputs/inputGR-LUAD-hg19.bed',
             target_genome_version = '',
             ascat = '', ascat_genome_version = '', chain = '', ploidy = '',
             drivers = '../TEST/results/tables/drivers/drivers-LUAD--hg19.csv',
             chr = '', amp = log2(4/2), gain = log2(2.5/2), loss = log2(1.5/2))
             

# drivers are optional
selected_chrGR <- NULL
if (!is.null(args$chr)) {
  selected_chrGR <- makeGRangesFromDataFrame(data.frame(chr = args$chr, 
                                                        start = 1, end = 2))
}

# Read patient inventory table ------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients)
patientsInv <- patientsInv[tumor_subtype %in% args$cancer_subtype]
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)

# Read segment wise copy number estimates from ASCAT --------------------------
ascat <- fread(args$ascat, header = T, stringsAsFactors = F, 
               select = c('patient', 'chr', 'startpos', 'endpos', 
                          'nMajor', 'nMinor'))
setnames(ascat, c('patient', 'startpos', 'endpos'),
         c('participant_id', 'start', 'end'))

# select only samples which are in mutation table
ascat <- ascat[participant_id %in% patientsInv$participant_id]
ascat[, chr := gsub('chr', '', chr)]

# add ID for each ASCAT region
ascat <- ascat[order(participant_id, chr, start, end)]
ascat[, cn_reg_id := .(1:.N), by = .(participant_id, chr)]
ascat[, cn_reg_id := apply(ascat[,.(participant_id, chr, cn_reg_id)], 1,
                           paste0, collapse = '-')]
ascat[, cn_reg_id := gsub(' ', '', cn_reg_id)]

ascat <- makeGRangesFromDataFrame(ascat, keep.extra.columns = T)
# perform liftover, if needed
if (args$ascat_genome_version != args$target_genome_version) {
  chain_to_use <- import.chain(args$chain)
  ascat <- liftOverGenomicRegion(ascat, chain_to_use, args$min.gapwidt, 
                                 args$min.width)
  rm(chain_to_use)
  gc()
}

# select chromosome
if (!is.null(selected_chrGR)) {
  seqlevelsStyle(selected_chrGR) <- seqlevelsStyle(ascat)[1]
  ascat <- ascat[seqnames(ascat) %in% seqlevels(selected_chrGR)]
  message('[', Sys.time(), '] Selected genomic regions residing on ',
          'chromosome ', args$chr, ' from ASCAT table')
}

gc()
message('[', Sys.time(), '] Read ', args$ascat)

# Read driver genes -----------------------------------------------------------
if (!is.null(args$drivers)) {
  drivers <- fread(args$drivers, header = T, stringsAsFactors = F, 
                   select = c('gr_id', 'gene_id', 'gene_name', 'FILTER', 'tier',
                              'is_known_cancer', 'known_in_tumor_subtype', 
                              'known_cancer_biotype'))
  drivers <- unique(drivers)
  drivers <- drivers[FILTER == 'PASS' & !is.na(tier)]
}

# Read all scanned for cancer driver genes genomic regions --------------------
message('[', Sys.time(), '] Started reading input genomic ranges file')
# GR = regions of interest
GR <- readBED12(args$genomic_regions)

# select only GR for driver genes, if they are given
if (!is.null(args$drivers)) {
  GR <- as.data.table(GR)
  GR <- merge(GR, drivers[,.(gr_id, gene_id, gene_name)],
              by = c('gr_id', 'gene_id', 'gene_name'))
  GR <- makeGRangesFromDataFrame(GR, keep.extra.columns = T)
  message('[', Sys.time(), '] Genomic regions were restricted to regions ',
          'which belong to identified driver genomic regions.')
}

# select chromosome, if requested
if (!is.null(selected_chrGR)) {
  seqlevelsStyle(selected_chrGR) <- seqlevelsStyle(GRs)[1]
  GRs <- GRs[seqnames(GRs) %in% seqlevels(selected_chrGR)]
  message('[', Sys.time(), '] Selected genomic regions residing on ',
          'chromosome ', args$chr)
}
message('[', Sys.time(), '] Finished reading input genomic ranges file')

# Intersect ASCAT copy number and scanned genomic regions ---------------------
# produce GRanges object which will contain coordinates of regions inside 
# scanned gene regions for which copy number estimation from ASCAT is available

message('[', Sys.time(), '] Started intersecting ASCAT copy number and ',
        'scanned genomic regions')

seqlevelsStyle(ascat) <- seqlevelsStyle(GRs)[1]
GRs_CN <- list()
all_participant <- unique(ascat$participant_id)
message('[', Sys.time(), '] Number of participants: ', length(all_participant))
# do it in a loop to preserve memory
for (participant_idx in 1:length(all_participant)) {
  participant <- all_participant[participant_idx]
  message('[', Sys.time(), '] ... ', 
          round(100 * participant_idx / length(all_participant), 2), '%')
  toAdd <- join_overlap_intersect_directed(GRs, 
                                           ascat[ascat$participant_id == participant])
  GRs_CN <- c(GRs_CN, toAdd)
  
  # freeing memory a bit
  ascat <- ascat[ascat$participant_id != participant] 
  rm(toAdd)
  gc()
} 
rm(GRs, ascat, all_participant) # clean memory
gc()

GRs_CN <- GRangesList(GRs_CN)
GRs_CN <- unlist(GRs_CN)
message('[', Sys.time(), '] Finished intersecting ASCAT copy number and ',
        'scanned genomic regions')

# Mark genomic regions which have CNA -----------------------------------------
# assign copy number status to each region
GRs_CN <- as.data.table(GRs_CN)
GRs_CN[, width := NULL]
# add ploidy
GRs_CN <- merge(GRs_CN, ploidyDT, by = 'participant_id', all.x = T)

# use log-defined cutoffs for amplifications, loss and gains
GRs_CN[, cn_type := character()]
GRs_CN[, total_cn := nMinor + nMajor]
GRs_CN[, log_to_comp := log2(total_cn/Ploidy)]
GRs_CN[log_to_comp < args$loss]$cn_type <- 'loss'
GRs_CN[log_to_comp > args$gain]$cn_type <- 'gain'
GRs_CN[log_to_comp > args$amp]$cn_type <- 'amp'
GRs_CN[, log_to_comp := NULL]

# Compute number of copies relative to ploidy ---------------------------------
GRs_CN[, n_copies := total_cn / Ploidy]
GRs_CN[, total_cn := NULL]
GRs_CN <- GRs_CN[,.(participant_id, seqnames, start, end, strand, name, 
                    cn_reg_id, nMajor, nMinor, cn_type, n_copies)]
gc()

# Read mutation maps to genomic regions ---------------------------------------
message('[', Sys.time(), '] Started reading ', args$muts_to_gr)
varsToGRmap <- fread(args$muts_to_gr, header = T, stringsAsFactors = F, 
                     select = c('participant_id', 'gr_name', 'key',
                                'var_class')) 
message('[', Sys.time(), '] Finished reading ', args$muts_to_gr)

