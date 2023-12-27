#!/usr/bin/env Rscript
# FILE: match_CN_and_scanned_genomic_regions.R --------------------------------
#
# DESCRIPTION: Annotates genomic regions scanned for de novo cancer driver 
# discovery with their copy number state.
#
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla match_CN_and_scanned_genomic_regions.R -h
#          to see the full list of options and their descriptions.
# EXAMPLE: 
# REQUIREMENTS: argparse, data.table, GenomicRanges, plyranges
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

subtypeHelp <- paste('A cancer subtype to select from patientsInv table. Only',
                     'mutations from patients with that cancer type will be',
                     'selected. In case an analysis of several cancer types',
                     'needed to be performed please run this script ',
                     'separetedly for each cancer type.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    help = subtypeHelp)

patientsHelp <- paste('Path to patientsInv table listing information about',
                      'patients, their cancer types and mutation files.',
                      'Minimal columns: participant_id, tumor_subtype,',
                      'participant_tumor_subtype, somatic_path,',
                      'somatic_genome, cohort_name')
parser$add_argument("-p", "--inventory_patients", required = T, 
                    type = 'character', help = patientsHelp)

grHelp <- paste('Part to file containing genomic ranges of region of',
                'interest, bed12 format')
parser$add_argument("-g", "--genomic_regions", required = T, 
                    type = 'character', help = grHelp)

targetVersionHelp <- paste('Target genome version')
parser$add_argument("-gv", "--target_genome_version", required = T,
                    type = 'character', help = targetVersionHelp)

chainHelp <- paste('Path to chain file in case genome version of ASCAT copy',
                   'number regions is not the same as --target_genome_version')
parser$add_argument("-l", "--chain", required = F, default = NULL,
                    type = 'character', help = chainHelp)

driversHelp <- paste('Path to file containing de-novo detected cancer drivers')
parser$add_argument("-d", "--drivers", required = F,  default = NULL,
                    type = 'character', help = driversHelp)

chrHelp <- paste('Chromosome to restrict matching of copy number to genomic ',
                 'regions.')
parser$add_argument("-chr", "--chr", required = F, default = NULL,
                    type = 'character', help = chrHelp)

ampHelp <- paste('Cut off (in log 2 scale) on copy number to determine',
                 'amplification. Default: log2(4/2) = 1')
parser$add_argument("-amp", "--amp", required = F, default = 1,
                    type = 'double', help = ampHelp)

gainHelp <- paste('Cut off (in log 2 scale) on copy number to determine',
                  'gain. The gain cut off should be < amplification cut off.',
                  'Default: log2(2.5/2) = 0.3219281.')
parser$add_argument("-gain", "--gain", required = F, default = 0.3219281,
                    type = 'double', help = gainHelp)

lossHelp <- paste('Cut off (in log 2 scale) on copy number to determine',
                  'loss. Default: log2(1.5/2) = -0.4150375')
parser$add_argument("-loss", "--loss", required = F, default = -0.4150375,
                    type = 'double', help = lossHelp)

minGapWhelp <- paste('Minimum length of a gap between lifted over regions',
                     'that will prevent them from being merged together.',
                     'Default = Inf (no reduce will be performed)')
parser$add_argument("-mgw", "--min_gapwidth", required = F, default = 'Inf',
                    type = 'character', help = minGapWhelp)

minWhelp <- paste('Minimum width of lifted over genomic regions.',
                  'Default = 0 (no regions will be filtered out).')
parser$add_argument("-mw", "--min_width", required = F, default = '0',
                    type = 'character', help = minWhelp)

outputHelp <- paste('Path to the output file')
parser$add_argument("-o", "--output", required = T, 
                    type = 'character', help = outputHelp)

args <- parser$parse_args()
args$min_gapwidt <- as.numeric(args$min_gapwidt)
args$min_width <- as.numeric(args$min_width)
# check_input_arguments_postproc(args, outputType = 'file')
if (args$amp <= args$gain) {
  stop('[', Sys.time(), '] Cut off on amplification should be > cut off on ',
       'gain.')
}

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
selected_chrGR <- NULL
if (!is.null(args$chr)) {
  selected_chrGR <- makeGRangesFromDataFrame(data.frame(chr = args$chr, 
                                                        start = 1, end = 2))
}
printArgs(args)

# Test arguments --------------------------------------------------------------
# args <- list(cancer_subtype = 'Panlung',
#              inventory_patients = 'data/inventory/inventory_patients.csv',
#              genomic_regions = 'completed_runs/2023-12-14/inputs/inputGR-Panlung-hg19.bed',
#              target_genome_version = 'hg19', chain = NULL, 
#              drivers = 'completed_runs/2023-12-14/results/tables/drivers/drivers-Panlung--hg19.csv',
#              chr = NULL, amp = log2(4/2), gain = log2(2.5/2), 
#              loss = log2(1.5/2), min_gapwidt = '', min_width = '', 
#              output = '')
# selected_chrGR <- NULL
# if (!is.null(args$chr)) {
#   selected_chrGR <- makeGRangesFromDataFrame(data.frame(chr = args$chr, 
#                                                         start = 1, end = 2))
# }

# Read patient inventory table ------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients)
patientsInv <- patientsInv[tumor_subtype %in% args$cancer_subtype]
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)
if (!'cn_segments_path' %in% colnames(patientsInv)) {
  stop('[', Sys.time(), '] ', args$inventory_patients, ' does not have ',
       'column cn_segments_path.')
}

# Read driver genes -----------------------------------------------------------
if (!is.null(args$drivers)) {
  message('[', Sys.time(), '] Started reading ', args$drivers)
  drivers <- fread(args$drivers, header = T, stringsAsFactors = F, 
                   select = c('gr_id', 'gene_id', 'gene_name', 'FILTER', 
                              'tier', 'is_known_cancer', 
                              'known_in_tumor_subtype', 
                              'known_cancer_biotype'))
  drivers <- unique(drivers)
  drivers <- drivers[FILTER == 'PASS' & !is.na(tier)]
  message('[', Sys.time(), '] Finished reading ', args$drivers)
  
  if (nrow(drivers) == 0) {
    stop('[', Sys.time(), '] no significant (FILTER is PASS and tier is ',
            'not NA) driver genes is found in ', args$drivers, ' table.')
  }
}

# Read segment wise copy number estimates -------------------------------------
message('[', Sys.time(), '] Started reading copy number segments')
cnSegs <- lapply(patientsInv$cn_segments_path, fread, header = T, 
                 stringsAsFactors = F, 
                 select = c('participant_id', 'chr', 'start', 'end', 
                            'nMajor', 'nMinor', 'Ploidy'))
cnSegs <- do.call(rbind, cnSegs)
setnames(cnSegs, 'Ploidy', 'ploidy')
cnSegs[, chr := gsub('chr', '', chr)]

cnSegs <- makeGRangesFromDataFrame(cnSegs, keep.extra.columns = T)
# perform liftover, if needed
if (unique(patientsInv$cn_segments_genome) != args$target_genome_version) {
  if (is.null(args$chain)) {
    stop('[', Sys.time(), '] target genome version (',
         args$target_genome_version, ') is not the same as genome version ',
         'of copy number segments (', unique(patientsInv$cn_segments_genome),
         '), but no chain file is provided.')
  }
  chain_to_use <- import.chain(args$chain)
  cnSegs <- liftOverGenomicRegion(cnSegs, chain_to_use, args$min_gapwidt, 
                                  args$min_width)
  rm(chain_to_use)
  gc()
}

# select chromosome
if (!is.null(selected_chrGR)) {
  seqlevelsStyle(selected_chrGR) <- seqlevelsStyle(cnSegs)[1]
  cnSegs <- cnSegs[seqnames(cnSegs) %in% seqlevels(selected_chrGR)]
  message('[', Sys.time(), '] Selected genomic regions residing on ',
          'chromosome ', args$chr, ' from copy number regions')
}

gc()
message('[', Sys.time(), '] Finished reading copy number segments')

# Read all scanned for cancer driver genes genomic regions --------------------
message('[', Sys.time(), '] Started reading ', args$genomic_regions)
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
  seqlevelsStyle(selected_chrGR) <- seqlevelsStyle(GR)[1]
  GR <- GR[seqnames(GR) %in% seqlevels(selected_chrGR)]
  message('[', Sys.time(), '] Selected genomic regions residing on ',
          'chromosome ', args$chr)
}
message('[', Sys.time(), '] Finished reading ', args$genomic_regions)

# Intersect copy number regions and scanned genomic regions -------------------
# produce GRanges object which will contain coordinates of regions inside 
# scanned gene regions for which copy number estimations are available
message('[', Sys.time(), '] Started intersecting copy number regions and ',
        'scanned genomic regions')

seqlevelsStyle(cnSegs) <- seqlevelsStyle(GR)[1]
cnSegs <- cnSegs[seqnames(cnSegs) %in% unique(seqnames(GR))]
GRs_CN <- list()
all_participants <- unique(cnSegs$participant_id)
message('[', Sys.time(), '] Number of participants: ', length(all_participants))
# do it in a loop to preserve memory
for (participant_idx in 1:length(all_participants)) {
  participant <- all_participants[participant_idx]
  message('[', Sys.time(), '] ... ', 
          round(100 * participant_idx / length(all_participants), 2), '%')
  toAdd <- join_overlap_intersect_directed(GR, 
                                           cnSegs[cnSegs$participant_id == participant])
  GRs_CN <- c(GRs_CN, toAdd)
  
  # freeing memory a bit
  cnSegs <- cnSegs[cnSegs$participant_id != participant] 
  rm(toAdd)
  gc()
} 
rm(GR, cnSegs, all_participants) # clean memory
gc()

GRs_CN <- GRangesList(GRs_CN)
GRs_CN <- unlist(GRs_CN)
message('[', Sys.time(), '] Finished intersecting copy number regions and ',
        'scanned genomic regions')

# Mark genomic regions which have CNA -----------------------------------------
# assign copy number status to each region
GRs_CN <- as.data.table(GRs_CN)
GRs_CN[, width := NULL]

# use log-defined cutoffs for amplifications, loss and gains
GRs_CN[, cn_type := character()]
GRs_CN[, total_cn := nMinor + nMajor]
GRs_CN[, log_to_comp := log2(total_cn/ploidy)]
GRs_CN[log_to_comp < args$loss]$cn_type <- 'loss'
GRs_CN[log_to_comp > args$gain]$cn_type <- 'gain'
GRs_CN[log_to_comp > args$amp]$cn_type <- 'amp'
GRs_CN[, log_to_comp := NULL]

# Compute number of copies relative to ploidy ---------------------------------
GRs_CN[, n_copies := total_cn / ploidy]
GRs_CN[, total_cn := NULL]
GRs_CN <- GRs_CN[,.(participant_id, seqnames, start, end, strand, gr_name, 
                    ploidy, nMajor, nMinor, cn_type, n_copies)]
gc()
 
# it can happen that a genomic region has multiple CNA inside. In such case
# we will handle cases like this by taking maximum deviating from neutrality 
# copy number score
GRs_CN[, is_max_dev_neutr := abs(n_copies - 1)]
GRs_CN[, is_max_dev_neutr := is_max_dev_neutr == max(is_max_dev_neutr), 
       by = .(participant_id, gr_name)]
GRs_CN[, n_copies_ge_maxDev := n_copies * is_max_dev_neutr]
# now n_copies_ge_maxDev is 0 in case there n_copies is not maximum deviating 
# from neutrality, so if we take max from n_copies_ge_maxDev it will fill 0s
# with proper values
GRs_CN[, n_copies_ge_maxDev := max(n_copies_ge_maxDev), 
       by = .(participant_id, gr_name)]
GRs_CN[, is_max_dev_neutr := NULL]

# Output to table -------------------------------------------------------------
write.table(GRs_CN, args$output, append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote copy number resolved genomic ranges to ',
        args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')