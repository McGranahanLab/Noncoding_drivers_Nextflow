#!/usr/bin/env Rscript
# FILE: calculate_n_driver_coding_and_noncoding_muts_per_patient.R ------------
#
# DESCRIPTION: An R script
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla calculate_n_driver_coding_and_noncoding_muts_per_patient.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: R v4.1.0, data.table
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  01.02.2023
# REVISION: 02.06.2024

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
Sys.sleep(sample(1:15, 1))
source(paste0(srcDir, '/custom_functions_postprocessing.R'))

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'calculate_n_driver_coding_and_noncoding_muts_per_patient.R')

subtypeHelp <- paste('A name of the tumor subtype in which drivers were',
                     'detected.')
parser$add_argument("-c", "--cancer_subtype", required = T, 
                    type = 'character', help = subtypeHelp)

analysisHelp <- paste('Path to inventory table containing details of the',
                      'future analysis to be conducted. Minimal columns:',
                      'tumor_subtype,', 'software,', 'gr_id,', 'gr_code,', 
                      'gr_file,', 'gr_upstr,', 'gr_downstr,', 'gr_genome,', 
                      'gr_excl_id,', 'gr_excl_code,', 'gr_excl_file,',
                      'gr_excl_upstr,', 'gr_excl_downstr,', 'gr_excl_genome,',
                      'blacklisted_codes.')
parser$add_argument("-a", "--inventory_analysis", required = T, 
                    type = 'character', help = analysisHelp)

patientsHelp <- paste('Path to table listing information about patients,',
                      'their cancer types and mutation files. Minimal',
                      'columns: participant_id, tumor_subtype,',
                      'participant_tumor_subtype, somatic_path,',
                      'somatic_genome, cohort_name.')
parser$add_argument("-p", "--inventory_patients", required = T, 
                    type = 'character', help = patientsHelp)

excludedHelp <- paste('Path to file(s) containing IDs of excluded patients.')
parser$add_argument("-ep", "--excluded_patients", required = F, 
                    type = 'character', nargs = '*', default = NULL, 
                    help = excludedHelp)

mutationsHelp <- paste('Path to detected driver mutations file in the ',
                       'cancer subtype. Columns is_driver, gr_id, gene_id,',
                       'gene_name, var_type, key, participant_id, n_total,',
                       'n_driverMut_low, n_driverMut_mle, n_driverMut_high,',
                       'n_tums_w_driver_mut are needed.')
parser$add_argument("-d", "--driver_mutations", required = T, 
                    type = 'character', default = NULL, help = mutationsHelp)

foldHelp <- paste0('Indicates, whether or not splice sites should be',
                   'considered as part of CDS for the counting.',
                   'Default: T.')
parser$add_argument("-f", "--fold_splicesites_in_coding", required = F,
                    default = 'T', choices = c('T', 'F'), type = 'character', 
                    help = foldHelp)

excludeCNVhelp <- paste0('Indicates, whether or not driver copy number',
                         'variantions should be excluded from counting.',
                         'Default: T.')
parser$add_argument("-e", "--exclude_cnv", required = F,
                    default = 'T', choices = c('T', 'F'), type = 'character', 
                    help = excludeCNVhelp)

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output file")

args <- parser$parse_args()
args$fold_splicesites_in_coding <- as.logical(args$fold_splicesites_in_coding)
args$exclude_cnv <- as.logical(args$exclude_cnv)
if (!is.null(args$excluded_patients)) {
  if (!all(sapply(args$excluded_patients, file.exists))) {
    stop('[', Sys.time(), '] some of the files submitted to ',
         '--excluded_patients do not exist.')
  }
}

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

CNV_VAR_TYPES <- c('amp', 'hom.del.')
LOW_CONFIDENCE <- c('only n. driver mutations')

# Test input arguments --------------------------------------------------------
# args <- list(cancer_subtype = 'Panlung',
#              inventory_analysis = 'data/inventory/inventory_analysis.csv',
#              inventory_patients = 'data/inventory/inventory_patients.csv',
#              excluded_patients = list("completed_runs/2023_12_25/inputs/hypermutated-Panlung.csv"), 
#              driver_mutations = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Panlung--hg19.csv",
#              exclude_cnv = T, fold_splicesites_in_coding = T, 
#              output = 'completed_runs/2023_12_25/results/tables/driver_mutations/driverMutationsCountsPerPatient-Panlung--hg19.csv')

# Read in patients inventory --------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients, 1)
patientsInv <- patientsInv[tumor_subtype %in% args$cancer_subtype]
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)

# Read in file with hypermutated patients IDs ---------------------------------
if (!is.null(args$excluded_patients)) {
  hypermutated <- lapply(args$excluded_patients, fread, header = T, 
                         stringsAsFactors = F,
                         select = c('tumor_subtype', 'participant_id'))
  message('[', Sys.time(), '] Read --excluded_patients: ', 
          paste(args$excluded_patients, collapse = ', '))
  hypermutated <- do.call(rbind, hypermutated)
  patientsInv <- patientsInv[!participant_id %in% hypermutated$participant_id]
  message('[', Sys.time(), '] Excluded hypermutated participants ids')
}

# Read in analysis inventory --------------------------------------------------
analysisInv <- readAnalysisInventory(args$inventory_analysis)
message('[', Sys.time(), '] Read --inventory_analysis: ', 
        args$inventory_analysis)
analysisInv <- analysisInv[tumor_subtype == args$cancer_subtype]

# Read driver mutations -------------------------------------------------------
colsToKeep <- c("is_driver", "gr_id", "gene_id", "gene_name", "key", 
                "var_type", "participant_id", "confidenceLvl",
                "prob_is_driver_mle")
driverMuts <- fread(args$driver_mutations, header = T, select = colsToKeep,
                    stringsAsFactors = F)
driverMuts <- unique(driverMuts)
driverMuts <- driverMuts[is_driver == T]
driverMuts[, is_driver := NULL]
message('[', Sys.time(), '] Read --driver_mutations: ', args$driver_mutations)

# Fold splice sites into coding regions, if requested -------------------------
# assign coding genomic regions - anything which contains CDS as gr_code
coding_gr_id <- unique(analysisInv[gr_code == 'CDS']$gr_id)
message('[', Sys.time(), '] Following genomic regions: ', 
        paste0(coding_gr_id, collapse = ', '), ', will be considered as ',
        'coding.')

if (args$fold_splicesites_in_coding) {
  message('[', Sys.time(), '] Will fold splice sites into corresponding ',
          'coding driver genetic elements.')
  # assign splice site genomic regions - anything which contains ss as gr_code
  ss_gr_id <- unique(analysisInv[gr_code == 'ss']$gr_id)
  if (length(coding_gr_id) != 0 & length(ss_gr_id) != 0) {
    message('[', Sys.time(), '] Following genomic regions: ', 
            paste0(ss_gr_id, collapse = ', '), ', will be considered as ',
            'containing splice sites.')
    driverMuts <- foldSplicSiteDriversIntoCodingDrivers(coding_gr_id, ss_gr_id,
                                                        driverMuts)
  } else {
    message('[', Sys.time(), '] Did not find either coding or splice site ',
            'regions. Can not fold splice site drivers into coding ones.')
  }
}

noncoding_gr_id <- setdiff(unique(analysisInv$gr_id), coding_gr_id)
msg <- paste('[', Sys.time(), '] Following genomic regions:', 
             paste0(noncoding_gr_id, collapse = ', '), ', will be considered',
             'as noncoding.')
if (args$fold_splicesites_in_coding) {
  msg <- paste(msg, "Driver genomic elements considered as containing splice",
               "sites will be treated as noncoding for genes which do not",
               "have coding part detected as driver.")
}
message(msg)

# Set driver genomic element type (coding vs noncoding) -----------------------
driverMuts[, gr_type := ifelse(gr_id %in% coding_gr_id, 'coding', 
                               'non-coding')]

# Count number of SNV/small indel driver mutations per patient ----------------
# count how many RESOLVED driver SNVs/small indels per patient
nDriverSNVresolved <- driverMuts[!var_type %in% CNV_VAR_TYPES]
nDriverSNVresolved <- nDriverSNVresolved[!confidenceLvl %in% LOW_CONFIDENCE & 
                                           confidenceLvl != 'passenger']
nDriverSNVresolved <- nDriverSNVresolved[,.(length(unique(key))), 
                                         by = .(participant_id, gr_type)]
setnames(nDriverSNVresolved, 'V1', 'n_resolved')
nDriverSNVresolved[, struct_type := 'SNVs/small indels']

# next, count how many unresolved driver mutations per patient we have
nDriverSNVunresolved <- driverMuts[!var_type %in% CNV_VAR_TYPES]
nDriverSNVunresolved <- nDriverSNVunresolved[confidenceLvl == LOW_CONFIDENCE]
nDriverSNVunresolved <- nDriverSNVunresolved[,.(participant_id, gr_id, 
                                                gene_name, gene_id, key,
                                                prob_is_driver_mle, gr_type)]
nDriverSNVunresolved <- nDriverSNVunresolved[!is.na(prob_is_driver_mle)]
nDriverSNVunresolved <- nDriverSNVunresolved[,.(sum(prob_is_driver_mle)),
                                             by = .(participant_id, gr_type)]
nDriverSNVunresolved[, V1 := round(V1)]
nDriverSNVunresolved <- nDriverSNVunresolved[V1 > 0]
setnames(nDriverSNVunresolved, 'V1', 'n_unresolved')
nDriverSNVunresolved[, struct_type := 'SNVs/small indels']

nDriverSNV <- merge(nDriverSNVresolved, nDriverSNVunresolved, 
                    by = c('participant_id', 'gr_type', 'struct_type'),
                    all = T)
nDriverSNV[is.na(n_resolved)]$n_resolved <- 0
nDriverSNV[is.na(n_unresolved)]$n_unresolved <- 0

# in case we do have CNA - handle them
# FIX CNA COUNTING! NEED cn_reg_id, put it in key

# Combine SNV/small indel and CNA driver mutations ----------------------------
nDriverMuts <- copy(nDriverSNV)
nDriverMuts[, N := n_resolved + n_unresolved]

templateDT <- expand.grid(c('coding', 'non-coding'),
                          unique(nDriverMuts$struct_type))
templateDT <- as.data.table(templateDT)
colnames(templateDT) <- c('gr_type', 'struct_type')

nDriverMuts <- split(nDriverMuts, by = 'participant_id')
nDriverMuts <- lapply(nDriverMuts,
                      function(x) merge(x, templateDT, all = T,
                                        by = c('gr_type', 'struct_type')))
nDriverMuts <- lapply(nDriverMuts,
                      function(x) x[, participant_id := participant_id[!is.na(participant_id)]])
nDriverMuts <- do.call(rbind, nDriverMuts)
nDriverMuts[is.na(nDriverMuts)] <- 0

# Add patients where no driver mutations were found --------------------------
patientsWithoutDrivers <- setdiff(unique(patientsInv$participant_id),
                                  unique(nDriverMuts$participant_id))

if (length(patientsWithoutDrivers) > 0) {
  patientsWithoutDrivers <- expand.grid(participant_id = patientsWithoutDrivers,
                                        gr_type = c('coding', 'non-coding'),
                                        struct_type = unique(nDriverMuts$struct_type))
  patientsWithoutDrivers <- as.data.table(patientsWithoutDrivers)
  patientsWithoutDrivers <- cbind(patientsWithoutDrivers,
                                  n_resolved = 0, n_unresolved = 0, N = 0)
  nDriverMuts <- rbind(nDriverMuts, patientsWithoutDrivers)
}

nDriverMuts <- nDriverMuts[order(participant_id, gr_type, struct_type)]

# Output ----------------------------------------------------------------------
write.table(nDriverMuts, args$output, append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote output to ', args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')