#!/usr/bin/env Rscript
# FILE: match_mutmulti_and_mutations.R ----------------------------------------
#
# DESCRIPTION: Annotates mutations in genomic regions scanned for de novo 
# cancer driver discovery with their mutation multiplicity.
#
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla match_CN_and_scanned_genomic_regions.R -h
#          to see the full list of options and their descriptions.
# EXAMPLE: 
# REQUIREMENTS: argparse, data.table, plyr
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  17.12.2020
# REVISION: 01.12.2023

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
suppressWarnings(suppressPackageStartupMessages(library(plyr)))

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'match_mutations_and_multiplicity.R')

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

mutstoGrHelp <- paste('Path to files containing information about mutations',
                      'mapping to genomic regions. Usually produced by',
                      'calculate_mutation_rates.R')
parser$add_argument("-m", "--muts_to_gr", required = T, 
                    type = 'character', help = mutstoGrHelp)

driversHelp <- paste('Path to file containing de-novo detected cancer drivers')
parser$add_argument("-d", "--drivers", required = F, 
                    type = 'character', help = driversHelp)

synClassHelp <- paste('Variant_Classification-s from MAF format which are',
                      'acceptable as markers of synonymous variants. If given',
                      'variants of those classes will be removed from',
                      'consideration. Suggested values: Silent')
parser$add_argument("-sc", "--synAcceptedClass", required = F, nargs = '+',
                    default = NULL, type = 'character', help = synClassHelp)

outputHelp <- paste('Path to the output file')
parser$add_argument("-o", "--output", required = T, 
                    type = 'character', help = outputHelp)

args <- parser$parse_args()
# check_input_arguments_postproc(args, outputType = 'file')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test inputs -----------------------------------------------------------------
# args <- list(cancer_subtype = 'Adenocarcinoma',
#              inventory_patients = 'data/inventory/inventory_patients.csv',
#              muts_to_gr = 'completed_runs/2023-12-14/results/mut_rates/mutMapToGR-Adenocarcinoma--hg19.csv',
#              drivers = 'completed_runs/2023-12-14/results/tables/drivers/drivers-Panlung--hg19.csv',
#              output = '')

# Read patient inventory table ------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients)
patientsInv <- patientsInv[tumor_subtype %in% args$cancer_subtype]
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)
if (!'mutmultiplicity_path' %in% colnames(patientsInv)) {
  stop('[', Sys.time(), '] ', args$inventory_patients, ' does not have ',
       'column mutmultiplicity_path')
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
    stop('[', Sys.time(), '] no significant (FILTER is PASS and tier is not ',
         'NA) driver genes is found in ', args$drivers, ' table.')
  }
}

# Read mutation map to genomic regions ----------------------------------------
varsToGRmap <- fread(args$muts_to_gr, header = T, stringsAsFactors = F, 
                     select = c('participant_id', 'gr_name', 'key', 'key_orig',
                                'var_class'))
message('[', Sys.time(), '] Read ', args$muts_to_gr)

# restrict to driver genes, if requested
if (!is.null(args$drivers)) {
  grNames <- apply(drivers[,.(gr_id, gene_id, gene_name)], 1, paste0, 
                   collapse = '--')
  varsToGRmap[, genomeVersion := gsub('--.*', '', gr_name)]
  varsToGRmap[, genomeVersion := paste0(genomeVersion, '--')]
  varsToGRmap[, gr_name_no_version := apply(varsToGRmap, 1,
                                            function(x) gsub(x['genomeVersion'],
                                                             '',
                                                             x['gr_name']))]
  varsToGRmap <- varsToGRmap[gr_name_no_version %in% grNames]
  varsToGRmap[, genomeVersion := NULL]
  varsToGRmap[, gr_name_no_version := NULL]
  message('[', Sys.time(), '] Mutations were restricted to mutations ',
          'which belong to identified driver genomic regions.')
}

# remove silent mutations, if requested
if (!is.null(args$synAcceptedClass)) {
  nbefore <- nrow(varsToGRmap)
  varsToGRmap <- varsToGRmap[!var_class %in% args$synAcceptedClass]
  nafter <- nrow(varsToGRmap)
  message('[', Sys.time(), '] Removed ', nbefore - nafter, ' mutations out ',
          'of ', nbefore, '(', 100*round((nbefore - nafter)/nbefore, 4), '%) ',
          'because they fell into ', 
          paste(args$synAcceptedClass, collapse = ', '), ' class(es).')
}

patientsInv <- patientsInv[participant_id %in% 
                             unique(varsToGRmap$participant_id)]

# Read mutation multiplicity files --------------------------------------------
message('[', Sys.time(), '] Started reading input mutation multiplicity ',
        'files')
n_files <- nrow(patientsInv)
mutMulti <- list()
for (i in 1:n_files) {
  mutMulti[[i]] <- readSomaticVars(patientsInv$mutmultiplicity_path[i])
  # because files in MAF format do contain Tumor_Sample_Barcode renamed to 
  # participant_id
  suppressWarnings(mutMulti[[i]][, participant_id := NULL])
  mutMulti[[i]] <- cbind(mutMulti[[i]],
                         participant_id = patientsInv$participant_id[i])
  mutMulti[[i]][, fileType := NULL]
  
  if ((100 * (i / n_files)) %% 2) {
    message('\t[', Sys.time(), '] Read: ', round(100 * (i / n_files), 2), '%')
  }
}
mutMulti <- as.data.table(do.call(rbind.fill, mutMulti))
message('[', Sys.time(), '] Finished reading input mutation multiplicity ',
        'files')

# Join mutation map table and mutation multiplicity table ---------------------
# in case liftover was performed, the key_orig should be used
if ('key_orig' %in% colnames(varsToGRmap)) {
  varsToGRmap[, key_for_merge := key_orig]
} else {
  varsToGRmap[, key_for_merge := key]
}

mutMulti[, key_for_merge := key]
mutMulti[, key_for_merge := gsub('^chr', '', key_for_merge)]
if (grepl('^chr', varsToGRmap$key_for_merge[1])) {
  mutMulti[, key_for_merge := paste0('chr', key_for_merge)]
}

varsToGRmap[, participant_id := as.character(participant_id)]
mutMulti[, participant_id := as.character(participant_id)]

message('[', Sys.time(), '] Started joining mutation map table and mutation ',
        'multiplicity table')
varsToGRmapMutMulti <- merge(varsToGRmap, 
                             mutMulti[,.(participant_id, key_for_merge,
                                         mut.multi)], all.x = T,
                             by = c('participant_id', 'key_for_merge'))
message('[', Sys.time(), '] Finished joining mutation map table and mutation ',
        'multiplicity table')

# Output to file --------------------------------------------------------------
varsToGRmapMutMulti[, key_for_merge := NULL]
write.table(varsToGRmapMutMulti, file = args$output, append = F, quote = F,
            sep = '\t', row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote multiplicity resolved mutation to ',
        args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')