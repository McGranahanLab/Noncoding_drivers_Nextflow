#!/usr/bin/env Rscript
# FILE: 1_check_patients_inventory.R ------------------------------------------
#
# DESCRIPTION: Checks if submitted by user patients inventory is valid.
#              Checks that 1) inventory has all needed columns participant_id, 
#                             tumor_subtype, somatic_path, somatic_genome, 
#                             cohort_name 
#                          2) all participant_id - tumor_subtype pairs are
#                             unique 
#                          3) tumor subtypes have > min_n_participants patients
#                          4) mutation files exist 
#                          5) genome version is the same for all variant files
#                          6) tumor_subtype do not contain - sign because it's
#                             going to be used in file naming 
#                          7) values in tumor_subtype column are not numbers
#             
# USAGE: Rscript --vanilla 1_check_patients_inventory.R \
#                --inventory_patients [path to your file]
#
# OPTIONS: --inventory_patients
#
# REQUIREMENTS: argparse, data.table
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  22.06.2023
# REVISION: 05.07.2023

box::use(./custom_functions[...])
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

# Functions -------------------------------------------------------------------
#' checkParticipantInventory
#' @description Cecks for validity inventory patient table.
#' @author Maria Litovchenko
#' @param pathToInventory path to file with participants inventory
#' @param cancer_subtype list of vector of cancer subtypes for which a subtype
#'        specific analysis should be performed
#' @param min_n_participants integer, minimal number of participants for cancer
#'        subtype specific analysis
#' @param cores integer, number of cores to use
#' @return void, if patient inventory is ok and terminates the script otherwise
#' @note checks that 1) inventory has all needed columns participant_id, 
#' tumor_subtype, participant_tumor_subtype, somatic_path, somatic_genome,
#' cohort_name 2) all participant_id - tumor_subtype pairs are unique 3) tumor
#' subtypes have > min_n_participants patients 4) mutation files exist 
#' 5) genome version is the same for all variant files 5) tumor_subtype do not
#' contain - sign because it's going to be used in file naming 6) values in
#' tumor_subtype column are not number
checkParticipantInventory <- function(pathToInventory, cancer_subtype = NULL, 
                                      min_n_participants = 1, cores = 1) {
  # read and check that all needed columns are present
  essenCols <- c('tumor_subtype', 'participant_tumor_subtype', 
                 'participant_id', 'somatic_path', 'somatic_genome',
                 'cohort_name')
  result <- checkHaveEssentialColumns(pathToInventory, essenCols,
                                      'participant', cores)
  
  # check, that values in tumor_subtype column are not numbers
  checkColumnNotNumber(result$tumor_subtype, 'tumor_subtype')
  
  # remove duplicated entries
  nbefore <- nrow(result)
  result <- result[!duplicated(result[,.(participant_id, tumor_subtype)])]
  nafter <- nrow(result)
  if (nafter < nbefore) {
    warning(paste0('[', Sys.time(), '] Removed ', nbefore - nafter, '(',
                   round(100 * (nbefore - nafter) / nbefore, 2), '%) entries ',
                   'from patient inventory file due to duplicated ',
                   'participant_id - tumor_subtype pairs'))
  }
  
  # determine if tumor type specific analysis was requested
  if (!is.null(cancer_subtype)) {
    nbefore <- nrow(result)
    result <- result[tumor_subtype %in% cancer_subtype]
    nafter <- nrow(result)
    message('[', Sys.time(), '] Selected ', nafter, ' out of ', nbefore, 
            'participants due request on ',
            paste(cancer_subtype, collapse = ','),' cancer type(s)')
    
  }
  
  # check that we have at least a minimum number of patients per tumor subtype
  n_particip_per_subtype <- result[,.(length(unique(participant_id))), 
                                   by = tumor_subtype]
  n_particip_per_subtype <- n_particip_per_subtype[V1 < min_n_participants]
  if (nrow(n_particip_per_subtype) != 0) {
    stop('[', Sys.time(), '] Less than ', min_n_participants, ' participants ',
         'found for ', 
         paste0(n_particip_per_subtype$tumor_subtype, collapse = ', '),
         ' in ', pathToInventory)
  }
  
  # check, that tumor_subtype does not contain - sign because it's going to 
  # be used in file naming
  if (any(grepl('-', result$tumor_subtype))) {
    stop('[', Sys.time(), '] values of tumor_subtype columns should not ',
         'contain - character as it will be used in future file names as ',
         'separator. Please rename. Offending values: ',
         paste(result$tumor_subtype[grepl('-', result$tumor_subtype)],
               collapse = ', '), '.')
  }
  
  # check that input mutation files exist
  fileExist <- file.exists(result$somatic_path)
  if (!all(fileExist)) {
    stop('[', Sys.time(), '] Files do not exist: ', 
         paste(result[!fileExist]$somatic_path, collapse = ',\n'))
  }
  
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = '1_check_patients_inventory.R')

inventoryHelp <- paste('Path to inventory table listing information about ',
                       'patients, their cancer types and mutation files. ',
                       'Minimal columns: participant_id, tumor_subtype, ',
                       'participant_tumor_subtype, somatic_path, ',
                       'somatic_genome, cohort_name')
parser$add_argument("-p", "--inventory_patients", required = T, 
                    type = 'character', help = inventoryHelp)

parser$add_argument("-a", "--inventory_analysis", required = T, 
                    type = 'character', help = 'a')

parser$add_argument("-b", "--inventory_blacklisted", required = F, 
                    type = 'character', help = 'b')

args <- parser$parse_args()

# Check the inventory----------------------------------------------------------
checkParticipantInventory(args$inventory_patients)
