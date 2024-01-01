#!/usr/bin/env Rscript
# FILE: custom_functions_postprocessing.R -------------------------------------
#
# DESCRIPTION: Custom created R functions which are used only in
#             
# USAGE: In R script, insert:
#        source('custom_functions_postprocessing.R')
#
# OPTIONS: None
#
# REQUIREMENTS: R v4.1.0
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  07.11.2023
# REVISION: 27.12.2023

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
suppressWarnings(suppressPackageStartupMessages(library(GenomeInfoDb)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(stats)))
suppressWarnings(suppressPackageStartupMessages(library(utils)))

# Parsing input arguments -----------------------------------------------------
#' check_input_arguments_preproc
#' @description Checks validity of submitted command line arguments. Terminates
#'              script if something is wrong
#' @author Maria Litovchenko
#' @param argsList list of submitted command line arguments
#' @param outputType string, type of output argument: folder or file
#' @return void
#' @export
check_input_arguments_preproc <- function(argsList, outputType = NULL) {
  lapply(argsList$inventory_patients, check_file_existence_msg)
  lapply(argsList$inventory_analysis, check_file_existence_msg)
  
  check_min_int_msg(argsList$min_n_participants, 'min_n_participants', 1)
  check_min_int_msg(argsList$min_depth, 'min_depth', 1)
  check_min_int_msg(argsList$min_tumor_vac, 'min_tumor_vac', 1)
  check_min_int_msg(argsList$min_tumor_vaf, 'min_tumor_vaf', 0)
  check_min_int_msg(argsList$max_germline_vaf, 'max_germline_vaf', 0)
  check_min_int_msg(argsList$max_germline_vac, 'max_germline_vac', 0)
  check_min_int_msg(argsList$max_n_vars, 'max_n_vars', 0)
  check_min_int_msg(argsList$min_reg_len, 'min_reg_len', 1)
  check_file_existence_msg(argsList$target_genome_path)
  check_file_existence_msg(argsList$chain)
  
  check_file_existence_msg(argsList$variants)
  check_file_existence_msg(argsList$genomic_regions)
  check_file_existence_msg(argsList$target_genome_chr_len)
  check_file_existence_msg(argsList$gene_name_synonyms)
  check_file_existence_msg(argsList$varanno_conversion_table)
  check_min_int_msg(argsList$bin_len, 'bin_len', 1)
  
  check_file_existence_msg(argsList$maf)
  check_file_existence_msg(argsList$bed)
  
  check_min_int_msg(argsList$cores, 'cores', 1)
  
  if (!is.null(argsList$output)) {
    if (!is.null(outputType)) {
      dirToCheck <- argsList$output
      if (outputType == 'file') {
        dirToCheck <- dirname(argsList$output)
      }
      if (!dir.exists(dirToCheck)) {
        dir.create(dirToCheck, recursive = T)
        message('[', Sys.time(), '] Created directory ', dirToCheck) 
      }
    } else {
      stop('[', Sys.time(), '] Please give outputType argument.')
    }
  }
}

# Functions : reading from dNdScv & NBR ---------------------------------------
#' readExpObsNmuts
#' @description Reads in result files from dNdScv & NBR into uniformal data 
#' table while extracting columns containing information about expected and 
#' observed number of mutations. 
#' @author Maria Litovchenko
#' @param filePath path to result file
#' @param infoStr named vector with at least items cancer_subtype, gr_id
#' @return data table with columns tumor_subtype, gr_id + n_syn, n_mis, n_non, 
#' n_spl, n_ind, mis_mle, tru_mle, mis_low, tru_low, mis_high, tru_high, 
#' ind_mle, ind_low, ind_high
readExpObsNmuts <- function(filePath, infoStr) {
  if (!file.exists(filePath)) {
    message('[', Sys.time(), '] File not found: ', filePath, '. Return empty ',
            'data table.')
    return(data.table())
  }
  
  if (!all(c('cancer_subtype', 'gr_id') %in% names(infoStr))) {
    stop('[', Sys.time(), '] readExpectedNmuts: infoStr vector should have ',
         'items under names cancer_subtype, gr_id')
  }
  infoStr <- infoStr[c('cancer_subtype', 'gr_id')]
  infoStr <- as.data.table(t(infoStr))
  setnames(infoStr, 'cancer_subtype', 'tumor_subtype')
  
  colsToGet <- c('gene_name', 'n_syn', 'n_mis', 'n_non', 'n_spl', 'n_ind', 
                 'mis_mle', 'tru_mle', 'mis_low', 'tru_low', 'mis_high', 
                 'tru_high', 'ind_mle', 'ind_low', 'ind_high', 'exp_syn', 
                 'exp_mis', 'exp_non', 'exp_spl',#colums for dndscv finish here
                 'region', 'obs_subs', 'obs_indels', 'obsexp_subs_mle', 
                 'obsexp_subs_low', 'obsexp_subs_high', 'obsexp_indels_mle',
                 'obsexp_indels_low', 'obsexp_indels_high', 'exp_subs',
                 'exp_indels')
  result <- suppressWarnings(fread(filePath, header = T, stringsAsFactors = F, 
                                   select = colsToGet))
  setnames(result, 
           c('gene_name', 'region', 'obs_subs', 'obs_indels', 
             'obsexp_subs_mle', 'obsexp_subs_low', 'obsexp_subs_high', 
             'obsexp_indels_mle', 'obsexp_indels_low', 'obsexp_indels_high'),
           c('gene_id', 'gene_id', 'n_subs', 'n_ind', 'subs_mle', 'subs_low', 
             'subs_high', 'ind_mle', 'ind_low', 'ind_high'), skip_absent = T)
  result <- cbind(result, infoStr)
  result
}

