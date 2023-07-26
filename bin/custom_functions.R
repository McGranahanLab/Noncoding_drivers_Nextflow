#!/usr/bin/env Rscript
# FILE: 0_functions.R ---------------------------------------------------------
#
# DESCRIPTION: All custom created R functions to perform de-novo cancer driver 
#              genetic element discovery.
#             
# USAGE: In R script, insert:
#        source('0_functions.R')
#
# OPTIONS: None
#
# REQUIREMENTS: R v4.1.0
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  22.06.2023
# REVISION: 05.07.2023

box::use(argparse[...])
box::use(data.table[...])

# GLOBAL ARGUMENTS ------------------------------------------------------------
#' @export
GR_CODES <- c('protein_coding', "3primeUTR", "5primeUTR", "CDS", "lincRNA", 
              "lincRNA_promoter", "lincRNA_ss", "miRNA", "misc_RNA", 
              "promoter", "rRNA", "snoRNA", "snRNA", "ss")

#' @export
SOFTWARE_GR_CODES <- list('activedriverwgs' = NULL, 'dndscv' = c('CDS'),
                          'mutpanning' = c('CDS'), 'chasmplus' = c('CDS'), 
                          'driverpower' = NULL, 'nbr' = NULL, 
                          'oncodrivefml' = NULL, 'oncodriveclustl' = NULL,
                          'digdriver' = NULL)
#' @export
acceptedChrNames <- c(c(1:24, 'X', 'Y'), paste0('chr', c(1:24, 'X', 'Y')))

# Parcing input arguments -----------------------------------------------------
#' check_file_existence_msg
#' @description Checks, whatever or not file exist and stops the execution if
#'              it is not. Prints an error.
#' @author Maria Litovchenko
#' @param filePath path to file
#' @return void
#' @export
check_file_existence_msg <- function(filePath) {
  if (!file.exists(filePath)) {
    stop('[', Sys.time(), '] File does not exist: ', filePath)
  }
}

#' check_min_int_msg
#' @description Checks, whatever numerical argument submitted via command line 
#'              to the script is not less than certain value. Stops and prints
#'              an error if not
#' @author Maria Litovchenko
#' @param value a value of submitted argument
#' @param value_name a name of submitted argument
#' @param min_value minimal value of submitted argument
#' @return void
#' @export
check_min_int_msg <- function(value, value_name, min_value) {
  if (value < min_value) {
    stop('[', Sys.time(), '] Minimal value for the argument --', value_name, 
         ' is ', min_value)
  }
}

#' check_input_arguments
#' @description Checks validity of submitted command line arguments. Terminates
#'              script if something is wrong
#' @author Maria Litovchenko
#' @param argsList list of submitted command line arguments
#' @param outputType string, type of output argument: folder or file
#' @return void
#' @export
check_input_arguments <- function(argsList, outputType = NULL) {
  if (!is.null(argsList$inventory_patients)) {
    lapply(argsList$inventory_patients, check_file_existence_msg)
  }
  if (!is.null(argsList$inventory_analysis)) {
    lapply(argsList$inventory_analysis, check_file_existence_msg)
  }
  if (!is.null(argsList$inventory_blacklisted)) {
    lapply(argsList$inventory_blacklisted, check_file_existence_msg)
  }
  if (!is.null(argsList$min_n_participants)) {
    check_min_int_msg(argsList$min_n_participants, 'min_n_participants', 1)
  }
  if (!is.null(argsList$min_depth)) {
    check_min_int_msg(argsList$min_depth, 'min_depth', 1)
  }
  if (!is.null(argsList$min_tumor_vac)) {
    check_min_int_msg(argsList$min_tumor_vac, 'min_tumor_vac', 1)
  }
  if (!is.null(argsList$min_tumor_vaf)) {
    check_min_int_msg(argsList$min_tumor_vaf, 'min_tumor_vaf', 0)
  }
  if (!is.null(argsList$max_germline_vaf)) {
    check_min_int_msg(argsList$max_germline_vaf, 'max_germline_vaf', 0)
  }
  if (!is.null(argsList$max_germline_vac)) {
    check_min_int_msg(argsList$max_germline_vac, 'max_germline_vac', 0)
  }
  if (!is.null(argsList$max_n_vars)) {
    check_min_int_msg(argsList$max_n_vars, 'max_n_vars', 0)
  }
  if (!is.null(argsList$min_reg_len)) {
    check_min_int_msg(argsList$min_reg_len, 'min_reg_len', 0)
  }
  if (!is.null(argsList$target_genome_path)) {
    check_file_existence_msg(argsList$target_genome_path)
  }
  if (!is.null(argsList$chain)) {
    check_file_existence_msg(argsList$chain)
  }
  
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
  
  if (!is.null(argsList$cores)) {
    check_min_int_msg(argsList$cores, 'cores', 1)
  }
}

#' printArgs
#' @description Prints submitted to the script arguments as messages.
#' @author Maria Litovchenko
#' @param argsList named list representing submitted to the script arguments
#' @return void
#' @export
printArgs <- function(argsList) {
  message("Submitted arguments:")
  for (argName in names(args)) {
    oneArg <- args[[argName]]
    if (length(oneArg) > 1 & !is.null(names(oneArg))) {
      msg <- paste(apply(data.table(names(oneArg), oneArg), 1, paste, 
                         collapse = ' - '), collapse = ',')
    } else {
      msg <- paste(oneArg, collapse = ', ')
    }
    message(argName, ':', msg)
  }
}

