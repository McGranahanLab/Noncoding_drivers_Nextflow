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

# Parcing input arguments -----------------------------------------------------
#' check_file_existance_msg
#' @description Checks, whatever or not file exist and stops the execution if
#'              it is not. Prints an error.
#' @author Maria Litovchenko
#' @param filePath path to file
#' @return void
#' @export
check_file_existance_msg <- function(filePath) {
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
#' @return void
#' @export
check_input_arguments <- function(argsList) {
  if (!is.null(args$inventoryP)) {
    lapply(args$inventoryP, check_file_existance_msg)
  }
  if (!is.null(args$inventory)) {
    lapply(args$inventory, check_file_existance_msg)
  }
  if (!is.null(args$blacklist_inventory)) {
    check_file_existance_msg(args$blacklist_inventory)
  }
  if (!is.null(args$chain)) {
    check_file_existance_msg(args$chain)
  }
  if (!is.null(args$min_n_participants)) {
    check_min_int_msg(args$min_n_participants, 'min_n_participants', 1)
  }
  if (!is.null(args$min_depth)) {
    check_min_int_msg(args$min_depth, 'min_depth', 1)
  }
  if (!is.null(args$min_tumor_vac)) {
    check_min_int_msg(args$min_tumor_vac, 'min_tumor_vac', 1)
  }
  if (!is.null(args$min_tumor_vaf)) {
    check_min_int_msg(args$min_tumor_vaf, 'min_tumor_vaf', 0)
  }
  if (!is.null(args$max_germline_vaf)) {
    check_min_int_msg(args$max_germline_vaf, 'max_germline_vaf', 0)
  }
  if (!is.null(args$max_germline_vac)) {
    check_min_int_msg(args$max_germline_vac, 'max_germline_vac', 0)
  }
  if (!is.null(args$max_n_vars)) {
    check_min_int_msg(args$max_n_vars, 'max_n_vars', 0)
  }
  if (!is.null(args$min_reg_len)) {
    check_min_int_msg(args$min_reg_len, 'min_reg_len', 0)
  }
  if (!is.null(args$target_genome_path)) {
    check_file_existance_msg(args$target_genome_path)
  }
  if (!is.null(args$output)) {
    if (!dir.exists(dirname(args$output))) {
      dir.create(dirname(args$output), recursive = T)
      message('[', Sys.time(), '] Created directory ', dirname(args$output)) 
    }
  }
  if (!is.null(args$cores)) {
    check_min_int_msg(args$cores, 'cores', 1)
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

# Checking inventories --------------------------------------------------------
#' checkHaveEssentialColumns
#' @description Checks that all columns needed to be in inventory table are 
#'              indeed present there.
#' @author Maria Litovchenko
#' @param filePath path to inventory
#' @param essentialCols string vector, containing essential column names
#' @param inventory_type string, type of inventory (patient, analysis, etc)
#' @param cores integer, number of cores to be used
#' @return data table, in case all columns are found and stops the script 
#'         otherwise.
#' @export
checkHaveEssentialColumns <- function(filePath, essentialCols, 
                                      inventory_type, cores = 1) {
  result <- suppressWarnings(fread(filePath, header = T, stringsAsFactors = F, 
                                   nThread = cores, select = essentialCols))
  if (!all(essentialCols %in% colnames(result))) {
    stop('[', Sys.time(), '] ', inventory_type, ' inventory must have ',
         'columns: ', paste(essentialCols, collapse = ', '))
  }
  result
}

#' checkColumnNotNumber
#' @description Checks that there are no numbers in the vector 
#' @author Maria Litovchenko
#' @param x vector
#' @param x_name name of that vector, if available
#' @return void, if no number in the vector was found and stops script 
#'         execution otherwise
#' @export
checkColumnNotNumber <- function(x, x_name = NULL) {
  x_to_int <- suppressWarnings(as.integer(x))
  
  if (any(!is.na(x_to_int))) {
    int_in_x <- x[!is.na(x_to_int)]
    stop('[', Sys.time(), '] Found a numerical value(s) in ', x_name, ': ',
         paste(int_in_x, collapse = ', '), '. Not allowed.')
  }
}