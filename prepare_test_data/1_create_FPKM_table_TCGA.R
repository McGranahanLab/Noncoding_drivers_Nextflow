#!/usr/bin/env Rscript
# FILE: create_FPKM_table_TCGA.R ----------------------------------------------
#
# DESCRIPTION: Creates two FPKM tables. First one with columns gene_id, 
# gene_name and one column per sample from TCGA data. Second one with columns
# gene_id, gene_name and one column per tumor subtype.
#
# USAGE: Run Rscript --vanilla create_FPKM_table_TCGA.R -h to see all options
# OPTIONS:
# EXAMPLE: 
# REQUIREMENTS: argparse, data.table
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  15.11.2023
# REVISION: 15.11.2023

# LIBRARIES -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

# FUNCTIONS -------------------------------------------------------------------
#' printArgs
#' @description Prints submitted to the script arguments as messages.
#' @author Maria Litovchenko
#' @param argsList named list representing submitted to the script arguments
#' @return void
#' @export
printArgs <- function(argsList) {
  message('[', Sys.time(), '] Submitted arguments: ')
  for (argName in names(argsList)) {
    oneArg <- argsList[[argName]]
    if (length(oneArg) > 1 & !is.null(names(oneArg))) {
      msg <- paste(apply(data.table(names(oneArg), oneArg), 1, paste, 
                         collapse = ' - '), collapse = ',')
    } else {
      msg <- paste(oneArg, collapse = ', ')
    }
    message('\t', argName, ':', msg)
  }
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'create_FPKM_table_TCGA.R')

sampleSheetHelp <- 'Path to TCGA sample sheet file'
parser$add_argument("-s", "--sample_sheet", required = T, 
                    type = 'character', help = sampleSheetHelp)

folderHelp <- paste('Path to folder where TCGA data were downloaded using',
                    'GDC Data Transfer Tool')
parser$add_argument("-f", "--folder", required = T, 
                    type = 'character', help = folderHelp)

subtypeHelp <- paste('Column name --sample_sheet table which contains tumor',
                     'subtype information. Default: "Project ID"')
parser$add_argument("-t", "--tumor_subtype_column", required = T, 
                    default = "Project ID", type = 'character', 
                    help = subtypeHelp)

aggrHelp <- paste('Boolean, indicating, whether or not an aggregate',
                  'across all tumor subtypes should be computed.')
parser$add_argument("-a", "--do_aggregate", required = F, default = 'F',
                    choices = c('T', 'F'), type = 'character', help = aggrHelp)

aggrNameHelp <- paste('Aggregated tumor subtype name, in case --do_aggregate',
                      'was requested.')
parser$add_argument("-an", "--aggregate_name", required = F, default = NULL,
                    type = 'character', help = aggrNameHelp)

outSampleHelp <- paste('Path to output table which will contain FPKM, one',
                       'column per sample')
parser$add_argument("-os", "--output_per_sample", required = T, 
                    type = 'character', help = outSampleHelp)

outSampleHelp <- paste('Path to output table which will contain FPKM, one',
                       'column per tumor subtype')
parser$add_argument("-ot", "--output_per_tumor_subtype", required = T, 
                    type = 'character', help = outSampleHelp)

args <- parser$parse_args()

args$do_aggregate <- as.logical(args$do_aggregate)
if (args$do_aggregate) {
  if (is.null(args$aggregate_name)) {
    stop('[', Sys.time(), '] --do_aggregate is set to T, but ',
         '--aggregate_name is not given.')
  } 
  if (args$aggregate_name == '') {
    stop('[', Sys.time(), '] --do_aggregate is set to T, but ',
         '--aggregate_name is an empty string.')
  }
}

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test inputs -----------------------------------------------------------------
# args <- list(sample_sheet = 'gdc_sample_sheet_test_cancer_drivers_pipeline.tsv',
#              folder = 'test_dataset/', 
#              data_type = "Gene Expression Quantification",
#              tumor_subtype_column = 'Project ID',
#              do_aggregate = T, aggregate_name = 'PANLUNG',
#              output_per_sample = '../data/assets_raw/TCGA_FPKM_per_sample.csv',
#              output_per_tumor_subtype = '../data/assets_raw/TCGA_FPKM_per_tumor_subtype.csv')

# Read sample_sheet  ----------------------------------------------------------
sample_sheet <- fread(args$sample_sheet, header = T, stringsAsFactors = F)
sample_sheet <- sample_sheet[`Data Type` == args$data_type]
setnames(sample_sheet, args$tumor_subtype_column, 'tumor_subtype')
message('[', Sys.time(), '] Read ', args$sample_sheet)
# create paths to corresponding files
sample_sheet[, filePath := paste0(args$folder, '/', `File ID`, '/', 
                                  `File Name`)]
# check that files exist
fileFound <- sapply(sample_sheet$filePath, file.exists)
if (!any(fileFound)) {
  stop('[', Sys.time(), '] File(s) not found: ', 
       paste0(sample_sheet$filePath[!fileFound], collapse = ', '))
}

# Read the expression tables --------------------------------------------------
message('[', Sys.time(), '] Started reading expression tables')
colsToKeep <- c('gene_id', 'gene_name', 'fpkm_unstranded')
exprTab <- lapply(unique(sample_sheet$tumor_subtype),
                  function(x)  apply(sample_sheet[tumor_subtype == x], 1, 
                                     function(y) fread(y['filePath'],
                                                       header = T,
                                                       stringsAsFactors = F, 
                                                       select = colsToKeep, 
                                                       col.names = c('gene_id',
                                                                     'gene_name',
                                                                     y['Sample ID']))))
exprTab <- lapply(exprTab, 
                  function(x) Reduce(function(...) merge(..., all = T, 
                                                         by = c('gene_id', 
                                                                'gene_name')), 
                                     x))
names(exprTab) <- unique(sample_sheet$tumor_subtype)
message('[', Sys.time(), '] Finished reading expression tables')

# Compute median expression per tumor subtype ---------------------------------
exprTab <- lapply(exprTab, 
                  function(x) setcolorder(x, c('gene_id', 'gene_name')))
medianExprTab <- lapply(exprTab, 
                        function(x) cbind(x[, 1:2], apply(x[, -2:-1], 1, 
                                                          median, na.rm = T)))
medianExprTab <- Reduce(function(...) merge(..., all = T, 
                                            by = c('gene_id', 'gene_name')), 
                        medianExprTab)
setnames(medianExprTab, setdiff(colnames(medianExprTab),
                                c('gene_id', 'gene_name')), names(exprTab))
message('[', Sys.time(), '] Computed median expression per tumor subtype')

# Perform aggregation, if requested -------------------------------------------
exprTab <- Reduce(function(...) merge(..., all = T, 
                                      by = c('gene_id', 'gene_name')), exprTab)
if (args$do_aggregate) {
  aggMedianExprTab <- cbind(exprTab[, 1:2], 
                            apply(exprTab[, -2:-1], 1, median, na.rm = T))
  setnames(aggMedianExprTab, setdiff(colnames(aggMedianExprTab),
                                     c('gene_id', 'gene_name')), 
           args$aggregate_name)
  medianExprTab <- merge(medianExprTab, aggMedianExprTab, all = T,
                         by = c('gene_id', 'gene_name'))
  message('[', Sys.time(), '] Computed aggregated across tumor subtypes ', 
          'median expression')
}

# Output to file --------------------------------------------------------------
write.table(exprTab, args$output_per_sample, append = F, quote = F, sep = '\t', 
            col.names = T, row.names = F)
message('[', Sys.time(), '] Wrote to ', args$output_per_sample)
write.table(medianExprTab, args$output_per_tumor_subtype, append = F, 
            quote = F, sep = '\t', col.names = T, row.names = F)
message('[', Sys.time(), '] Wrote to ', args$output_per_tumor_subtype)

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')