#!/usr/bin/env Rscript
# FILE: create_FPKM_table_TCGA.R ----------------------------------------------
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

sampleSheetHelp <- paste('Path to TCGA sample sheet file')
parser$add_argument("-s", "--sample_sheet", required = T, 
                    type = 'character', help = sampleSheetHelp)

manifestHelp <- paste('Path to TCGA manifest file')
parser$add_argument("-m", "--manifest", required = T, 
                    type = 'character', help = manifestHelp)

projectHelp <- paste('Projects IDs to select, i.e. TCGA-LUAD, TCGA-LUSC')
parser$add_argument("-p", "--project_id", required = F, default = NULL,
                    type = 'character', help = projectHelp, nargs = '*')

outputHelp <- paste('Path to output resulting table')
parser$add_argument("-o", "--output", required = T, 
                    type = 'character', help = outputHelp)
args <- parser$parse_args()

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test inputs -----------------------------------------------------------------
# args <- list(sample_sheet = '/Users/maria/Downloads/gdc_sample_sheet.2023-11-15.tsv',
#              manifest = '/Users/maria/Downloads/gdc_manifest_20231115_125954.txt',
#              project_id = list('TCGA-LUAD', 'TCGA-LUSC'),
#              sample_type = list("Primary Tumor"), 
#              data_type = list("Allele-specific Copy Number Segment",
#                               "Masked Somatic Mutation",
#                               "Gene Expression Quantification"),
#              n_samples = 100,
#              output_sample_sheet = 'gdc_sample_sheet_test_cancer_drivers_pipeline.tsv',
#              output_manifest = 'gdc_manifest_test_cancer_drivers_pipeline.txt')

# Read in & filter sample sheet -----------------------------------------------
sample_sheet <- fread(args$sample_sheet, header = T, stringsAsFactors = F)
if (!is.null(args$project_id)) {
  sample_sheet <- sample_sheet[`Project ID` %in% args$project_id]
}
sample_sheet <- do.call(rbind, apply(sample_sheet, 1, splitRow, ', '))
if (!is.null(args$sample_type)) {
  sample_sheet <- sample_sheet[`Sample Type` %in% args$sample_type]
}
sample_sheet <- sample_sheet[`Data Type` %in% args$data_type]

# Select cases (tumors) for which all data_type are present -------------------
n_files <- sample_sheet[,.(length(unique(`Data Type`))), by = `Case ID`]
n_files <- n_files[V1 == length(args$data_type)]
sample_sheet <- sample_sheet[`Case ID` %in% n_files$`Case ID`]

# Read in manifest ------------------------------------------------------------
manifest <- fread(args$manifest, header = T, stringsAsFactors = F)
manifest <- manifest[id %in% sample_sheet$`File ID`]
setkey(sample_sheet, 'File ID')
manifest <- cbind(manifest, sample_sheet[manifest$id][, c('Case ID',
                                                          'Data Category'),
                                                      with = F])
manifest <- manifest[order(-size)]

# Select n_samples cases (tumors) which have the biggest sizes ----------------
selected_cases <- unique(manifest$`Case ID`)[1:args$n_samples]

manifest <- manifest[`Case ID` %in% selected_cases]
manifest <- split(manifest, manifest$`Data Category`)
manifest <- lapply(manifest, 
                   function(x) x[order(-size)][,.SD[1], by = `Case ID`])
manifest <- do.call(rbind, manifest)

sample_sheet <- sample_sheet[`File ID` %in% manifest$id]

# Output to file --------------------------------------------------------------
write.table(sample_sheet, args$output_sample_sheet, append = F, quote = F,
            sep = '\t', col.names = T, row.names = F)
manifest[, `Case ID` := NULL]
manifest[, `Data Category` := NULL]
write.table(manifest, args$output_manifest, append = F, quote = F,
            sep = '\t', col.names = T, row.names = F)

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')
