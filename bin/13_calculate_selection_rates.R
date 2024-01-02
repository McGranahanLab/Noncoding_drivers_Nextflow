#!/usr/bin/env Rscript
# FILE: calculate_selection_rates.R -------------------------------------------
#
# DESCRIPTION: Calculates selection rates of cancer driver genomic elements 
# detected in the tumor subtypes of uniform histological origin.
#
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla calculate_selection_rates.R -h
#          to see the full list of options and their descriptions.
# EXAMPLE: 
# REQUIREMENTS: data.table, plyr
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  17.12.2020
# REVISION: 26.12.2023

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
Sys.sleep(sample(1:10, 1))
source(paste0(srcDir, '/custom_functions_postprocessing.R'))

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))

# Functions : calculate selection rates ---------------------------------------
#' calculateGlobalMisTru
#' @description computes global dN/dS ratios from all genes except gene of 
#' interest (to normalise the differences for the genes being tested in 
#' function compare_dNdScv_selection)
#' @author Maria Litovchenko & Ariana Huebner
#' @param raw_results_DT
#' @param goi gene_id of interest
#' @param goi_gr_id genomic region of regions of interest
#' @return data table with columns tumor_subtype, gr_id, gene_id, n_mis/n_sub,
#' n_non, n_spl/n_ind, exp_mis/exp_sub, exp_non, exp_spl/exp_indels, 
#' wmis_global/wsub_global, wtru_global/wind_global
calculateGlobalMisTru <- function(raw_results_DT, goi, goi_gr_id) {
  result <- data.table(gene_id = goi, gr_id = goi_gr_id) 
  if ('exp_non' %in% colnames(raw_results_DT)) { #dNdScv
    wmis <- sum(raw_results_DT[gene_id != goi]$n_mis)
    wmis <- wmis / sum(raw_results_DT[gene_id != goi]$exp_mis)
    wtru <- sum(rowSums(raw_results_DT[gene_id != goi][,.(n_non, n_spl)]))
    wtru <- wtru / sum(rowSums(raw_results_DT[gene_id != goi][,.(exp_non, 
                                                                 exp_spl)]))
    result <- cbind(result, wmis_global = wmis, wtru_global = wtru)
  } else { # NBR
    wsub <- sum(raw_results_DT[gene_id != goi]$n_subs)
    wsub <- wsub / sum(raw_results_DT[gene_id != goi]$exp_subs)
    wind <- sum(raw_results_DT[gene_id != goi]$n_indels)
    wind <- wind / sum(raw_results_DT[gene_id != goi]$exp_indels)
    result <- cbind(result, wsub_global = wsub, wind_global = wind)
  }
  colsToKeep <- intersect(c('tumor_subtype', 'gene_id',
                            'n_mis', 'n_non', 'n_spl', 'n_subs', 'n_ind',
                            'exp_mis', 'exp_non', 'exp_spl', 'exp_subs', 
                            'exp_indels'), colnames(raw_results_DT))
  result <- merge(raw_results_DT[, colsToKeep, with = F], result,
                  by = 'gene_id')
  result
}

#' calculateGlobalMisTruList
#' @description Applies function calculateGlobalMisTru to list of data tables
#' containing raw results of driver discovery by dNdScv/NBR on one genomic 
#' region and one tumor subtype
#' @author Maria Litovchenko & Ariana Huebner
#' @param goi gene_id of interest
#' @param goi_gr_id genomic region of regions of interest
#' @param raw_results_list 
#' @return data table with columns tumor_subtype, gr_id, gene_id, n_mis, 
#' exp_mis, n_non, exp_non, n_spl, exp_spl, wmis_global, wtru_global
calculateGlobalMisTruList <- function(goi, goi_gr_id, raw_results_list) {
  raw_results_goi <- lapply(raw_results_list, 
                            function(x) x[gr_id == goi_gr_id])
  raw_results_goi <- raw_results_goi[sapply(raw_results_goi, 
                                            function(x) nrow(x) != 0)]
  raw_results_goi <- raw_results_goi[sapply(raw_results_goi, 
                                            function(x) goi %in% x$gene_id)]
  result <- lapply(raw_results_goi, calculateGlobalMisTru, goi, goi_gr_id)
  result <- as.data.table(do.call(rbind.fill, result))
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'calculate_selection_rates.R')

driversHelp <- paste('Path to file(s) containing de-novo detected cancer drivers')
parser$add_argument("-d", "--drivers", required = T, nargs = '+', 
                    type = 'character', help = driversHelp)

subtypeHelp <- paste('A cancer subtype to select from patientsInv table. Only',
                     'mutations from patients with that cancer type will be',
                     'selected. In case an analysis of several cancer types',
                     'needed to be performed please run this script ',
                     'separetedly for each cancer type.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    help = subtypeHelp)

grIDhelp <- paste0('Genomic regions IDs corresponding to files submitted in',
                   '--run_results argument. Each genomic region ID should be',
                   'mentioned only once.')
parser$add_argument("-g", "--gr_id", required = T, type = 'character',
                    help = grIDhelp, nargs = '+')

resHelp <- paste('A list of file paths to results of software runs on current',
                 'cancer subtype genomic region pair. The order of paths',
                 'should match the order of software names in --software',
                 'argument')
parser$add_argument("-r", "--run_results", required = T, type = 'character',
                    nargs = '+', help = resHelp)

outputHelp <- paste('Path to the output file')
parser$add_argument("-o", "--output", required = T, 
                    type = 'character', help = outputHelp)

args <- parser$parse_args()
# check_input_arguments_postproc(args, outputType = 'file')
names(args$run_result) <- unlist(args$gr_id)
# check that each each gr_id is only one

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test inputs -----------------------------------------------------------------
# submit only drivers for single origin tumor subtypes
# args <- list(drivers = list('completed_runs/2023-12-14/results/tables/drivers/drivers-Adenocarcinoma--hg19.csv',
#                             'completed_runs/2023-12-14/results/tables/drivers/drivers-Squamous_cell--hg19.csv'),
#              cancer_subtype = 'Adenocarcinoma',
#              gr_id = list('CDS', '3primeUTR', '5primeUTR', 'enhancer',
#                           'lincRNA', 'lincRNA_promoter', 'lincRNA_ss',
#                           'miRNA', 'promoter', 'shortRNA', 'ss'),
#              run_result = list('completed_runs/2023_12_25/results/dndscv/dndscvResults-Adenocarcinoma-CDS-hg19.csv', 
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-3primeUTR-hg19.csv',
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-5primeUTR-hg19.csv',
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-enhancer-hg19.csv',
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-lincRNA-hg19.csv',
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-lincRNA_promoter-hg19.csv',
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-lincRNA_ss-hg19.csv',
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-miRNA-hg19.csv',
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-promoter-hg19.csv',
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-shortRNA-hg19.csv',
#                                'completed_runs/2023-12-14/results/nbr/nbrResults-Adenocarcinoma-ss-hg19.csv'),
#              output = 'selectionRates-Adenocarcinoma--hg19.csv')

# Read driver genes -----------------------------------------------------------
message('[', Sys.time(), '] Started reading ', paste0(args$drivers, 
                                                      collapse = ', '))
drivers <- lapply(args$drivers, fread, header = T, stringsAsFactors = F, 
                  select = c('tumor_subtype', 'gr_id', 'gene_id', 'gene_name',
                             'tier', 'FILTER', 'is_known_cancer', 
                             'known_in_tumor_subtype'))
drivers <- do.call(rbind, drivers)
drivers <- unique(drivers)
drivers <- drivers[FILTER == 'PASS' & !is.na(tier)]
drivers[, FILTER := NULL]
message('[', Sys.time(), '] Finished reading ',
        paste0(args$drivers, collapse = ', '))

if (nrow(drivers) == 0) {
  stop('[', Sys.time(), '] no significant (FILTER is PASS and tier is not ',
       'NA) driver genes is found in ', args$drivers, ' table.')
}

# Read in observed and estimated number of driver mutations - dNdScv, NBR -----
message('[', Sys.time(), '] Started reading ',
        paste(args$run_result, collapse = ','))
mleDriveMuts <- lapply(names(args$run_result), 
                       function(x) readExpObsNmuts(args$run_result[[x]],
                                                   c(cancer_subtype = args$cancer_subtype,
                                                     gr_id = x)))
mleDriveMuts <- mleDriveMuts[sapply(mleDriveMuts, function(x) nrow(x) != 0)]
names(mleDriveMuts) <- names(args$run_result)
message('[', Sys.time(), '] Finished reading ',
        paste(args$run_result, collapse = ','))

# Compare dNdScv/NBR selection rates for each tumor subtype pair --------------
globalMisTru <- apply(unique(drivers[,.(gr_id, gene_id)]), 1, 
                      function(x) calculateGlobalMisTruList(x['gene_id'],
                                                            x['gr_id'],
                                                            mleDriveMuts))
globalMisTru <- as.data.table(do.call(rbind.fill, globalMisTru))
setnames(globalMisTru, 'tumor_subtype', 'scored_in_tumor_subtype')
globalMisTru <- merge(drivers, globalMisTru, by = c('gr_id', 'gene_id'), 
                      all = T)

# Write to output -------------------------------------------------------------
write.table(globalMisTru, args$output, append = F, quote = F, sep = '\t',
            row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote selection rates in ', args$cancer_subtype,
        'of detected in histologically uniform cancer subtypes driver ',
        'genomic elements to ', args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')