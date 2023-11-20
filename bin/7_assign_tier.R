#!/usr/bin/env Rscript
# FILE: 7_assign_tier.R -------------------------------------------------------
#
# DESCRIPTION: A script to assign tiers to raw de-novo detected cancer driver
# genes. Tiers are assigned based on adjusted for multiple testing combined 
# p-values (i.e. Brown method) for all scanned genomic regions for a given 
# tumor subtype (one tumor subtype).
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
# CREATED:  17.12.2020
# REVISION: 17.11.2023

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

# FUNCTIONS: assign tier ------------------------------------------------------
#' is_tier
#' @author Maria Litovchenko
#' @description Checks if all criteria for tier assigment are satisfied or not
#' @param tierDefVect
#' @param combAdjPdt data table with results of de-novo driver gene detection.
#' It must have columns containing raw p-values for all used de-novo driver 
#' callers, FDR-adjusted p-values for all used de-novo driver callers, column 
#' with combined p-values and a column indicating if gene is known cancer gene
#' or not.
#' @param rawPcols pattern which column names containing raw p-values for each
#' method should follow.
#' @param fdrPcols pattern which column names containing FDR corrected 
#' p-values for each method should follow.
#' @param mergedPcol a column name carrying FDR adjusted combined p-values
#' @param knownCancerGeneCol column name which indicates whatever or not gene 
#' is a known cancer gene.
#' @return vector of T or F indicating for every row in combAdjPdt if all 
#' criteria for that gene to be driver gene are satisfied or not.
is_tier <- function(tierDefVect, combAdjPdt, rawPcols, fdrPcols, mergedPcol,
                    knownCancerGeneCol) {
  cutoffs <- suppressWarnings(as.numeric(tierDefVect))
  names(cutoffs) <- names(tierDefVect)
  
  # 1. check that raw p-values of individual methods are < cutoff and that 
  # number of those methods is at least nIndivFDRsoft_sign
  passRaw <- combAdjPdt[, rawPcols, with = F] < cutoffs['indivRaw_cutoff']
  passRaw <- rowSums(passRaw, na.rm = T)
  passRaw <- passRaw >= cutoffs["nIndivRawSoft_sign"]
  
  # 2. check that FDR corrected p-values of individual methods are < cutoff
  # and that number of those methods is at least nIndivFDRsoft_sign
  passFDR <- combAdjPdt[, fdrPcols, with = F] < cutoffs["indivFDR_cutoff"]
  passFDR <- rowSums(passFDR, na.rm = T)
  passFDR <- passFDR >= cutoffs["nIndivFDRsoft_sign"]
  
  # 3. check that merged FDR corrected combined p-value is < mergedFDR_cutoff
  passMerged <- unlist(combAdjPdt[, mergedPcol, with = F])
  passMerged <- passMerged < cutoffs["mergedFDR_cutoff"]
  
  result <- passRaw & passFDR & passMerged
  
  # restrict to known cancer genes, if required
  if (as.logical(tierDefVect["restrictToKnownCancer"])) {
    result[which(!unlist(combAdjPdt[, knownCancerGeneCol, with = F]))] <- F
  }
  result
}

#' assignTier
#' @description Assigns tier of significance for gene based on tier defining 
#' table.
#' @author Maria Litovchenko 
#' @param combAdjPdt data table with results of de-novo driver gene detection.
#' It must have columns containing raw p-values for all used de-novo driver 
#' callers, FDR-adjusted p-values for all used de-novo driver callers, column 
#' with combined p-values and a column indicating if gene is known cancer gene
#' or not.
#' @param tierDefDT data table which defines tiers based on cutoff on 
#' individual FDR methods, merged FDR methods, etc. Should have columns tier,
#' indivRaw_cutoff, nIndivRawSoft_sign, indivFDR_cutoff, nIndivFDRsoft_sign,
#' restrictToKnownCancer and mergedFDR_cutoff
#' @param rawColPtrn pattern which column names containing raw p-values for
#' each method should follow.
#' @param fdrColPtrn pattern which column names containing FDR corrected 
#' p-values for each method should follow.
#' @param mergedPcol a column name carrying FDR adjusted combined p-values
#' @param knownCancerGeneCol column name which indicates whatever or not gene 
#' is a known cancer gene.
#' @return data table with tier column added.
assignTier <- function(combAdjPdt, tierDefDT, rawColPtrn = '[.]raw_p$', 
                       fdrColPtrn = '[.]bh_p$', mergedPcol = 'brown.bh_p', 
                       knownCancerGeneCol = 'is_known_cancer') {
  # columns with raw p-values per driver calling method
  rawCols <- grep(rawColPtrn, colnames(combAdjPdt), value = T)
  # columns with FDR-adjusted p-values per driver calling method
  fdrCols <- grep(fdrColPtrn, colnames(combAdjPdt), value = T)
  
  tierStatus <- apply(tierDefDT, 1, is_tier, combAdjPdt, rawCols, fdrCols, 
                      mergedPcol, knownCancerGeneCol)
  tierStatus <- as.data.table(tierStatus)
  # if a gene has tier of a higher score, it should not have a tier of a lower 
  # one
  tierStatus[, tier := rowSums(tierStatus)]
  tierStatus[tier != 0]$tier <- apply(tierStatus[tier != 0], 1,
                                      function(x) which(x == T)[1])
  tierStatus[tier == 0]$tier <- NA
  
  result <- cbind(combAdjPdt, tier = tierStatus$tier)
  result
}

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))

# Test inputs -----------------------------------------------------------------
args <- list(combined_p_values_tables = list('../TEST/results/tables/combined_p_values/combinedP_LUAD-CDS_hg19.csv',
                                             '../TEST/results/tables/combined_p_values/combinedP_LUAD-5primeUTR_hg19.csv'),
             inventory_tier = '../data/inventory/inventory_tier_definition.csv',
             output = '../TEST/results/tables/rawDrivers_LUAD--_hg19.csv')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Read inventory defining tiers -----------------------------------------------
tierDefine_inventory <- fread(args$inventory_tier, header = T, 
                              stringsAsFactors = F)

# Read tables with combined p values ------------------------------------------
combinedPs <- lapply(args$combined_p_values_tables, fread, header = T, 
                     stringsAsFactors = F)
combinedPs <- do.call(rbind.fill, combinedPs)
combinedPs <- as.data.table(combinedPs)
  
# Adjust FDR each driver calling method & for combined p-values ---------------
# get software names
software <- gsub(".raw_p", "", grep('[.]raw_p$', colnames(combinedPs),
                                    value = T))
# get column names of columns containing p values
pValCols <- grep('[.]comb_p$|[.]raw_p$', colnames(combinedPs), value = T)

# perform multiple testing correction
combinedPsAdj <- apply(combinedPs[, pValCols, with = F], 2, p.adjust, 
                       method = 'BH')
combinedPsAdj <- as.data.table(combinedPsAdj)
colnames(combinedPsAdj) <- gsub('.comb_p|.raw_p', '.bh_p', 
                                colnames(combinedPsAdj))
combinedPsAdj <- cbind(combinedPs, combinedPsAdj) 
rm(combinedPs)

# Assign tier based on all p-value combining methods --------------------------
pValCombMethods <- gsub('.comb_p', '', grep('.comb_p', pValCols, value = T))

for (pCombMet in pValCombMethods) {
  combinedPsAdj <- assignTier(combinedPsAdj, tierDefine_inventory,
                              rawColPtrn = paste0(software, '.raw_p$',
                                                  collapse = '|'),
                              fdrColPtrn = paste0(software, '.bh_p$', 
                                                  collapse = '|'),
                              mergedPcol = paste0(pCombMet, '.bh_p'), 
                              knownCancerGeneCol = 'is_known_cancer')
  setnames(combinedPsAdj, 'tier', paste0(pCombMet, '.tier'))
  message('[', Sys.time(), '] Assigned tier based on ', pCombMet)
}

# [SAVE] tiered results as table ----------------------------------------------
write.table(combinedPsAdj, args$output, append = F, quote = F,  sep = '\t', 
            row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote raw drivers table to ', args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')