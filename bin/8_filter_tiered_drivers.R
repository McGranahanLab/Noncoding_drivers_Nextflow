#!/usr/bin/env Rscript
# FILE: filter_tiered_drivers.R -----------------------------------------------
#
# DESCRIPTION: Applies filters on minimal number of successfully run software,
# minimal number of mutations and/or patients, maximum local mutation rate /
# synonymous mutation rate, maximum quantile of genomic region length to tiered
# drivers.
#
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla filter_tiered_drivers.R -h
#          to see the full list of options and their descriptions.
# EXAMPLE: 
# REQUIREMENTS: argparse, data.table
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

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'filter_tiered_drivers.R')

analysisHelp <- paste('Path to inventory table containing details of the',
                      'future analysis to be conducted. Minimal columns:',
                      'tumor_subtype,', 'software,', 'gr_id,', 'gr_code,', 
                      'gr_file,', 'gr_upstr,', 'gr_downstr,', 'gr_genome,', 
                      'gr_excl_id,', 'gr_excl_code,', 'gr_excl_file,',
                      'gr_excl_upstr,', 'gr_excl_downstr,', 'gr_excl_genome,',
                      'blacklisted_codes.')
parser$add_argument("-a", "--inventory_analysis", required = T, 
                    type = 'character', help = analysisHelp)

subtypeHelp <- paste('A cancer subtype to select from patientsInv table. Only',
                     'mutations from patients with that cancer type will be',
                     'selected. In case an analysis of several cancer types',
                     'needed to be performed please run this script ',
                     'separetedly for each cancer type.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    help = subtypeHelp)

tieredPhelp <- paste('Path to table containing tiered genomic regions,',
                     'result of assign_tier.R script')
parser$add_argument("-t", "--tiered_pvals", required = T, type = 'character',
                    help = tieredPhelp)

minSoftCodHelp <- paste('Minimal number of de-novo driver detecting software',
                        'successfully run on coding genomic regions')
parser$add_argument("-msc", "--min_n_soft_cod", required = T, type = 'integer',
                    help = minSoftCodHelp, default = 3)

minSoftNChelp <- paste('Minimal number of de-novo driver detecting software',
                       'successfully run on noncoding genomic regions')
parser$add_argument("-msnc", "--min_n_soft_noncod", required = T,
                    type = 'integer', default = 2, help = minSoftNChelp)

minMutsHelp <- paste('Minimal number of somatic mutations detected in',
                     'genomic region.')
parser$add_argument("-mm", "--min_n_muts", required = F, type = 'integer',
                    default = 0, help = minMutsHelp)

minPatientsHelp <- paste('Minimal number of patients with a mutation in',
                         'genomic region.')
parser$add_argument("-mp", "--min_n_patients", required = F, type = 'integer',
                    default = 0, help = minPatientsHelp)

maxLocalMutRateQhelp <- paste('Maximum local mutation rate (value between',
                              '0 and 1)')
parser$add_argument("-mlmrq", "--max_local_mut_rate_q", required = F, 
                    type = 'double', default = 1, help = maxLocalMutRateQhelp)

maxGrMutRateQhelp <- paste('Maximum genomic region - specific mutation rate',
                           'quantile (value between 0 and 1)')
parser$add_argument("-mgmrq", "--max_gr_mut_rate_q", required = F, 
                    type = 'double', default = 1, help = maxGrMutRateQhelp)

maxSybGrMutRateQhelp <- paste('Maximum genomic region - specific mutation ',
                              'rate quantile (value between 0 and 1) ',
                              'computed based on synonymous mutations only. ',
                              'This filter will only be applied to coding ',
                              'genomic regions.')
parser$add_argument("-mgsmrq", "--max_gr_syn_mut_rate_q", required = F, 
                    type = 'double', default = 1, help = maxSybGrMutRateQhelp)

maxGrLenQhelp <- paste('Maximum quantile of genomic region length',
                       '(value between 0 and 1)')
parser$add_argument("-mglq", "--max_gr_len_q", required = F, type = 'double',
                    default = 1, help = maxGrLenQhelp)

removeOlfactoryHelp <- paste('Wherther or not genomic regions associated to ',
                             'olfactory genes whould be removed.')
parser$add_argument("-of", "--remove_olfactory", required = F, 
                    type = 'character', default = 'T', choices = c('T', 'F'),
                    help = removeOlfactoryHelp)

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output file")

args <- parser$parse_args()
args$remove_olfactory <- as.logical(args$remove_olfactory)
# check_input_arguments_postproc(args, outputType = 'file')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test arguments --------------------------------------------------------------
# args <- list(inventory_analysis = '../data/inventory/inventory_analysis_tcga_short.csv',
#              cancer_subtype = 'LUAD',
#              tiered_pvals = '../TEST/results/tables/tiered_drivers/assignedTier-LUAD--hg19.csv',
#              min_n_soft_cod = 3,
#              min_n_soft_noncod = 2, min_n_muts = 3, min_n_patients = 3,
#              max_local_mut_rate_q = 1, max_gr_syn_mut_rate_q = 0.99, 
#              max_gr_len_q = 0.99, remove_olfactory = T, 
#              output = 'drivers-LUAD--hg19.csv')

# Read in analysis inventory --------------------------------------------------
analysisInv <- readAnalysisInventory(args$inventory_analysis)
message('[', Sys.time(), '] Read --inventory_analysis: ', 
        args$inventory_analysis)
analysisInv <- analysisInv[tumor_subtype == args$cancer_subtype]
# assign coding genomic regions - anything which contains CDS as gr_code
coding_gr_id <- unique(analysisInv[gr_code == 'CDS']$gr_id)
rm(analysisInv)
gc()

# Read in tiered p-values -----------------------------------------------------
tiered_pvals <- fread(args$tiered_pvals, header = T, stringsAsFactors = F)
message('[', Sys.time(), '] Read --tiered_pvals ', args$tiered_pvals)

# Initialize FILTER -----------------------------------------------------------
tiered_pvals[, FILTER := 'PASS']

# Mark genomic elements not scanned with enough of software -------------------
if (args$min_n_soft_cod > 0 | args$min_n_soft_noncod > 0) {
  # get columns with raw P values
  rawPcols <- grep('[.]raw_p$', colnames(tiered_pvals), value = T)
  tiered_pvals[, nSoft := apply(tiered_pvals[, rawPcols, with = F], 1, 
                                function(x) sum(!is.na(x)))]
  tiered_pvals[, minNsoft := ifelse(gr_id %in% coding_gr_id, 
                                    args$min_n_soft_cod, 
                                    args$min_n_soft_noncod)]
  pass <- tiered_pvals$nSoft >= tiered_pvals$minNsoft
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER,
                                       '; low N. software')
  tiered_pvals[, nSoft := NULL]
  tiered_pvals[, minNsoft := NULL]
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
          'regions which have too low number of successfully run software.')
  rm(pass)
}

# Mark genomic elements with not enough of N mutations or patients ----------
if (args$min_n_muts > 0) {
  pass <- !is.na(tiered_pvals$nMuts_total) & 
    tiered_pvals$nMuts_total >= args$min_n_muts
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER, 
                                       '; low N. mutations')
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
          'regions which have too low number of mutations.')
  rm(pass)
}

if (args$min_n_patients > 0) {
  pass <- !is.na(tiered_pvals$nParts_total) & 
    tiered_pvals$nParts_total >= args$min_n_patients
    
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER,
                                       '; low N. participants') 
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
          'regions which have too low number of mutated participants.')
  rm(pass)
}

# Mark genes with high local mutational rate  -------------------------------
if ('meanMutRateQuant_local' %in% colnames(tiered_pvals)) {
  if (args$max_local_mut_rate_q < 1) {
    pass <- tiered_pvals$meanMutRateQuant_local <= 100 * args$max_local_mut_rate_q |
      is.na(tiered_pvals$meanMutRateQuant_local)
    tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER,
                                     '; high local mut.rate')
    
    message('[', Sys.time(), '] Found ', sum(!pass), '(',
            round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
            'regions which have too high local mutation rate.')
    rm(pass)
  }
} else {
  if (!'meanMutRateQuant_local' %in% colnames(tiered_pvals)) {
    message('[', Sys.time(), '] Filtering based on mean mutation rate ',
            'quantile was not performed, because column ',
            'meanMutRateQuant_local was not found.')
  }
}

# Mark by genomic regions specific mutational rate/synonymous mut.rate --------
if (args$max_gr_mut_rate_q < 1 | args$max_gr_syn_mut_rate_q < 1) {
  # set a flag to determine, if filtering by synonymous mutations should be
  # performed
  is_cod <- tiered_pvals$gr_id %in% coding_gr_id
  doSyn <- F
  if (args$max_gr_syn_mut_rate_q < 1 & 
      'synMeanMutRateQuant' %in% colnames(tiered_pvals)) {
    doSyn <- T
  }
  if (args$max_gr_syn_mut_rate_q < 1 & !doSyn) {
    message('[', Sys.time(), '] Filtering by synonymous mutation rate in ',
            'coding regions will not be performed because column ',
            'synMeanMutRateQuant was not found. Will use instead column ',
            'meanMutRateQuant.')
  }
  
  # Assign dummy column values which will be replaced with real ones in case
  # max_gr_mut_rate_q or max_gr_syn_mut_rate_q is < 1
  tiered_pvals[, meanMutRateQuant_cutOff := 110]
  tiered_pvals[, meanMutRateQuant_toFltr := -1]
  if (args$max_gr_mut_rate_q < 1) {
    tiered_pvals[, meanMutRateQuant_toFltr := meanMutRateQuant]
    tiered_pvals[, meanMutRateQuant_cutOff := 100 * args$max_gr_mut_rate_q]
  }
  if (doSyn) {
    tiered_pvals[is_cod]$meanMutRateQuant_toFltr <- tiered_pvals[is_cod]$synMeanMutRateQuant
    tiered_pvals[is_cod]$meanMutRateQuant_cutOff <- 100 * args$max_gr_syn_mut_rate_q
  }
  
  pass <- tiered_pvals$meanMutRateQuant_toFltr <= tiered_pvals$meanMutRateQuant_cutOff |
    is.na(tiered_pvals$meanMutRateQuant_toFltr)
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER,
                                       '; high genomic region mut.rate')
  tiered_pvals[, meanMutRateQuant_toFltr := NULL]
  tiered_pvals[, meanMutRateQuant_cutOff := NULL]
  if (doSyn) {
    tiered_pvals[is_cod]$FILTER <- gsub('; high genomic region mut.rate',
                                        '; high syn genomic region mut.rate', 
                                        tiered_pvals[is_cod]$FILTER)
    rm(is_cod)
  }
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
          'regions which have too high genomic region - specific ',
          'mutation rate.')
  if (doSyn) {
    message('[', Sys.time(), '] Synonymous mutation rate was used in the ',
            'filter above for coding regions')
  }
  rm(pass)
}

# Mark olfactory genes --------------------------------------------------------
if (args$remove_olfactory & !'is_olfactory' %in% colnames(tiered_pvals)) {
  message('[', Sys.time(), '] Filtering of olfactory genes will not be ',
          'performed because column is_olfactory was not found.')
}
if (args$remove_olfactory & 'is_olfactory' %in% colnames(tiered_pvals)) {
  pass <- !tiered_pvals$is_olfactory
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER,
                                       '; olfactory')
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
          'regions which are linked to olfactory genes.')
  rm(pass)
}

# Mark by quantile of length --------------------------------------------------
if (args$max_gr_len_q < 1) {
  pass <- is.na(tiered_pvals$gr_lenQuant) |
          tiered_pvals$gr_lenQuant <= 100 * args$max_gr_len_q
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER, 
                                       '; long gene')
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
          'regions which are too long.')
  rm(pass)
}

# Mark out gene based on being not expressed ----------------------------------
exprCols <- grep('expr_in', colnames(tiered_pvals), value = T)
if (length(exprCols) > 0) {
  for (eCol in exprCols) {
    eVals <- unlist(tiered_pvals[, eCol, with = F])
    pass <- eVals | is.na(eVals)
    filterName <- paste0('; not expr. in ', gsub('expr_in_', '', eCol))
    tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER,
                                         filterName) 
    
    message('[', Sys.time(), '] Found ', sum(!pass), '(',
            round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
            'regions which are not expressed in ', gsub('expr_in_', '', eCol),
            '.')
    rm(pass)
  }
} 

# Remove PASS at the start of FILTER for entries which didn't pass ------------
tiered_pvals[, FILTER := gsub('^PASS; ', '', FILTER)]

# Print summary table of filtering --------------------------------------------
filterStats <- as.data.table(table(tiered_pvals[FILTER != 'PASS']$FILTER))
# add number of instances with assigned tier which was filtered out
if (nrow(tiered_pvals[!is.na(tier) & FILTER != 'PASS']) > 0) {
  filterStats <- merge(filterStats, 
                       as.data.table(table(tiered_pvals[!is.na(tier) & 
                                                          FILTER != 'PASS']$FILTER)),
                       by = 'V1', all.x = T)
  filterStats[is.na(N.y)]$N.y <- 0
} else {
  filterStats[, N.y := 0]
}
setnames(filterStats, c('V1', 'N', 'N.x', 'N.y'),
         c('Filter', 'N genomic regions', 'N genomic regions', 
           'N genomic regions with assigned tier'), skip_absent = T)
filterStats <- filterStats[order(-`N genomic regions`)]
filterStats[, `% genomic regions` := round(100 * `N genomic regions` / 
                                             nrow(tiered_pvals), 1)]
message(paste0(capture.output(knitr::kable(filterStats, format = "markdown")), 
               collapse = '\n'))

# Save filtered results as a table --------------------------------------------
colOrderToPrint <- c('tumor_subtype', 'gr_id', 'gene_id', 'gene_name', 
                     'is_known_cancer', 'FILTER', 'participant_tumor_subtype', 
                     'nMuts', 'nMutsUniq', 'nParts', 'nMuts_total', 
                     'nMutsUniq_total', 'nParts_total', 
                     'meanMutRate', 'meanMutRateQuant',
                     sort(grep('.raw_p$', colnames(tiered_pvals), 
                               value = T)),
                     sort(grep('.comb_p$', colnames(tiered_pvals), 
                               value = T)),
                     sort(grep('.bh_p$', colnames(tiered_pvals), 
                               value = T)),
                     sort(grep('.tier$', colnames(tiered_pvals), value = T)))
setcolorder(tiered_pvals, colOrderToPrint)
write.table(tiered_pvals, args$output, append = F, quote = F,  sep = '\t', 
            row.names = F, col.names = T)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')