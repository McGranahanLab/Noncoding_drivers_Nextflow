#!/usr/bin/env Rscript
# FILE: filter_tiered_drivers.R -----------------------------------------------
#
# DESCRIPTION: Applies filters on minimal number of successfully run software,
# minimal number of mutations and/or patients, maximum local mutation rate /
# synonymous mutation rate, maximum quantile of genomic region length to tiered
# drivers.
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

useSynHelp <- paste('Whether or not quantile of synonumous mutation rate',
                    'should be used as genomic region - specific mutation',
                    'rate for coding genomic regions.')
parser$add_argument("-syn", "--use_synonumous_mut_rate", required = F, 
                    type = 'character', default = 'F', help = useSynHelp,
                    choices = c('T', 'F'))

maxGrMutRateQhelp <- paste('Maximum genomic region - specific mutation rate',
                           'quantile (value between 0 and 1)')
parser$add_argument("-mgmrq", "--max_gr_mut_rate_q", required = F, 
                    type = 'double', default = 1, help = maxGrMutRateQhelp)

maxGrLenQhelp <- paste('Maximum quantile of genomic region length',
                       '(value between 0 and 1)')
parser$add_argument("-mglq", "--max_gr_len_q", required = F, type = 'double',
                    default = 1, help = maxGrLenQhelp)

removeOlfactoryHelp <- paste('Maximum quantile of genomic region length',
                             '(value between 0 and 1)')
parser$add_argument("-o", "--remove_olfactory", required = F, type = 'character',
                    default = 1, help = removeOlfactoryHelp)


as.logical(args$use_synonumous_mut_rate)
as.logical(args$remove_olfactory)

codingGRid = list('CDS')
min_n_soft_cod = 3
min_n_soft_noncod = 2
min_n_muts = 3
min_n_patients = 3
max_local_mut_rate_q = 0.99
use_synonumous_mut_rate = T
max_gr_mut_rate_q = 0.99
max_gr_len_q = 0.99

# Test arguments --------------------------------------------------------------
args <- list(inventory_analysis = '../data/inventory/inventory_analysis_tcga_short.csv',
             cancer_subtype = 'LUAD',
             tiered_pvals = '../TEST/results/tables/tiered_drivers/assignedTier-LUAD--hg19.csv',
             min_n_soft_cod = 3,
             min_n_soft_noncod = 2, min_n_muts = 3, min_n_patients = 3,
             max_local_mut_rate_q = 0.99, use_synonumous_mut_rate = T,
             max_gr_mut_rate_q = 0.99, max_gr_len_q = 0.99, 
             remove_olfactory = T)

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
    pass <- tiered_pvals$meanMutRateQuant < 100 * args$max_local_mut_rate_q &
      !is.na(tiered_pvals$meanMutRateQuant_local)
    tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER,
                                     '; high local mut.rate')
    
    message('[', Sys.time(), '] Found ', sum(!pass), '(',
            round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
            'regions which have too low number of mutated participants.')
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
if (args$use_synonumous_mut_rate & 
    !'synMeanMutRateQuant' %in% colnames(tiered_pvals)) {
  message('[', Sys.time(), '] ',
}

if (args$use_synonumous_mut_rate) {
  
  if (!c('synMeanMutRateQuant') %in% colnames(tiered_pvals)) {
    msg <- 'Filtering based on mean genome region '
    if (args$use_synonumous_mut_rate) {
      msg <- paste(msg, 'SYNONYMOUS')
    }
    msg <- paste(msg, 'mutation rate quantile was not performed, because ',
                 'column ', colToUse, ' was not found.')
    message('[', Sys.time(), '] ', msg)
  }
  
  
  
  # coding regions
  grMutRate <- unlist(tiered_pvals[, 'synMeanMutRateQuant', with = F])
  # is.na(grMutRate) because a gene can be mutated, but have 0 synonymous muts
  pass <- (tiered_pvals$gr_id %in% coding_gr_id) &
    (is.na(grMutRate) | grMutRate < 100 * args$max_gr_mut_rate_q)
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER, 
                                       '; high syn gr mut.rate')
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / nrow(tiered_pvals[gr_id %in% coding_gr_id]), 
                2), '%) genomic regions which have too high synonymous ',
          'mutation rate.')
  rm(pass, grMutRate)
  
  # noncoding regions
  grMutRate <- unlist(tiered_pvals[, 'meanMutRateQuant', with = F])
  pass <- (!tiered_pvals$gr_id %in% coding_gr_id)  &
    (!is.na(grMutRate) & grMutRate < 100 * args$max_gr_mut_rate_q)
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER, 
                                       '; high genomic region mut.rate')
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / 
                  nrow(tiered_pvals[!gr_id %in% coding_gr_id]), 2), 
          '%) genomic regions which have too high gemonic ',
          'mutation rate.')
  rm(pass, grMutRate)
  
} else {
  
}

colToUse <- 'meanMutRateQuant'
filterName <-  '; high genomic region mut.rate'
if (args$use_synonumous_mut_rate) {
  colToUse <- 'synMeanMutRateQuant'
  filterName <-  '; high syn gr mut.rate'
}

if (colToUse %in% colnames(tiered_pvals)) {
  if (args$max_gr_mut_rate_q < 1) {
    
    grMutRate <- unlist(tiered_pvals[, colToUse, with = F])
    pass <- grMutRate < !is.na(grMutRate) & 100 * args$max_gr_mut_rate_q 
    tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER, 
                                         filterName)
    
    message('[', Sys.time(), '] Found ', sum(!pass), '(',
            round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
            'regions which have too', gsub(';', '', filterName))
    rm(pass)
  }
} else {
  if (!colToUse %in% colnames(tiered_pvals)) {
    msg <- 'Filtering based on mean genome region '
    if (args$use_synonumous_mut_rate) {
      msg <- paste(msg, 'SYNONYMOUS')
    }
    msg <- paste(msg, 'mutation rate quantile was not performed, because ',
                 'column ', colToUse, ' was not found.')
    message('[', Sys.time(), '] ', msg)
  }
}

# Mark olfactory genes --------------------------------------------------------
if (args$remove_olfactory & 'is_olfactory' %in% colnames(tiered_pvals)) {
  pass <- !tiered_pvals$is_olfactory
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER,
                                       '; olfactory')
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
          'regions which are olfactory.')
  rm(pass)
}

# Mark by quantile of length --------------------------------------------------
if (args$max_gr_len_q < 1) {
  pass <- !is.na(tiered_pvals$gr_lenQuant) &
          tiered_pvals$gr_lenQuant < 100 * args$max_gr_len_q
  tiered_pvals[!pass]$FILTER <- paste0(tiered_pvals[!pass]$FILTER, 
                                       '; long gene')
  
  message('[', Sys.time(), '] Found ', sum(!pass), '(',
          round(100 * sum(!pass) / nrow(tiered_pvals), 2), '%) genomic ',
          'regions which too long.')
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
            'regions which are not expressed in ', eCol, '.')
    rm(pass)
  }
} 

# Remove PASS at the start of FILTER for entries which didn't pass ------------
tiered_pvals[, FILTER := gsub('^PASS; ', '', FILTER)]

# Print summary table of filtering --------------------------------------------

# [SAVE] Tiered results as a table --------------------------------------------
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