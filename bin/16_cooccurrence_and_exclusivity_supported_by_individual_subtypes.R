#!/usr/bin/env Rscript
# FILE: cooccurrence_and_exclusivity_supported_by_individual_subtypes.R -------
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
# CREATED:  17.12.2021
# REVISION: 28.01.2024

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

# Functions -------------------------------------------------------------------
#' check_significance_in_uniform_subtype
#' @description Checks if a gene pair was found significant (based on the raw
#' p-value) in the runs of DISCOVER made on histologically uniform tumor 
#' subtypes.
#' @author Maria Litovchenko
#' @param gene_pair named vector with members gr_id_1, gene_id_1, gr_id_2,
#'                  gene_id_2, comparison, mode
#' @param discoverUniformSubtypes data table with DISCOVER runs performed on
#' histologically uniform tumor subtypes
#' @param p_val_max numeric, maximum p-value which can be achieved by the gene
#'        pair in the individual subtypes
#' @return logical, TRUE if there was a tumor subtype which have p-value for 
#'         that gene pair < p_val_max
check_significance_in_uniform_subtype <- function(gene_pair,
                                                  discoverUniformSubtypes, 
                                                  p_val_max = 0.1) {
  gene_pair <- as.data.table(t(gene_pair))
  gene_pair <- gene_pair[, c('gr_id_1', 'gene_id_1', 'gr_id_2', 'gene_id_2', 
                             'comparison', 'mode'), with = F]
  gene_pair_r <- copy(gene_pair)
  setnames(gene_pair_r,
           c('gr_id_1', 'gene_id_1', 'gr_id_2', 'gene_id_2'),
           c('gr_id_2', 'gene_id_2', 'gr_id_1', 'gene_id_1'))
  gene_pair <- rbind(gene_pair, gene_pair_r)
  rm(gene_pair_r)
  # retrieve comparisons made in the individual tumor subtypes
  indivDiscover <- merge(discoverUniformSubtypes, gene_pair, 
                         by = c('gr_id_1', 'gene_id_1', 'gr_id_2', 'gene_id_2', 
                                'comparison', 'mode'))
  result <- F
  if (nrow(indivDiscover) != 0) {
    if (min(indivDiscover$p.value_subtypeBiasRemoved, 
            na.rm = T) <= p_val_max) {
      result <- T
    }
  }
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'cooccurrence_and_exclusivity_supported_by_individual_subtypes.R')

subtypeHelp <- paste('A histologically composite subtype, i.e. Panlung.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    help = subtypeHelp)

discoverHelp <- 'Path to DISCOVER result for the composite subtype'
parser$add_argument("-d", "--discover_composite_subtype", required = T, 
                    type = 'character', help = discoverHelp)

unifSubtypeHelp <- paste('List of histologically uniform tumor subtypes from',
                         'which the composite subtype contains')
parser$add_argument("-cu", "--uniformal_tumor_subtypes", required = T,
                    nargs = '+', type = 'character', help = unifSubtypeHelp)

discoverUnifHelp <- paste('Paths to DISCOVER result files for the',
                          'histologically uniform subtypes')
parser$add_argument("-du", "--discover_uniformal_subtypes", nargs = '+',
                    required = T, type = 'character', help = discoverUnifHelp)

pvalHelp <- paste('Cut off on raw p-value for significance of mutual ',
                  'incompatibility or co-occurence in histologically ',
                  'uniform subtypes')
parser$add_argument("-pval", "--p_max_unif_subtype", required = T, 
                    type = 'character', help = pvalHelp)

outputHelp <- paste('Path to the output file')
parser$add_argument("-o", "--output", required = T, 
                    type = 'character', help = outputHelp)

args <- parser$parse_args()
args$p_adj_max_comp_subtype <- as.numeric(args$p_adj_max_comp_subtype)
args$pval <- as.numeric(args$pval)
if (length(args$uniformal_tumor_subtypes) != 
    length(args$discover_uniformal_subtypes)) {
  stop('[', Sys.time(), '] Please submit the same number of uniformal ',
       'tumor subtypes and DISCOVER result files for uniformal tumor ',
       'subtypes.')
}

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test arguments --------------------------------------------------------------
# args <- list(cancer_subtype = 'Panlung', 
#              discover_composite_subtype = '', 
#              uniformal_tumor_subtypes = '',
#              discover_uniformal_subtypes = '',
#              p_max_unif_subtype = "",
#              output = '')

# Read DISCOVER results for composite subtype ---------------------------------
discover <- fread(args$discover_composite_subtype, header = T, 
                  stringsAsFactors = F)
message('[', Sys.time(), '] Read DISCOVER results for composite subtype ',
        args$cancer_subtype, ": ", args$discover_composite_subtype)

# Read DISCOVER results for histologically uniform subtypes -------------------
discoverUnifSubtypes <- lapply(args$discover_uniformal_subtypes, 
                               fread, header = T, stringsAsFactors = F)
message('[', Sys.time(), '] Read DISCOVER results for uniformal subtypes ',
        paste0(args$uniformal_tumor_subtypes, collapse = ', '),
        ": ", paste0(args$discover_uniformal_subtypes, collapse = ', '))

discoverUnifSubtypes <- do.call(rbind, discoverUnifSubtypes)
discoverUnifSubtypes <- discoverUnifSubtypes[p.value_subtypeBiasRemoved <= args$p_max_unif_subtype]
discoverUnifSubtypes <- discoverUnifSubtypes[,.(gene_id_1, gr_id_1, 
                                                gene_id_2, gr_id_2)]
discoverUnifSubtypes <- unique(discoverUnifSubtypes)

# Check that linked genomic elements are supported by signal from histologically uniform subtypes ----
discover[, supported_by_individual_subtypes := F]
if (nrow(discoverUnifSubtypes) > 0) {
  discoverUnifSubtypes[, pair_id := apply(discoverUnifSubtypes[,.(gene_id_1,
                                                                  gr_id_1, 
                                                                  gene_id_2, 
                                                                  gr_id_2)], 1,
                                          paste, collapse = "--")]
  discoverUnifSubtypes[, pair_id_rev := apply(discoverUnifSubtypes[,.(gene_id_1,
                                                                      gr_id_1, 
                                                                      gene_id_2,
                                                                      gr_id_2)],
                                              1, paste, collapse = "--")]
  
  discover[, pair_id := apply(discover[,.(gene_id_1, gr_id_1, 
                                          gene_id_2, gr_id_2)], 1,
                              paste, collapse = "--")]
  supported <- discover$pair_id %in% discoverUnifSubtypes$pair_id |
    discover$pair_id %in% discoverUnifSubtypes$pair_id_rev
  discover[supported]$supported_by_individual_subtypes <- T
}

message('[', Sys.time(), '] Updated DISCOVER results for ', 
        args$cancer_subtype, ' taking into account support or lack of ',
        'thereof from uniformal cancer subtypes')

# Write table with results of co-occurrence/exclusivity to file ---------------
write.table(discover, args$output, append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote results of DISCOVER (co-occurrence & ',
        'exclusivity analysis) run to ', args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')