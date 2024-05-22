#!/usr/bin/env Rscript
# FILE:cooccurrence_and_exclusivity_composite_subtypes_postprocessing.R -------
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
# REVISION: 08.05.2024

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
#' @description Checks if ONE gene pair was found significant (based on the raw
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
    if (min(indivDiscover$p.value, na.rm = T) <= p_val_max) {
      result <- T
    }
  }
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'cooccurrence_and_exclusivity_composite_subtypes_postprocessing.R')

subtypeHelp <- paste('A histologically composite subtype, i.e. Panlung.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    help = subtypeHelp)

discoverHelp <- 'Path to DISCOVER result for the composite subtype'
parser$add_argument("-d", "--discover_composite_subtype", required = T, 
                    type = 'character', help = discoverHelp)

discoverUnifHelp <- paste('Paths to DISCOVER result files for the',
                          'histologically uniform subtypes')
parser$add_argument("-du", "--discover_uniformal_subtypes", nargs = '+',
                    required = T, type = 'character', help = discoverUnifHelp)

pvalHelp <- paste('Cut off on raw p-value for significance of mutual ',
                  'incompatibility or co-occurence in histologically ',
                  'uniform subtypes')
parser$add_argument("-pval", "--p_max_unif_subtype", required = T, 
                    type = 'character', help = pvalHelp)

specHelp <- paste('Path to file with subtype specificity of all detected ',
                  'driver genomic elements. Must have columns: gr_id, ',
                  'gene_id, gene_name, tumor_subtype_1, tumor_subtype_2, ',
                  'specificity_mode, tumor_subtype_spec')
parser$add_argument("-spec", '--subtype_specificity', required = F, 
                    default = NULL, type = 'character', help = specHelp)

specModeHelp <- paste("Degree of drivers's subtype specificity considered to",
                      "be enough to cause mutual incompatibility caused by",
                      "tumor subtype specificity. One or several of ",
                      "'specific' or 'preferential'.")
parser$add_argument("-spec_mode", '--specificity_mode', required = F, 
                    choices = c('specific', 'preferential'), nargs = '+',
                    default = NULL, type = 'character', help = specModeHelp)

outputHelp <- paste('Path to the output file')
parser$add_argument("-o", "--output", required = T, 
                    type = 'character', help = outputHelp)

args <- parser$parse_args()
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
#              discover_composite_subtype = 'completed_runs/2023_12_25/results/discover/discoverResults-Panlung--hg19.csv', 
#              discover_uniformal_subtypes = list('completed_runs/2023_12_25/results/discover/discoverResults-Adenocarcinoma--hg19.csv',
#                                                 'completed_runs/2023_12_25/results/discover/discoverResults-Squamous_cell--hg19.csv'),
#              p_max_unif_subtype = 0.05,
#              subtype_specificity = 'completed_runs/2023_12_25/results/tables/subtype_specificity/subtypeSpecificity---hg19.csv',
#              specificity_mode = 'specific',
#              output = '')

# Read DISCOVER results of composite subtype ----------------------------------
discover <- fread(args$discover_composite_subtype, header = T, 
                  stringsAsFactors = F,
                  select = c('tumor_subtype', 'comparison', 
                             'mode', 'gr_id_1', 'gene_id_1', 
                             'gene_name_1', 'gr_id_2', 
                             'gene_id_2', 'gene_name_2', 
                             'p.value', 'p.adj'))
message('[', Sys.time(), '] Read DISCOVER results of composite subtype ',
        args$cancer_subtype, ": ", args$discover_composite_subtype)

# Read DISCOVER results for histologically uniform subtypes -------------------
discoverUnifSubtypes <- lapply(args$discover_uniformal_subtypes, 
                               fread, header = T, stringsAsFactors = F,
                               select = c('tumor_subtype', 'comparison', 
                                          'mode', 'gr_id_1', 'gene_id_1', 
                                          'gene_name_1', 'gr_id_2', 
                                          'gene_id_2', 'gene_name_2', 
                                          'p.value'))
message('[', Sys.time(), '] Read DISCOVER results for uniformal subtypes ',
        ": ", paste0(args$discover_uniformal_subtypes, collapse = ', '))

discoverUnifSubtypes <- do.call(rbind, discoverUnifSubtypes)
uniformal_tumor_subtypes <- unique(discoverUnifSubtypes$tumor_subtype)
discoverUnifSubtypes <- discoverUnifSubtypes[p.value <= args$p_max_unif_subtype]
discoverUnifSubtypes <- discoverUnifSubtypes[,.(gr_id_1, gene_id_1, gene_name_1, 
                                                gr_id_2, gene_id_2, gene_name_2)]
discoverUnifSubtypes <- unique(discoverUnifSubtypes)

# Check that signal is supported by signal from histologically uniform subtypes ----
if (nrow(discoverUnifSubtypes) > 0) {
  discoverUnifSubtypesRev <- copy(discoverUnifSubtypes)
  setnames(discoverUnifSubtypes, 
           c('gr_id_1', 'gene_id_1', 'gene_name_1', 
             'gr_id_2', 'gene_id_2', 'gene_name_2'),
           c('gr_id_2', 'gene_id_2', 'gene_name_2', 
             'gr_id_1', 'gene_id_1', 'gene_name_1'))
  discoverUnifSubtypes <- rbind(discoverUnifSubtypes, discoverUnifSubtypesRev)
  rm(discoverUnifSubtypesRev)  
  discoverUnifSubtypes[, supported_by_individual_subtypes := T]
  
  discover <- merge(discover, discoverUnifSubtypes, all.x = T,
                    by = c('gr_id_1', 'gene_id_1', 'gene_name_1', 
                           'gr_id_2', 'gene_id_2', 'gene_name_2'))
  discover[is.na(supported_by_individual_subtypes)]$supported_by_individual_subtypes <- F
} else {
  discover[, supported_by_individual_subtypes := logical()]
}

message('[', Sys.time(), '] Updated DISCOVER results for ', 
        args$cancer_subtype, ' by adding column ',
        'supported_by_individual_subtypes')

# Read tumor subtype specificity of drivers -----------------------------------
subtypeSpecificity <- data.table()
if (!is.null(args$subtype_specificity)) {
  subtypeSpecificity <- fread(args$subtype_specificity, header = T, 
                              stringsAsFactors = F,
                              select = c("tumor_subtype_1", "tumor_subtype_2",
                                         "gr_id", "gene_id", "gene_name",
                                         "tumor_subtype_spec",
                                         "specificity_mode"))
  subtypeSpecificity <- unique(subtypeSpecificity)
  message('[', Sys.time(), '] Read ', args$subtype_specificity)
  
  # select only histologically uniformal subtypes which composite subtype is
  # composed of
  subtypeSpecificity <- subtypeSpecificity[tumor_subtype_1 %in% uniformal_tumor_subtypes &
                                             tumor_subtype_2 %in% uniformal_tumor_subtypes]
  subtypeSpecificity <- subtypeSpecificity[specificity_mode %in%
                                             args$specificity_mode]
  subtypeSpecificity <- subtypeSpecificity[,.(gr_id, gene_id, gene_name,
                                              tumor_subtype_spec)]
  subtypeSpecificity <- unique(subtypeSpecificity)
  
  if (nrow(subtypeSpecificity) == 0) {
    message('[', Sys.time(), '] No tumor subtype specific driver genomic ',
            'regions is found.')
  } else {
    idCols <- c('gr_id', 'gene_id', 'gene_name')
    subtypeSpecificity <- split(subtypeSpecificity, drop = T, by = idCols)
    subtypeSpecificity <- lapply(subtypeSpecificity, 
                                 function(x) cbind(unique(x[, idCols, 
                                                            with = F]),
                                                   tumor_subtype_spec = list(x$tumor_subtype_spec)))
    subtypeSpecificity <- do.call(rbind, subtypeSpecificity)
  }
} else {
  message('[', Sys.time(), '] --subtype_specificity is not given. Will not ',
          'be able to account for tumor subtype specificity of driver ',
          'genomic regions in search for mutual co-occurrence/exclusivity.')
}

# Mark driver pairs which can be incompatible due to their subtype specs ------
if (nrow(subtypeSpecificity) != 0) {
  setnames(subtypeSpecificity, 
           c('gr_id', 'gene_id', 'gene_name', 'tumor_subtype_spec'), 
           c('gr_id_1', 'gene_id_1', 'gene_name_1', 'tumor_subtype_spec_1'))
  discover <- merge(discover, subtypeSpecificity, all.x= T,
                    by = c('gr_id_1', 'gene_id_1', 'gene_name_1'))
  setnames(subtypeSpecificity, 
           c('gr_id_1', 'gene_id_1', 'gene_name_1', 'tumor_subtype_spec_1'), 
           c('gr_id_2', 'gene_id_2', 'gene_name_2', 'tumor_subtype_spec_2'))
  discover <- merge(discover, subtypeSpecificity, all.x= T,
                    by = c('gr_id_2', 'gene_id_2', 'gene_name_2'))
  incompatibleDueToSubtype <- apply(discover[,.(tumor_subtype_spec_1, 
                                                tumor_subtype_spec_2)], 1,
                                    function(x) length(intersect(x[1][[1]],
                                                                 x[2][[1]])) == 0 &
                                      !is.na(x[1][[1]]) & !is.na(x[2][[1]]))
  incompatibleDueToSubtype <- sapply(incompatibleDueToSubtype, 
                                     function(x) if(length(x) == 0) {
                                       return(F)
                                     } else {
                                       return(x)
                                     })
  discover[, incompatible_due_subtype_specificity := incompatibleDueToSubtype]
  discover[, tumor_subtype_spec_1 := sapply(tumor_subtype_spec_1, paste,
                                           collapse = ',')]
  discover[tumor_subtype_spec_1 == '']$tumor_subtype_spec_1 <- NA
  discover[, tumor_subtype_spec_2 := sapply(tumor_subtype_spec_2, paste,
                                           collapse = ',')]
  discover[tumor_subtype_spec_2 == '']$tumor_subtype_spec_2 <- NA
  message('[', Sys.time(), '] Updated DISCOVER results for ', 
          args$cancer_subtype, ' by adding column ',
          'incompatible_due_subtype_specificity.')
} else {
  discover[, incompatible_due_subtype_specificity := logical()]
}

# Write table with results of co-occurrence/exclusivity to file ---------------
write.table(discover, args$output, append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote results of DISCOVER (co-occurrence & ',
        'exclusivity analysis) run to ', args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')