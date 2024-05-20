#!/usr/bin/env Rscript
# FILE: cooccurrence_and_exclusivity.R ----------------------------------------
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
# REVISION: 19.12.2023

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
suppressWarnings(suppressPackageStartupMessages(library(discover)))

# Functions: co-occurrence / exclusivity  -------------------------------------
#' getBinaryMatrixOfOccurrence
#' @description Produces binary matrix of occurrence of mutation in samples
#' @author Maria Litovchenko
#' @param varsDT data table which every row of which is a mutation to sample
#' match
#' @param columns vector of length 2 - column names to be used to create a 
#' binary matrix of occurrence, i.e. gene_name and participant_id
#' @return matrix with columns[1] in rows and columns[2] in columns 
#' @note each gene if mutated several times will be counted just 1 time
getBinaryMatrixOfOccurrence <- function(varsDT, columns = c('gene_name',
                                                            'participant_id')){
  if (length(columns) != 2) {
    stop('[', Sys.time(), '] Length of columns should be 2')
  }
  
  result <- varsDT[, columns, with = F]
  result <- result[!duplicated(result)]
  result[, presence := 1]
  result <- dcast(result, as.formula(paste(columns[1], '~', columns[2])),
                  value.var = "presence")
  result[is.na(result)] <- 0
  
  result <- as.data.frame(result)
  rownames(result) <- as.vector(result[, 1])
  result <- result[, -1]
  result <- as.matrix(result)
  result
}

#' coocurOrExclWithDiscover
#' @description Runs test on mutual exclusivity or co-occurrence with use of 
#' discover package.
#' @author Maria Litovchenko
#' @param driverDT data table with information about mutations in driver 
#' genomic elements in all patients in the cohort.
#' @param minNpatients minimal number of mutated patients
#' @param ... arguments to pass to getBinaryMatrixOfOccurrence function
#' @return data table with columns Gene_1, Gene_2, mode, p.value, p_adj, 
#' plotLab
coocurOrExclWithDiscover <- function(driverDT, minNpatients = 10, ...) {
  result <- suppressWarnings(data.table(Gene_1 = character(), 
                                        Gene_2 = character(),
                                        mode = character(),
                                        p.value = numeric()))
  
  if (nrow(driverDT) > 2) {
    matrOcc <- getBinaryMatrixOfOccurrence(driverDT, ...)
    
    if (nrow(matrOcc) <= 3) {
      return(result)
    }
    
    events <- discover.matrix(matrOcc)
    subset <- rowSums(matrOcc) >= minNpatients
    
    for (mode in c('exclusivity', 'co-occurrence')) {
      alt <- ifelse(mode == 'exclusivity', 'less', 'greater')
      modeTested <- pairwise.discover.test(events[subset, ],
                                           alternative = alt)
      modeTested <- suppressWarnings(melt(as.data.table(modeTested$p.values, 
                                                        keep.rownames = T)))
      colnames(modeTested) <- c('Var1', 'Var2', 'p.value')
      modeTested[Var1 == Var2]$p.value <- NA
      
      modeTested[, mode := mode]
      setnames(modeTested, c('Var1', 'Var2'), c('Gene_1', 'Gene_2'))
      result <- rbind(result, modeTested)
      rm(modeTested)
    }
  }
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'cooccurrence_and_exclusivity.R')

subtypeHelp <- paste('A cancer subtype to select from patientsInv table. Only',
                     'mutations from patients with that cancer type will be',
                     'selected. In case an analysis of several cancer types',
                     'needed to be performed please run this script ',
                     'separetedly for each cancer type.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    help = subtypeHelp)

analysisHelp <- paste('Path to inventory table containing details of the',
                      'future analysis to be conducted. Minimal columns:',
                      'tumor_subtype,', 'software,', 'gr_id,', 'gr_code,', 
                      'gr_file,', 'gr_upstr,', 'gr_downstr,', 'gr_genome,', 
                      'gr_excl_id,', 'gr_excl_code,', 'gr_excl_file,',
                      'gr_excl_upstr,', 'gr_excl_downstr,', 'gr_excl_genome,',
                      'blacklisted_codes.')
parser$add_argument("-a", "--inventory_analysis", required = T, 
                    type = 'character', help = analysisHelp)

driversHelp <- paste('Path to file containing de-novo detected cancer drivers')
parser$add_argument("-d", "--drivers", required = F,
                    type = 'character', help = driversHelp)

mutstoGrHelp <- paste('Path to files containing information about mutations',
                      'mapping to genomic regions. Usually produced by',
                      'calculate_mutation_rates.R')
parser$add_argument("-m", "--muts_to_gr", required = T, 
                    type = 'character', help = mutstoGrHelp)

synClassHelp <- paste('Variant_Classification-s from MAF format which are',
                      'acceptable as markers of synonymous variants. If given',
                      'variants of those classes will be removed from',
                      'consideration. Suggested values: Silent')
parser$add_argument("-sc", "--synAcceptedClass", required = F, nargs = '+',
                    default = NULL, type = 'character', help = synClassHelp)

foldHelp <- paste0('Indicates, whether or not splice sites should be',
                   'considered as part of CDS for DISCOVER analysis.',
                   'Default: T.')
parser$add_argument("-f", "--fold_splicesites_in_coding", required = F,
                    default = 'T', choices = c('T', 'F'), type = 'character', 
                    help = foldHelp)

minPatientsHelp <- paste('Minimal number of mutated patients required to run',
                         'DISCOVER.')
parser$add_argument("-minP", '--min_patients_discover', required = F, 
                    default = NULL, type = 'integer', help = minPatientsHelp)

outputHelp <- paste('Path to the output file')
parser$add_argument("-o", "--output", required = T, 
                    type = 'character', help = outputHelp)

args <- parser$parse_args()
args$fold_splicesites_in_coding <- as.logical(args$fold_splicesites_in_coding)
# check_input_arguments_postproc(args, outputType = 'file')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

emptyTable <- data.table(tumor_subtype = character(), 
                         comparison = character(), mode = character(), 
                         gr_id_1 = character(), gene_id_1 = character(),
                         gene_name_1 = character(), gr_id_2 = character(), 
                         gene_id_2 = character(), gene_name_2 = character(),
                         p.value = numeric(), p.adj = numeric())

# Test arguments --------------------------------------------------------------
# args <- list(cancer_subtype = 'Adenocarcinoma',
#              inventory_analysis = 'data/inventory/inventory_analysis.csv',
#              drivers = 'completed_runs/2023_12_25/results/tables/drivers/drivers-Adenocarcinoma--hg19.csv',
#              muts_to_gr = 'completed_runs/2023_12_25/results/mut_rates/mutMapToGR-Adenocarcinoma--hg19.csv',
#              synAcceptedClass = 'Silent', 
#              fold_splicesites_in_coding = T, min_patients_discover = 10, 
#              output = '')

# Read driver genes -----------------------------------------------------------
message('[', Sys.time(), '] Started reading ', args$drivers)
drivers <- fread(args$drivers, header = T, stringsAsFactors = F, 
                 select = c('gr_id', 'gene_id', 'gene_name', 'FILTER', 'tier'))
drivers <- unique(drivers)
drivers <- drivers[FILTER == 'PASS' & !is.na(tier)]
drivers <- drivers[,.(gr_id, gene_id, gene_name)]
message('[', Sys.time(), '] Finished reading ', args$drivers)
if (nrow(drivers) == 0) {
  stop('[', Sys.time(), '] no significant (FILTER is PASS and tier is not ',
       'NA) driver genes is found in ', args$drivers, ' table.')
}

# Read mutation map to genomic regions ----------------------------------------
varsToGRmap <- fread(args$muts_to_gr, header = T, stringsAsFactors = F, 
                     select = c('participant_id', 'gr_id', 'gene_id', 
                                'gene_name', 'key', 'var_class'))
message('[', Sys.time(), '] Read ', args$muts_to_gr)

# restrict to driver genes
varsToGRmap <- merge(varsToGRmap, drivers, 
                     by = c('gr_id', 'gene_id', 'gene_name'))
message('[', Sys.time(), '] Mutations were restricted to mutations ',
        'which belong to identified driver genomic regions.')

# remove silent mutations, if requested
if (!is.null(args$synAcceptedClass)) {
  nbefore <- nrow(varsToGRmap)
  varsToGRmap <- varsToGRmap[!var_class %in% args$synAcceptedClass]
  nafter <- nrow(varsToGRmap)
  message('[', Sys.time(), '] Removed ', nbefore - nafter, ' mutations out ',
          'of ', nbefore, '(', 100*round((nbefore - nafter)/nbefore, 4), '%) ',
          'because they fell into ', 
          paste(args$synAcceptedClass, collapse = ', '), ' class(es).')
}

# Read in analysis inventory --------------------------------------------------
analysisInv <- readAnalysisInventory(args$inventory_analysis)
message('[', Sys.time(), '] Read --inventory_analysis: ', 
        args$inventory_analysis)
analysisInv <- analysisInv[tumor_subtype == args$cancer_subtype]

# Fold splice sites into coding regions, if requested -------------------------
if (args$fold_splicesites_in_coding) {
  message('[', Sys.time(), '] Will fold splice sites into corresponding ',
          'coding driver genetic elements.')
  
  # assign coding genomic regions - anything which contains CDS as gr_code
  coding_gr_id <- unique(analysisInv[gr_code == 'CDS']$gr_id)
  message('[', Sys.time(), '] Following genomic regions: ', 
          paste0(coding_gr_id, collapse = ', '), ', will be considered as ',
          'coding.')
  # assign splice site genomic regions - anything which contains ss as gr_code
  ss_gr_id <- unique(analysisInv[gr_code == 'ss']$gr_id)
  message('[', Sys.time(), '] Following genomic regions: ', 
          paste0(ss_gr_id, collapse = ', '), ', will be considered as ',
          'containing splice sites.')
  
  if (length(coding_gr_id) != 0 & length(ss_gr_id) != 0) {
    driverMuts <- foldSplicSiteDriversIntoCodingDrivers(coding_gr_id, ss_gr_id,
                                                        driverMuts)
    message('[', Sys.time(), "] Driver genomic elements considered as ",
            "containing splice sites will be treated as noncoding for genes ",
            "which do not have coding part detected as driver.")
  } else {
    message('[', Sys.time(), '] Did not find either coding or splice site ',
            'regions. Can not fold splice site drivers into coding ones.')
  }
}

# Filter drivers by number of mutated tumors ----------------------------------
# we can't use nParts_total column in drivers data table to perform filtering
# by number of mutated tumors due to folding of splice sites into CDS. This 
# changes number of tumors and we have to re-calculate it.
nTumors <- varsToGRmap[,.(length(unique(participant_id))), 
                       by = .(gr_id, gene_id, gene_name)]
if (any(nTumors$V1 < args$min_patients_discover)) {
  nbefore <- nrow(nTumors)
  nafter <- nrow(nTumors[V1 >= args$min_patients_discover])
  
  if (nbefore > nafter) {
    message('[', Sys.time(), '] Excluded ', nbefore - nafter, '(',
            round(100*(nbefore - nafter)/nbefore, 2), '%), driver genomic ',
            'elements from consideration because number of mutated ',
            'tumors was < ', args$min_patients_discover, 
            '. Removed elements: ',
            paste(apply(nTumors[V1 < args$min_patients_discover], 1, 
                        function(x) paste(x['gr_id'], 'of', x['gene_name'])), 
                  collapse = ', '))
    nTumors <- nTumors[V1 >= args$min_patients_discover]
    varsToGRmap <- merge(varsToGRmap, nTumors,
                         by = c('gr_id', 'gene_id', 'gene_name'))
    varsToGRmap[, V1 := NULL]
  }
  
  if (nrow(varsToGRmap) == 0) {
    write.table(emptyTable, args$output, append = F, quote = F, sep = '\t', 
                col.names = T, row.names = F)
    message('[', Sys.time(), '] No driver genomic regions passed requested ',
            'thresholds. Analysis will not be performed.')
    stop_quietly()
  }
}

# Co-occurrence & exclusivity of coding/noncoding/all drivers -----------------
# assemble list of gr_id(s) which correspond to coding, noncoding and all 
# regions.
coding_gr_id <- unique(analysisInv[gr_code == 'CDS']$gr_id)
noncoding_gr_id <- setdiff(unique(analysisInv$gr_id), coding_gr_id)
gr_ids_list <- list(coding = coding_gr_id, noncoding = noncoding_gr_id,
                    all = c(coding_gr_id, noncoding_gr_id))
# create a unique name containing both gr_id and gene_name/id
varsToGRmap[, gr_name := apply(varsToGRmap[,.(gr_id, gene_id, gene_name)], 1,
                               paste0, collapse = '--')]
varsToGRmap[, gr_name := gsub(' ', '', gr_name)]

if (length(unique(varsToGRmap$gr_name)) == 1) {
  write.table(emptyTable, args$output, append = F, quote = F, sep = '\t', 
              col.names = T, row.names = F)
  message('[', Sys.time(), '] Only one driver genomic region passed all the ',
          'filters. Can not perform DISCOVER run.')
  stop_quietly()
}

# perform co-occurrence & exclusivity analysis
discover <- lapply(gr_ids_list, 
                   function(x) coocurOrExclWithDiscover(varsToGRmap[gr_id %in% x],
                                                        args$min_patients_discover,
                                                        c('gr_name', 
                                                          'participant_id')))
discover <- lapply(names(gr_ids_list),  
                   function(x) discover[[x]][, comparison := x])
discover <- do.call(rbind, discover)
if (nrow(discover) == 0) {
  write.table(emptyTable, args$output, append = F, quote = F, sep = '\t', 
              col.names = T, row.names = F)
  message('[', Sys.time(), '] DISCOVER run returned back empty.')
  stop_quietly()
}

# get back the gene and region names
discover[, Gene_1 := as.character(Gene_1)]
discover[, Gene_2 := as.character(Gene_2)]
discover <- cbind(do.call(rbind, 
                          lapply(discover$Gene_1,
                                 function(x) unlist(strsplit(x, '--')))),
                  do.call(rbind, 
                          lapply(discover$Gene_2,
                                 function(x) unlist(strsplit(x, '--')))),
                  discover)
colnames(discover)[1:6] <- c('gr_id_1', 'gene_id_1', 'gene_name_1',
                             'gr_id_2', 'gene_id_2', 'gene_name_2')
discover[, tumor_subtype := args$cancer_subtype]
discover[, Gene_1 := NULL]
discover[, Gene_2 := NULL]
setcolorder(discover, c('tumor_subtype', 'comparison', 'mode', 
                        'gr_id_1', 'gene_id_1', 'gene_name_1',
                        'gr_id_2', 'gene_id_2', 'gene_name_2', 'p.value'))

# Perform multiple testing correction -----------------------------------------
discover[, p.adj := p.adjust(p.value, 'BH'), by = 'mode']

# Write tables with results of co-occurrence/exclusivity to files -------------
write.table(discover, args$output, append = F, quote = F, sep = '\t', 
            col.names = T, row.names = F)
message('[', Sys.time(), '] Wrote results of DISCOVER (co-occurrence & ',
        'exclusivity analysis) run to ', args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')