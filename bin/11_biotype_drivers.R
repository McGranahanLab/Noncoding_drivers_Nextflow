#!/usr/bin/env Rscript
# FILE: biotype_drivers.R -----------------------------------------------------
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
# CREATED:  17.12.2020
# REVISION: 18.12.2023

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
#' assignNcopiesRelToPloidyGroups
#' @description 
#' @author Maria Litovchenko
#' @param varsDT data table with computed number of mutant copies relative to 
#' ploidy. Must have columns: participant_id, gr_name and a column with 
#' computed number of mutant copies relative to ploidy. Name of that column
#' should be given to nCopiesCol.
#' @param nCopiesCol name of the column with computed number of mutant copies
#' relative to ploidy. 
#' @return data table with columns gr_name, val_group, perc_group, n_total,
#' where val_group is the assigned group, perc_group is the % of mutated 
#' patients in that genomic range having an alteration falling within the 
#' group, n_total is the total number of mutated patients in that genomic range
assignNcopiesRelToPloidyGroups <- function(varsDT, nCopiesCol) {
  # assign a group to each value of number of mutant copies relative to ploidy
  valGroupLvl <- c('[0;0.125)' = 0, '[0.125;0.375)' = 0.125,
                   '[0.375;0.625)' = 0.375, '[0.625;0.875)' = 0.625,
                   '[0.875;1)' = 0.875, '[1;1.25)' = 1, '[1.25;1.75)' = 1.25,
                   '[1.75;2.25)' = 1.75, '[2.25;3]' = 2.25, '[3;4)' = 3, 
                   '[4;5)' = 4, '[5;6)' = 5, '6+' = 6)
  result <- copy(varsDT)
  result[, val_group := names(valGroupLvl)[1]]
  nCopies <- unlist(varsDT[, nCopiesCol, with = F])
  for (i in 1:length(valGroupLvl)) {
    result[nCopies >= valGroupLvl[i]]$val_group <- names(valGroupLvl)[i]
  } 
  
  # calculate % patients per genomic region falling within each group
  result[, perc_group := length(unique(participant_id)), 
         by = .(gr_name, val_group)]
  result[, n_total := length(unique(participant_id)), by = gr_name]
  result[, perc_group := perc_group/n_total]
  result <- unique(result[,.(gr_name, val_group, perc_group, n_total)])
  
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'biotype_drivers.R')

subtypeHelp <- paste('A cancer subtype to select from patientsInv table. Only',
                     'mutations from patients with that cancer type will be',
                     'selected. In case an analysis of several cancer types',
                     'needed to be performed please run this script ',
                     'separetedly for each cancer type.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    help = subtypeHelp)

cnHelp <- paste('Path to a file with copy number information for a scanned',
                'genomic regions. Columns participant_id, gr_name, nMinor,',
                'nMajor and cn_type are needed.')
parser$add_argument("-cn", "--cn", required = T, type = 'character',
                    help = cnHelp)

mutHelp <- paste('Path to a file with a mutations (SNV and small indels)',
                 'information for a scanned genomic regions. Columns',
                 'participant_id, gr_name, key and mut.multi are needed.')
parser$add_argument("-m", "--mutmult", required = T, type = 'character',
                    help = mutHelp)

driversHelp <- paste('Path to file containing de-novo detected cancer drivers')
parser$add_argument("-d", "--drivers", required = F,
                    type = 'character', help = driversHelp)

minNmutHelp <- paste('Minimal number of patients in which a mutation (SNV or',
                     'small indel) in a genomic element should be found so',
                     'that mutations would be considered for biotyping.')
parser$add_argument("-minMut", "--min_n_patient_mut", required = F,
                    type = 'integer', default = 0, help = minNmutHelp)

minNcnaHelp <- paste('Minimal number of patients in which a amplification/',
                     'gain/loss in a genomic element should be found so',
                     'that copy number would be considered for biotyping.')
parser$add_argument("-minCNA", "--min_n_patient_cna", required = F,
                    type = 'integer', default = 0, help = minNcnaHelp)

weakTSGhelp <- paste('A munimal cut off on percentage of biallelically',
                     'inactivated patients from the total number of patients',
                     'mutated (SNV and small indels) in a genomic region for',
                     'the region to be assigned as weak tumor suppressor gene',
                     '(TSG). In reverse, for a weak oncogene (OG) the ',
                     'percentage should be below this number. Should range',
                     'between 0 and 1.')
parser$add_argument("-wTSG", "--weak_tsg", required = T, type = 'double', 
                    default = 0.33, help = weakTSGhelp)

tsgHelp <- paste('A munimal cut off on percentage of biallelically',
                 'inactivated patients from the total number of patients',
                 'mutated (SNV and small indels) in a genomic region for',
                 'the region to be assigned as tumor suppressor gene (TSG).',
                 'In reverse, for an oncogene (OG) the percentage should be',
                 'below this number. Should range between 0 and 1.')
parser$add_argument("-TSG", "--tsg", required = T, type = 'double',
                    default = 0.50, help = tsgHelp)

weakOGhelp <- paste('A munimal cut off on percentage of patients with an',
                    'amplification from the total number of patients',
                    'with a copy number variatons in a genomic region for',
                    'the region to be assigned as a weak oncogene (OG).',
                    'In reverse, for a weak tumor supressor gene (TSG) the',
                    'percentage should be below this number. Should range',
                    'between 0 and 1.')
parser$add_argument("-wOG", "--weak_og", required = T, type = 'double',
                    default = 0.33, help = weakOGhelp)

ogHelp <- paste('A munimal cut off on percentage of patients with an',
                'amplification from the total number of patients',
                'with a copy number variatons in a genomic region for',
                'the region to be assigned as an oncogene (OG). In reverse,',
                'for a tumor supressor gene (TSG) the percentage should be',
                'below this number. Should range between 0 and 1.')
parser$add_argument("-OG", "--og", required = T, type = 'double',
                    default = 0.50, help = ogHelp)

outputHelp <- paste('Path to the output file')
parser$add_argument("-o", "--output", required = T, 
                    type = 'character', help = outputHelp)

args <- parser$parse_args()
# check_input_arguments_postproc(args, outputType = 'file')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test arguments --------------------------------------------------------------
# args <- list(cancer_subtype = 'Panlung',
#              cn = 'completed_runs/2023-12-14/results/tables/drivers_cn_multiplicity_annotated/driversAnnotatedWitCN-Panlung--hg19.csv', 
#              mutmult = 'completed_runs/2023-12-14/results/tables/drivers_cn_multiplicity_annotated/driversAnnotatedWitMultipl-Panlung--hg19.csv',
#              drivers = 'completed_runs/06_12_2023/results/tables/drivers/drivers-Panlung--hg19.csv',
#              min_n_patient_mut = 5, min_n_patient_cna = 10, 
#              weak_tsg = 0.33, tsg = 0.5, weak_og = 0.33, og = 0.5,
#              output = '')

# Read driver genes -----------------------------------------------------------
if (!is.null(args$drivers)) {
  message('[', Sys.time(), '] Started reading ', args$drivers)
  drivers <- fread(args$drivers, header = T, stringsAsFactors = F, 
                   select = c('gr_id', 'gene_id', 'gene_name', 'FILTER', 
                              'tier', 'is_known_cancer', 
                              'known_in_tumor_subtype', 
                              'known_cancer_biotype'))
  drivers <- unique(drivers)
  drivers <- drivers[FILTER == 'PASS' & !is.na(tier)]
  message('[', Sys.time(), '] Finished reading ', args$drivers)
  if (nrow(drivers) == 0) {
    stop('[', Sys.time(), '] no significant (FILTER is PASS and tier is not ',
         'NA) driver genes is found in ', args$drivers, ' table.')
  }
  drivers[, gr_name_no_version := apply(drivers[,.(gr_id, gene_id, gene_name)],
                                        1, paste0, collapse = '--')]
}

# Read copy number resolved scanned genomic regions ---------------------------
cnDT <- fread(args$cn, header = T, stringsAsFactors = F)
cnDT[, participant_id := as.character(participant_id)]
message('[', Sys.time(), '] Read ', args$cn)

# restrict to driver genes, if requested
if (!is.null(args$drivers)) {
  cnDT[, genomeVersion := gsub('--.*', '', gr_name)]
  cnDT[, genomeVersion := paste0(genomeVersion, '--')]
  cnDT[, gr_name_no_version := apply(cnDT, 1,
                                     function(x) gsub(x['genomeVersion'], '',
                                                      x['gr_name']))]
  cnDT <- cnDT[gr_name_no_version %in% unique(drivers$gr_name_no_version)]
  cnDT[, genomeVersion := NULL]
  cnDT[, gr_name_no_version := NULL]
  message('[', Sys.time(), '] Copy number regions were restricted to ones ',
          'which belong to identified driver genomic regions.')
}

# Read multiplicity of mutations within scanned genomic regions ---------------
mutmultDT <- fread(args$mutmult, header = T, stringsAsFactors = F)
mutmultDT[, participant_id := as.character(participant_id)]
message('[', Sys.time(), '] Read ', args$mutmult)

# restrict to driver genes, if requested
if (!is.null(args$drivers)) {
  mutmultDT[, genomeVersion := gsub('--.*', '', gr_name)]
  mutmultDT[, genomeVersion := paste0(genomeVersion, '--')]
  mutmultDT[, gr_name_no_version := apply(mutmultDT, 1,
                                          function(x) gsub(x['genomeVersion'], 
                                                           '', x['gr_name']))]
  mutmultDT <- mutmultDT[gr_name_no_version %in% 
                           unique(drivers$gr_name_no_version)]
  mutmultDT[, genomeVersion := NULL]
  mutmultDT[, gr_name_no_version := NULL]
  message('[', Sys.time(), '] Mutations were restricted to mutations ',
          'which belong to identified driver genomic regions.')
}

# Compute mutant fraction for SNVs --------------------------------------------
mutmultDT <- merge(mutmultDT, 
                   unique(cnDT[,.(participant_id, gr_name, nMinor, nMajor)]),
                   all.x = T, by = c('participant_id', 'gr_name'))
mutmultDT[, mut_fract := mut.multi/(nMinor + nMajor)]
mutmultDT[mut_fract > 1]$mut_fract <- 1
message('[', Sys.time(), '] Computed mutant allele fraction for genomic ',
        'regions containing mutation.')

# inform about mutations for which mut_fract is NA
nbefore <- nrow(mutmultDT)
nafter <- nrow(mutmultDT[!is.na(mut_fract)])
if (nbefore > nafter) {
  message('[', Sys.time(), '] Removed ', nbefore - nafter, '(', 
          100*round((nbefore - nafter)/nbefore, 4), 
          '%) mutations with NA mutant fraction.')
  mutmultDT <- mutmultDT[!is.na(mut_fract)]
}

# Assign a group of mutant fraction for SNVs per gene -------------------------
mutFractDT <- mutmultDT[,.(participant_id, gr_name, key, mut_fract)]
# in case gene has multiple mutations - average
mutFractDT[, mut_fract := mean(mut_fract),by = .(participant_id, gr_name)]
mutFractDT <- unique(mutFractDT[,.(participant_id, gr_name, mut_fract)])

# assign groups mutant fraction groups
mutFractDT <- assignNcopiesRelToPloidyGroups(mutFractDT, 'mut_fract')
mutFractDT[val_group == '[0.875;1)']$val_group <- '[0.875;1]'
mutFractDT[val_group == '[1;1.25)']$val_group <- '[0.875;1]'
mutFractDT[, perc_group := sum(perc_group),
           by = .(gr_name, val_group, n_total)]
mutFractDT <- unique(mutFractDT)

# remove genes which are mutated in < args$min_n_patient_mut
nbefore <- length(unique(mutFractDT$gr_name))
mutFractDT <- mutFractDT[n_total >= args$min_n_patient_mut]
nafter <- length(unique(mutFractDT$gr_name))
if (nbefore > nafter) {
  message('[', Sys.time(), '] Removed ', nbefore - nafter, '(', 
          100*round((nbefore - nafter)/nbefore, 4), '%) genomic elements ',
          'from SNV table because they were mutated in < ',
          args$min_n_patient_mut, ' patients.')
}

# Assign a group of N.copies rel. to ploidy for genes w/o SNV but with CNA ----
# Patient - genomic elements pairs with mutations are considered  
# separately (even if they also have CN not overlapping the mutation) with 
# the priority given to mutation. Therefore, we can exclude genomic elements
# from consideration if they overlap with mutations
cnDT <- merge(cnDT, 
              cbind(unique(mutmultDT[,.(participant_id, gr_name)]), 
                    'has_snv' = T),
              by = c('participant_id', 'gr_name'), all.x = T)
cnDT <- cnDT[is.na(has_snv)]
cnDT[, has_snv := NULL]

# select patient - genomic element pairs which have at least one CNA 
cnDT <- merge(cnDT, unique(cnDT[!is.na(cn_type)][,.(participant_id, gr_name)]),
              by = c('participant_id', 'gr_name'))

# assign groups mutant fraction groups
cnFractDT <- assignNcopiesRelToPloidyGroups(cnDT, 'n_copies_ge_maxDev')

# select genomic elements with CNA found in > min n patients
nbefore <- length(unique(cnFractDT$gr_name))
cnFractDT <- cnFractDT[n_total >= args$min_n_patient_cna]
nafter <- length(unique(cnFractDT$gr_name))
if (nbefore > nafter) {
  message('[', Sys.time(), '] Removed ', nbefore - nafter, '(', 
          100*round((nbefore - nafter)/nbefore, 4), '%) genomic elements ',
          'from CN table because they were mutated in < ',
          args$min_n_patient_cna, ' patients.')
}

# Assign biotype --------------------------------------------------------------
if (nrow(mutFractDT) == 0 | nrow(cnFractDT) == 0 | 
    length(intersect(cnDT$gr_name, mutFractDT$gr_name)) == 0) {
  emptyTable <- data.table(gr_id = character(), gene_id = character(), 
                           gene_name = character(), val_group = character(),
                           perc_gr = numeric(), n_total = integer(),
                           alt_type = character(), percInact = numeric(),
                           percAmp = numeric(), inferred_biotype = character(),
                           leading = character())
  write.table(emptyTable, file = args$output, sep = "\t", row.names = F, 
              col.names = T, append = F, quote = F)
  message('[', Sys.time(), '] Biotyping can not be performed due to absence',
          'of at least one genomic genomic region which would satisfy all ',
          'filters.')
  stop_quietly()
}

# for mutations, calculate % patients (from total number of patients mutated in
# a gene) with biallelic inactivation for every gene
mutFractDT[, percInact := sum(perc_group[val_group == '[0.875;1]']), 
           by = gr_name]
# for CNA calculate the opposite - percentage of amplifications
delGroups <- c('[0;0.125)', '[0.125;0.375)', '[0.375;0.625)', '[0.625;0.875)',
               '[0.875;1)', '[1;1.25)')
cnFractDT[, percAmp:= sum(perc_group[!val_group %in% delGroups], na.rm = T),
          by = gr_name]

# assign biotype based on % of biallelic inactivation and % of amplifications
infBiotype <- merge(unique(mutFractDT[,.(gr_name, percInact)]), 
                    unique(cnFractDT[,.(gr_name, percAmp)]), by = 'gr_name',
                    all = T)
infBiotype[, inferred_biotype := character()]
infBiotype[percInact >= args$weak_tsg & 
             percAmp < args$weak_og]$inferred_biotype <- 'weak TSG'
infBiotype[percInact >= args$tsg & percAmp < args$og]$inferred_biotype <- 'TSG'
infBiotype[percInact < args$weak_tsg & 
             percAmp >= args$weak_og]$inferred_biotype <- 'weak OG'
infBiotype[percInact < args$tsg & percAmp >= args$og]$inferred_biotype <- 'OG'

# determine the leading process: CNA or mutations (mostly needed for plotting)
infBiotype[, leading := character()]
infBiotype[percInact > percAmp | is.na(percAmp)]$leading <- 'Mut'
infBiotype[percInact < percAmp | is.na(percInact)]$leading <- 'CN'

# Output to file --------------------------------------------------------------
infBiotype_out <- merge(cbind(mutFractDT[,.(gr_name, val_group, perc_group,
                                            n_total)], alt_type = 'Mut'),
                        infBiotype, by = 'gr_name')
infBiotype_out <- rbind(infBiotype_out, 
                        merge(cbind(cnFractDT[,.(gr_name, val_group, 
                                                 perc_group, n_total)], 
                                    alt_type = 'CN'), infBiotype, 
                              by = 'gr_name'))
# parse gr_name for future convenience
grNameParsed <- do.call(rbind, lapply(lapply(unique(infBiotype_out$gr_name),
                                             strsplit, '--'), unlist))
colnames(grNameParsed) <- c('target_genome_version', 'gr_id', 'gene_id',
                            'gene_name')
grNameParsed <- as.data.table(grNameParsed)
grNameParsed[, target_genome_version := NULL]
grNameParsed[, gr_name := unique(infBiotype_out$gr_name)]
infBiotype_out <- merge(grNameParsed, infBiotype_out, by = 'gr_name')
infBiotype_out[, gr_name := NULL]

# in case driver genes were given and annotated with known cancer biotype
# status
if (!is.null(args$drivers)) {
  if (any(c('is_known_cancer', 'known_cancer_biotype') %in% 
          colnames(drivers))) {
    colsToGet <- intersect(c('gr_id', 'gene_id', 'gene_name',
                             'is_known_cancer', 'known_cancer_biotype'),
                           colnames(drivers))
    infBiotype_out <- merge(infBiotype_out,
                            unique(drivers[, colsToGet, with = F]),
                            by = c('gr_id', 'gene_id', 'gene_name'), all.x = T)
    
    message('[', Sys.time(), '] Information about known cancer biotype was ',
            'added.')
  }
}

write.table(infBiotype_out, file = args$output, append = F, quote = F,
            sep = '\t', row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote inferred cancer biotype to ', args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')