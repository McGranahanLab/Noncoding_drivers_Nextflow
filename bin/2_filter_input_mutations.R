#!/usr/bin/env Rscript
# FILE: 2_filter_input_mutations.R ---------------------------------------
#
# DESCRIPTION: Creates input mutation files for the requested de-novo cancer 
# driver genes calling software (i.e. chasm+, digdriver, dndscv, mutpanning, 
# nbr, oncodrivefml) from submitted patient and analysis inventories. Can
# filter mutations by tumor and germline MAF, VAC, depth, black&white regions.
# Also can liftover mutation to the target genome version. Will also output one
# MAF file containing all mutation in all patients for a tumor subtype.
#
#
# USAGE: Rscript --vanilla 2_filter_input_mutations.R \
#                --inventory_patients [path to your file] \
#                --inventory_analysis [path to your file] \
#                --inventory_blacklisted [path to your file] \ # optional
#                --cancer_subtype [cancer subtype of interest] \
#                --min_depth [min coverage of tumor/germline] \
#                --min_tumor_vac [min tumor VAC] \
#                --min_tumor_vaf [min tumor VAF] \
#                --max_germline_vaf [max germline VAC] \
#                --max_germline_vac [max tumor VAC] \
#                --max_n_vars [max n variants per patient] \
#                --target_genome_path [path to fasta file of target genome] \
#                --target_genome_version hg19 \
#                --chain [path to chain file for liftover] \
#                --output [path to file to write files to] \
#                --cores [number of cores]
#
# OPTIONS: Run 
#          Rscript --vanilla 2_filter_input_mutations.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: 
# BUGS: --
# NOTES: It may seem excessive to filter mutations by location in black&white
# regions because genomic regions of interest will exclude those regions 
# anyway. Yet, in some tools, i.e. dNdScv we will not be able to completely 
# remove regions overlapping black&white regions as that will disrupt internal
# gene structure. Therefore, it is best to filter mutations by location in 
# black&white regions as well.
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  03.11.2020
# REVISION: 26.07.2023

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
source(paste0(srcDir, '/custom_functions_preprocessing.R'))

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges)))
suppressWarnings(suppressPackageStartupMessages(library(maftools)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
options(scipen = 999)

# Functions: lifting over mutations -------------------------------------------
#' liftOverVariants
#' @description IF NESSECARY, lifts over coordinates from original genome given by 
#' inDT to targetGenome with using chain chainObj
#' @author Maria Litovchenko
#' @param inDT data table with genomic regions. Columns chr, start, end, key, 
#' somatic_genome, ref, var are required. Key column should contain a unique 
#' identifier of a genomic region from the table, i.e. number of keys should be 
#' the same as number of rows in inDT. origGenome column should contain the 
#' same value for all rows.
#' @param targetGenome target genome version
#' @param chainObj chain to move genome coordinates from original genome 
#' version to the target one
#' @return lifted over genomic coordinates with added columns key_orig,
#' chr_orig, start_orig, end_orig, somatic_genome_orig
liftOverVariants <- function(inDT, targetGenome, chainObj) {
  if (unique(inDT$somatic_genome) != targetGenome) {
    message('[', Sys.time(), '] A liftover from ', unique(inDT$somatic_genome),
            ' to ', targetGenome, ' for ', unique(inDT$participant_id),
            ' will be performed.')
    
    # perform lift over
    inGR <- makeGRangesFromDataFrame(inDT, keep.extra.columns = T)
    seqlevelsStyle(inGR) <- seqlevelsStyle(chainObj)
    result <- liftOver(inGR, chainObj)
    result <- unlist(result)
    rm(inGR) # free memory 
    
    # inform about not lifted variants
    notLifted <- inDT[!key %in% result$key]
    notLiftedPerc <- round(100 * length(notLifted) / nrow(inDT), 2)
    message('[', Sys.time(), '] ', length(notLifted), '(', notLiftedPerc,
            '%) variants were not lifted over for ', 
            unique(inDT$participant_id))
    
    # inform if lifted variants are not uniquely mapped to the new genome
    if (max(table(result$key)) > 1) {
      nonUniqMap <- table(result$key)
      nonUniqMap <- nonUniqMap[nonUniqMap > 1]
      nonUniqMap <- data.table(ID = names(nonUniqMap), nVars = nonUniqMap)
      nonUniqMap <- apply(nonUniqMap, 1, paste, collapse = ': ')
      nonUniqMap <- paste(nonUniqMap, collapse = ', ')
      message('[', Sys.time(), '] lift over for ', unique(inDT$participant_id),
              ' had variants mapped non-uniquely: ', nonUniqMap)
    }
    
    # update key column
    keyUPD <- data.table(chr = as.character(seqnames(result)),
                         pos = start(result), ref = result$ref,
                         var = result$var)
    result$keyUpd <-  gsub(' ', '', apply(keyUPD, 1, paste, collapse = ':'))
    colnames(mcols(result)) <- gsub('somatic_genome', 'somatic_genome_orig',
                                    colnames(mcols(result)))
    result$somatic_genome <- targetGenome
    result <- as.data.table(result)
    result[, strand := NULL]
    result[, width := NULL]
    setnames(result, 'seqnames', 'chr')
  } else {
    result <- copy(inDT)
    
    # update key column
    keyUPD <- result[,.(chr, start, ref, var)]
    result[, keyUpd :=  gsub(' ', '', apply(keyUPD, 1, paste, collapse = ':'))]
    colnames(result) <- gsub('somatic_genome', 'somatic_genome_orig',
                             colnames(result))
    result$somatic_genome <- targetGenome
  }
  
  # add columns with original chr, start and end columns
  origPos <- inDT[,.(key, chr, start, end)]
  setnames(origPos, c('chr', 'start', 'end'), 
           c('chr_orig', 'start_orig', 'end_orig'))
  result <- merge(result, origPos, by = 'key', all.x = T)
  setnames(result, c('key', 'keyUpd'), c('key_orig', 'key'))
  setcolorder(result, colnames(inDT))
  
  result
}

# Functions: filtering mutations by black& white regions ----------------------
#' filterVarsByBlackWhiteList
#' @description 
#' @author Maria Litovchenko
#' @param varsDT data table with variants. Columns chr, start, end are 
#' essential
#' @param bwFile path to file with black- or white- listed genomic regions.
#' @param bwName string, a code name for black- or white- listed genomic 
#' regions file
#' @param chrStyle character, one of NCBI or UCSC which determine chromosome
#' naming style (1 or chr1 respectively). Final result will have this 
#' chromosome naming style.
#' @param cores number of cores
#' @param roiGR GRanges, regions of interest. If given (not NULL), only regions
#' of interest will be read from the file.
#' @return a filtered data table with variants
filterVarsByBlackWhiteList <- function(varsDT, bwFile, chrStyle, bwName, 
                                       bwFileType, cores = 1) {
  if (chrStyle != 'NCBI' & chrStyle != 'UCSC') {
    stop('[', Sys.time(), '] chrStyle should be one of NCBI or UCSC. Current ',
         'value: ', chrStyle, '.')
  }
  
  message('[', Sys.time(), '] Started filtering variants by ', bwName)
  
  # read black- or white- listed files
  annoGR <- readBWregions(bwFile, chrStyle, cores)
  # convert variants to Granges
  result <- makeGRangesFromDataFrame(varsDT, keep.extra.columns = T)
  seqlevelsStyle(result) <- chrStyle
  before <- nrow(varsDT)
  # find overlaps between black- or white- listed regions and variants
  annoOvrl <- as.data.table(findOverlaps(result, annoGR, ignore.strand = F))
  colnames(annoOvrl) <- c('varIdx', 'annoIdx')
  # filtering
  if (bwFileType == 'white') {
    result <- varsDT[unique(annoOvrl$varIdx)]
  } else {
    result <- varsDT[setdiff(1:nrow(varsDT), unique(annoOvrl$varIdx))]
  }
  after <- nrow(result)
  if (after != before) {
    message('[', Sys.time(), '] Removed ', before - after, '(',
            round(100 * (before - after) / before, 2), '%) mutations because',
            ' of filtering by ', bwName)
  }
  message('[', Sys.time(), '] Finished filtering variants by ', bwName)
  
  result
}

# Functions: exporting --------------------------------------------------------
#' mutTabToMAF
#' @description Converts variant data table to MAF format
#' @author Maria Litovchenko
#' @param varDT variant data table, essential columns: key, chr, start, end,
#' ref, var, Func.refGene, Gene.refGene, GeneDetail.refGene,
#' ExonicFunc.refGene, AAChange.refGene, participant_id, somatic_genome
#' @return data table in MAF format
mutTabToMAF <- function(varDT) {
  result <- copy(varDT)
  setnames(result, c('chr', 'start', 'end', 'ref', 'var'),
           c('Chr', 'Start', 'End', 'Ref', 'Alt'))
  
  tmpFile <- do.call(paste0, replicate(20, sample(LETTERS, 1, T), F))
  write.table(result, tmpFile, append = F, quote = F, sep = '\t',
              col.names = T, row.names = F)
  result <- annovarToMaf(tmpFile, Center = 'NONAME', tsbCol = 'participant_id',
                         refBuild = unique(varDT$somatic_genome),
                         table = 'refGene')
  result <- as.data.table(result)
  if (nrow(result[is.na(Variant_Classification)]) != 0) {
    message('[', Sys.time(), '] Found ', 
            nrow(result[is.na(Variant_Classification)]), ' mutations with NA ',
            'in Variant_Classification. Will replace it with values from ',
            'ExonicFunc.refGene.')
    result[is.na(Variant_Classification)]$Variant_Classification <-
      result[is.na(Variant_Classification)]$ExonicFunc.refGene
  }
  file.remove(tmpFile)
  
  result[, sample_id := NULL]
  result$Strand <- '+'
  result[, Tumor_Seq_Allele1 := Tumor_Seq_Allele2]
  
  # Variant_Type reported by annovarToMaf is not correct in cases of small 
  # indels. Correct it.
  result[, Variant_Type := struct_type]
  
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'create_input_mutation_files.R')

patientsHelp <- paste('Path to patientsInv table listing information about',
                      'patients, their cancer types and mutation files.',
                      'Minimal columns: participant_id, tumor_subtype,',
                      'participant_tumor_subtype, somatic_path,',
                      'somatic_genome, cohort_name')
parser$add_argument("-p", "--inventory_patients", required = T, 
                    type = 'character', help = patientsHelp)

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

blackListHelp <- paste('Path to patientsInv table containing details of the',
                       'black&white lists to use. Minimal columns:',
                       'list_name,', 'file_path,', 'file_genome,', 
                       'file_type,', 'score_column,', 'min_value,',
                       'max_value.')
parser$add_argument("-b", "--inventory_blacklisted", required = F, 
                    type = 'character', default = NULL, help = blackListHelp)

parser$add_argument('-d', '--min_depth', required = F, type = 'integer',
                    default = 30, 
                    help = paste("Minimal mutations' depth of coverage ,",
                                 'both tumor and germline. Default: 30'))

parser$add_argument('-t_vac', '--min_tumor_vac', required = F, 
                    type = 'integer', default = 10, 
                    help = paste('Minimal variant allele count in tumor. ',
                                 'Default: 10'))
parser$add_argument('-t_vaf', '--min_tumor_vaf', required = F, 
                    type = 'double', default = 5.0, 
                    help = paste('Minimal variant allele frequency in tumor. ',
                                 'Default: 5, meaning 5 percent'))
parser$add_argument('-g_vaf', '--max_germline_vaf', required = F, 
                    type = 'double', default = 1.0, 
                    help = paste('Maximum variant allele frequency in ',
                                 'germline. Default: 1, meaning 1 percent'))
parser$add_argument('-g_vac', '--max_germline_vac', required = F, 
                    type = 'integer', default = 5, 
                    help = paste('Maximum variant allele count in germline. ',
                                 'Default: 5'))

maxNvarsHelp <- paste("Maximum number of variants per patient. Patients with", 
                      "variants exceeding this number will be removed.")
parser$add_argument("-v", "--max_n_vars", required = F, default = 90000,
                    type = 'integer', help = maxNvarsHelp)

targetGenomePathHelp <- paste('Path to the fasta file, genome version of',
                              'which should be same as ',
                              '--target_genome_version.')
parser$add_argument("-f", "--target_genome_path", required = T, 
                    type = 'character', help = targetGenomePathHelp)

targetGenomeHelp <- paste("Genome version, i.e. hg19, to which input variants",
                          "files for software should be brought.",
                          "Default: hg19.")
parser$add_argument("-g", "--target_genome_version", required = F,
                    default = 'hg19', type = 'character',
                    help = targetGenomeHelp)

targetGenomeChrHelp <- paste('Path to the tab-separated file containing ', 
                             'chromosomal lengths of the target genome. ',
                             'Must have 3 columns: chr, start(1) and length ',
                             'of the chromosome. No header.')
parser$add_argument("-cl", "--target_genome_chr_len", required = F,
                    default = NULL, type = 'character', 
                    help = targetGenomeChrHelp)

chainHelp <- paste('Path to chain file in case genome version of mutations is',
                   'not the same as --target_genome_version')
parser$add_argument("-l", "--chain", required = F, default = NULL,
                    type = 'character', help = chainHelp)

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output file")

coresHelp <- 'How many cores the script should use. Default: 1.'
parser$add_argument("-n", "--cores", required = F, type = 'integer', 
                    default = 1, help = coresHelp)

args <- parser$parse_args()
check_input_arguments_preproc(args, outputType = 'file')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test input args -------------------------------------------------------------
# args <- list(inventory_patients = '../data/inventory/inventory_patients_tcga.csv',
#              inventory_analysis = '../data/inventory/inventory_analysis_tcga.csv',
#              inventory_blacklisted = '../data/inventory/inventory_blacklist_tcga.csv', 
#              cancer_subtype = 'LUSC', min_depth = 30, 
#              min_tumor_vac = 10, max_germline_vac = 5,
#              min_tumor_vaf = 5.0, max_germline_vaf = 1.0, 
#              max_n_vars = 90000,
#              target_genome_version = 'hg19',
#              target_genome_path = '../data/assets/reference_genome/Homo_sapiens_assembly19.fasta',
#              target_genome_chr_len = '../data/assets/reference_genome/Homo_sapiens_assembly19.chrom.sizes.bed',
#              chain = '../data/assets/reference_genome/hg38ToHg19.over.chain',
#              output = 'test.maf', cores = 2)

# Read inventories ------------------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients, args$cores)
# check that if somatic_genome in patientsInv is not the same as target one
# chain file is submitted
if (any(!unique(patientsInv$somatic_genome) %in% args$target_genome_version) &
    is.null(args$chain)) {
  stop('[', Sys.time(), '] genome version of some mutation tables is not ',
       'the same as --target_genome_version, but no chain file is provided')
}
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)

analysisInv <- readAnalysisInventory(args$inventory_analysis, args$cores)
message('[', Sys.time(), '] Read --inventory_analysis: ', 
        args$inventory_analysis)

if (!is.null(args$inventory_blacklisted)) {
  bwInv <- readBlacklistInventory(args$inventory_blacklisted, args$cores)
  if (is.null(bwInv)) {
    args$inventory_blacklisted <- NULL
    print(paste0('[', Sys.time(), '] Black & white lists inventory file ',
                 '(--inventory_blacklisted) was empty. Proceeding without ',
                 'it.'))
  } else {
    message('[', Sys.time(), '] Read --inventory_blacklisted: ', 
            args$inventory_blacklisted)
  }
}

# select cancer subtype
patientsInv <- patientsInv[tumor_subtype %in% args$cancer_subtype]
analysisInv <- analysisInv[tumor_subtype %in% args$cancer_subtype]

# select only black&white lists which will be used in the future
if (!is.null(args$inventory_blacklisted)) {
  uniqBWcodes <- unlist(unique(analysisInv$blacklisted_codes))
  if (!any(is.na(uniqBWcodes))) {
    if (any(uniqBWcodes != '')) {
      bwInv <- bwInv[list_name %in% uniqBWcodes]
    } else {
      bwInv <- NULL
    }
  } else {
    bwInv <- NULL
  }
}

# Read mutations files --------------------------------------------------------
message('[', Sys.time(), '] Started reading input mutation files')
n_files <- nrow(patientsInv)
allVars <- list()
for (i in 1:n_files) {
  allVars[[i]] <- readSomaticVars(patientsInv$somatic_path[i], args$cores)
  # because files in MAF format do contain Tumor_Sample_Barcode renamed to 
  # participant_id
  suppressWarnings(allVars[[i]][, participant_id := NULL])
  allVars[[i]] <- cbind(allVars[[i]], 
                        patientsInv[,.(participant_id, somatic_genome, 
                                       cohort_name)][i, ])
  
  if ((100 * (i / n_files)) %% 2) {
    message('\t[', Sys.time(), '] Read: ', round(100 * (i / n_files), 2), '%')
  }
}

allVars <- as.data.table(do.call(rbind.fill, allVars))
message('[', Sys.time(), '] Finished reading input mutation files')

# Bring mutations' chromosome notation to the one of reference genome ---------
# determine, if chromosomal names in fasta file have 'chr' in front of them or
# not. This will help us to output variants with chromosomal names in the
# proper format matching reference genome
outChrStyle <- getSeqlevelsStyle(args$target_genome_path)

allVars[, chr := gsub('chr', '', chr)] # this covers NCBI
if (outChrStyle == 'UCSC') {
  allVars[, chr := paste0('chr', chr)] 
}

suppressMessages(suppressWarnings(gc()))

# Filter variants: noncanonical chr, duplicated, N var/ref, same alt/ref ------
# 1) Non canonical chr
before <- nrow(allVars)
allVars <- allVars[chr %in% acceptedChrNames]
after <- nrow(allVars)
if (after != before) {
  message('[', Sys.time(), '] Removed ', before - after, ' out of ', before, 
          ' (', round(100 * (before - after)/ before), '%) mutations because ',
          'of non canonical chromosome')
}

# 2) duplicated mutations
before <- nrow(allVars)
allVars <- allVars[!duplicated(allVars)]
after <- nrow(allVars)
if (after != before) {
  message('[', Sys.time(), '] Removed ', before - after, ' out of ', before, 
          ' (', round(100 * (before - after)/ before), '%) mutations because ',
          'of duplication')
}

# 3) N as ref or alt allele
before <- nrow(allVars)
allVars <- allVars[ref != 'N' & var != 'N']
after <- nrow(allVars)
if (after != before) {
  message('[', Sys.time(), '] Removed ', before - after, ' out of ',  before, 
          ' (', round(100 * (before - after)/ before), '%) mutations because ',
          'of N in ref/alt')
}

# 4) equal alt and ref
before <- nrow(allVars)
allVars <- allVars[ref != var]
after <- nrow(allVars)
if (after != before) {
  message('[', Sys.time(), '] Removed ', before - after, ' out of ', before, 
          ' (',  round(100 * (before - after)/before), '%) mutations because ',
          'of ref = alt')
}

suppressMessages(suppressWarnings(gc()))

# Filter variants : VAF, ALT counts, etc --------------------------------------
# 1) minimum depth of coverage for mutations
if (!is.null(args$min_depth)) {
  if ('t_depth' %in% colnames(allVars)) {
    # inform about number of entries and participants with NA in t_depth
    nNA <- nrow(allVars[is.na(t_depth)])
    if (nNA > 0) {
      nParticipantsNA <- unique(allVars[is.na(t_depth)]$participant_id)
      nParticipantsNA <- paste(nParticipantsNA, collapse = '; ')
      message('[', Sys.time(), '] Found ', nNA, ' entries with NA in t_depth ',
              'column in participants: ', nParticipantsNA, '. ',
              'The entries will NOT be filtered out. If it is undesirable ',
              'behaviour, consider adding column(s) t_depth to the ',
              'corresponding mutational tables and set values to 0.')
    }
    
    before <- nrow(allVars)
    allVars <- allVars[t_depth >= args$min_depth | is.na(t_depth)]
    after <- nrow(allVars)
    if (after != before) {
      message('[', Sys.time(), '] Removed ', before - after, ' out of ', 
              before, ' (', round(100 * (before - after)/ before), '%) ',
              'mutations because of low (< ', args$min_depth, ') tumor depth.')
    }
  } else {
    warning('[', Sys.time(), '] an argument is given to --min_depth, but ',
            'column t_depth is not found in the mutation table. Filtering ',
            'will not be performed.')
  }
  if ('n_depth' %in% colnames(allVars)) {
    # inform about number of entries and participants with NA in n_depth
    nNA <- nrow(allVars[is.na(n_depth)])
    if (nNA > 0) {
      nParticipantsNA <- unique(allVars[is.na(n_depth)]$participant_id)
      nParticipantsNA <- paste(nParticipantsNA, collapse = '; ')
      message('[', Sys.time(), '] Found ', nNA, ' entries with NA in n_depth ',
              'column in participants: ', nParticipantsNA, '. ',
              'The entries will NOT be filtered out. If it is undesirable ',
              'behaviour, consider adding column(s) n_depth to the ',
              'corresponding mutational tables and set values to 0.')
    }
    
    before <- nrow(allVars)
    allVars <- allVars[n_depth >= args$min_depth | is.na(n_depth)]
    after <- nrow(allVars)
    if (after != before) {
      message('[', Sys.time(), '] Removed ', before - after, ' out of ', 
              before, ' (', round(100 * (before - after)/ before), '%) ',
              'mutations because of low (< ', args$min_depth, ') tumor depth.')
    }
  } else {
    warning('[', Sys.time(), '] an argument is given to --min_depth, but ',
            'column n_depth is not found in the mutation table. Filtering ',
            'will not be performed.')
  }
}

# 2) minimum VAF for tumor mutations
if (!is.null(args$min_tumor_vaf)) {
  if ('t_alt_count' %in% colnames(allVars) & 
      't_depth' %in% colnames(allVars)) {
    # inform about number of entries and participants with NA in t_alt_count or 
    # t_depth
    nNA <- nrow(allVars[is.na(t_alt_count) | is.na(t_depth)])
    if (nNA > 0) {
      nParticipantsNA <- allVars[is.na(t_alt_count) | is.na(t_depth)]
      nParticipantsNA <- unique(nParticipantsNA$participant_id)
      nParticipantsNA <- paste(nParticipantsNA, collapse = '; ')
      message('[', Sys.time(), '] Found ', nNA, ' entries with NA in t_depth ',
              't_alt_count or t_depth columns in participants: ', 
              nParticipantsNA, '. The entries will NOT be filtered out. If ',
              'it is undesirable behaviour, consider adding column(s) ',
              't_depth and t_alt_count to the corresponding mutational ',
              'tables and set values to 0.')
    }
    
    before <- nrow(allVars)
    allVars <- allVars[100 * t_alt_count / t_depth  >= args$min_tumor_vaf |
                         is.na(t_alt_count) | is.na(t_depth)]
    after <- nrow(allVars)
    if (after != before) {
      message('[', Sys.time(), '] Removed ', before - after, ' out of ',
              before, ' (', round(100 * (before - after)/ before),
              '%) mutations because of low (< ', args$min_tumor_vaf,
              ') tumor VAF')
    }
  } else {
    if ('t_maxVAF' %in% colnames(allVars)) {
      # inform about number of entries and participants with NA in t_maxVAF 
      nNA <- nrow(allVars[is.na(t_maxVAF)])
      if (nNA > 0) {
        nParticipantsNA <- unique(allVars[is.na(t_maxVAF)]$participant_id)
        nParticipantsNA <- paste(nParticipantsNA, collapse = '; ')
        message('[', Sys.time(), '] Found ', nNA, ' entries with NA in ',
                't_maxVAF column in participants: ',  nParticipantsNA, 
                '. The entries will NOT be filtered out. If it ',
                'is undesirable behaviour, consider adding column t_maxVAF ',
                'to the corresponding mutational tables and set up proper ',
                'values.')
      }
      
      before <- nrow(allVars)
      allVars <- allVars[t_maxVAF  >= args$min_tumor_vaf | is.na(t_maxVAF)]
      after <- nrow(allVars)
      if (after != before) {
        message('[', Sys.time(), '] Removed ', before - after, ' out of ', 
                before, ' (', round(100 * (before - after)/ before),  
                '%) because low (< ', args$min_tumor_vaf, ') tumor VAF')
      }
    } else {
      warning('[', Sys.time(), '] an argument is given to --min_tumor_vaf, ',
              'but columns t_alt_count and t_depth or t_maxVAF are not found ',
              'in the mutation table. Filtering will not be performed.')
    }
  }
}

# 3) minimum VAC for tumor mutations
if (!is.null(args$min_tumor_vac)) {
  if ('t_alt_count' %in% colnames(allVars)) {
    # inform about number of entries and participants with NA in t_alt_count
    nNA <- nrow(allVars[is.na(t_alt_count)])
    if (nNA > 0) {
      nParticipantsNA <- unique(allVars[is.na(t_alt_count)]$participant_id)
      nParticipantsNA <- paste(nParticipantsNA, collapse = '; ')
      message('[', Sys.time(), '] Found ', nNA, ' entries with NA in ',
              't_alt_count column in participants: ', nParticipantsNA, '. ',
              'The entries will NOT be filtered out. If it is undesirable ',
              'behaviour, consider adding column(s) t_alt_count to the ',
              'corresponding mutational tables and set values to 0.')
    }
    
    before <- nrow(allVars)
    allVars <- allVars[t_alt_count >= args$min_tumor_vac | is.na(t_alt_count)]
    after <- nrow(allVars)
    if (after != before) {
      message('[', Sys.time(), '] Removed ', before - after, ' out of ',
              before, ' (', round(100 * (before - after)/ before),
              '%) mutations because of low (< ', args$min_tumor_vac,
              ') tumor VAC')
    }
  } else {
    warning('[', Sys.time(), '] an argument is given to --min_tumor_vac, but ',
            'column t_alt_count is not found in the mutation table. ',
            'Filtering will not be performed.')
  }
}

# 4) maximum germline VAF 
if(!is.null(args$max_germline_vaf)) {
  if ('n_alt_count' %in% colnames(allVars) & 
      'n_depth' %in% colnames(allVars)) {
    # inform about number of entries and participants with NA in n_alt_count or 
    # n_depth
    nNA <- nrow(allVars[is.na(n_alt_count) | is.na(n_depth)])
    if (nNA > 0) {
      nParticipantsNA <- allVars[is.na(n_alt_count) | is.na(n_depth)]
      nParticipantsNA <- unique(nParticipantsNA$participant_id)
      nParticipantsNA <- paste(nParticipantsNA, collapse = '; ')
      message('[', Sys.time(), '] Found ', nNA, ' entries with NA in ',
              'n_alt_count or n_depth columns in participants: ', 
              nParticipantsNA, '. The entries will NOT be filtered out. If ',
              'it is undesirable behaviour, consider adding column(s) ',
              'n_alt_count and n_depth to the corresponding mutational ',
              'tables and set up proper values.')
    }
    
    before <- nrow(allVars)
    allVars <- allVars[100 * n_alt_count / n_depth  <= args$max_germline_vaf | 
                         is.na(n_alt_count) | is.na(n_depth)]
    after <- nrow(allVars)
    if (after != before) {
      message('[', Sys.time(), '] Removed ', before - after, ' out of ', 
              before, ' (', round(100 * (before - after)/ before), '%) ',
              'mutations because of high (> ', args$max_germline_vaf, 
              ') germline VAF')
    }
  } else {
    warning('[', Sys.time(), '] an argument is given to --max_germline_vaf, ',
            'but columns n_alt_count and n_depth are not found in the ',
            'mutation table. Filtering will not be performed.')
  }
}

# 5) maximum germline VAC
if (!is.null(args$max_germline_vac)) {
  if ('n_alt_count' %in% colnames(allVars)) {
    # inform about number of entries and participants with NA in n_alt_count
    nNA <- nrow(allVars[is.na(n_alt_count)])
    if (nNA > 0) {
      nParticipantsNA <- unique(allVars[is.na(n_alt_count)]$participant_id)
      nParticipantsNA <- paste(nParticipantsNA, collapse = '; ')
      message('[', Sys.time(), '] Found ', nNA, ' entries with NA in ',
              'n_alt_count column in participants: ', nParticipantsNA, '. ',
              'The entries will NOT be filtered out. If it is undesirable ',
              'behaviour, consider adding column(s) n_alt_count to the ',
              'corresponding mutational tables and set values to 0.')
    }
    
    before <- nrow(allVars)
    allVars <- allVars[n_alt_count <= args$max_germline_vac | 
                         is.na(n_alt_count)]
    after <- nrow(allVars)
    if (after != before) {
      message('[', Sys.time(), '] Removed ', before - after, ' out of ',
              before, ' (', round(100 * (before - after)/before),
              '%) mutations because of high (> ', args$max_germline_vac,
              ') germline VAC')
    }
  } else {
    warning('[', Sys.time(), '] an argument is given to --max_germline_vac, ',
            'but column n_alt_count is not found in the mutation table. ',
            'Filtering will not be performed.')
  }
}

suppressMessages(suppressWarnings(gc()))

# Liftover variants to target genome ------------------------------------------
if (any(unique(allVars$somatic_genome) != args$target_genome_version)) {
  # read in chain file
  chain <- import.chain(args$chain)
  
  allVars <- split(allVars, allVars$participant_id)
  message('[', Sys.time(), '] Started liftover')
  allVars <- mclapply(allVars, liftOverVariants, args$target_genome_version, 
                      chain, mc.cores = args$cores)
  allVars <- do.call(rbind, allVars)
  message('[', Sys.time(), '] Finished liftover')
}

suppressMessages(suppressWarnings(gc()))

# Filter variants: location outside target genome chromosomes -----------------
if (!is.null(args$target_genome_chr_len)) {
  message('[', Sys.time(), '] Lengths of chromosomes are given, will check ',
          'that mutations do not surpass them.')
  nbefore <- nrow(allVars)
  allVars <- filterVarsByBlackWhiteList(varsDT = allVars,
                                        bwFile = args$target_genome_chr_len,
                                        chrStyle = outChrStyle,
                                        bwName = 'chromosome length',
                                        bwFileType = 'white',
                                        cores = args$cores)
  nafter <- nrow(allVars)
  if (nbefore == nafter) {
    message('[', Sys.time(), '] All mutations are within chromosomal ',
            'boundaries')
  } else {
    message('[', Sys.time(), '] WARNING: Removed', before - after, '(',
            round(100 * (before - after) / before, 2), '%) mutations because',
            ' they were outside chromosomal boundaries')
  }
}

# Filter variants: location in black-/white- listed regions -------------------
if ('blacklisted_codes' %in% colnames(analysisInv)) {
  if (!is.null(bwInv)) {
    for (i in 1:nrow(bwInv)) {
      allVars <- filterVarsByBlackWhiteList(varsDT = allVars,
                                            bwFile = bwInv$file_path[i],
                                            chrStyle = outChrStyle,
                                            bwName = bwInv$list_name[i],
                                            bwFileType = bwInv$file_type[i],
                                            cores = args$cores)
      suppressMessages(suppressWarnings(gc()))
    }
  }
}

suppressMessages(suppressWarnings(gc()))

# Filter patients: number of variants does not exceed max_n_vars --------------
varsPerParticip <- allVars[,.N, by = participant_id]

if (any(varsPerParticip$N > args$max_n_vars)) {
  hypermutated <- varsPerParticip[N > args$max_n_vars]$participant_id
  hypermutated <- patientsInv[participant_id %in% hypermutated]
  
  # output table with hypermutated samples
  write.table(hypermutated, append = F, sep = '\t', row.names = F, quote = F,
              col.names = T,
              file = paste0('hypermutated-', args$cancer_subtype, '.csv'))
  
}
before <- nrow(varsPerParticip)
removedPatients <- varsPerParticip[N > args$max_n_vars]$participant_id
varsPerParticip <- varsPerParticip[N <= args$max_n_vars]
after <- nrow(varsPerParticip)
if (after != before) {
  message('[', Sys.time(), '] Removed ', before - after, '(',
          round(100 * (before - after) / before, 2), '%) patients because',
          ' number of mutations exceeded ', args$max_n_vars, '. Removed ',
          'patients: ', paste(removedPatients, collapse = ', '))
}
allVars <- allVars[participant_id %in% varsPerParticip$participant_id]

suppressMessages(suppressWarnings(gc()))

# Save annotated cancer subtype mutation table as MAF file --------------------
# bring allVars chromosome notation to the one of reference genome, yes, again
allVars[, chr := gsub('chr', '', chr)] # this covers NCBI
if (outChrStyle == 'UCSC') {
  allVars[, chr := paste0('chr', chr)] 
}

patientsInv[, participant_id := as.character(participant_id)]
allVars[, participant_id := as.character(participant_id)]
setkey(patientsInv, participant_id)
allVars[, participant_tumor_subtype := patientsInv[allVars$participant_id]$participant_tumor_subtype]
allVarsMAF <- mutTabToMAF(allVars)

# sort by chromosome & start
allVarsMAF[, Chromosome := factor(Chromosome, orderChromosomes(Chromosome))]
allVarsMAF <- allVarsMAF[order(Chromosome, Start_Position)]

write.table(allVarsMAF, args$output, append = F, sep = '\t', row.names = F,
            col.names = T, quote = F)

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')