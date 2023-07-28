#!/usr/bin/env Rscript
# FILE: create_input_mutation_files.R -----------------------------------------
#
# DESCRIPTION:
#
# USAGE: 
#
# OPTIONS: 
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

box::use(./custom_functions[...])

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(maftools))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(plyr))
options(scipen = 999)

# FUNCTIONS: reading original input files with mutations in ANNOVAR format ----
#' getVarStructTypeFromAnnovar
#' @description Gets structural type of a variant, i.e. SNP, MNP, INS, DEL.
#' @author Maria Litovchenko
#' @param dt data table with variants, columns ref and var are needed
#' @return data table with added columns struct_type and mut_len which give
#' info about variants structural type (INS, DEL, SNP, MNP) and length of
#' mutation
getVarStructTypeFromAnnovar <- function(dt) {
  result <- copy(dt)
  result[, struct_type := 'NA']
  
  result[ref %in% c('A', 'T', 'G', 'C') & 
           var %in% c('A', 'T', 'G', 'C')]$struct_type <- 'SNP'
  result[ref == '-']$struct_type <- 'INS'
  result[!ref %in% c('A', 'T', 'G', 'C', '-')]$struct_type <- 'MNP'
  result[grepl('^-.*', var)]$struct_type <- 'DEL'
  
  if(nrow(result[struct_type == 'NA']) != 0) {
    unknown <- result[struct_type == 'NA']
    stop('[', Sys.time(), '] Unknown variant type:\n',
         paste0(paste(colnames(unknown), collapse = '\t'), '\n',
                paste(unlist(unknown), collapse = '\t')))
  }
  
  result[, mut_len := 1]
  result[struct_type=='MNP']$mut_len <- apply(result[struct_type == 'MNP'], 1,
                                              function(x) max(nchar(x['ref']),
                                                              nchar(x['var'])))
  result[struct_type=='INS']$mut_len <- nchar(result[struct_type=='INS']$var)
  result[struct_type=='DEL']$mut_len <- nchar(result[struct_type=='DEL']$var)-1
  
  result
}

#' transformAllelesFromAnnovarToNorm
#' @description Translates mutation table from AnnoVar format written reference
#' and alternatives alleles to usual ones.
#' @author Maria Litovchenko
#' @param varDT data table with mutations. Have to have columns start, end, 
#' ref, var, struct_type
#' @return data table with updated columns start, end, ref and var. The row 
#' order will be the same.
transformAllelesFromAnnovarToNorm <- function(varDT) {
  result <- copy(varDT)
  result[, idx := 1:nrow(varDT)] # to keep the same order
  
  # ref allele should not have values which would follow pattern ^-.*, unless 
  # it is equal to -
  unexpectRef <- grepl('^-.*', varDT$ref) & varDT$ref != '-'
  if (sum(unexpectRef) != 0) {
    unexpectRef <- which(grepl('^-.*', varDT$ref) & varDT$ref != '-')
    message('[', Sys.time(), '] Found unexpected reference allele. Offending ',
            'line: ', paste(varDT[unexpectRef[0]], collapse = '\t'))
  }
  
  # process SNPs. Final position should be start + 1
  result[struct_type == 'SNP']$end <- result[struct_type == 'SNP']$start + 1
  # MNPs: adjust end so it corresponds to the end of mnp
  result[struct_type == 'MNP']$end <- result[struct_type == 'MNP']$start + 
    nchar(result[struct_type == 'MNP']$ref)
  # process insertions: in Annovar their position is the actual one - 1. 
  # Checked manually.
  result[struct_type == 'INS']$start <- result[struct_type == 'INS']$start + 1
  # from digdriver: final position of the interval should be the start 
  # position + 1
  result[struct_type == 'INS']$end <- result[struct_type == 'INS']$start + 1
  # deletions: transform alleles into their usual shape with - as alternative 
  # allele and stop as start + number of deleted bases + 1
  dels <- result[struct_type == 'DEL']
  if (nrow(dels) != 0) {
    message('[', Sys.time(), '] Found deletions written in annovar format, ',
            'i.e. ', dels$var[1], '. Will reformat it to usual shape with - ',
            'as alternative.')
    dels[, ref := gsub('^-', '', var)]
    dels[, var := '-']
    dels[, start := start + 1]
    # despite that in DIG it's written to add 1, in their mutation files 1 is
    # not added.
    dels[, end := start + nchar(ref)]
    result <- rbind(result[struct_type != 'DEL'], dels)
  }
  
  result <- result[order(idx)]
  if (!identical(result$idx, 1:nrow(result))) {
    stop('[', Sys.time(), '] Extra or missing entries in variant tables')
  }
  result[, idx := NULL]
  result
}

readAnnovarMutFile <- function(filePath, cores = 1) {
  colsToSelect <- c('key', 'chr', 'start', 'stop', 'ref', 'var', 
                    'Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 
                    'GeneDetail.refGene', 'AAChange.refGene', 'Use.For.Plots',
                    'Use.For.Plots.Indel')
  result <- suppressWarnings(fread(filePath, header = T, select = colsToSelect, 
                                   stringsAsFactors = F, nThread = cores))
  setnames(result, 'stop', 'end')
  
  if (any(c(23, '23', 'chr23') %in% result$chr)) {
    message('[', Sys.time(), '] Found that chr X is encoded as 23. Changed it',
            ' to X.')
    result[, chr := gsub('23', 'X', chr)]
    result[, chr := gsub('24', 'Y', chr)]
    result[, chr := gsub('25', 'M', chr)]
  }
  
  # assign structural variant type
  result <- getVarStructTypeFromAnnovar(result)
  result <- transformAllelesFromAnnovarToNorm(result)
  suppressWarnings(result[, stop := NULL])
  
  result
}

# FUNCTIONS: reading original imput files with mutations in MAF format --------
#' getVarStructType
#' @description Gets structural type of a variant, i.e. SNP, MNP, INS, DEL.
#' @author Maria Litovchenko
#' @param dt data table with variants, columns ref and var are needed
#' @return data table with added columns struct_type and mut_len which give
#' info about variants structural type (INS, DEL, SNP, MNP) and length of
#' mutation
#' @note only data tables NOT IN ANNOVAR FORMAT are processed
getVarStructType <- function(dt) {
  result <- copy(dt)
  result[, struct_type := 'NA']
  
  result[ref %in% c('A', 'T', 'G', 'C') & 
           var %in% c('A', 'T', 'G', 'C')]$struct_type <- 'SNP'
  result[ref == '-']$struct_type <- 'INS'
  result[var == '-']$struct_type <- 'DEL'
  result[(nchar(ref) > 1 & nchar(var) > 1) | 
           (nchar(ref) == 1 & ref != '-' & nchar(var) > 1) |
           (nchar(var) == 1 & var != '-' & nchar(ref) > 1)]$struct_type <- 'MNP'
  
  if(nrow(result[struct_type == 'NA']) != 0) {
    unknown <- result[struct_type == 'NA']
    stop('[', Sys.time(), '] Unknown variant type:\n',
         paste0(paste(colnames(unknown), collapse = '\t'), '\n',
                paste(unlist(unknown), collapse = '\t')))
  }
  
  result[, mut_len := 1]
  result[struct_type=='MNP']$mut_len <- apply(result[struct_type == 'MNP'], 1,
                                              function(x) max(nchar(x['ref']),
                                                              nchar(x['var'])))
  result[struct_type=='INS']$mut_len <- nchar(result[struct_type=='INS']$var)
  result[struct_type=='DEL']$mut_len <- nchar(result[struct_type=='DEL']$var)
  
  result
}

readMafMutFile <- function(filePath, cores = 1) {
  colsToSelect <- c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position',
                    'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 
                    'Gene', 'Variant_Classification', 'Amino_acids',
                    't_depth', 't_ref_count', 't_alt_count', 'n_depth', 
                    'n_ref_count', 'n_alt_count', 't_maxVAF')
  result <- suppressWarnings(fread(filePath, header = T, select = colsToSelect, 
                                   stringsAsFactors = F, nThread = cores))
  setnames(result, c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 
                     'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
                     'Gene', 'Variant_Classification', 'Amino_acids'),
           c('participant_id', 'chr', 'start', 'end', 'ref', 'var',
             'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene'))
  
  # bring ExonicFunc.refGene to annovar format
  codingExonicRg <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
                      "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
                      "Silent", "Translation_Start_Site", "Nonstop_Mutation", 
                      "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
                      "Unknown")
  result[, Func.refGene := '']
  result[ExonicFunc.refGene %in% codingExonicRg]$Func.refGene <- 'exonic'
  result[ExonicFunc.refGene == "3'UTR"]$Func.refGene <- 'UTR3'
  result[ExonicFunc.refGene == "5'UTR"]$Func.refGene <- 'UTR5'
  result[ExonicFunc.refGene == "3'Flank"]$Func.refGene <- 'downstream'
  result[ExonicFunc.refGene == "5'Flank"]$Func.refGene <- 'upstream'
  result[ExonicFunc.refGene == "IGR"]$Func.refGene <- 'intergenic'
  result[ExonicFunc.refGene == "Intron"]$Func.refGene <- 'intronic'
  result[ExonicFunc.refGene == "RNA"]$Func.refGene <- 'ncRNA_exonic'
  result[ExonicFunc.refGene == "Splice_Site"]$Func.refGene <- 'splicing'
  result[, GeneDetail.refGene := NA]
  
  if (any(c(23, '23', 'chr23') %in% result$chr)) {
    message('[', Sys.time(), '] Found that chr X is encoded as 23. Changed it',
            ' to X.')
    result[, chr := gsub('23', 'X', chr)]
    result[, chr := gsub('24', 'Y', chr)]
    result[, chr := gsub('25', 'M', chr)]
  }
  
  # add key
  result[, key := apply(result[,.(chr, start, ref, var)], 1,
                        paste, collapse = ':')]
  result[, key := gsub(' ', '', key)]
  
  # assign structural variant type
  result <- getVarStructType(result)
  
  result
}

# FUNCTIONS: reading original input files with mutations ----------------------
#' getFileType
#' @description Retrieves file type (annovar or maf) from the file's header
#' @author Maria Litovchenko
#' @param inPath path to file with somatic variants, AnnoVar format
#' @return string with file type, one of annovar or maf
getFileType <- function(inPath) {
  fileHeader <- colnames(fread(inPath, nrows = 1))
  result <- 'annovar'
  if (any(grepl('Tumor_Seq_Allele', fileHeader))) {
    result <- 'maf'
  }
  result
}

#' readSomaticVars
#' @description Reads in file with variants annotated by AnnoVar or MAF 
#' formatted file  to data table
#' @author Maria Litovchenko
#' @param filePath path to file with somatic variants, AnnoVar format
#' @param cores number of cores to use (used in fread only)
#' @return data.table with columns participant_id (maf only), chr, start, end, 
#' ref, var, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, struct_type,
#' mut_len. Optional columns: t_depth, t_ref_count, t_alt_count, n_depth,
#' n_ref_count, n_alt_count, t_maxVAF.
readSomaticVars <- function(filePath, cores = 1) {
  fileType <- getFileType(filePath)
  
  if (fileType == 'annovar') {
    result <- readAnnovarMutFile(filePath, cores)
  } else {
    result <- readMafMutFile(filePath, cores)
  }
  
  result[, fileType := fileType]
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
                    default = NULL, help = subtypeHelp)

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
                    default = '', type = 'character',
                    help = targetGenomePathHelp)

targetGenomeHelp <- paste("Genome version, i.e. hg19, to which input variants",
                          "files for software should be brought.",
                          "Default: hg19.")
parser$add_argument("-g", "--target_genome_version", required = F,
                    default = 'hg19', type = 'character',
                    help = targetGenomeHelp)

chainHelp <- paste('Path to chain file in case genome version of mutations is',
                   'not the same as --target_genome_version')
parser$add_argument("-l", "--chain", required = F, 
                    type = 'character', help = chainHelp)

softwareHelp <- paste("Format of the mutational file output (name of the",
                      "software the file will be used with). One or several",
                      "from the list: chasmplus, digdriver, dndscv,",
                      "mutpanning, nbr, oncodrivefml, maf. If maf is given,",
                      'only one MAF file containing all variants will be',
                      'created')
parser$add_argument("-s", "--software", required = T, type = 'character',
                    nargs = "+", help = softwareHelp,
                    choices = c('chasmplus', 'digdriver', 'dndscv', 
                                'mutpanning', 'nbr', 'oncodrivefml', 'maf'))

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output folder")

coresHelp <- 'How many cores the script should use. Default: 1.'
parser$add_argument("-n", "--cores", required = F, type = 'integer', 
                    default = 1, help = coresHelp)

# Test input args -------------------------------------------------------------
args <- list(inventory_patients = '../data/inventory/inventory_patients_tcga.csv',
             inventory_analysis = '../data/inventory/inventory_analysis_tcga.csv',
             blacklist_inventory = '../data/inventory/inventory_blacklist_tcga.csv', 
             cancer_subtype = 'LUSC', min_depth = 30, 
             min_tumor_vac = 10, max_germline_vac = 5,
             min_tumor_vaf = 5.0, max_germline_vaf = 1.0, 
             max_n_vars = 90000,
             target_genome_version = 'hg19',
             target_genome_path = '../data/assets/reference_genome/hg19.fa',
             chain = '../data/assets/reference_genome/hg38ToHg19.over.chain', 
             software = list('chasmplus', 'digdriver', 'dndscv',
                             'mutpanning', 'nbr', 'oncodrivefml'),
             output = 'test', cores = 2)

# Olala ----------
args <- parser$parse_args()
check_input_arguments(args, outputType = 'folder')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# READ inventories ------------------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients, args$cores)
# check that if somatic_genome in patientsInv is not the same as target one
# chain file is submitted
if (any(!unique(patientsInv$somatic_genome) %in% args$target_genome_version) &
    is.null(args$chain)) {
  stop('[', Sys.time(), '] genome version of some mutation tables is not ',
       'the same as --target_genome_version, but no chain file is provided')
}

analysisInv <- readAnalysisInventory(args$inventory_analysis, args$cores)
if (!is.null(args$blacklist_inventory)) {
  bwInv <- readBlacklistInventory(args$blacklist_inventory, args$cores)
}

# select cancer subtype
patientsInv <- patientsInv[tumor_subtype %in% args$cancer_subtype]
analysisInv <- analysisInv[tumor_subtype %in% args$cancer_subtype]

message('[', Sys.time(), '] Read inventories')

# READ mutations files --------------------------------------------------------
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

# FILTER variants: noncanonical chr, duplicated, N var/ref, same alt/ref ------
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

# FILTER variants : VAF, ALT counts, etc --------------------------------------
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

# 5) Maximum germline VAC
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

# LIFTOVER variants to target genome ------------------------------------------
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

# FILTER variants: location in black-/white- listed regions -------------------
if (!is.null(bwInv)) {
  for (i in 1:nrow(bwInv)) {
    allVars <- filterVarsByBlackWhiteList(varsDT = allVars,
                                          bwFile = bwInv$file_path[i],
                                          chrStyle = outChrStyle,
                                          bwName = bwInv$list_name[i],
                                          bwFileType = bwInv$file_type[i],
                                          bwScoreCol = bwInv$score_column[i],
                                          bwScoreMin = bwInv$min_value[i],
                                          bwScoreMax = bwInv$max_value[i],
                                          cores = args$cores)
  }
}

# FILTER patients: number of variants does not exceed max_n_vars --------------
varsPerParticip <- allVars[,.N, by = participant_id]
before <- nrow(varsPerParticip)
removedPatients <- varsPerParticip[N > args$max_n_vars]$participant_id
varsPerParticip <- varsPerParticip[N <= args$max_n_vars]
after <- nrow(varsPerParticip)
if (after != before) {
  message('[', Sys.time(), '] Removed ', before - after, '(',
          round(100 * (before - after) / before, 2), '%) patients becasuse',
          ' number of mutations exceeded', args$max_n_vars, '. Removed ',
          'patients: ', paste(removedPatients, collapse = ', '))
}
allVars <- allVars[participant_id %in% varsPerParticip$participant_id]

# SAVE annotated cohort data as MAF file & as table ---------------------------
# bring allVars chromosome notation to the one of reference genome, yes, again
allVars[, chr := gsub('chr', '', chr)] # this covers NCBI
if (outChrStyle == 'UCSC') {
  allVars[, chr := paste0('chr', chr)] 
}

patientsInv[, participant_id := as.character(participant_id)]
allVars[, participant_id := as.character(participant_id)]
setkey(patientsInv, participant_id)
# do merge instead
allVars[, patient_tumor_subtype := patientsInv[allVars$participant_id]$tumor_subtype]
allVarsMAF <- mutTabToMAF(allVars)

write.table(allVarsMAF, append = F, sep = '\t', row.names = F, col.names = T,
            file = paste0(args$output, '/inputMutations-', args$cancer_subtype, 
                          '-', args$target_genome_version, '.maf'), quote = F)

# [OUTPUT] to files in driver-calling input format ----------------------------
# bring allVars chromosome notation to the one of reference genome, yes, again
allVars[, chr := gsub('chr', '', chr)] # this covers NCBI
if (outChrStyle == 'UCSC') {
  allVars[, chr := paste0('chr', chr)] 
}

# output variants for all formats
for (software in setdiff('maf', args$software)) {
  message('[', Sys.time(), '] Outputting variants for ', software)
  outfile <- paste0(args$output, '/' , software, '/variants/', software, 
                    '-inputMutations-', args$cancer_subtype, '-',
                    args$target_genome_version, '.csv')
  
  writeVarsToFile(allVars, software, outfile)
  if (software == 'mutpanning') {
    # mutpanning also requires txt file with samples, similar to patientsInv
    mutpanInvent <- patientsInv[, c('participant_id', 'participant_id',
                                    'tumor_subtype')]
    colnames(mutpanInvent) <- c('ID',	'Sample', 'Subtype')
    mutpanInvent[, ConfidenceLevel := 'custom']
    mutpanInvent[, Study := 'custom']
    mutpanInvent[, Cohort := 'custom']
    
    # output to file
    outfile <- paste0(args$output, '/', software, '/variants/',
                      software, '-patientsInv-', args$cancer_subtype, '-',
                      args$target_genome_version, '.csv')
    write.table(mutpanInvent, append = F, quote = F, row.names = F, outfile,
                col.names = T, sep = '\t')
  }
}

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')