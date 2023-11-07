#!/usr/bin/env Rscript
# FILE: 1_check_inventories.R -------------------------------------------------
#
# DESCRIPTION: Checks if submitted by user patients', analysis and black&white
# lists inventories are valid.
# Checks performed for patients' inventory:
# 1) inventory has all needed columns participant_id, tumor_subtype, 
#    participant_tumor_subtype, somatic_path, somatic_genome, cohort_name and 
#    columns tumor_subtype, participant_tumor_subtype, somatic_genome and 
#    cohort_name are not numbers
# 2) participant_id is linked to one and only one participant_tumor_subtype
# 3) all participant_id - tumor_subtype pairs are unique 
# 4) tumor_subtype do not contain - sign because it's going to be used in 
#    file naming
# 5) mutation files exist 
# 6) genome version is the same for all variant files 
#
# Checks performed for black & white lists inventory:
# 1) all needed columns are present (list_name, file_path, file_genome, 
#    file_type, score_column, min_value and max_value) are present and columns
#    list_name, file_type, file_genome and score_column are not numbers
# 2) each list_name appears once and only once
# 3) each file_path appears once and only once
# 4) files listed in file_path exist
# 5) codes in file_type are "white" or "black" 
# 6) min_value and max_value are given or 
# 7) black& white lists are on the same genome version as target genome version
#
# Checks performed for analysis inventory:
# 1) all needed columns are present and columns tumor_subtype, gr_id, gr_code,
#    gr_genome, gr_excl_id, gr_excl_code, gr_excl_genome, blacklisted_codes
#    are not numbers
# 2) files exist
# 3) software are one of the accepted ones
# 4) software can handle gr_code, i.e dndscv can only process CDS
# 5) gr_code and gr_excl_code do not contain - sign because it's going to be 
#    used in file naming.
# 6) in gr_file and gr_excl_file only gtf or bed files are supplied
# 7) if gr_code / gr_excl_code is given, that the file is provided as well. 
# 8) gr_code is acceptable, if gr_file and gr_excl_file are gtf 
# 9) for promoters and splice sites upstream and downstream are not 0 at the
#    same time if gtf file is given 
# 10) within same gr_id, all lines have same value of blacklisted codes 
# 11) within same tumor_subtype, all lines have same value of blacklisted 
#     codes
# 12) there are <= 2 genome versions, including target genome version
# 13) files in gr_file and gr_excl_file always have the same genome version 
#     assigned to the same file
# 14) gr_id is defined by the same files and upstream/downstream combinations
#     across tumor_subtypes and software
# 15) check, that if dndscv is requested, a gtf file is given
# 16) check, that if digdriver is requested, a gtf file is given, one per 
#     tumor
# 17) If union_percentage/intersect_percentage column are given, check that 
#     it is NA for rows with gr_code equals to CDS, and values in it in range
#     of 0 to 100 for the other rows              
#
# USAGE: Rscript --vanilla 1_check_patients_inventory.R \
#                --inventory_patients [path to your file] \
#                --inventory_analysis [path to your file] \
#                --inventory_blacklisted [path to your file] \ # optional
#                --min_n_participants 9 \ # optional
#                --target_genome_version hg19
#
# OPTIONS: Run 
#          Rscript --vanilla 1_check_patients_inventory.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: argparse, data.table
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  22.06.2023
# REVISION: 26.07.2023

box::use(./custom_functions_preprocessing[...])
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

# Functions -------------------------------------------------------------------
#' checkHaveEssentialColumns
#' @description Checks that all columns needed to be in inventory table are 
#'              indeed present there.
#' @author Maria Litovchenko
#' @param filePath path to inventory
#' @param essentialCols string vector, containing essential column names
#' @param inventory_type string, type of inventory (patient, analysis, etc)
#' @param cores integer, number of cores to be used
#' @return data table, in case all columns are found and stops the script 
#'         otherwise.
#' @export
checkHaveEssentialColumns <- function(filePath, essentialCols, 
                                      inventory_type, cores = 1) {
  result <- suppressWarnings(fread(filePath, header = T, stringsAsFactors = F, 
                                   nThread = cores))
  if (!all(essentialCols %in% colnames(result))) {
    stop('[', Sys.time(), '] ', inventory_type, ' inventory must have ',
         'columns: ', paste(essentialCols, collapse = ', '), '. Absent ',
         'columns: ', paste(setdiff(essentialCols, colnames(result)),
                            collapse = ', '))
  }
  result
}

#' checkColumnNotNumber
#' @description Checks that there are no numbers in the vector 
#' @author Maria Litovchenko
#' @param x vector
#' @param x_name name of that vector, if available
#' @return void, if no number in the vector was found and stops script 
#'         execution otherwise
#' @export
checkColumnNotNumber <- function(x, x_name = NULL) {
  x_to_int <- suppressWarnings(as.integer(x))
  
  if (any(!is.na(x_to_int))) {
    int_in_x <- x[!is.na(x_to_int)]
    stop('[', Sys.time(), '] Found a numerical value(s) in ', x_name, ': ',
         paste(int_in_x, collapse = ', '), '. Not allowed.')
  }
}

#' checkParticipantInventory
#' @description Checks for validity inventory patient table. Checks that:
#' 1) inventory has all needed columns participant_id, tumor_subtype, 
#'    participant_tumor_subtype, somatic_path, somatic_genome, cohort_name and 
#'    columns tumor_subtype, participant_tumor_subtype, somatic_genome and 
#'    cohort_name are not numbers
#' 2) participant_id is linked to one and only one participant_tumor_subtype
#' 3) all participant_id - tumor_subtype pairs are unique 
#' 4) tumor_subtype do not contain - sign because it's going to be used in 
#'    file naming
#' 5) mutation files exist 
#' 6) genome version is the same for all variant files 
#' @author Maria Litovchenko
#' @param inventoryPath path to file with participants inventory
#' @param cores integer, number of cores to use
#' @return data table with the inventory
checkParticipantInventory <- function(inventoryPath, cores = 1) {
  # 1) read and check that all needed columns are present
  essenCols <- c('tumor_subtype', 'participant_tumor_subtype', 
                 'participant_id', 'somatic_path', 'somatic_genome',
                 'cohort_name')
  result <- checkHaveEssentialColumns(inventoryPath, essenCols,
                                      'participant', cores)
  
  # check, that values in tumor_subtype column are not numbers
  checkColumnNotNumber(result$cohort_name, 'cohort_name')
  checkColumnNotNumber(result$tumor_subtype, 'tumor_subtype')
  checkColumnNotNumber(result$somatic_genome, 'somatic_genome')
  checkColumnNotNumber(result$participant_tumor_subtype, 
                       'participant_tumor_subtype')
  
  # 2) participant_id is linked to one and only one participant_tumor_subtype
  n_tumorSubt_by_samp <- result[,(length(unique(participant_tumor_subtype))),
                                by = participant_id]
  if (any(n_tumorSubt_by_samp$V1 > 1)) {
    stop('[', Sys.time(), '] participant_id(s)', 
         paste(n_tumorSubt_by_samp[V1 > 1]$participant_id, collapse = ','),
         ' have > 1 unique participant_tumor_subtype assigned to each of ',
         'them.')
  }
  
  # 3) all participant_id - tumor_subtype pairs are unique
  nbefore <- nrow(result)
  result <- result[!duplicated(result[,.(participant_id, tumor_subtype)])]
  nafter <- nrow(result)
  if (nafter < nbefore) {
    stop('[', Sys.time(), '] Detected ', nbefore - nafter, '(',
         round(100 * (nbefore - nafter) / nbefore, 2), '%) entries in the ',
         'patient inventory where participant_id - tumor_subtype pairs are ',
         'duplicated')
  }
  
  # 4) check, that tumor_subtype does not contain - sign because it's going to 
  # be used in file naming
  if (any(grepl('-', result$tumor_subtype))) {
    stop('[', Sys.time(), '] values of tumor_subtype columns should not ',
         'contain - character as it will be used in future file names as ',
         'separator. Please rename. Offending values: ',
         paste(result$tumor_subtype[grepl('-', result$tumor_subtype)],
               collapse = ', '), '.')
  }
  
  # 5) check that input mutation files exist
  fileExist <- file.exists(result$somatic_path)
  if (!all(fileExist)) {
    stop('[', Sys.time(), '] Files do not exist: ', 
         paste(result[!fileExist]$somatic_path, collapse = ',\n'))
  }
  
  # 6) genome version is the same for all variant files 
  n_genomes <- length(unique(result$somatic_genome))
  if (n_genomes > 1) {
    stop('[', Sys.time(), '] all mutation files should have the same genome ',
         'version.')
  }
  
  result
}

#' checkAnalysisInventory
#' @description Checks that in analysis inventory : 
#' 1) all needed columns are present and columns tumor_subtype, gr_id, gr_code,
#'    gr_genome, gr_excl_id, gr_excl_code, gr_excl_genome, blacklisted_codes
#'    are not numbers
#' 2) files exist
#' 3) software are one of the accepted ones
#' 4) software can handle gr_code, i.e dndscv can only process CDS
#' 5) gr_code and gr_excl_code do not contain - sign because it's going to be 
#'    used in file naming.
#' 6) in gr_file and gr_excl_file only gtf or bed files are supplied
#' 7) if gr_code / gr_excl_code is given, that the file is provided as well. 
#' 8) gr_code is acceptable, if gr_file and gr_excl_file are gtf 
#' 9) for promoters and splice sites upstream and downstream are not 0 at the
#'    same time if gtf file is given 
#' 10) within same gr_id, all lines have same value of blacklisted codes
#' 11) within same tumor_subtype, all lines have same value of blacklisted 
#'     codes
#' 12) there are <= 2 genome versions, including target genome version
#' 13) files in gr_file and gr_excl_file always have the same genome version 
#'     assigned to the same file
#' 14) gr_id is defined by the same files and upstream/downstream combinations
#'     across tumor_subtypes and software
#' 15) check, that if dndscv is requested, a gtf file is given
#' 16) check, that if digdriver is requested, a gtf file is given, one per 
#'     tumor
#' 17) If union_percentage/intersect_percentage column are given, check that 
#'     it is NA for rows with gr_code equals to CDS, and values in it in range
#'     of 0 to 100 for the other rows 
#' @author Maria Litovchenko 
#' @param inventoryPath analysis inventory path
#' @param acceptedRegCodes character vector: accepted region codes in case gtf
#' file is given.
#' @param acceptedSoftware named list with names = software names and values 
#' genomic regions codes which that software can handle. If software can handle
#' any regions, put NULL.
#' @param targetGenomeVersion character, genome version, in which final files
#' should be.
#' @param inventoryType string, one of "analysis" or "background_regions"
#' @param cores integer, number of cores to use
#' @return data table with the inventory
checkAnalysisInventory <- function(inventoryPath, acceptedRegCodes, 
                                   acceptedSoftware, targetGenomeVersion,
                                   cores = 1) {
  
  # 1) all needed columns are present and columns tumor_subtype, gr_id, gr_code,
  #    gr_genome, gr_excl_id, gr_excl_code, gr_excl_genome, blacklisted_codes
  #    are not numbers
  essenCols <- c('tumor_subtype', 'software', 'gr_id', 'gr_code', 'gr_file',
                 'gr_upstr', 'gr_downstr', 'gr_genome', 'gr_excl_id',
                 'gr_excl_code', 'gr_excl_file', 'gr_excl_upstr',
                 'gr_excl_downstr', 'gr_excl_genome', 'blacklisted_codes')
  result <- checkHaveEssentialColumns(inventoryPath, essenCols,
                                      'analysis', cores)
  # check, that values in tumor_subtype column are not numbers
  checkColumnNotNumber(result$gr_id, 'gr_id')
  checkColumnNotNumber(result$gr_code, 'gr_code')
  checkColumnNotNumber(result$gr_genome, 'gr_genome')
  checkColumnNotNumber(result$gr_excl_id, 'gr_excl_id')
  checkColumnNotNumber(result$gr_excl_code, 'gr_excl_code')
  checkColumnNotNumber(result$tumor_subtype, 'tumor_subtype')
  checkColumnNotNumber(result$gr_excl_genome, 'gr_excl_genome')
  checkColumnNotNumber(result$blacklisted_codes, 'blacklisted_codes') 
  
  result[is.na(gr_upstr)]$gr_upstr <- 0
  result[is.na(gr_downstr)]$gr_downstr <- 0
  result[is.na(gr_excl_upstr)]$gr_excl_upstr <- 0
  result[is.na(gr_excl_downstr)]$gr_excl_downstr <- 0
  result[result == ''] <- NA # replace empty values with NA
  
  # 2) check, that files exist
  # data table with just all submitted files and their extensions.
  result[, gr_file_ext := gsub('.*[.]', '', gsub('.gz$', '', gr_file))]
  result[, gr_excl_file_ext := gsub('.*[.]', '', gsub('.gz$', '', 
                                                      gr_excl_file))]
  files <- rbind(setnames(result[,.(gr_file, gr_file_ext)], 
                          c('gr_file', 'gr_file_ext'), c('filePath', 'ext')),
                 setnames(result[,.(gr_excl_file, gr_excl_file_ext)],
                          c('gr_excl_file', 'gr_excl_file_ext'),
                          c('filePath', 'ext')))
  files <- unique(files[!is.na(filePath)])
  notExist <- unique(files[!file.exists(filePath)])
  if (nrow(notExist) != 0) {
    stop('[', Sys.time(), '] Files do not exist: ', 
         paste(notExist$filePath, collapse = ', '))
  }
  
  # 3) check, that software are one of the accepted ones
  if (any(!result$software %in% names(acceptedSoftware))) {
    stop('[', Sys.time(), '] software ', 
         paste(unique(result[!software %in% 
                               names(acceptedSoftware)]$software), 
               collapse = ', '), ' are not supported.')
  }
  
  # 4) check, that software can handle gr_code, i.e dndscv can only process CDS
  setkey(result, 'software')
  for (soft in unique(result$software)) {
    isOk <- ifelse(is.null(acceptedSoftware[[soft]]), T, 
                   all(result[soft]$gr_code %in% acceptedSoftware[[soft]]))
    if (!isOk) {
      stop('[', Sys.time(), '] ', soft, ' does not support ', 
           paste(result[soft][!gr_code %in% acceptedSoftware[[soft]]],
                 collapse = ', '), ' regions.')
    }
  }
  
  # 5) check, that gr_code and gr_excl_code do not contain - sign because it's
  # going to be used in file naming
  if (any(grepl('-', result$gr_code) | grepl('-', result$gr_excl_code))) {
    stop('[', Sys.time(), '] values of gr_code and gr_excl_code columns ',
         'should not contain - character as it will be used in future file ', 
         'names as separator. Please rename.')
  }
  
  # 6) check that in gr_file and gr_excl_file only gtf or bed files are 
  # supplied
  notSupported <- files[!ext %in% c('gtf', 'bed')]
  if (nrow(notSupported) != 0) {
    stop('[', Sys.time(), '] Only files with bed.gz, bed, gtf.gz and gtf are ',
         'accepted in gr_file and gr_excl_file columns. Offending values: ',
         paste(notSupported$filePath, collapse = ', '))
  }
  
  # 7) check that if gr_code / gr_excl_code is given, that the file is provided
  # as well
  fileNotGiven <- result[(!is.na(gr_code) & is.na(gr_file)) | 
                           (!is.na(gr_excl_code) & is.na(gr_excl_file))]
  fileNotGiven <- unique(fileNotGiven)
  if (nrow(fileNotGiven) != 0) {
    stop('[', Sys.time(), '] target and/or excluded regions (',
         paste(c(fileNotGiven$gr_code, fileNotGiven$gr_excl_file),
               collapse = ', '), ') are not given the file to be read from.')
  }
  
  # 8) check, that gr_code is acceptable, if gr_file and gr_excl_file are gtf
  if (nrow(result[gr_file_ext == 'gtf']) != 0) {
    resultGTF <- result[gr_file_ext == 'gtf']
    if (any(!resultGTF$gr_code %in% acceptedRegCodes)) {
      stop('[', Sys.time(), '] Regions codes ', 
           paste0(unique(resultGTF[!gr_code %in% acceptedRegCodes]$gr_code), 
                  collapse = ', '), ' are not acceptable if gtf file is given',
           ' in gr_file. Please use custom bed file, if needed.')
    }
    rm(resultGTF)
  }
  if (nrow(result[gr_excl_file_ext == 'gtf']) != 0) {
    resultGTF <- result[gr_excl_file_ext == 'gtf']
    if (any(!resultGTF$gr_excl_code %in% acceptedRegCodes)) {
      stop('[', Sys.time(), '] Regions codes ', 
           paste0(unique(resultGTF[!gr_excl_code %in% 
                                     acceptedRegCodes]$gr_excl_code), 
                  collapse = ', '), ' are not acceptable if gtf file is given',
           ' in gr_file. Please use custom bed file, if needed.')
    }
    rm(resultGTF)
  }
  
  # 9) check, that for promoters and splice sites upstream and downstream are 
  #    not 0 at the same time
  resultPromSS <- result[gr_code %in% c("lincRNA_promoter", "lincRNA_ss",
                                        "promoter", "ss") & 
                           gr_file_ext == 'gtf']
  resultPromSS <- unique(resultPromSS)
  if (nrow(resultPromSS) != 0) {
    resultPromSS_both0 <- resultPromSS[gr_upstr == 0 & gr_downstr == 0]
    if (nrow(resultPromSS_both0) != 0) {
      stop('[', Sys.time(), '] gr_upstr and gr_downstr can not be 0 at ',
           'the same time for regions codes lincRNA_promoter, lincRNA_ss, ',
           'promoter, ss and GFT file in gr_file_ext. Offending values:\n',
           paste(apply(resultPromSS_both0, 1, paste, collapse = '\t'), 
                 collapse = '\t'))
    }
  }
  resultPromSS <- result[gr_excl_code %in% c("lincRNA_promoter", "lincRNA_ss",
                                             "promoter", "ss")& 
                           gr_excl_file_ext == 'gtf']
  resultPromSS <- unique(resultPromSS)
  if (nrow(resultPromSS) != 0) {
    resultPromSS_both0 <- resultPromSS[gr_excl_upstr == 0 & 
                                         gr_excl_downstr == 0]
    if (nrow(resultPromSS_both0) != 0) {
      stop('[', Sys.time(), '] gr_excl_upstr and gr_excl_downstr can not be ',
           '0 at the same time for regions codes lincRNA_promoter, ',
           'lincRNA_ss, promoter, ss and GFT file in gr_file_ext. Offending ',
           'values:\n', paste(apply(resultPromSS_both0, 1, paste,
                                    collapse = '\t'), collapse = '\t'))
    }
  }
  
  # 10) check, that within same gr_id, all lines have same value of blacklisted 
  #     codes
  # convert blacklisted_codes to lists
  blackListsDT <- data.table(blacklisted_codes = unique(result$blacklisted_codes))
  blackListsDT <- blackListsDT[complete.cases(blackListsDT)]
  if (nrow(blackListsDT) > 0) {
    blackListsDT[, blacklisted_codes_parsed := strsplit(blacklisted_codes,';')]
    colsOrder <- colnames(result)
    result <- merge(result, blackListsDT, by = 'blacklisted_codes', all.x = T)
    result[, blacklisted_codes := NULL]
    setnames(result, 'blacklisted_codes_parsed', 'blacklisted_codes')
    setcolorder(result, colsOrder)
    rm(colsOrder, blackListsDT)
  }
  result <- split(result, result$gr_id)
  sameBlackList <- sapply(result, 
                          function(x) length(unique(x$blacklisted_codes)) == 1)
  if (!all(sameBlackList)) {
    stop('[', Sys.time(), '] Rows with the same gr_id should have the same ',
         'blacklisted_codes. Offending gr_id: ',
         paste(names(sameBlackList[!sameBlackList]), collapse = ', '), '.')
  }
  result <- do.call(rbind, result)
  
  # 11) within same tumor_subtype, all lines have same value of blacklisted 
  #     codes
  result <- split(result, result$tumor_subtype)
  sameBlackList <- sapply(result, 
                          function(x) length(unique(x$blacklisted_codes)) == 1)
  if (!all(sameBlackList)) {
    stop('[', Sys.time(), '] Rows with the same tumor_subtype should have ',
         'the same blacklisted_codes. Offending tumor_subtype: ',
         paste(names(sameBlackList[!sameBlackList]), collapse = ', '), '.')
  }
  result <- do.call(rbind, result)
  
  
  # 12) check, that there are <= 2 genome versions, including target genome
  #     version
  genomeVers <- c(result$gr_genome, result$gr_excl_genome, targetGenomeVersion)
  genomeVers <- unique(genomeVers)
  genomeVers <- genomeVers[!is.na(genomeVers)]
  if (length(unique(genomeVers)) > 2) {
    stop('[', Sys.time(), '] No more than 2 genome versions, including ',
         '--target_genome_version are allowed in the analysis. Current genome ',
         'versions: ', paste(genomeVers, collapse = ', '), '.')
  }
  
  # 13) check, that files in gr_file and gr_excl_file always have the same 
  #     genome version assigned to the same file.
  fileToGenomeVers <- rbind(setNames(result[,.(gr_file, gr_genome)],
                                     c('filePath', 'genomeVersion')),
                            setNames(result[,.(gr_excl_file, gr_excl_genome)],
                                     c('filePath', 'genomeVersion')))
  fileToGenomeVers <- fileToGenomeVers[!duplicated(fileToGenomeVers)]
  fileToGenomeVers <- fileToGenomeVers[!is.na(filePath)]
  if (sum(complete.cases(fileToGenomeVers)) != nrow(fileToGenomeVers)) {
    stop('[', Sys.time(), '] Some files in gr_file and/or gr_excl_file ',
         'columns do not have genome version assigned in gr_genome and/or ',
         'gr_excl_genome columns.')
  }
  fileToGenomeVers <- fileToGenomeVers[,.(length(unique(genomeVersion))),
                                       by = filePath]
  if (max(fileToGenomeVers$V1) > 1) {
    stop('[', Sys.time(), '] File(s) ', 
         paste(fileToGenomeVers[V1 > 1]$filePath, collapse = ', '),
         'are given different genome version on different row. That is not ',
         'correct.')
  }
  
  # 14) check that gr_id is defined by the same files and upstream/downsteam 
  #     combinations across tumor_subtypes and softwares genome version
  grIdUniq <- lapply(split(result[, !colnames(result) %in% 
                                    c('restrictedTest', 'union_percentage',
                                      'intersect_percentage'),
                                  with = F], by = 'gr_id', drop = T),
                     split, by = c('tumor_subtype', 'software'), drop = T)
  grIdUniq <- lapply(grIdUniq, 
                     function(x) lapply(x, 
                                        function(y) y[, setdiff(colnames(y),
                                                                c('tumor_subtype',
                                                                  'software')),
                                                      with = F]))
  grIdUniq <- sapply(grIdUniq, function(x) length(unique(x)))
  if (any(grIdUniq != 1)) {
    stop('[', Sys.time(), '] gr_id(s) ', paste(names(grIdUniq[grIdUniq != 1]),
                                               collapse = ', '),
         ' are defined by different combinations of files/upstream/down',
         'stream values across tumor_subtype & software columns. This is not ',
         'permitted in order to perform fair comparisons. Please rename ',
         'gr_ids.')
  }
  
  # 15) check, that if dndscv is requested, a gtf file is given
  if ('dndscv' %in% result$software & 
      any(!result[software == 'dndscv']$gr_file_ext %in% 'gtf')) {
    stop('[', Sys.time(), '] Creation of regions input for dndscv from bed ',
         'files is not supported. Please submit gtf')
  }
  
  # 16) check, that if digdriver is requested, a gtf file is given, one per 
  #     tumor
  if ('digdriver' %in% result$software) {
    digGTFcheck <- result[software == 'digdriver']
    digGTFcheck <- digGTFcheck[,.(tumor_subtype, gr_code, gr_file_ext)]
    digGTFcheck <- unique(digGTFcheck)
    digNotOk <- digGTFcheck[,.(sum(gr_code == 'CDS') > 0,
                               sum(gr_file_ext == 'gtf') > 0),
                            by = .(tumor_subtype)]
    digNotOk <- digNotOk[V1 == F | V2 == F]
    if (nrow(digNotOk) != 0) {
      stop('[', Sys.time(), '] GTF file containing genome annotation is ',
           'essential for run of digdriver, even if you are only ',
           'interested in noncoding regions. Please add CDS (as gr_code) ',
           'with the correponding GTF file for each of the following tumor ',
           'subtypes: ', paste(digNotOk$tumor_subtype, collapse = ', '),
           ' to the analysis inventory')
    }
  }
  
  result[, gr_file_ext := NULL]
  result[, gr_excl_file_ext := NULL]
  
  # 17) If union_percentage/intersect_percentage column are given, check that 
  #     it is NA for rows with gr_code equals to CDS, and values in it in range
  #     of 0 to 100 for the other rows 
  if (!'union_percentage' %in% colnames(result)) {
    result[, union_percentage := NA]
  } else {
    if (nrow(result[gr_code == 'CDS' & !is.na(union_percentage)]) != 0) {
      stop('[', Sys.time(), '] Found entries in analysis table with gr_code ',
           'set to CDS and union_percentage not NA. Coding regions can not ',
           'be subjected to union/intersection due to presence of internal ',
           'structure (codons).')
    }
    if (min(result$union_percentage, na.rm = T) <= 0 | 
        max(result$union_percentage, na.rm = T) > 100) {
      stop('[', Sys.time(), '] Values in column union_percentage should ',
           'range from > 0 to <= 100')
    }
    n_uniq_union_perc <- result[,.(length(unique(union_percentage))), 
                                by = gr_id]
    n_uniq_union_perc <- min(n_uniq_union_perc$V1)
    if (n_uniq_union_perc > 1) {
      stop('[', Sys.time(), '] Found several union_percentage for one gr_id ',
           '. This is not allowed.')
    }
  }
  if (!'intersect_percentage' %in% colnames(result)) {
    result[, intersect_percentage := NA]
  } else {
    if (nrow(result[gr_code == 'CDS' & !is.na(intersect_percentage)]) != 0) {
      stop('[', Sys.time(), '] Found entries in analysis table with gr_code ',
           'set to CDS and intersect_percentage not NA. Coding regions can ',
           'not be subjected to union/intersection due to presence of ',
           'internal structure (codons).')
    }
    if (min(result$intersect_percentage, na.rm = T) <= 0 | 
        max(result$intersect_percentage, na.rm = T) > 100) {
      stop('[', Sys.time(), '] Values in column intersect_percentage should ',
           'range from > 0 to <= 100')
    }
    n_uniq_inter_perc <- result[,.(length(unique(intersect_percentage))), 
                                by = gr_id]
    n_uniq_inter_perc <- min(n_uniq_inter_perc$V1)
    if (n_uniq_inter_perc > 1) {
      stop('[', Sys.time(), '] Found several intersect_percentage for one ',
           'gr_id. This is not allowed.')
    }
  }
  result
}

#' checkBlacklistInventory
#' @description Checks for validity inventory table with black and white listed
#' regions. Checks that:
#' 1) all needed columns are present (list_name, file_path, file_genome, 
#'    file_type, score_column, min_value and max_value) are present and columns
#'    list_name, file_type, file_genome and score_column are not numbers
#' 2) each list_name appears once and only once
#' 3) each file_path appears once and only once
#' 4) files listed in file_path exist
#' 5) codes in file_type are "white" or "black" 
#' 6) min_value and max_value are given or 
#' 7) black& white lists are on the same genome version as target 
#'    genome version.
#' @author Maria Litovchenko
#' @param inventoryPath black&white lists inventory analysis path
#' @param targetGenomeVersion character, genome version, in which final files
#' should be.
#' @return data table with black&white lists inventory or NULL in case file was
#' empty
checkBlacklistInventory <- function(inventoryPath, targetGenomeVersion, 
                                    cores = 1) {
  if (file.size(inventoryPath) == 0) {
    return(NULL)
  }
  
  # 1) read and check that all needed columns are present
  essenCols <- c('list_name', 'file_path', 'file_genome', 'file_type', 
                 'score_column', 'min_value', 'max_value')
  result <- checkHaveEssentialColumns(inventoryPath, essenCols, 
                                      'black&white lists', cores)
  
  # check, that values in list_name, file_genome, file_type and score_column 
  # columns are not numbers
  checkColumnNotNumber(result$list_name, 'list_name')
  checkColumnNotNumber(result$file_type, 'file_type')
  checkColumnNotNumber(result$file_genome, 'file_genome')
  checkColumnNotNumber(result$score_column, 'score_column')
  
  # 2) each file_path appears once and only once
  n_appearences <- result[,.N, by = file_path]
  if (any(n_appearences$V1 > 1)) {
    stop('[', Sys.time(), '] ', 
         paste(n_appearences[V1 > 1]$file_path, collapse = ','), ' values ',
         'appear more than once in file_path column in the black&white lists ',
         'inventory')
  }
  
  # 3) each list_name appears once and only once
  n_appearences <- result[,.N, by = list_name]
  if (any(n_appearences$V1 > 1)) {
    stop('[', Sys.time(), '] ', 
         paste(n_appearences[V1 > 1]$list_name, collapse = ','), ' values ',
         'appear more than once in list_name column in the black&white lists ',
         'inventory')
  }
  
  # 4) files listed in file_path exist
  fileExist <- file.exists(result$file_path)
  if (!all(fileExist)) {
    stop('[', Sys.time(), '] Files do not exist: ', 
         paste(result[!fileExist]$file_path, collapse = ',\n'))
  }
  
  # 5) codes in file_type are "white" or "black"
  if (!all(result$file_type %in% c('white', 'black'))) {
    stop('[', Sys.time(), '] file_type in --blacklist_inventory should be ',
         'white or black', )
  }
  
  # 6) in case column score_column is not NA, min_value and max_value are given
  result[, min_value := as.numeric(min_value)]
  result[, min_value := as.numeric(max_value)]
  if (any(is.na(result[!is.na(score_column)]$min_value)) |  
      any(is.na(result[!is.na(score_column)]$max_value))) {
    stop('[', Sys.time(), '] min_value and max_value should be numeric ',
         '(not NA) in black&white lists inventory table if value in ',
         'score_column is not NA')
  }
  
  # 7) check, black& white lists are on the same genome version as target 
  #    genome version.
  genomeVers <- c(result$file_genome, targetGenomeVersion)
  genomeVers <- unique(genomeVers)
  genomeVers <- genomeVers[!is.na(genomeVers)]
  if (length(unique(genomeVers)) > 1) {
    stop('[', Sys.time(), '] Please submit files in black&white lists ',
         'inventories on the same genome version as --target_genome_version ',
         '(', targetGenomeVersion, '). Lift over of black&white lists can be ',
         'time & memory consuming and therefore is not implemented in this ',
         'pipeline.')
  }
  
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = '1_check_inventories.R')

patientsHelp <- paste('Path to inventory table listing information about',
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

blackListHelp <- paste('Path to inventory table containing details of the',
                       'black&white lists to use. Minimal columns:',
                       'list_name,', 'file_path,', 'file_genome,', 
                       'file_type,', 'score_column,', 'min_value,',
                       'max_value.')
parser$add_argument("-b", "--inventory_blacklisted", 
                    required = F, type = 'character', default = NULL,
                    help = blackListHelp)

parser$add_argument('-m', '--min_n_participants', 
                    required = F, type = 'integer', default = 1, 
                    help = 'Minimal number of participants per tumor subtype')

parser$add_argument("-g", "--target_genome_version", 
                    required = T, type = 'character', 
                    help = 'Target genome version, i.e. hg19')

args <- parser$parse_args()

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')

# check arguments are correct
check_input_arguments_preproc(args)
printArgs(args)

# Check the inventory----------------------------------------------------------
patientsInv <- checkParticipantInventory(args$inventory_patients)
print(paste0('[', Sys.time(), '] Patients inventory is OK'))

if (!is.null(args$inventory_blacklisted)) {
  bwInv <- checkBlacklistInventory(args$inventory_blacklisted, 
                                   targetGenomeVersion = args$target_genome_version)
  if (is.null(bwInv)) {
    args$inventory_blacklisted <- NULL
    print(paste0('[', Sys.time(), '] Black & white lists inventory file was ',
                 'empty. Proceeding without it.'))
  } else {
    print(paste0('[', Sys.time(), '] Black & white lists inventory is OK'))
  }
}

analysisInv <- checkAnalysisInventory(args$inventory_analysis,
                                      acceptedRegCodes = GR_CODES, 
                                      acceptedSoftware = SOFTWARE_GR_CODES, 
                                      targetGenomeVersion = args$target_genome_version)

# check, that if tumor subtype is requested for analysis there are samples for
# it in patientsInv
if(!all(unique(analysisInv$tumor_subtype) %in% 
        unique(patientsInv$tumor_subtype))) {
  noSamplesTumorSubtype <- setdiff(unique(analysisInv$tumor_subtype), 
                                     unique(patientsInv$tumor_subtype))
  stop('[', Sys.time(), '] Following tumor subtypes were requested for ',
       'analysis but do not have any corresponding to them samples: ',
       paste(noSamplesTumorSubtype, collapse = ', '), '.')
}

# check that if tumor subtype is requested for analysis it has a minimum 
# number of samples
n_samp_by_subtype <- patientsInv[,.(length(unique(participant_id))), 
                                 by = tumor_subtype]
if (any(n_samp_by_subtype$V1 < args$min_n_participants)) {
  stop('[', Sys.time(), '] Following tumor subtypes have < ',
       args$min_n_participants, ' number of participants: ',
       paste(n_samp_by_subtype[V1 < args$min_n_participants]$tumor_subtype,
             ', '), '.')
}

# check that analysis inventory does not contain black/white list codes which 
# are not listed in black/white lists inventory
if (!is.null(args$inventory_blacklisted)) {
  blackListCodesAnalysis <- unique(analysisInv$blacklisted_codes)
  blackListCodesAnalysis <- lapply(blackListCodesAnalysis, 
                                   function(x) strsplit(x, ';'))
  blackListCodesAnalysis <- unique(unlist(blackListCodesAnalysis))
  if(!all(blackListCodesAnalysis %in% bwInv$list_name)) {
    stop('[', Sys.time(), '] black&white list codes ', 
         paste(setdiff(blackListCodesAnalysis, bwInv$list_name), 
               collapse =  ', '), 'are requested in the analysis inventory ',
         'but not given in the black&white lists inventory.')
  }
}

# check, that across all inventories there are just 2 different genome versions
genomesVersions <- c(patientsInv$somatic_genome, analysisInv$gr_genome,
                     analysisInv$gr_excl_genome)
if (!is.null(args$inventory_blacklisted)) {
  genomesVersions <- c(genomesVersions, bwInv$file_genome)
}
genomesVersions <- unique(genomesVersions)
genomesVersions <- genomesVersions[!is.na(genomesVersions)]
if (length(genomesVersions) > 2) {
  stop('[', Sys.time(), '] > 2 different genome versions are detected across ',
       'all inventory tables. Current genome versions: ',
       paste0(genomesVersions, collapse = ', '), '. Please restrict it to ',
       'the maximum of 2 different genome versions.')
}

print(paste0('[', Sys.time(), '] Analysis inventory is OK'))

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')
