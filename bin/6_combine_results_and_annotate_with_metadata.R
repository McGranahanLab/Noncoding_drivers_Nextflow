#!/usr/bin/env Rscript
# FILE: combine_results_and_annotate_with_metadata.R --------------------------
#
# DESCRIPTION: A script to merge results of all software into one table (but 
#              not yet combine p-values with Brown method!) and annotate them 
#              with number of patients, number of mutations, local mutation 
#              rate, cds length, GTEx and TCGA expression
#
# USAGE: run in RStudio
# OPTIONS:
# EXAMPLE: 
# REQUIREMENTS: R 4.0.2
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  10.12.2020
# REVISION: 19.01.2022

suppressMessages(library(data.table)) 
suppressMessages(library(GenomicRanges))
suppressMessages(library(maftools))
suppressMessages(library(plyr))
suppressMessages(library(plyranges))
suppressMessages(library(poolr))
suppressMessages(library(reshape2))
suppressMessages(library(rtracklayer))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))

# [FUNCTIONS] Reading data ----------------------------------------------------
# Functions to read in R various files

#' printArgs
#' @description Prints submitted to the script arguments as messages.
#' @author Maria Litovchenko
#' @param argsList named list representing submitted to the script arguments
#' @return void
printArgs <- function(argsList) {
  message("Submitted arguments:")
  for (argName in names(args)) {
    oneArg <- args[[argName]]
    if (length(oneArg) > 1 & !is.null(names(oneArg))) {
      msg <- paste(apply(data.table(names(oneArg), oneArg), 1, paste, 
                         collapse = ' - '), collapse = ',')
    } else {
      msg <- paste(oneArg, collapse = ', ')
    }
    message(argName, ':', msg)
  }
}


#' readAndCheckAnalysisInventory
#' @description Reads into data table analysis inventory and checks that: 
#' 1) in gr_file and gr_excl_file only gtf or bed files are supplied 2) files 
#' exist 3) gr_code is acceptable, if gr_file and gr_excl_file are gtf 
#' 4) software are one of the accepted ones 5) software can handle gr_code,
#' i.e dndscv can only process CDS 6) for promoters and splice sites upstream 
#' and downstream are not 0 at the same time if gtf file is given 7) if 
#' gr_code / gr_excl_code is given, that the file is provided as well. 
#' 8) gr_code and gr_excl_code do not contain - sign because it's going to be 
#' used in file naming. 9) within same gr_id, all lines have same value of 
#' blacklisted codes 10) there are <= 2 genome versions, including target 
#' 11) files in gr_file and gr_excl_file always have the same genome version 
#' assigned to the same file. 12) gr_id is defined by the same files and 
#' upstream/downsteam combinations across tumor_subtypes and softwares
#' genome version
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
#' @return data table, analysis inventory.
readAndCheckAnalysisInventory <- function(inventoryPath, acceptedRegCodes, 
                                          acceptedSoftware, 
                                          targetGenomeVersion, 
                                          inventoryType = 'analysis') {
  if (!inventoryType %in% c("analysis", "background_regions")) {
    stop('[', Sys.time(), '] inventoryType should be analysis or ',
         'background_regions')
  }
  colsToRead <- c('tumor_subtype', 'restrictedTest', 'software', 'gr_id', 
                  'gr_code', 'gr_file', 'gr_upstr', 'gr_downstr', 'gr_genome',
                  'gr_excl_id', 'gr_excl_code', 'gr_excl_file', 
                  'gr_excl_upstr', 'gr_excl_downstr', 'gr_excl_genome', 
                  'blacklisted_codes', 'union_percentage', 
                  'intersect_percentage')
  result <- suppressWarnings(fread(inventoryPath, header = T, sep = '\t', 
                                   stringsAsFactors = F,
                                   select = colsToRead))
  
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
  
  # check that in gr_file and gr_excl_file only gtf or bed files are supplied
  result[, gr_file_ext := gsub('.*[.]', '', gsub('.gz$', '', gr_file))]
  result[, gr_excl_file_ext := gsub('.*[.]', '', gsub('.gz$', '', 
                                                      gr_excl_file))]
  # data table with just all submitted files and their extensions.
  files <- rbind(setnames(result[,.(gr_file, gr_file_ext)], 
                          c('gr_file', 'gr_file_ext'), c('filePath', 'ext')),
                 setnames(result[,.(gr_excl_file, gr_excl_file_ext)],
                          c('gr_excl_file', 'gr_excl_file_ext'),
                          c('filePath', 'ext')))
  files <- unique(files[!is.na(filePath)])
  notSupported <- files[!ext %in% c('gtf', 'bed')]
  if (nrow(notSupported) != 0) {
    stop('[', Sys.time(), '] Only files with bed.gz, bed, gtf.gz and gtf are ',
         'accepted in gr_file and gr_excl_file columns. Offending values: ',
         paste(notSupported$filePath, collapse = ', '))
  }
  
  # check, that files exist
  notExist <- unique(files[!file.exists(filePath)])
  if (nrow(notExist) != 0) {
    stop('[', Sys.time(), '] Files do not exist: ', 
         paste(notExist$filePath, collapse = ', '))
  }
  
  # check that if gr_code / gr_excl_code is given, that the file is provided
  # as well
  fileNotGiven <- result[(!is.na(gr_code) & is.na(gr_file)) | 
                           (!is.na(gr_excl_code) & is.na(gr_excl_file))]
  fileNotGiven <- unique(fileNotGiven)
  if (nrow(fileNotGiven) != 0) {
    stop('[', Sys.time(), '] target and/or excluded regions (',
         paste(c(fileNotGiven$gr_code, fileNotGiven$gr_excl_file),
               collapse = ', '), ') are not given the file to be read from.')
  }
  
  # check, that gr_code is acceptable, if gr_file and gr_excl_file are gtf
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
  
  if (inventoryType == "analysis") {
    # check, that software are one of the accepted ones
    if (any(!result$software %in% names(acceptedSoftware))) {
      stop('[', Sys.time(), '] software ', 
           paste(unique(result[!software %in% 
                                 names(acceptedSoftware)]$software), 
                 collapse = ', '), ' are not supported.')
    }
    
    # check, that software can handle gr_code, i.e dndscv can only process CDS
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
    
    # check, that if dndscv is requested, a gtf file is given
    if ('dndscv' %in% result$software & 
        any(!result['dndscv']$gr_file_ext %in% 'gtf')) {
      stop('[', Sys.time(), '] Creation of regions input for dndscv from bed ',
           'files is not supported. Please submit gtf')
    }
    
    # check, that if digdriver is requested, a gtf file is given, one per tumor
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
  }
  
  # check, that for promoters and splice sites upstream and downstream are not
  # 0 at the same time
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
  
  # check, that gr_code and gr_excl_code do not contain - sign because it's
  # going to be used in file naming
  if (any(grepl('-', result$gr_code) | grepl('-', result$gr_excl_code))) {
    stop('[', Sys.time(), '] values of gr_code and gr_excl_code columns ',
         'should not contain - character as it will be used in future file ', 
         'names as separator. Please rename.')
  }
  
  # check, that within same gr_id, all lines have same value of blacklisted 
  # codes
  result <- split(result, result$gr_id)
  sameBlackList <- sapply(result, 
                          function(x) length(unique(x$blacklisted_codes)) == 1)
  if (!all(sameBlackList)) {
    stop('[', Sys.time(), '] Rows with the same gr_id should have the same ',
         'blacklisted_codes. Offending gr_id: ',
         paste(names(sameBlackList[!sameBlackList]), collapse = ', '), '.')
  }
  result <- do.call(rbind, result)
  
  # check, that there are <= 2 genome versions, including target genome version
  genomeVers <- c(result$gr_genome, result$gr_excl_genome, targetGenomeVersion)
  genomeVers <- unique(genomeVers)
  genomeVers <- genomeVers[!is.na(genomeVers)]
  if (length(unique(genomeVers)) > 2) {
    stop('[', Sys.time(), '] No more than 2 genome versions, including ',
         '--target_genome_version are allowed in the analysis. Current genome ',
         'versions: ', paste(genomeVers, collapse = ', '), '.')
  }
  
  # check, that files in gr_file and gr_excl_file always have the same genome
  # version assigned to the same file.
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
  
  # check that gr_id is defined by the same files and upstream/downsteam 
  # combinations across tumor_subtypes and softwares genome version
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
  
  # convert blacklisted_codes to lists
  blackListsDT <- data.table(blacklisted_codes = unique(result$blacklisted_codes))
  blackListsDT <- blackListsDT[complete.cases(blackListsDT)]
  if (nrow(blackListsDT) > 0) {
    blackListsDT[, blacklisted_codes_parsed := strsplit(blacklisted_codes,',')]
    colsOrder <- colnames(result)
    result <- merge(result, blackListsDT, by = 'blacklisted_codes', all.x = T)
    result[, blacklisted_codes := NULL]
    setnames(result, 'blacklisted_codes_parsed', 'blacklisted_codes')
    setcolorder(result, colsOrder)
    rm(colsOrder, blackListsDT)
  }
  
  result[, is_coding := F]
  result[gr_id %in% result[gr_code == 'CDS']$gr_id]$is_coding <- T
  message('[', Sys.time(), '] Following gr_id will be considered as coding: ',
          paste(unique(result[is_coding == T]$gr_id), collapse = ', '))
  
  result[, gr_file_ext := NULL]
  result[, gr_excl_file_ext := NULL]
  
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
  if (!'restrictedTest' %in% colnames(result)) {
    result[, restrictedTest := F]
  }
  
  result
}

#' parseBED12regName
#' @description Parces name string(s) from BED12 files containing information
#' about genomic regions of interest. Usually those files are used with 
#' driverpower.
#' @author Maria Litovchenko
#' @param nameStr vector of strings
#' @return data table with number of rows = length of nameStr with columns 
#'         target_genome_version, gr_id, gene_id, gene_name
parseBED12regName <- function(nameStr, sepStr = '--') {
  result <- lapply(nameStr, strsplit, sepStr)
  result <- do.call(rbind, lapply(result, function(x) x[[1]]))
  result <- as.data.table(result)
  result[, name := nameStr]
  colnames(result) <- c('target_genome_version', 'gr_id', 'gene_id',
                        'gene_name', 'name')
  result
}

#' readSoftwareResults
#' @description Reads in result files from dndscv, driverpower, mutpanning, 
#' nbr, oncodrivefml into uniformal data table.
#' @author Maria Litovchenko
#' @param filePath path to result file
#' @param softName name of the software filePath corresponds to. One of dndscv,
#' driverpower, mutpanning, nbr, oncodrivefml
#' @param infoStr named vector with at least items tumor_subtype, gr_id
#' @return data table with columns software, tumor_subtype, gr_id
readSoftwareResults <- function(filePath, softName, infoStr) {
  if (!file.exists(filePath)) {
    message('[', Sys.time(), '] File not found: ', filePath, '. Return empty ',
            'data table.')
    return(data.table())
  }
  if (!softName %in% c('activedriverwgs', 'dndscv', 'driverpower',
                       'mutpanning', 'nbr', 'oncodrivefml', 'oncodriveclustl',
                       'digdriver')) {
    stop('[', Sys.time(), '] readSoftwareResults: softName should be one of ',
         'activedriverwgs, dndscv, driverpower, mutpanning, nbr, ',
         'oncodrivefml, oncodriveclustl, digdriver. Current value: ', softName)
  }
  
  if (!all(names(c('tumor_subtype', 'gr_id') %in% infoStr))) {
    stop('[', Sys.time(), '] readSoftwareResults: infoStr vector should have ',
         'items under names tumor_subtype, gr_id')
  }
  infoStr <- infoStr[c('tumor_subtype', 'gr_id')]
  infoStr <- as.data.table(t(infoStr))
  
  message('[', Sys.time(), '] Reading: ', filePath)
  
  if (softName == 'activedriverwgs') {
    result <- fread(filePath, header = T, stringsAsFactors = F, 
                    select = c('id', 'pp_element'))
    setnames(result, c('id', 'pp_element'), c('gene_id', 'raw_p'), 
             skip_absent = T)
    result <- cbind(result, infoStr)
  }
  
  if (softName %in% c('driverpower', 'digdriver')) {
    result <- fread(filePath, header = T, stringsAsFactors = F)
    
    colsToRead <- c('binID', 'CADD_p')
    if (softName == 'digdriver') {
      colsToRead <- c('ELT', 'PVAL_MUT_BURDEN')
      if (!'PVAL_SNV_BURDEN' %in% colnames(result) &
          !'PVAL_INDEL_BURDEN' %in% colnames(result) ) {
        stop('[', Sys.time(), '] DigDriver output lacks all 3 columns: ',
             'PVAL_SNV_BURDEN, PVAL_INDEL_BURDEN, PVAL_MUT_BURDEN') 
      }
      if (!'PVAL_MUT_BURDEN' %in% colnames(result) & 
          'PVAL_SNV_BURDEN' %in% colnames(result)) {
        colsToRead <- c('ELT', 'PVAL_SNV_BURDEN')
      }
      if (!'PVAL_MUT_BURDEN' %in% colnames(result) & 
          'PVAL_INDEL_BURDEN' %in% colnames(result)) {
        colsToRead <- c('ELT', 'PVAL_INDEL_BURDEN')
      }
    }
    
    result <- result[, colnames(result) %in% colsToRead, with = F]
    setnames(result, colsToRead, c('binID', 'raw_p'), skip_absent = T)
    result <- cbind(result, parseBED12regName(result$binID))
    result <- result[, !colnames(result) %in% c('binID', 'name', 
                                                'target_genome_version'),
                     with = F]
    setnames(result, 'gr_id', 'grID')
    result <- cbind(result, infoStr)
    result <- result[grID == gr_id]
    result[, grID := NULL]
  }
  
  if (softName == 'dndscv') {
    result <- fread(filePath, header = T, stringsAsFactors = F)
    if (!'pglobal_cv' %in% colnames(result)) {
      message('[', Sys.time(),  '] pglobal_cv is not found as a columns in ',
              'dNdScv output. Most likely that there were no indels. Using ',
              'pallsubs_cv instead. pallsubs_cv will be modified according ',
              'to formula: ',
              'pglobal_cv := 1 - pchisq(-2 * (log(pallsubs_cv)), df = 4) ',
              'see dNdScv code as reference')
      result[, pglobal_cv := 1 - pchisq(-2 * (log(pallsubs_cv)), df = 4)]
    }
    result <- result[,.(gene_name, pglobal_cv)]
    setnames(result,  'pglobal_cv', 'raw_p')
    result <- cbind(result, infoStr)
  }
  
  if (softName == 'mutpanning') {
    result <- fread(filePath, header = T, stringsAsFactors = F,
                    select = c('Name', 'Significance'))
    setnames(result, c('Name', 'Significance'), c('gene_name', 'raw_p'),
             skip_absent = T)
    result <- cbind(result, infoStr)
  }
  
  if (softName == 'nbr') {
    result <- fread(filePath, header = T, stringsAsFactors = F,
                    select = c('region', 'pval_subs_CV', 'pval_indels_CV',
                               'pval_both_CV'))
    result[, pval_subs_CV := as.double(pval_subs_CV)]
    result[, pval_indels_CV := as.double(pval_indels_CV)]
    if (!'pval_both_CV' %in% colnames(result)) {
      result[, pval_both_CV := NA]
    } else {
      result[, pval_both_CV := as.double(pval_both_CV)]
    }
    result[is.na(pval_both_CV)]$pval_both_CV <- result[is.na(pval_both_CV)]$pval_subs_CV
    result[is.na(pval_both_CV)]$pval_both_CV <- result[is.na(pval_both_CV)]$pval_indels_CV
    result <- result[,.(region, pval_both_CV)]
    setnames(result, c('region', 'pval_both_CV'), c('gene_id', 'raw_p'),
             skip_absent = T)
    if (nrow(result[is.na(raw_p)]) != 0) {
      message('[', Sys.time(), '] Found ', nrow(result[is.na(raw_p)]),
              ' entries in ', filePath, ' with all pval_subs_CV, ',
              'pval_indels_CV, pval_both_CV equal to NA. Removed them.')
      result <- result[!is.na(raw_p)]
    }
    result <- cbind(result, infoStr)
  }
  
  if (softName == 'oncodrivefml') {
    result <- fread(filePath, header = T, stringsAsFactors = F,
                    select = c('GENE_ID', 'SYMBOL', 'P_VALUE'))
    setnames(result, c('GENE_ID', 'SYMBOL', 'P_VALUE'), 
             c('gene_id', 'gene_name', 'raw_p'), skip_absent = T)
    result <- cbind(result, infoStr)
  }
  
  if (softName == 'oncodriveclustl') {
    result <- fread(filePath, header = T, stringsAsFactors = F,
                    select = c('ENSID', 'SYMBOL', 'P_ANALYTICAL'))
    setnames(result, c('ENSID', 'SYMBOL', 'P_ANALYTICAL'),
             c('gene_id', 'gene_name', 'raw_p'), skip_absent = T)
    result <- cbind(result, infoStr)
  }
  
  result[, software := softName]
  finalCols <- c('software', 'tumor_subtype', 'gr_id', 'gene_id', 'gene_name',
                 'raw_p')
  absentCols <- setdiff(finalCols, colnames(result))
  if (length(absentCols) > 0) {
    result <- cbind(result, 
                    setNames(data.table(matrix(nrow = nrow(result), 
                                               ncol = length(absentCols))),
                             absentCols))
  }
  result <- setcolorder(result, finalCols)
  result
}

#' getGeneIDtoGeneNameMap
#' @description Retrieved gene_id (ensembl_gene_id) tp gene_name (symbol) map
#' from bed12 files
#' @author Maria Litovchenko
#' @param bed12Files list of paths to bed12 files containing genomic ranges
#' @return data table with columns gene_id and gene_name
getGeneIDtoGeneNameMap <- function(bed12Files) {
  result <- lapply(bed12Files, fread, select = 4)
  result <- unique(unlist(result))
  result <- parseBED12regName(result)
  result <- result[,.(gene_id, gene_name)]
  result <- result[!duplicated(result)]
  
  if (max(result[,.N, by = gene_id]$N) > 1) {
    message('[', Sys.time(), '] gene_id to gene_name relationship is not ',
            'unique')
  }
  
  result
}

#' extractGeneToIDmap
#' @description Extracts a data table with columns gene_name, gene_id, 
#' gene_biotype, transcript_id, transcript_biotype (if available) from GRanges
#' extracted from gtf. transcript_biotype is not available for earlier versions
#' of ensemble gtf, i.e. v75 doesn't have it.
#' @author Maria Litovchenko
#' @param gtfGR GRanges object, result of import command applied to gtf file
#' @param tx_biotypes allowed transcript biotypes
#' @return data table with columns gene_name, gene_id, gene_biotype, 
#' transcript_id, transcript_biotype (if available). If tx_biotypes is not 
#' NULL, selection on transcript_biotype/gene_biotype will be performed.
extractGeneToIDmap <- function(gtfGR, tx_biotypes = NULL) {
  result <- as.data.table(mcols(gtfGR))
  result <- result[, intersect(c('gene_name', 'gene_id', 'gene_biotype',
                                 'transcript_id', 'transcript_biotype'),
                               colnames(result)), with = F]
  if ('transcript_biotype' %in% colnames(result)) {
    result <- result[complete.cases(transcript_biotype)]
  } else {
    result <- result[complete.cases(gene_biotype)]
  }
  result <- result[!duplicated(result)]
  
  if (!is.null(tx_biotypes)) {
    if ('transcript_biotype' %in% colnames(result)) {
      result <- result[transcript_biotype %in% tx_biotypes]
    } else {
      result <- result[gene_biotype %in% tx_biotypes]
    }
  }
  
  result
}

#' liftOverGenomicRegion
#' @description Lifts over coordinates from original genome to target genome
#' using chain chainObj
#' @author Maria Litovchenko
#' @param inGR input genomic regions. mcols are not used
#' @param chainObj chain to move genome coordinates from original genome 
#' version to the target one
#' @param min.gapwidt minimum length of a gap between lifted over regions that
#'                    will prevent them from being merged together. 
#'                    Default = Inf (no reduce will be performed).
#' @param min.width minimum width of lifted over genomic regions. Will only be
#'                  used if min.gapwidt is given and is not infinite. 
#'                  Default = 0 (no regions will be filtered out). 
#' @return lifted over genomic coordinates
liftOverGenomicRegion <- function(inGR, chainObj, min.gapwidt = Inf, 
                                  min.width = 0) {
  message('[', Sys.time(), '] Started liftover')
  if (class(inGR) == 'CompressedGRangesList') {
    inGRunlist <- unlist(inGR)
  } else {
    inGRunlist <- inGR
  }
  # change chromosomal format, if needed
  seqlevelsStyle(inGRunlist) <- seqlevelsStyle(chainObj)[1]
  
  # perform liftover
  result <- liftOver(inGRunlist, chainObj)
  # inform about not lifted regions
  notLifted <- unlist(lapply(width(result), 
                             function(x) identical(x, integer(0))))
  notLiftedPerc <- round(100 * sum(notLifted) / length(inGRunlist), 2)
  message('[', Sys.time(), '] ', sum(notLifted), ' (', notLiftedPerc,
          '%) genomic regions were not lifted over.')
  # inform if a region was split into 2
  splitted <- unlist(lapply(width(result), function(x) length(x) > 1))
  splittedPerc <- round(100 * sum(splitted) / length(splitted), 2)
  message('[', Sys.time(), '] ', sum(splitted), ' (', splittedPerc,
          '%) genomic regions were splitted into 2 and more.')
  
  # merge regions, if min.gapwidt is not infinite
  if (!is.infinite(min.gapwidt)) {
    message('[', Sys.time(), '] Reduction of lifted over genomic regions ',
            'located closer than ', min.gapwidt, 'bp will be performed.')
    result <- reduce(result, min.gapwidt = min.gapwidt)
    # reduce removed all mcols, add them back
    not_empty <- which(!notLifted)
    lo_length <- elementNROWS(result)[!notLifted]
    idx_inOrigin <- rep(not_empty, lo_length)
    result <- unlist(result)
    mcols(result) <- mcols(inGRunlist)[idx_inOrigin, ]
    
    # remove small regions, min.width is given
    if (min.width > 0) {
      message('[', Sys.time(), '] min.width is given. Summary of region ',
              'lengths before removal of regions with width < ', min.width)
      width_sum <- round(summary(width(result)), 0)
      message(paste(paste0(names(width_sum), collapse = '\t'), '\n', 
                    paste0(width_sum, collapse = '\t')))
      result <- result[width(result) >= min.width]
      message('[', Sys.time(), '] Summary of region lengths after removal of ',
              'regions with width < ', min.width)
      width_sum <- round(summary(width(result)), 0)
      message(paste(paste0(names(width_sum), collapse = '\t'), '\n', 
                    paste0(width_sum, collapse = '\t')))
    }
  } else {
    result <- unlist(result)
  }
  
  message('[', Sys.time(), '] Ratio of total region length after and ',
          'before liftover: ', round(sum(width(result)) / sum(width(inGR)), 3))
  message('[', Sys.time(), "] Ratio of regions' numbers after and before ",
          'liftover: ', round(length(result) / length(inGR), 3))
  
  # set style of chromosomal names to the same as input data
  seqlevelsStyle(result) <- seqlevelsStyle(inGR)
  
  message('[', Sys.time(), '] Finished liftover')
  result
}

#' gtfToTxDB
#' @description Reads GTF file into list containing TxDB object and a table 
#' with gene names, IDs, biotypes, as well as transcript ones.
#' @author Maria Litovchenko
#' @param gtfPath path to GTF file
#' @param acceptedChrCodes
#' @param doLiftOver
#' @param chainObj
#' @return list with 3 items: 1) type, string = gtf 2) geneMap: data table with
#' gene names, IDs, biotypes, as well as transcript ones 3) gtfTxdb txdb object
gtfToTxDB <- function(gtfPath, acceptedChrCodes = NULL, doLiftOver = F, 
                      chainObj = NULL) {  
  if (doLiftOver & is.null(chainObj)) {
    stop('[', Sys.time(), '] LiftOver is requested, but chain object for ',
         'lifover is NULL')
  }
  
  message('[', Sys.time(), '] Started reading ', gtfPath)
  result <- import(gtfPath)
  
  if (!is.null(acceptedChrCodes)) {
    result <- result[seqnames(result) %in% acceptedChrCodes]
    seqlevels(result, pruning.mode = "coarse") <- as.character(sort(unique(seqnames(result))))
    result <- sort(result)
  }
  
  if (any(c(23, '23', 'chr23') %in% seqlevels(result))) {
    message('[', Sys.time(), '] Found that chr X is coded with 23. Changed it',
            ' to X.')
    result <- as.data.table(result)
    result[, seqnames := gsub('23', 'X', seqnames)]
    result[, seqnames := gsub('24', 'Y', seqnames)]
    result[, seqnames := gsub('25', 'M', seqnames)]
    result <- makeGRangesFromDataFrame(result, keep.extra.columns = T)
  }
  seqlevelsStyle(result) <- "NCBI"
  message('[', Sys.time(), '] Finished reading ', gtfPath)
  
  if (doLiftOver & !is.null(chainObj)) {
    result <- liftOverGenomicRegion(inGR = result, chainObj)
  }
  
  # get a table with gene names, IDs, biotypes, as well as transcript ones
  geneMap <- extractGeneToIDmap(result)
  # convert to txdb & put to result list
  result <- list(type = 'gtf', 'geneMap' = geneMap, 
                 'gtfTxdb' = makeTxDbFromGRanges(result))
  result
}

# [FUNCTIONS] Annotation ------------------------------------------------------
#' fillInGeneIDs
#' @description Fills in missing gene_id -s based on gene_name or gene_name
#' synonyms and a map from gene_id to gene_name.
#' @author Maria Litovchenko
#' @param DT data table with columns gene_id and gene_name
#' @param id_To_Name_Map data table, map from gene_id to gene_name, columns 
#' gene_name and gene_id are essential
#' @param name_syns data table showing list of synonyms for gene_name. Should 
#' have columns gene_name and idx. Gene names with the same idx are considered
#' to be synonyms.
#' @return DT with filled in gene_id, if it was possible to do so.
fillInGeneIDs <- function(DT, id_To_Name_Map, name_syns = NULL) {
  missIDs <- DT[is.na(gene_id)]
  setkey(id_To_Name_Map, gene_name)
  
  # case 1: if gene_name is present in the idToNameMap, then retrieve gene_id 
  # using idToNameMap
  fixByName <- missIDs[gene_name %in% id_To_Name_Map$gene_name]$gene_name
  fixByName <- unique(fixByName)
  fixByName <- id_To_Name_Map[fixByName]
  # 1 to 1 match between gene_name and gene_id
  fixByName <- fixByName[order(gene_name, gene_id)]
  fixByName <- fixByName[,.SD[1], by = gene_name] 
  # fixByName now has 2 columns: gene_name and gene_id
  
  # case 2: if gene_name is NOT present in the id_To_Name_Map, try to find 
  # alternative synonymous gene_name which is present in id_To_Name_Map
  if (!is.null(name_syns)) {
    fixBySyn <- missIDs[!gene_name %in% id_To_Name_Map$gene_name]$gene_name
    fixBySyn <- unique(fixBySyn)
    fixBySyn <- name_syns[gene_name %in% fixBySyn]
    # record which gene name out of all synonyms was in the dt in order to 
    # trace it back
    setnames(fixBySyn, 'gene_name', 'gene_name_orig')
    # create a data table with original gene name in gene_name_orig column
    # and all possible synonyms in the gene_name column
    fixBySyn <- merge(name_syns[idx %in% fixBySyn$idx], fixBySyn, by = 'idx')
    fixBySyn[, idx := NULL]
    fixBySyn <- fixBySyn[gene_name != gene_name_orig]
    # restrict only to synonyms present in id_To_Name_Map
    fixBySyn <- fixBySyn[gene_name %in% id_To_Name_Map$gene_name]
    # 1 to 1 match between gene_name_orig and gene_name
    fixBySyn <- fixBySyn[order(gene_name_orig, gene_name)]
    fixBySyn <- fixBySyn[,.SD[1], by = gene_name_orig]
    # add gene_id
    fixBySyn <- merge(fixBySyn, id_To_Name_Map, by = 'gene_name')
    # here, some genes can be problematic as several gene_name_orig can be
    # matched to one gene_name - gene_id pair. To make it 1 to 1 match we 
    # could either select the most significant one or select a first one.
    # Let's go with the first one.
    fixBySyn <- fixBySyn[,.SD[1], by = gene_name_orig]
  }
  
  # to manage fixBySyn and fixByName together:
  fixByName[, gene_name_orig := gene_name]
  fixByNameOrSyn <- rbind(fixByName, fixBySyn)
  # now, there still can be a situation then several gene_name_orig can be
  # matched to one gene_name - gene_id pair, i.e. then first gene_name_orig
  # comes from fixByName and another one from fixBySyn. To make it 1 to 1 match
  # let's select first one. Otherwise we'll have different results for the same
  # gene from the same tool.
  fixByNameOrSyn <- fixByNameOrSyn[order(gene_id, gene_name_orig, gene_name)]
  fixByNameOrSyn <- fixByNameOrSyn[,.SD[1], by = gene_id]
  
  notFixed <- missIDs[!gene_name %in% fixByNameOrSyn$gene_name_orig]
  setnames(fixByNameOrSyn, c('gene_name_orig', 'gene_name', 'gene_id'),
           c('gene_name', 'gene_name_upd', 'gene_id_upd'))
  fixed <- merge(missIDs, fixByNameOrSyn, by = 'gene_name')
  fixed[, gene_id := NULL]
  fixed[, gene_name := NULL]
  setnames(fixed, c('gene_name_upd', 'gene_id_upd'), c('gene_name', 'gene_id'))
  
  result <- as.data.table(rbind(DT[!is.na(gene_id)], fixed, notFixed))
  result
}

#' fillInGeneName
#' @description Fills in missing gene_names -s based on gene_id a map from 
#' gene_id to gene_name.
#' @author Maria Litovchenko
#' @param DT data table with columns gene_id and gene_name
#' @param id_To_Name_Map data table, map from gene_id to gene_name, columns 
#' gene_name and gene_id are essential
#' @return DT with filled in gene_name, if it was possible to do so.
fillInGeneName <- function(DT, id_To_Name_Map) {
  missNames <- DT[is.na(gene_name)]
  setkey(id_To_Name_Map, gene_id)
  
  # gene_id is present in the id_To_Name_Map, retrieve gene_name using it
  fixByID <- missNames[gene_id %in% id_To_Name_Map$gene_id]$gene_id
  fixByID <- unique(fixByID)
  fixByID <- id_To_Name_Map[fixByID]
  # 1 to 1 match between gene_name and gene_id
  fixByID <- fixByID[order(gene_id, gene_name)]
  fixByID <- fixByID[,.SD[1], by = gene_id]
  setkey(fixByID, 'gene_id')
  # fixByID now has 2 columns: gene_name and gene_id
  
  notFixed <- missNames[!gene_id %in% fixByID$gene_id]
  fixed <- missNames[gene_id %in% fixByID$gene_id]
  fixed[, gene_name := fixByID[fixed$gene_id]$gene_name]
  result <- rbind(DT[!is.na(gene_name)], fixed, notFixed)
  result
}

# [FUNCTIONS] Other -----------------------------------------------------------
#' getUniqQuantileBreaks
#' @description Creates data table containing unique breaks and their labels
#' based on data vector
#' @param x numeric vector
#' @param nQuants number of quantiles
#' @return data table with columns topBound (use it as breaks) and quant (use 
#' it as labels)
getUniqQuantileBreaks <- function(x, nQuants = 100) {
  # add quantile of local mutation rate as well
  quantCuts <- data.table(quant = seq(0, 1, length.out = 100), 
                          topBound = quantile(x, seq(0, 1,length.out = 100)))
  quantCuts <- quantCuts[,.(quant = max(quant)), by = topBound]
  quantCuts[topBound == max(topBound)]$topBound <- Inf
  quantCuts[, quant := round(100 * quant, 2)]
  quantCuts
}

# Test inputs -----------------------------------------------------------------
regionCodes <- c('protein_coding', "3primeUTR", "5primeUTR", "CDS", "lincRNA", 
                 "lincRNA_promoter", "lincRNA_ss", "miRNA", "misc_RNA", 
                 "promoter", "rRNA", "snoRNA", "snRNA", "ss")
softRegCodes <- list('activedriverwgs' = NULL, 'dndscv' = c('CDS'), 
                     'mutpanning' = c('CDS'), 'chasmplus' = c('CDS'), 
                     'driverpower' = NULL, 'nbr' = NULL, 'oncodrivefml' = NULL,
                     'digdriver' = NULL)
acceptedChrNames <- c(c(1:24, 'X', 'Y'), paste0('chr', c(1:24, 'X', 'Y')))

args <- list(analysis_inventory = 'inventory/inventory_analysis_2023-03-21_local.csv',
             patients_inventory = 'inventory/inventory_patients_2023-03-21_local.csv',
             removed_patients = c('213000584', '213000692', '215003222', 
                                  '216000685', '216000960', '216001081',
                                  '221003178', '221004448'), 
             inputs_dir = 'inputs/2023-03-21/', 
             results_dir = 'results/2023-03-21/',
             target_genome_version = 'hg19',
             muts_to_gr_dir = 'results/2023-03-21/mutation_rate/',
             muts_to_gr_ptrn = 'mutRate-.*-hg19-mutMapToGR.csv',
             var_cat_ptrn = 'mutRate-.*-hg19-varCatEnrich.csv',
             gene_name_synonyms =  '_assets/additional_db_processed_2023-03-21/hgnc_complete_set_2022-07-01_proc.csv',
             known_cancer_genes = '_assets/additional_db_processed_2023-03-21/cgc_intogen_mc3_knownCancerGenes_lungGenesOnlyCOSMIC.csv',
             known_db_to_use = list('CGC'),
             olfactory_genes = '_assets/additional_db_raw/olfactory_barnes_2020.csv',
             gtex = '_assets/additional_db_processed_2023-03-21/GTEx_expression_in_tumorsubtypes.csv',
             tcga = '_assets/additional_db_processed_2023-03-21/TCGA_expression_in_tumorsubtypes.csv',
             output = 'results/2023-03-21_synRemoved/tables/raw_drivers_CGConly.csv')

if (!file.exists(args$analysis_inventory)) {
  stop('[',  Sys.time(), '] File does not exist: ', args$analysis_inventory)
}
if (!dir.exists(args$inputs_dir)) {
  stop('[',  Sys.time(), '] Folder does not exist: ', args$inputs_dir)
}
if (!dir.exists(args$results_dir)) {
  stop('[',  Sys.time(), '] Folder does not exist: ', args$results_dir)
}
if (!is.null(args$muts_to_gr_dir)) {
  if (!dir.exists(args$muts_to_gr_dir)) {
    stop('[',  Sys.time(), '] Folder does not exist: ', args$muts_to_gr_dir)
  }
}
if (!is.null(args$muts_to_gr_dir) & is.null(args$muts_to_gr_ptrn)) {
  stop('[',  Sys.time(), '] --muts_to_gr_ptrn is needed then ',
       '--muts_to_gr_dir is given.')
}
if (is.null(args$muts_to_gr_dir) & !is.null(args$muts_to_gr_ptrn)) {
  stop('[',  Sys.time(), '] --muts_to_gr_dir is needed then ',
       '--muts_to_gr_ptrn is given.')
}
if (!is.null(args$var_cat_ptrn) & is.null(args$muts_to_gr_dir)) {
  stop('[',  Sys.time(), '] --muts_to_gr_dir is needed then ',
       '--var_cat_ptrn is given.')
}
if (!file.exists(args$gene_name_synonyms)) {
  stop('[',  Sys.time(), '] File does not exist: ', args$gene_name_synonyms)
}
if (!file.exists(args$known_cancer_genes)) {
  stop('[',  Sys.time(), '] File does not exist: ', args$known_cancer_genes)
}
if (!is.null(args$olfactory_genes) & !file.exists(args$olfactory_genes)) {
  stop('[',  Sys.time(), '] File does not exist: ', args$olfactory_genes)
}
if (!is.null(args$gtex) & !file.exists(args$gtex)) {
  stop('[',  Sys.time(), '] File does not exist: ', args$gtex)
}
if (!is.null(args$tcga) & !file.exists(args$tcga)) {
  stop('[',  Sys.time(), '] File does not exist: ', args$tcga)
}
if (!dir.exists(dirname(args$output))) {
  dir.create(dirname(args$output), recursive = T)
  message('[', Sys.time(), '] Created diretory: ', dirname(args$output))
}

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# [READ] in inventory table and assign path to results ------------------------
patients_inventory <- fread(args$patients_inventory, header = T, 
                            stringsAsFactors = F)
patients_inventory <- patients_inventory[!participant_id %in% 
                                           args$removed_patients]
cohortSizeDT <- patients_inventory[,.(length(unique(participant_id))),
                                   by = tumor_subtype]
setnames(cohortSizeDT, 'V1', 'cohort_size')
cohortSizeDT <- rbind(cohortSizeDT, 
                      data.table(tumor_subtype = 'PANCAN',
                                 cohort_size = sum(cohortSizeDT[tumor_subtype != 'MET_PANCAN']$cohort_size)))

# read in inventory, check validity of the arguments
analisysInv <- readAndCheckAnalysisInventory(args$analysis_inventory, 
                                             acceptedRegCodes = regionCodes,
                                             acceptedSoftware = softRegCodes,
                                             targetGenomeVersion = args$target_genome_version)
analisysInv <- analisysInv[software != 'chasmplus']
# PATCH
analisysInv <- analisysInv[!(software == 'nbr' & gr_code == 'CDS')]
analisysInv[, restricted := T]

# get paths to gtf files to compute CDS length
gtfPaths <- analisysInv[gr_code == 'CDS' & grepl('gtf$|gtf.gz$', gr_file)]
gtfPaths <- unique(gtfPaths$gr_file)
if (length(gtfPaths) == 0) {
  stop('[', Sys.time(), '] No GTF files found in ', args$analysis_inventory,
       '. GTF files are needed for gene length calculation to perform ',
       'filtering baed on CDS length quantile (--max_cds_len_q)')
}

# restrict only to the fields we need
analisysInv <- analisysInv[,.(tumor_subtype, software, gr_id, restricted)]
analisysInv <- analisysInv[!duplicated(analisysInv)]

# add path to bed12 file with all scanned regions. This will also help us
# to get ensembl id to gene name matches because they are encoded in binID
analisysInv[, gr_bed12 := apply(analisysInv, 1, 
                                function(x) paste('inputGR', 
                                                  x['tumor_subtype'], 
                                                  args$target_genome_version, 
                                                  sep = '-'))]
analisysInv[, gr_bed12 := paste0(args$inputs_dir, '/', gr_bed12, '.bed')]

# add path to MAF file with all scanned mutations
analisysInv[, maf := apply(analisysInv, 1, 
                           function(x) paste('inputMutations', 
                                             x['tumor_subtype'], 
                                             args$target_genome_version, 
                                             sep = '-'))]
analisysInv[, maf := paste0(args$inputs_dir, '/', maf, '.maf')]

# add path to the results file:
analisysInv[, result_file := apply(analisysInv, 1, 
                                   function(x) paste(x['software'],
                                                     'results', 
                                                     x['tumor_subtype'], 
                                                     ifelse(x['software'] == 'driverpower',
                                                            '-', x['gr_id']), 
                                                     args$target_genome_version, 
                                                     sep = '-'))]
analisysInv[, result_file := paste0(args$results_dir, '/', software, '/',
                                    result_file, '.csv')]
# inform that certain files were not found
analisysInv[, fileFound := sapply(analisysInv$result_file, file.exists)]
analisysInv <- analisysInv[order(software)]
if (sum(!analisysInv$fileFound) != 0) {
  message('[', Sys.time(), '] Result files were not found: ', 
          paste(unique(analisysInv[fileFound == F]$result_file),
                collapse = '\n'))
}
# check, if there is a software for which we do not have results at all &
# exclude it
softNotRun <- analisysInv[,.(sum(fileFound)/.N), by = software][V1 == 0]
if (nrow(softNotRun) != 0) {
  message('[', Sys.time(), '] No results found for software(s): ', 
          paste(softNotRun$software, collapse = ', '), '. They will be ',
          'excluded from the analysis.')
  analisysInv <- analisysInv[!software %in% softNotRun$software]
}
analisysInv[, fileFound := NULL]
if (nrow(analisysInv) == 0) {
  stop('[', Sys.time(), '] Nothing to analyse!')
}

# [READ] in raw results of software runs --------------------------------------
rawResults <- apply(analisysInv, 1,
                    function(x) readSoftwareResults(x['result_file'], 
                                                    x['software'], x))
rawResults <- do.call(rbind, rawResults)

# Unify gene_id and gene_name across various software -------------------------
# create gene_id to gene_name map from bed12
GENE_ID_TO_NAME <- getGeneIDtoGeneNameMap(unique(analisysInv$gr_bed12))
GENE_NAME_SYNS <- NULL
# read data table with synonyms, if available 
if (!is.null(args$gene_name_synonyms)) {
  GENE_NAME_SYNS <- fread(args$gene_name_synonyms, header = T, 
                          stringsAsFactors = F)
}
rawResults <- fillInGeneIDs(rawResults, GENE_ID_TO_NAME, GENE_NAME_SYNS)
rawResults <- fillInGeneName(rawResults, GENE_ID_TO_NAME)

nbefore <- nrow(rawResults)
rawResults <- rawResults[!is.na(gene_id)]
nafter <- nrow(rawResults)
message('[', Sys.time(), '] Removed ', nbefore - nafter, ' entries out of ',
        nbefore, '(', 100*round((nbefore - nafter)/nbefore, 4), '%) ',
        'from the raw results table due to the absence of gene_id')

nbefore <- nrow(rawResults)
rawResults <- rawResults[!is.na(gene_name)]
nafter <- nrow(rawResults)
message('[', Sys.time(), '] Removed ', nbefore - nafter, ' entries out of ',
        nbefore, '(', 100*round((nbefore - nafter)/nbefore, 4), '%) ',
        'from the raw results table due to the absence of gene_name')

# [READ] in known cancer genes ------------------------------------------------
known_cancer <- fread(args$known_cancer_genes, header = T, 
                      stringsAsFactors = F)

# Add known cancer gene status ------------------------------------------------
if (!is.null(args$known_db_to_use)) {
  genes_in_db <- sapply(known_cancer$db_name, 
                        function(x) any(unlist(strsplit(x, ', ')) %in% 
                                          args$known_db_to_use))
  known_cancer_selected <- known_cancer[genes_in_db]
  
  rawResults[, is_known_cancer := gene_id %in% known_cancer_selected$gene_id]
  rawResults[, is_driver_lung := gene_id %in% 
               known_cancer_selected[is_driver_lung == T]$gene_id]
  rm(known_cancer_selected)
} else {
  rawResults[, is_known_cancer := gene_id %in% known_cancer$gene_id]
  rawResults[, is_driver_lung := gene_id %in% 
               known_cancer[is_driver_lung == T]$gene_id]
}
# add column in_db indicating in which data bases gene was already listed
rawResults <- merge(rawResults, unique(known_cancer[,.(gene_id, db_name)]),
                    by = 'gene_id', all.x = T)
setnames(rawResults, 'db_name', 'in_db')

# [ANNOTATE] with CDS length --------------------------------------------------
message('[', Sys.time(), '] Found following GTF files: ', 
        paste(gtfPaths, collapse = ', '), ' will use them to extract gene ',
        'lengths and compute quantile lengths')
gtf <- gtfToTxDB(gtfPaths, doLiftOver = F, chainObj = NULL)
gtf$geneIdToName <- gtf$geneMap[, intersect(c('gene_name', 'gene_id', 
                                              'gene_biotype'),
                                            colnames(gtf$geneMap)), with = F]
gtf$geneIdToName <- gtf$geneIdToName[!duplicated(gtf$geneIdToName)]

cdsLength <- as.data.table(transcriptLengths(gtf$gtfTxdb, with.cds_len = T))
# count CDS length per gene as the longest CDS length (basically the longest 
# transcript )
cdsLength[, maxCDSlen := max(cds_len), by = gene_id]
cdsLength <- cdsLength[maxCDSlen == cds_len & cds_len != 0]
cdsLength <- cdsLength[,.(gene_id, cds_len)]
cdsLength <- cdsLength[!duplicated(cdsLength)]
cdsLength <- merge(cdsLength, 
                   unique(gtf$geneIdToName[,.(gene_id, gene_name)]),
                   by = 'gene_id', all.x = T)

quantBreaks <- getUniqQuantileBreaks(cdsLength$cds_len)
cdsLength[, cdsLenQuant := cut(cds_len, breaks = c(-1, quantBreaks$topBound),
                               labels = quantBreaks$quant)]
cdsLength[, cdsLenQuant := as.character(cdsLength$cdsLenQuant)]
cdsLength[, cdsLenQuant := as.numeric(cdsLength$cdsLenQuant)]
cdsLength[which(cdsLength$cds_len == 0)]$cdsLenQuant <- NA
setnames(cdsLength, 'cds_len', 'cdsLen')

rawResults <- merge(rawResults, cdsLength[,.(gene_id, cdsLen, cdsLenQuant)], 
                    by = 'gene_id', all.x = T)
rm(cdsLength, gtf)

# [ANNOTATE] with olfactory gene status ---------------------------------------
if (!is.null(args$olfactory_genes)) {
  olfactory <- fread(args$olfactory_genes, header = T, sep = '\t', 
                     stringsAsFactors = F, select = c('gene_id', 'gene_name'))
  # choose, which column to use for merging
  nOvrlName <- sum(olfactory$gene_name %in% rawResults$gene_name)
  nOvrlId <- sum(olfactory$gene_id %in% rawResults$gene_id)
  mergeKey <- ifelse(nOvrlId >= nOvrlName, 'gene_id', 'gene_name')
  olfactory <- unlist(olfactory[, c(mergeKey), with = F])
  rawResults[, isOlfact := unlist(rawResults[, mergeKey, with = F]) %in%
               olfactory]
}

# [READ] mutations to genomic regions -----------------------------------------
if (!is.null(args$muts_to_gr_dir)) {
  varsToGRfiles <- list.files(args$muts_to_gr_dir, full.names = T,
                              pattern = args$muts_to_gr_ptrn)
  
  message('[', Sys.time(), '] Started reading mutation to genome region maps')
  # do in a loop for better memory efficiency. We'll read each file and check
  # if tumor_subtype in the file is the same as one of rawResults one
  varsToGRmap <- data.table()
  acceptedTumorTypes <- unique(rawResults$tumor_subtype)
  for (i in 1:length(varsToGRfiles)) {
    message('[', Sys.time(), '] Processing ', varsToGRfiles[i])
    mutGRmap_i <- fread(varsToGRfiles[i], header = T, stringsAsFactors = F)
    if (!'tumor_subtype' %in% colnames(mutGRmap_i)) {
      stop('[', Sys.time(), '] tumor_subtype column is not in ',  
           varsToGRfiles[i])
    }
    varsToGRmap <- rbind(varsToGRmap, 
                         mutGRmap_i[tumor_subtype %in% acceptedTumorTypes])
    rm(mutGRmap_i)
  }
  # add pos (position) to varsToGRmap so that we can count number of unique 
  # positions
  varsToGRmap[, pos := gsub(':[ATGC].*', '', key)]
  varsToGRmap[, pos := gsub('.*:', '', pos)]
  message('[', Sys.time(), '] Finished reading mutation to genome region maps')
  
  if (nrow(varsToGRmap) == 0){
    stop('[', Sys.time(), '] did not read any entries for mean mutation ',
         'rates')
  }
  
  # PATCH. Need to insert into arguments if silent should be removed
  varsToGRmap <- varsToGRmap[gr_id != 'CDS' | 
                               (gr_id == 'CDS' & var_class != 'Silent')]
  message('[', Sys.time(), '] removed silent mutations from consideration ',
          'for CDS regions')
}

# [ANNOTATE] with number of patients, mutations, loc. mut. rate to results-----
if (!is.null(args$muts_to_gr_dir)) {
  # count number of mutations (nMuts), unique mutation sites (nMutsUniq) and 
  # number of participants with mutations per genomic region PER PATIENT TUMOR 
  # SUBTYPE. Counting per patient tumor subtype is done to preserve information
  # regarding mutations distribution in the composite (i.e. created out of 
  # several independent tumor types) tumor_subtype. We can count number of 
  # mutations as .N and not as length(unique(key)) because mutations can occur  
  # at the same spot in different patients. 
  mutStatsPerGR <- varsToGRmap[,.(nMuts = .N, 
                                  nMutsUniq = length(unique(pos)),
                                  nParts = length(unique(participant_id))),
                               by = .(gr_id, gene_id, gene_name, tumor_subtype,
                                      patient_tumor_subtype)]
  mutStatsPerGR <- mutStatsPerGR[!duplicated(mutStatsPerGR)]
  
  # nParts_total will hold a total number of participants with mutation in that
  # region regardless of patient_tumor_subtype
  mutStatsPerGR_total <- varsToGRmap[,.(nMuts_total = .N, 
                                        nMutsUniq_total = length(unique(pos)),
                                        nParts_total = length(unique(participant_id))),
                                     by = .(gr_id, gene_id, gene_name, 
                                            tumor_subtype)]
  mutStatsPerGR_total <- mutStatsPerGR_total[!duplicated(mutStatsPerGR_total)]
  mutStatsPerGR <- merge(mutStatsPerGR, mutStatsPerGR_total, all = T,
                         by = c('gr_id', 'gene_id', 'gene_name', 
                                'tumor_subtype'))
  rm(mutStatsPerGR_total)
  
  # add cohort size
  mutStatsPerGR <- merge(mutStatsPerGR, cohortSizeDT, by = 'tumor_subtype', 
                         all.x = T)
  
  # add info about loc. mut. rate
  locMutRateDT <- varsToGRmap[, c('gr_id', 'gene_id', 'gene_name', 
                                  'tumor_subtype', 
                                  grep('meanMutRate|MeanMutRate', 
                                       colnames(varsToGRmap), value = T)),
                              with = F]
  locMutRateDT <- locMutRateDT[!duplicated(locMutRateDT)]
  mutStatsPerGR <- merge(mutStatsPerGR, locMutRateDT,
                         by = c('gr_id', 'gene_id', 'gene_name', 
                                'tumor_subtype'), all.x = T)
  rm(locMutRateDT)
  
  rawResults <- merge(rawResults, mutStatsPerGR, all.x = T, 
                      allow.cartesian = T,
                      by = c('gr_id', 'gene_id','gene_name', 'tumor_subtype'))
  
  
  mismatched <- rawResults[is.na(patient_tumor_subtype) & raw_p < 0.05]
  mismatched <- unique(mismatched[,.(tumor_subtype, gr_id, gene_name, 
                                     software)])
  message('[', Sys.time(), '] Found ', nrow(mismatched), '(',
          round(100 * nrow(mismatched) / 
                  nrow(unique(rawResults[,.(tumor_subtype, gr_id, gene_name, 
                                            software)])), 4), '%) ', 
          'tumor subtype-gene - genomic region - software combinations where ',
          'raw p value is significant, but no mutations were assigned to ',
          'that region. This can be cause by slight misalignment of ',
          'annotations between softwares and genomic ranges. dNdScv and ',
          'DIGdriver are prone to that due to RefCDS creation and inability ',
          'to remove positions affected by low mappability, etc. MutPanning ',
          'is prone to that due do genome annotation being fixed inside ',
          'the software. Overview of distribution across software: ',
          paste0(capture.output(knitr::kable(mismatched[,.N, by = software], 
                                             format = "markdown")), 
                 collapse = '\n'), 
          '\n\nOverview of distribution across genomic regions: ',
          paste0(capture.output(knitr::kable(mismatched[,.N, 
                                                        by = gr_id][order(-N)], 
                                             format = "markdown")), 
                 collapse = '\n'))
  
  rm(mutStatsPerGR, varsToGRmap)
}

# [ANNOTATE] with GTEx expression status (True/False/NA) ----------------------
if (!is.null(args$gtex)) {
  message('[', Sys.time(), '] Started reading GTEx expression table')
  gtex <- fread(args$gtex, header = T, stringsAsFactors = F)
  message('[', Sys.time(), '] Finished reading GTEx expression table')
  
  rawResults <- split(rawResults, by = 'tumor_subtype')
  for (tumorSubtype in names(rawResults)) {
    if (tumorSubtype %in% colnames(gtex)) {
      message('[', Sys.time(), '] Found GTEx expression data for ', 
              tumorSubtype, ', annotating.')
      expressedInGTEx <- unlist(gtex[, tumorSubtype, with = F])
      expressedInGTEx <- gtex[expressedInGTEx != 0][,.(gene_name, gene_id)]
      rawResults[[tumorSubtype]][, exprInGTEx := F]
      rawResults[[tumorSubtype]][!gene_name %in% gtex$gene_name & 
                                   !gene_id %in% gtex$gene_id]$exprInGTEx <- NA
      rawResults[[tumorSubtype]][gene_name %in% 
                                   expressedInGTEx$gene_name]$exprInGTEx <- T
      rawResults[[tumorSubtype]][gene_id %in% 
                                   expressedInGTEx$gene_id]$exprInGTEx <- T
    } else {
      message('[', Sys.time(), '] Did not find GTEx expression data for ', 
              tumorSubtype)
      rawResults[[tumorSubtype]][, exprInGTEx := as.logical(NA)]
    }
  }
  rawResults <- do.call(rbind, rawResults)
}

# [ANNOTATE] with TCGA expression status (True/False/NA) ----------------------
if (!is.null(args$tcga)) {
  message('[', Sys.time(), '] Started reading TCGA expression table')
  tcga <- fread(args$tcga, header = T, stringsAsFactors = F)
  message('[', Sys.time(), '] Started reading TCGA expression table')
  
  rawResults <- split(rawResults, by = 'tumor_subtype')
  for (tumorSubtype in names(rawResults)) {
    if (tumorSubtype %in% colnames(tcga)) {
      message('[', Sys.time(), '] Found TCGA expression data for ', 
              tumorSubtype, ', annotating.')
      expressedInTCGA <- unlist(tcga[, tumorSubtype, with = F])
      expressedInTCGA <- tcga[expressedInTCGA != 0][,.(gene_name, gene_id)]
      rawResults[[tumorSubtype]][, exprInTCGA := F]
      rawResults[[tumorSubtype]][!gene_name %in% tcga$gene_name & 
                                   !gene_id %in% tcga$gene_id]$exprInTCGA <- NA
      rawResults[[tumorSubtype]][gene_name %in% 
                                   expressedInTCGA$gene_name]$exprInTCGA <- T
      rawResults[[tumorSubtype]][gene_id %in% 
                                   expressedInTCGA$gene_id]$exprInTCGA <- T
    } else {
      message('[', Sys.time(), '] Did not find GTEx expression data for ', 
              tumorSubtype)
      rawResults[[tumorSubtype]][, exprInTCGA := as.logical(NA)]
    }
  }
  rawResults <- do.call(rbind, rawResults)
}
# [SAVE] Raw results of software run as table ---------------------------------
write.table(rawResults, args$output, append = F, quote = F,  sep = '\t', 
            row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote raw drivers table to ', args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')