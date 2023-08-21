#!/usr/bin/env Rscript
# FILE: 0_functions.R ---------------------------------------------------------
#
# DESCRIPTION: All custom created R functions to perform de-novo cancer driver 
#              genetic element discovery.
#             
# USAGE: In R script, insert:
#        source('0_functions.R')
#
# OPTIONS: None
#
# REQUIREMENTS: R v4.1.0
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  22.06.2023
# REVISION: 05.07.2023

box::use(argparse[...])
box::use(data.table[...])
box::use(GenomeInfoDb[...])
box::use(GenomicRanges[...])
box::use(rtracklayer[...])


# GLOBAL ARGUMENTS ------------------------------------------------------------
#' @export
GR_CODES <- c('protein_coding', "3primeUTR", "5primeUTR", "CDS", "lincRNA", 
              "lincRNA_promoter", "lincRNA_ss", "miRNA", "misc_RNA", 
              "promoter", "rRNA", "snoRNA", "snRNA", "ss")

#' @export
SOFTWARE_GR_CODES <- list('activedriverwgs' = NULL, 'dndscv' = c('CDS'),
                          'mutpanning' = c('CDS'), 'chasmplus' = c('CDS'), 
                          'driverpower' = NULL, 'nbr' = NULL, 
                          'oncodrivefml' = NULL, 'oncodriveclustl' = NULL,
                          'digdriver' = NULL)
#' @export
acceptedChrNames <- c(c(1:24, 'X', 'Y'), paste0('chr', c(1:24, 'X', 'Y')))

# Parcing input arguments -----------------------------------------------------
#' check_file_existence_msg
#' @description Checks, whatever or not file exist and stops the execution if
#'              it is not. Prints an error.
#' @author Maria Litovchenko
#' @param filePath path to file
#' @return void
#' @export
check_file_existence_msg <- function(filePath) {
  if (!file.exists(filePath)) {
    stop('[', Sys.time(), '] File does not exist: ', filePath)
  }
}

#' check_min_int_msg
#' @description Checks, whatever numerical argument submitted via command line 
#'              to the script is not less than certain value. Stops and prints
#'              an error if not
#' @author Maria Litovchenko
#' @param value a value of submitted argument
#' @param value_name a name of submitted argument
#' @param min_value minimal value of submitted argument
#' @return void
#' @export
check_min_int_msg <- function(value, value_name, min_value) {
  if (value < min_value) {
    stop('[', Sys.time(), '] Minimal value for the argument --', value_name, 
         ' is ', min_value)
  }
}

#' check_input_arguments
#' @description Checks validity of submitted command line arguments. Terminates
#'              script if something is wrong
#' @author Maria Litovchenko
#' @param argsList list of submitted command line arguments
#' @param outputType string, type of output argument: folder or file
#' @return void
#' @export
check_input_arguments <- function(argsList, outputType = NULL) {
  if (!is.null(argsList$inventory_patients)) {
    lapply(argsList$inventory_patients, check_file_existence_msg)
  }
  if (!is.null(argsList$inventory_analysis)) {
    lapply(argsList$inventory_analysis, check_file_existence_msg)
  }
  if (!is.null(argsList$inventory_blacklisted)) {
    lapply(argsList$inventory_blacklisted, check_file_existence_msg)
  }
  if (!is.null(argsList$min_n_participants)) {
    check_min_int_msg(argsList$min_n_participants, 'min_n_participants', 1)
  }
  if (!is.null(argsList$min_depth)) {
    check_min_int_msg(argsList$min_depth, 'min_depth', 1)
  }
  if (!is.null(argsList$min_tumor_vac)) {
    check_min_int_msg(argsList$min_tumor_vac, 'min_tumor_vac', 1)
  }
  if (!is.null(argsList$min_tumor_vaf)) {
    check_min_int_msg(argsList$min_tumor_vaf, 'min_tumor_vaf', 0)
  }
  if (!is.null(argsList$max_germline_vaf)) {
    check_min_int_msg(argsList$max_germline_vaf, 'max_germline_vaf', 0)
  }
  if (!is.null(argsList$max_germline_vac)) {
    check_min_int_msg(argsList$max_germline_vac, 'max_germline_vac', 0)
  }
  if (!is.null(argsList$max_n_vars)) {
    check_min_int_msg(argsList$max_n_vars, 'max_n_vars', 0)
  }
  if (!is.null(argsList$min_reg_len)) {
    check_min_int_msg(argsList$min_reg_len, 'min_reg_len', 0)
  }
  if (!is.null(argsList$target_genome_path)) {
    check_file_existence_msg(argsList$target_genome_path)
  }
  if (!is.null(argsList$chain)) {
    check_file_existence_msg(argsList$chain)
  }
  
  if (!is.null(argsList$output)) {
    if (!is.null(outputType)) {
      dirToCheck <- argsList$output
      if (outputType == 'file') {
        dirToCheck <- dirname(argsList$output)
      }
      if (!dir.exists(dirToCheck)) {
        dir.create(dirToCheck, recursive = T)
        message('[', Sys.time(), '] Created directory ', dirToCheck) 
      }
    } else {
      stop('[', Sys.time(), '] Please give outputType argument.')
    }
  }
  
  if (!is.null(argsList$cores)) {
    check_min_int_msg(argsList$cores, 'cores', 1)
  }
}

#' printArgs
#' @description Prints submitted to the script arguments as messages.
#' @author Maria Litovchenko
#' @param argsList named list representing submitted to the script arguments
#' @return void
#' @export
printArgs <- function(argsList) {
  message('[', Sys.time(), '] Submitted arguments: ')
  for (argName in names(argsList)) {
    oneArg <- argsList[[argName]]
    if (length(oneArg) > 1 & !is.null(names(oneArg))) {
      msg <- paste(apply(data.table(names(oneArg), oneArg), 1, paste, 
                         collapse = ' - '), collapse = ',')
    } else {
      msg <- paste(oneArg, collapse = ', ')
    }
    message('\t', argName, ':', msg)
  }
}

# Reading inventories ---------------------------------------------------------
#' readParticipantInventory
#' @description Reads participant inventory into data table
#' @author Maria Litovchenko
#' @param inventoryPath path to file with participants inventory
#' @param cores integer, number of cores to use
#' @return data table with the inventory
#' @export
readParticipantInventory <- function(inventoryPath, cores = 1) {
  result <- fread(inventoryPath, sep = ',', header = T, stringsAsFactors = F, 
                  select = c('tumor_subtype', 'participant_id', 
                             'participant_tumor_subtype', 'somatic_path', 
                             'somatic_genome', 'cohort_name'),
                  nThread = cores)
  result[, tumor_subtype := as.character(tumor_subtype)]
  result[, participant_tumor_subtype := as.character(participant_tumor_subtype)]
  result[, participant_id := as.character(participant_id)]
  result[, cohort_name := as.character(cohort_name)]
  
  result
}

#' readAnalysisInventory
#' @description Reads analysis inventory into data table
#' @author Maria Litovchenko
#' @param inventoryPath path to file with analysis inventory
#' @param cores integer, number of cores to use
#' @return data table with the inventory
#' @export
readAnalysisInventory <- function(inventoryPath, cores = 1) {
  essenCols <- c('tumor_subtype', 'software', 'gr_id', 'gr_code', 'gr_file',
                 'gr_upstr', 'gr_downstr', 'gr_genome', 'gr_excl_id',
                 'gr_excl_code', 'gr_excl_file', 'gr_excl_upstr',
                 'gr_excl_downstr', 'gr_excl_genome', 'blacklisted_codes')
  result <- fread(inventoryPath, sep = ',', header = T, stringsAsFactors = F,
                  select = essenCols, nThread = cores)
  
  result[, tumor_subtype := as.character(tumor_subtype)]
  result[, gr_id := as.character(gr_id)]
  result[, gr_code := as.character(gr_code)]
  result[, gr_upstr := as.integer(gr_upstr)]
  result[, gr_downstr := as.integer(gr_downstr)]
  result[, gr_genome := as.character(gr_genome)]
  result[, gr_excl_id := as.character(gr_excl_id)]
  result[, gr_excl_code := as.character(gr_excl_code)]
  result[, gr_excl_upstr := as.integer(gr_excl_upstr)]
  result[, gr_excl_downstr := as.integer(gr_excl_downstr)]
  result[, gr_excl_genome := as.character(gr_excl_genome)]
  result[, blacklisted_codes := as.character(blacklisted_codes)]
  blackListsDT <- data.table(blacklisted_codes = unique(result$blacklisted_codes))
  blackListsDT <- blackListsDT[!is.na(blacklisted_codes)]
  if (nrow(blackListsDT) > 0) {
    blackListsDT[, blacklisted_codes_parsed := strsplit(blacklisted_codes,';')]
    colsOrder <- colnames(result)
    result <- merge(result, blackListsDT, by = 'blacklisted_codes', all.x = T)
    result[, blacklisted_codes := NULL]
    setnames(result, 'blacklisted_codes_parsed', 'blacklisted_codes')
    setcolorder(result, colsOrder)
    rm(colsOrder, blackListsDT)
  }
  
  result
}

#' readBlacklistInventory
#' @description Reads black&white lists inventory into data table
#' @author Maria Litovchenko
#' @param inventoryPath black&white lists inventory analysis path
#' @param cores integer, number of cores to use
#' @return data table with black&white lists inventory
#' @export
readBlacklistInventory <- function(inventoryPath, cores = 1) {
  # 1) read and check that all needed columns are present
  essenCols <- c('list_name', 'file_path', 'file_genome', 'file_type', 
                 'score_column', 'min_value', 'max_value')
  result <- fread(inventoryPath, sep = ',', header = T, stringsAsFactors = F,
                  select = essenCols, nThread = cores)
  
  result[, list_name := as.character(list_name)]
  result[, file_genome := as.character(file_genome)]
  result[, file_type := as.character(file_type)]
  result[, score_column := as.character(score_column)]
  result
}

# Reading black&white regions lists -------------------------------------------
#' readMappabilityTracks
#' @description  Reads mappability tracks to data table. Accepts only bigwig
#' (bw), bed and gtf files. Bed files can be gziped.
#' @author Maria Litovchenko
#' @param filePath path to the file.
#' @param cores number of cores to use (used in fread only)
#' @param roiGR GRanges, regions of interest. If given (not NULL), only regions
#' of interest will be read from the file.
#' @return data table.
#' @export
readMappabilityTracks <- function(filePath, cores = 1, roiGR = NULL) {
  fileExt <- gsub('.gz$', '', tolower(filePath))
  fileExt <- gsub('.*[.]', '', fileExt)
  
  if (!grepl('bed|bigwig|bw|gtf', fileExt)) {
    msg <- paste('Unrecognized mappability file format:', filePath,
                 'only bigwig (bw), bed and gtf formats are acceptable.',
                 'Please provide files with bw or bed in file extension.')
    stop(msg)
  }
  
  if (is.null(roiGR)) {
    mapTrack <- import(filePath)
  } else {
    mapTrack <- import(filePath, which = roiGR)
  }
  
  if (fileExt == 'gtf') {
    message('[', Sys.time(), '] GTF file as black-/white- list file was ',
            'submitted. All entries in the file will be used.')
    mapTrack <- disjoin(reduce(mapTrack))
  }
  
  mapTrack <- as.data.table(mapTrack)
  setnames(mapTrack, 'seqnames', 'chr')
  
  mapTrack
}

#' filterBWregions
#' @description Filters black- or white- listed regions according to the 
#' requested range of values in bwScoreCol column.
#' @author Maria Litovchenko
#' @param bwGRs GRanges object representation of black and white regions
#' @param bwScoreCol string, name of the column containing scores on which 
#' regions should be filtred.
#' @param bwScoreMin numeric, minimum value of a score
#' @param bwScoreMax numeric, maximum value of a score
#' @return filtered Granges
#' @export
filterBWregions <- function(bwGRs, bwScoreCol = NA, bwScoreMin = NA, 
                            bwScoreMax = NA) {
  if (!is.na(bwScoreCol)) {
    result <- bwGRs[as.vector(bwGRs[, bwScoreCol, with = F] >= bwScoreMin), ]
    result <- result[as.vector(result[, bwScoreCol, with = F] <= bwScoreMax), ]
  } else {
    result <- copy(bwGRs)
  }
  result
}

#' readInAndFilterBWregions
#' @description Reads in and filters, if required, file with black- or white-
#' listed genomic regions.
#' @author Maria Litovchenko
#' @param bwFile path to file with black- or white- listed genomic regions.
#' @param chrStyle character, one of NCBI or UCSC which determine chromosome
#' naming style (1 or chr1 respectively). Final result will have this 
#' chromosome naming style.
#' @param bwScoreCol string, name of the column containing scores on which 
#' regions should be filtred.
#' @param bwScoreMin numeric, minimum value of a score
#' @param bwScoreMax numeric, maximum value of a score
#' @param roiGR GRanges, regions of interest. If given (not NULL), only regions
#' of interest will be read from the file.
#' @return GRanges object
#' @export
readInAndFilterBWregions <- function(bwFile, chrStyle, bwScoreCol = NA, 
                                     bwScoreMin = NA, bwScoreMax = NA, 
                                     cores = 1, roiGR = NULL) {
  if (chrStyle != 'NCBI' & chrStyle != 'UCSC') {
    stop('[', Sys.time(), '] chrStyle should be one of NCBI or UCSC.')
  }
  
  message('[', Sys.time(), '] Started reading ', bwFile)
  if (!is.null(roiGR)) {
    nbefore <- 0
    result <- data.table()
    for (chrID in sort(unique(as.character(seqnames(roiGR))))) {
      message('[', Sys.time(), ']\t processing chromosome ', chrID)
      oneChr <- readMappabilityTracks(bwFile, cores, 
                                      roiGR[seqnames(roiGR) == chrID])
      nbefore <- nbefore + length(oneChr)
      result <- rbind(result, filterBWregions(oneChr, bwScoreCol, 
                                              bwScoreMin, bwScoreMax))
    }
  } else {
    result <- readMappabilityTracks(bwFile, cores)
    nbefore <- length(result)
    result <- filterBWregions(result, bwScoreCol, bwScoreMin, bwScoreMax)
  }
  result <- as.data.table(result)
  if (!is.na(bwScoreCol)) {
    nafter <- length(result)
    message('[', Sys.time(), '] Removed ', nbefore - nafter, '(',
            round(100 * (nbefore - nafter) / nbefore, 2), '%) entries from ', 
            bwFile, ' due to not passing cut offs on ', bwScoreCol)
    
  }
  message('[', Sys.time(), '] Finished reading ', bwFile)
  
  result <- makeGRangesFromDataFrame(result)
  seqlevelsStyle(result) <- chrStyle
  
  result
}

# Misc -----------------------------------------------------------------------
#' getSeqlevelsStyle
#' @description Detemines a seqlevelstyle (aka chromosome naming style) from
#' the reference genome
#' @author Maria Litovchenko
#' @param fastaPath path to fasta file
#' @return UCSC, in case naming format is chr1, and NCBI otherwise
#' @export
getSeqlevelsStyle <- function(fastaPath) {
  # read just the first line
  result <- fread(fastaPath, nrows = 1, header = F)$V1
  result <- ifelse(grepl('>chr', result), 'UCSC', 'NCBI')
  result
}