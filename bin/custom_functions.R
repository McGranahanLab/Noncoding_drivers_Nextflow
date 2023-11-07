#!/usr/bin/env Rscript
# FILE: custom_functions.R ----------------------------------------------------
#
# DESCRIPTION: Custom created R functions which are used in both data 
#              preprocessing before de-novo cancer driver genetic element 
#              discovery and de-novo discovery post-processing.
#             
# USAGE: In R script, insert:
#        source('custom_functions.R')
#
# OPTIONS: None
#
# REQUIREMENTS: R v4.1.0
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  07.11.2023
# REVISION: 07.11.2023

box::use(argparse[...])
box::use(data.table[...])
box::use(GenomeInfoDb[...])
box::use(GenomicRanges[...])
box::use(rtracklayer[...])
box::use(utils[...])

# GLOBAL ARGUMENTS ------------------------------------------------------------
#' @export
GR_CODES <- c("promoter", "5primeUTR", "CDS", "ss", "3primeUTR", 
              "lincRNA_promoter", "lincRNA", "lincRNA_ss",
              "miRNA", "misc_RNA", "rRNA", "snoRNA", "snRNA")

#' @export
SOFTWARE_GR_CODES <- list('activedriverwgs' = NULL, 'dndscv' = c('CDS'),
                          'mutpanning' = c('CDS'), 'chasmplus' = c('CDS'), 
                          'driverpower' = NULL, 'nbr' = NULL, 
                          'oncodrivefml' = NULL, 'oncodriveclustl' = NULL,
                          'digdriver' = NULL)
#' @export
acceptedChrNames <- c(c(1:24, 'X', 'Y'), paste0('chr', c(1:24, 'X', 'Y')))

# Parsing input arguments -----------------------------------------------------
#' check_file_existence_msg
#' @description Checks, whatever or not file exist and stops the execution if
#'              it is not. Prints an error.
#' @author Maria Litovchenko
#' @param filePath path to file
#' @return void
#' @export
check_file_existence_msg <- function(filePath) {
  if (!is.null(filePath)) {
    if (!file.exists(filePath)) {
      stop('[', Sys.time(), '] File does not exist: ', filePath)
    }
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
  if (!is.null(value)) {
    if (value < min_value) {
      stop('[', Sys.time(), '] Minimal value for the argument --', value_name, 
           ' is ', min_value)
    }
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
                 'gr_excl_downstr', 'gr_excl_genome', 'blacklisted_codes',
                 'union_percentage', 'intersect_percentage')
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
  if (!'union_percentage' %in% colnames(result)) {
    result[, union_percentage := NA]
  }
  if (!'intersect_percentage' %in% colnames(result)) {
    result[, intersect_percentage := NA]
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
  if (file.size(inventoryPath) == 0) {
    return(NULL)
  }
  
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

# Reading bed12 files ---------------------------------------------------------
#' parseBED12regName
#' @description Parces name string(s) from BED12 files containing information
#' about genomic regions of interest. Usually those files are used with 
#' driverpower.
#' @author Maria Litovchenko
#' @param nameStr vector of strings
#' @return data table with number of rows = length of nameStr with columns 
#'         target_genome_version, gr_id, gene_id, gene_name, name
#' @export
parseBED12regName <- function(nameStr, sepStr = '--') {
  result <- lapply(nameStr, strsplit, sepStr)
  result <- do.call(rbind, lapply(result, function(x) x[[1]]))
  result <- as.data.table(result)
  result[, name := nameStr]
  colnames(result) <- c('target_genome_version', 'gr_id', 'gene_id',
                        'gene_name', 'name')
  result
}

#' readBED12
#' @description Reads BED12 file into GRanges
#' @author Maria Litovchenko
#' @param bedPath path to bed12 file
#' @return GRanges with mcols target_genome_version, gr_id, gene_id, gene_name,
#'         gr_name
#' @export
readBED12 <- function(bedPath) {
  message('[', Sys.time(), '] Started reading input genomic ranges file')
  result <- import(bedPath)
  result <- unlist(blocks(result))
  result_parsed <- parseBED12regName(unique(names(result)))
  
  setnames(result_parsed, 'name', 'gr_name')
  setkey(result_parsed, gr_name)
  mcols(result) <- result_parsed[names(result)]
  rm(result_parsed)
  
  suppressMessages(suppressWarnings(gc()))
  
  message('[', Sys.time(), '] Finished reading input genomic ranges file')
  message('[', Sys.time(), '] rtracklayer import function while reading ',
          'bed12 assumes that it is 0-based. Therefore it adds 1 to all ',
          'coordinates. We will correct it.')
  start(result) <- start(result) - 1
  end(result) <- end(result) - 1
  result
}

# Liftover --------------------------------------------------------------------
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
#' @param min.width minimum width of lifted over genomic regions.
#'                  Default = 0 (no regions will be filtered out). 
#' @return lifted over genomic coordinates
#' @export
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
  } else {
    result <- unlist(result)
  }
  
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
  
  message('[', Sys.time(), '] Ratio of total region length after and ',
          'before liftover: ', round(sum(width(result)) / sum(width(inGR)), 3))
  message('[', Sys.time(), "] Ratio of regions' numbers after and before ",
          'liftover: ', round(length(result) / length(inGR), 3))
  
  # set style of chromosomal names to the same as input data
  seqlevelsStyle(result) <- seqlevelsStyle(inGR)
  
  message('[', Sys.time(), '] Finished liftover')
  result
}

#' checkLiftOverByTr
#' @description Checks results of liftover. For every transcript checks that it
#'              was lifted to one chromosome and one strand. If it's not the
#'              case will remove them.
#' @author Maria Litovchenko
#' @param origGR CompressedGRangesList/GRangesList/GRanges from which liftover
#'               was performed.
#' @param loGR GRangesList of lifted regions, result of liftOverGenomicRegion
#'             function
#' @return GRanges
#' @export
checkLiftOverByTr <- function(origGR, loGR) {
  if (class(origGR) == 'CompressedGRangesList') {
    origGRunlist <- unlist(origGR)
  } else {
    origGRunlist <- origGR
  }
  
  # identify transcripts, where there is a non-unique combination of chromosome
  # and strand after liftover
  lo_tr_chrStrand <- data.table(chr = as.character(seqnames(loGR)),
                                strand = as.character(strand(loGR)),
                                transcript_id = loGR$transcript_id)
  lo_tr_chrStrand <- unique(lo_tr_chrStrand)
  weirdLO <- lo_tr_chrStrand[,.N, by = transcript_id]
  weirdLO <- weirdLO[N > 1]
  if (nrow(weirdLO) != 0) {
    weirdLO[, N := NULL]
    
    # add biotype, if available 
    if ('transcript_biotype' %in% colnames(mcols(loGR))) {
      biotype_to_id <- mcols(loGR)[, c('transcript_id', 'transcript_biotype')]
      biotype_to_id <- as.data.table(unique(biotype_to_id))
      weirdLO <- merge(weirdLO, biotype_to_id, by = 'transcript_id', all.x = T)
    } else {
      weirdLO[, transcript_biotype := 'unknown']
    }
    
    # inform user about situation
    msgTab <- weirdLO[,.(length(unique(transcript_id))), 
                      by = 'transcript_biotype']
    setnames(msgTab, 'V1', 'nWeird')
    # ... if biotype is available, calculate number of transcripts in the 
    #     original regions per biotype
    if ('transcript_biotype' %in% colnames(mcols(origGRunlist))) {
      orig_biotype_count <- mcols(origGRunlist)[, c('transcript_id',
                                                    'transcript_biotype')]
      orig_biotype_count <- unique(as.data.table(orig_biotype_count))
      orig_biotype_count <- orig_biotype_count[,.(length(unique(transcript_id))),
                                               by = 'transcript_biotype']
      setnames(orig_biotype_count, 'V1', 'N')
      orig_biotype_count[, Ntotal := length(unique(origGRunlist$transcript_id))]
    } else {
      orig_biotype_count <- data.table(transcript_biotype = 'unknown',
                                       N = length(unique(origGRunlist$transcript_id)))
      orig_biotype_count[, Ntotal := N]
    }
    msgTab <- merge(msgTab, orig_biotype_count, by = 'transcript_biotype',
                    all.x = T)
    msgTab[is.na(msgTab)] <- 0
    msgTab[, nWeirdPerc := round(100 * nWeird/N, 2)]
    msgTab[, nWeirdPercTotal := round(100 * nWeird/Ntotal, 2)]
    msgTab <- msgTab[,.(transcript_biotype, nWeird, 
                        nWeirdPerc, nWeirdPercTotal)]
    msgTab <- msgTab[order(-nWeirdPerc)]
    colnames(msgTab) <- c('transcript_biotype', 'N. trascripts', 
                          '% trascripts, that biotype',
                          '% trascripts, all biotypes')
    message('[', Sys.time(), '] Some transcripts were lifted over to > one ',
            'chromosome and/or > one strand. During creation of TxDB it will ',
            'lead to the transcripts being discarded. Here is an overview of ',
            'numbers:')
    message(paste0(capture.output(knitr::kable(msgTab, format = "markdown")),
                   collapse = '\n'))
    message('[', Sys.time(), '] These transcripts will be removed. ')
    
    result <- loGR[!loGR$transcript_id %in% weirdLO$transcript_id]
  } else {
    message('[', Sys.time(), '] No transcript was lifted over to different ',
            'chromosome/strand')
    result <- copy(loGR)
  }
  
  result
}