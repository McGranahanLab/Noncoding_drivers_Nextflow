#!/usr/bin/env Rscript
# FILE: custom_functions_preprocessing.R --------------------------------------
#
# DESCRIPTION: Custom created R functions to perform data preprocessing before
#              de-novo cancer driver genetic element discovery.
#             
# USAGE: In R script, insert:
#        source('custom_functions_preprocessing.R')
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
# REVISION: 07.11.2023

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
suppressWarnings(suppressPackageStartupMessages(library(GenomeInfoDb)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(stats)))
suppressWarnings(suppressPackageStartupMessages(library(utils)))

# Parsing input arguments -----------------------------------------------------
#' check_input_arguments_preproc
#' @description Checks validity of submitted command line arguments. Terminates
#'              script if something is wrong
#' @author Maria Litovchenko
#' @param argsList list of submitted command line arguments
#' @param outputType string, type of output argument: folder or file
#' @return void
#' @export
check_input_arguments_preproc <- function(argsList, outputType = NULL) {
  lapply(argsList$inventory_patients, check_file_existence_msg)
  lapply(argsList$inventory_analysis, check_file_existence_msg)
  lapply(argsList$inventory_blacklisted, check_file_existence_msg)
  check_min_int_msg(argsList$min_n_participants, 'min_n_participants', 1)
  check_min_int_msg(argsList$min_depth, 'min_depth', 1)
  check_min_int_msg(argsList$min_tumor_vac, 'min_tumor_vac', 1)
  check_min_int_msg(argsList$min_tumor_vaf, 'min_tumor_vaf', 0)
  check_min_int_msg(argsList$max_germline_vaf, 'max_germline_vaf', 0)
  check_min_int_msg(argsList$max_germline_vac, 'max_germline_vac', 0)
  check_min_int_msg(argsList$max_n_vars, 'max_n_vars', 0)
  check_min_int_msg(argsList$min_reg_len, 'min_reg_len', 1)
  check_file_existence_msg(argsList$target_genome_path)
  check_file_existence_msg(argsList$chain)
  
  check_file_existence_msg(argsList$variants)
  check_file_existence_msg(argsList$genomic_regions)
  check_file_existence_msg(argsList$target_genome_chr_len)
  check_file_existence_msg(argsList$gene_name_synonyms)
  check_file_existence_msg(argsList$varanno_conversion_table)
  
  check_file_existence_msg(argsList$maf)
  check_file_existence_msg(argsList$bed)
  
  check_min_int_msg(argsList$cores, 'cores', 1)
  
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

# Gene names synonyms ---------------------------------------------------------
#' getGeneSymbolSynonyms
#' @description Retrieves from the table of gene names synonyms all synonumous 
#'              names for a given gene.
#' @author Maria Litovchenko
#' @param geneSymbol vector of gene symbols
#' @param geneNameSynsDT data table with information about gene name synonyms.
#'                       Should have columns idx and gene_name. Synonymous gene
#'                       names will have the same idx.
#' @return data.table with columns gene_name_orig, containing submitted as 
#'         input geneSymbol and gene_name_syn - retrieved synonymous names.
#' @export
getGeneSymbolSynonyms <- function(geneSymbols, geneNameSynsDT) {
  geneSymbols <- unique(geneSymbols)
  
  result <- data.table(gene_name_orig = geneSymbols, 
                       gene_name_syn = as.character(NA))
  
  if (any(geneSymbols %in% geneNameSynsDT$gene_name)) {
    setkey(geneNameSynsDT, 'gene_name')
    result <- geneNameSynsDT[idx %in% geneNameSynsDT[geneSymbols]$idx]
    setkey(result, 'gene_name')
    foundSyn <- intersect(geneSymbols, result$gene_name)
    # unfortunately, geneSymbols can contain synonymous values
    result <- lapply(foundSyn, 
                     function(x) cbind(gene_name_orig = x,
                                       result[idx %in% result[x]$idx]))
    result <- do.call(rbind, result)
    result <- result[gene_name_orig != gene_name]
    setnames(result, 'gene_name', 'gene_name_syn')
    result <- result[,.(gene_name_orig, gene_name_syn)]
    
    foundSyn <- length(intersect(geneSymbols, result$gene_name_orig))
    message('[', Sys.time(), '] Found synonymous gene names for ',
            foundSyn, '(', round(100 * foundSyn/length(geneSymbols)), 
            '%) genes.')
  }
  result
}

# Bining into quantiles -------------------------------------------------------
#' getUniqQuantileBreaks
#' @description Creates data table containing unique breaks and their labels
#' based on data vector
#' @param x numeric vector
#' @param nQuants number of quantiles
#' @return data table with columns topBound (use it as breaks) and quant (use 
#' it as labels)
#' @export
getUniqQuantileBreaks <- function(x, nQuants = 1001) {
  quantCuts <- data.table(quant = seq(0, 1, length.out = nQuants), 
                          topBound = quantile(x, seq(0, 1,
                                                     length.out = nQuants)))
  quantCuts <- quantCuts[,.(quant = max(quant)), by = topBound]
  quantCuts[topBound == max(topBound)]$topBound <- Inf
  quantCuts[, quant := round(100 * quant, 1)]
  quantCuts
}

#' assignQuantileBreaks
#' @description Assigns quantile in which a numeric value falls in the 
#'              distribution
#' @author Maria Litovchenko
#' @param inDT input data table
#' @param colName name of the column containing numerical values. Cells with 0
#'        will not be taken into account
#' @return input data table with extra column named colName+Quant containing 
#'         quantile.
#' @export
assignQuantileBreaks <- function(inDT, colName) {
  setnames(inDT, colName, 'coi')
  
  quantBreaks <- getUniqQuantileBreaks(inDT[coi != 0]$coi)
  result <- copy(inDT)
  result[, coiQuant := cut(coi, breaks = c(-1, quantBreaks$topBound),
                           labels = quantBreaks$quant)]
  result[, coiQuant := as.numeric(as.character(coiQuant))]
  result[coi == 0]$coiQuant <- NA
  
  setnames(result, c('coi', 'coiQuant'), c(colName, paste0(colName, 'Quant')))
  result
}

# Misc -----------------------------------------------------------------------
#' getSeqlevelsStyle
#' @description Determines a seqlevelstyle (aka chromosome naming style) from
#' the reference genome
#' @author Maria Litovchenko
#' @param fastaPath path to fasta file
#' @return UCSC, in case naming format is chr1, and NCBI otherwise
#' @export
getSeqlevelsStyle <- function(fastaPath) {
  # read just the first line
  result <- readLines(fastaPath, n = 1)
  result <- ifelse(grepl('>chr', result), 'UCSC', 'NCBI')
  result
}

#' orderChromosomes
#' @description Orders chromosomal names
#' @author Maria Litovchenko
#' @param chrs vector of chromosomes
#' @return vector of unique chromosomal names, ordered (first numeric 
#' chromosomes, i.e. 1-22, and then 'character' chromosomes like X and Y)
#' @export
orderChromosomes <- function(chrs) {
  result <- data.table(chr = unique(chrs))
  result[, chrInt := gsub('^chr', '', chr)]
  suppressWarnings(result[, chrInt := as.integer(chrInt)])
  result[chr == 'X']$chrInt <- 23
  result[chr == 'Y']$chrInt <- 24
  result <- result[order(chrInt, chr)]
  result$chr
}