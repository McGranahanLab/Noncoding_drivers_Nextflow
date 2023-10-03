#!/usr/bin/env Rscript
# FILE: create_genomic_regions.R ----------------------------------------------
#
# DESCRIPTION: a script for creating input genomic regions for various cancer
#              driver programs.
#
# USAGE: Rscript --vanilla
#
# OPTIONS:
#
# EXAMPLE:
#
#
# REQUIREMENTS:
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  08.10.2020
# REVISION: 29.09.2023


box::use(./custom_functions[...])

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dndscv)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicFeatures)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges)))
suppressWarnings(suppressPackageStartupMessages(library(plyranges)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(strex)))
options(scipen = 999)

# FUNCTIONS: conversion to txdb -----------------------------------------------
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
  result <- result[, intersect(c('gene_name', 'gene_id', 
                                 'gene_type', 'gene_biotype',
                                 'transcript_id', 
                                 'transcript_type', 'transcript_biotype'),
                               colnames(result)), with = F]
  
  # gene_type and gene_biotype were included above in order to access gene 
  # biotype. In some GTF files the column name is gene_type and in the others
  # it is gene_biotype. However, they never should be both in the same file.
  if ('gene_type' %in% colnames(result) & 
      'gene_biotype' %in% colnames(result)) {
    stop('[', Sys.time(), '] Malformed GFT: both columns gene_type and ',
         'gene_biotype are found.')
  }
  if ('transcript_type' %in% colnames(result) & 
      'transcript_biotype' %in% colnames(result)) {
    stop('[', Sys.time(), '] Malformed GFT: both columns transcript_type and ',
         'transcript_biotype are found.')
  }
  setnames(result, c('transcript_type', 'gene_type'), 
           c('transcript_biotype', 'gene_biotype'), skip_absent = T)
  
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

#' gtfToTxDB
#' @description Reads GTF file into list containing TxDB object and a table 
#' with gene names, IDs, biotypes, as well as transcript ones.
#' @author Maria Litovchenko
#' @param gtfPath path to GTF file
#' @param chrStyle character, one of NCBI or UCSC which determine chromosome
#' naming style (1 or chr1 respectively). Final result will have this 
#' chromosome naming style.
#' @param acceptedChrCodes vector with accepted chromosomal code.
#' @param doLiftOver boolean, indicates, whatever a liftover from one genome
#' version to another should be performed
#' @param chain chain file for liftover. Leave NULL if liftover is FALSE
#' @return list with 3 items: 1) type, string = gtf 2) geneMap: data table with
#' gene names, IDs, biotypes, as well as transcript ones 3) gtfTxdb txdb object
gtfToTxDB <- function(gtfPath, chrStyle, acceptedChrCodes = NULL,
                      doLiftOver = F, chainObj = NULL) {  
  if (doLiftOver & is.null(chainObj)) {
    stop('[', Sys.time(), '] LiftOver is requested, but chain object for ',
         'lifover is NULL')
  }
  
  message('[', Sys.time(), '] Started reading ', gtfPath)
  result <- import(gtfPath)
  
  if (!is.null(acceptedChrCodes)) {
    newSeqLevels <- intersect(as.character(unique(seqnames(sort(result)))),
                              acceptedChrCodes)
    result <- keepSeqlevels(result, value = as.factor(newSeqLevels),
                            pruning.mode = "coarse")
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
  message('[', Sys.time(), '] Finished reading ', gtfPath)
  
  if (doLiftOver & !is.null(chainObj)) {
    resultLO <- liftOverGenomicRegion(inGR = result, chainObj)
    if ('transcript_id' %in% colnames(mcols(result))) {
      result <- checkLiftOverByTr(result, resultLO)
    } else {
      result <- unlist(resultLO)
    }
  }
  seqlevelsStyle(result) <- intersect(c('UCSC', 'NCBI'), chrStyle)
  
  # get a table with gene names, IDs, biotypes, as well as transcript ones
  geneMap <- extractGeneToIDmap(result)
  # convert to txdb & put to result list
  result <- list(type = 'gtf', 'geneMap' = geneMap, 
                 'gtfTxdb' = makeTxDbFromGRanges(result))
  result
}

#' bedToGR
#' @description Reads in the data table a bed file with additional genomic
#' data base.
#' @author Maria Litovchenko
#' @param bedPath path to bed file. It should have header. Minimal columns:
#' chr, start, end, strand, gene_id, gene_name, rCode
#' @param acceptedChrCodes vector with accepted chromosomal code.
#' @param chrStyle character, one of NCBI or UCSC which determine chromosome
#' naming style (1 or chr1 respectively). Final result will have this 
#' chromosome naming style.
#' @param doLiftOver boolean, indicates, whatever a liftover from one genome
#' version to another should be performed
#' @param chain chain file for liftover. Leave NULL if liftover is FALSE
#' @return GRanges object.
bedToGR <- function(bedPath, chrStyle, acceptedChrCodes = NULL, doLiftOver = F, 
                    chainObj = NULL) {
  if (doLiftOver & is.null(chainObj)) {
    stop('[', Sys.time(), '] LiftOver is requested, but chain object for ',
         'lifover is NULL')
  }
  
  message('[', Sys.time(), '] Started reading ', bedPath)
  result <- fread(bedPath, header = T, stringsAsFactors = F)
  minCols <- c('chr', 'start', 'end', 'strand', 'gene_id', 'gene_name', 
               'rCode')
  if (any(!minCols %in% colnames(result))) {
    stop('[', Sys.time(), '] Column(s) ', 
         paste(minCols[!minCols %in% colnames(result)], collapse = ', '), ' ',
         'are not found in ', bedPath)
  }
  message('[', Sys.time(), '] Finished reading ', bedPath)
  
  if (!is.null(acceptedChrCodes)) {
    result <- result[chr %in% acceptedChrCodes]
  }
  
  if (any(c(23, '23', 'chr23') %in% result$chr)) {
    message('[', Sys.time(), '] Found that chr X is coded with 23. Changed it',
            ' to X.')
    result[, chr := gsub('23', 'X', chr)]
    result[, chr := gsub('24', 'Y', chr)]
    result[, chr := gsub('25', 'M', chr)]
    result <- makeGRangesFromDataFrame(result, keep.extra.columns = T)
  }
  result[, chr := gsub('chr', '', chr)]
  
  result <- makeGRangesFromDataFrame(result, keep.extra.columns = T)
  
  if (doLiftOver & !is.null(chainObj)) {
    resultLO <- liftOverGenomicRegion(inGR = result, chainObj)
    if ('transcript_id' %in% colnames(mcols(result))) {
      result <- checkLiftOverByTr(result, resultLO)
    } else {
      result <- unlist(resultLO)
    }
  }
  seqlevelsStyle(result) <- intersect(c('UCSC', 'NCBI'), chrStyle)
  
  result
}

# FUNCTIONS: region extraction from txdb --------------------------------------
#' addRegsCount
#' @description Adds an index within a gene for every region (exon) of a gene,
#' as well as adds column showing total exon count.
#' @author Maria Litovchenko
#' @param gr GRanges object, gene_id column in mcols is needed.
#' @return GRanges object with added columns n_regs and reg_idx, which contain
#' number of regions per gene and region's index within the gene respectively.
addRegsCount <- function(gr) {
  result <- sort(gr)
  
  mcolsUpd <- as.data.table(mcols(result))
  mcolsUpd[, strand := as.character(strand(result))]
  # to return mcols in the exact order as before
  mcolsUpd[, idx := 1:nrow(mcolsUpd)]
  # number of regions
  mcolsUpd[, n_regs := .N, by = gene_id]
  
  # have to split by strand because for forward strand 1st exon is the most 
  # left one, and for the reverse strand - the most right one
  mcolsUpd_F <- mcolsUpd[strand == '+']
  mcolsUpd_R <- mcolsUpd[strand == '-']
  mcolsUpd_F[, reg_idx := 1:.N, by = gene_id]  # Forward
  mcolsUpd_R[, reg_idx := .N:1, by = gene_id]  # Reverse
  mcolsUpd <- rbind(mcolsUpd_F, mcolsUpd_R, mcolsUpd[!strand %in% c('+', '-')],
                    fill = T)
  mcolsUpd <- as.data.table(mcolsUpd)
  mcolsUpd <- mcolsUpd[order(idx)]
  mcolsUpd[, idx := NULL]
  mcolsUpd[, strand := NULL]
  mcols(result) <- mcolsUpd
  result
}

#' extract_ss_txdb
#' @description Extracts splice sites coordinates from txdb object
#' @param txdbObj - txdb, result of makeTxDbFromGRanges or similar function
#' @param upstr length of intronic region upstream of acceptor site
#' @param downstr length of intronic region downstream of donor site.
#' @return GRanges
extract_ss_txdb <- function(txdbObj, upstr, downstr) {
  # sometimes it can happen that gene has a UTR split between 2 exons. In such
  # case intronsByTranscript will also give the intron sandwitched by UTR 
  # pieces. We don't want that, we only want splice sites located betwen 2 
  # protein coding exons. That's why we're going to extract gaps between CDS.
  intronsGR <- intronsByTranscript(txdbObj, use.names = T)
  # remove genes without introns
  intronsGR <- intronsGR[unique(names(unlist(intronsGR)))]
  intronsGR <- unlist(intronsGR)
  
  # extract UTRs
  utrGR <- c(threeUTRsByTranscript(txdbObj, use.names = T), 
             fiveUTRsByTranscript(txdbObj, use.names = T))
  utrGR <- unlist(utrGR)
  utrGR <- utrGR[which(names(utrGR) %in% 
                         intersect(names(utrGR), names(intronsGR)))]
  
  # add 1 base from each side of UTR so they could be overlapped with introns
  utrGR <- utrGR %>% anchor_5p() %>% mutate(width = width + 1)
  utrGR <- utrGR %>% anchor_3p() %>% mutate(width = width + 1)
  
  # find introns between UTRs
  if (length(utrGR) != 0) {
    ovrl <- as.data.table(findOverlaps(intronsGR, utrGR))
    colnames(ovrl) <- c('intronIdx', 'utrIdx')
    ovrl[, intronID := names(intronsGR)[intronIdx]]
    ovrl[, utrID := names(utrGR)[utrIdx]]
    ovrl <- ovrl[intronID == utrID]
    # remove introns between UTRs
    intronsGR <- intronsGR[setdiff(1:length(intronsGR), 
                                   unique(ovrl$intronIdx))]
  }
  
  donor <- intronsGR %>% anchor_5p() %>% mutate(width = downstr)
  acceptor <- intronsGR %>% anchor_3p() %>% mutate(width = upstr)
  result <- c(donor, acceptor)
  result <- unique(result)
  result$transcript_id <- names(result)
  result
}

#' extract_ss_Granges
#' @description Extracts splice sites from GRaneges object
#' @author Maria Litovchenko
#' @param gr GRanges object, gene_id column in mcols is needed.
#' @param upstr length of intronic region upstream of acceptor site
#' @param downstr length of intronic region downstream of donor site.
#' @return Granges object with extracted splice sites
#' @note use this function if and only if txdb version of annotation is not 
#' available.
extract_ss_Granges <- function(gr, upstr, downstr) {
  if (upstr == 0 & downstr == 0) {
    stop('[', Sys.time(), '] Attempt to extract splice sites with 0 down',
         'stream and upstream length.')
  }
  # check on strand
  if (any(strand(gr) == '*') & downstr != upstr) {
    message('[', Sys.time(), '] It is not possible to extract splice sites ',
            'in case strand is not set and requested downstream and ',
            'upstream lengths differ. Will remove regions without + or - ',
            'strand.')
    nbefore <- length(gr)
    result <- gr[strand(gr) %in% c('+', '-')]
    nafter <- length(result)
    message('[', Sys.time(), '] Removed ', nbefore - nafter, '(',
            round(100 * (nbefore - nafter) / nbefore, 2), '%) entries ', 
            'because they did not have assigned strand ')
  } else {
    result <- copy(gr)
  }
  
  # add exon count and maximum number of exons. It's faster than though the
  # loop of GRanges
  result <- addRegsCount(result)
  # restrict to genes with > 1 exon
  result <- result[result$n_regs > 1]
  
  # result$reg_idx != result$n_regs because there is no splice  
  # site after the last exon
  resultDown <- result[result$reg_idx != result$n_regs]
  resultDown <- shift_downstream(resultDown, downstr + 1) %>% 
    anchor_3p() %>% mutate(width = downstr + 1)
  # result$reg_idx != 1 because there is no splice site before 1st exon
  resultUp <- result[result$reg_idx != 1]
  resultUp <- shift_upstream(resultUp, upstr + 1) %>% anchor_5p() %>% 
    mutate(width = upstr + 1)
  result <- sort(c(resultDown, resultUp))
  
  mcols(result) <- mcols(result)[, !colnames(mcols(result)) %in%
                                   c('n_regs', 'reg_idx')]
  result
}

#' extract_promoters_Granges
#' @description Extracts promoters from GRaneges object
#' @author Maria Litovchenko
#' @param gr GRanges object, gene_id column in mcols is needed.
#' @param upstr length to take upstream of TSS.
#' @param downstr length to take downstream of TSS.
#' @note use this function if and only if txdb version of annotation is not 
#' available.
extract_promoters_Granges <- function(gr, upstr, downstr) {
  if (upstr == 0 & downstr == 0) {
    stop('[', Sys.time(), '] Attempt to extract promoters with 0 downstream ',
         'and upstream length.')
  }
  
  if (any(strand(gr) == '*')) {
    message('[', Sys.time(), '] It is not possible to extract promoters ',
            'in case strand is not set. Will remove regions without + or - ',
            'strand.')
    nbefore <- length(gr)
    result <- gr[strand(gr) %in% c('+', '-')]
    nafter <- length(result)
    message('[', Sys.time(), '] Removed ', nbefore - nafter, '(',
            round(100 * (nbefore - nafter) / nbefore, 2), '%) entries ', 
            'because they did not have assigned strand ')
  } else {
    result <- copy(gr)
  }
  
  # add exon count and maximum number of exons. It's faster than though the
  # loop of GRanges
  result <- addRegsCount(result)
  # restrict only to the first exon
  result <- result[result$reg_idx == 1]
  result <- result %>% anchor_5p() %>% mutate(width = 1) %>% 
    anchor_5p() %>% stretch(extend = downstr) %>% anchor_3p() %>% 
    stretch(extend = upstr)
  
  result <- sort(result)
  mcols(result) <- mcols(result)[, !colnames(mcols(result)) %in%
                                   c('n_regs', 'reg_idx')]
  result
}

#' extractRegionsByCode_txdb
#' @description Extracts coordinates of desired type of genomic regions (i.e. 
#' CDS, 3primeUTR, etc) from txdb object
#' @author Maria Litovchenko
#' @param rCode region code, one of 3primeUTR, 5primeUTR, CDS, lincRNA, 
#' lincRNA_promoter, lincRNA_ss, miRNA, misc_RNA, promoter, rRNA, snoRNA, 
#' snRNA, ss, upstream, downstream
#' @param txdbObj txdb object 
#' @param geneBiotypeDT a data table with columns gene_name, gene_id, 
#' gene_biotype, transcript_id, transcript_biotype (if available). 
#' transcript_biotype is not available for earlier versions of ensemble gtf, 
#' i.e. v75 doesn't have it. Usually would be extracted from GRanges object, 
#' preceeding creation of txdbObj
#' @param upstr upstream length to extract. For splice sites it will be length
#' of intronic region upstream of acceptor site. For promoters it will be 
#' length upstream TSS.
#' @param downstr dowstream length to extract. For splice sites it would be 
#' length of intronic region downstream of donor site. For promoters it will be 
#' length downstream TSS.
#' @param reduceDisjoin boolean, indicates, whatever regions needed to be 
#' reduced and disjoined on by gene level
#' @return list of genomic ranges. One item = one genes. Extra columns: 
#' gene_id, gene_name, rCode
#' @note rCode is going to be matched against transcript_biotype (if available)
#' and gene_biotype otherwise. For CDS, ss, 3primeUTR, 5primeUTR and promoter 
#' transcript_biotype/gene_biotype is set to protein_coding. transcript_biotype
#' is not available for earlier versions of ensemble gtf, i.e. v75 doesn't have
#' it.
extractRegionsByCode_txdb <- function(rCode, txdbObj, geneBiotypeDT, upstr = 0,
                                      downstr = 0, reduceDisjoin = T) {
  upstr <- as.integer(upstr)
  downstr <- as.integer(downstr)
  message('[', Sys.time(), '] Extracting ', rCode, ' upstream length: ', 
          upstr, ', downstream length: ', downstr)
  
  # extract regions of interest
  result <- switch(rCode, 
                   "promoter" = promoters(txdbObj, upstream = upstr, 
                                          downstream = downstr, use.names = T),
                   "5primeUTR" = fiveUTRsByTranscript(txdbObj, use.names = T),
                   "CDS" = cdsBy(txdbObj, by = 'tx', use.names = T),  # cdsBy returns without 3 and 5 UTR
                   'ss' = extract_ss_txdb(txdbObj, upstr, downstr),
                   "3primeUTR" = threeUTRsByTranscript(txdbObj, use.names = T),
                   'lincRNA_promoter' = promoters(txdbObj, upstream = upstr, 
                                                  downstream = downstr,
                                                  use.names = T),
                   'lincRNA' = exonsBy(txdbObj, 'tx', use.names = T),
                   'lincRNA_ss' = extract_ss_txdb(txdbObj, upstr, downstr),
                   'miRNA' = exonsBy(txdbObj, 'tx', use.names = T),
                   'misc_RNA' = exonsBy(txdbObj, 'tx', use.names = T),
                   'rRNA' = exonsBy(txdbObj, 'tx', use.names = T),
                   'snoRNA' = exonsBy(txdbObj, 'tx', use.names = T),
                   'snRNA' = exonsBy(txdbObj, 'tx', use.names = T))
  
  if (class(result) == 'CompressedGRangesList') {
    result <- unlist(result)
  }
  if (length(result) == 0) {
    message('[', Sys.time(), '] 0 regions is extracted! Returning empty ',
            'GRanges')
    return(result)
  }
  
  # select proper transcript biotype
  trBioType <- switch(rCode, "promoter" = 'protein_coding',
                      "5primeUTR" = 'protein_coding', "CDS" = 'protein_coding', 
                      'ss' = 'protein_coding', "3primeUTR" = 'protein_coding', 
                      'lincRNA_promoter' = 'lincRNA', 'lincRNA' = 'lincRNA', 
                      'lincRNA_ss' = 'lincRNA', 'miRNA' = 'miRNA', 
                      'misc_RNA' = 'misc_RNA', 'rRNA' = 'rRNA', 
                      'snoRNA' = 'snoRNA', 'snRNA' = 'snRNA')
  if ('transcript_biotype' %in% colnames(geneBiotypeDT)) {
    mcolsResult <- geneBiotypeDT[transcript_biotype %in% trBioType]
  } else {
    mcolsResult <- geneBiotypeDT[gene_biotype %in% trBioType]
  }
  result <- result[names(result) %in% mcolsResult$transcript_id]
  
  if (length(result) == 0) {
    message('[', Sys.time(), '] 0 regions is extracted! Returning empty ',
            'GRanges')
    return(result)
  }
  
  # add upstream and downstream bases, if requested (except ss and promoter)
  if (!rCode %in% c('ss', 'promoter', 'lincRNA_promoter', 'lincRNA_ss',
                    'upstream', 'downstream') & 
      (upstr != 0 | downstr != 0)) {
    upstrRegs <- NULL
    downstrRegs <- NULL
    if (upstr != 0) {
      upstrRegs <- shift_upstream(result, upstr) %>% anchor_5p() %>% 
        mutate(width = upstr)
    }
    if (downstr != 0) {
      downstrRegs <- shift_downstream(result, downstr) %>% anchor_3p() %>% 
        mutate(width = downstr)
    }
    result <- c(upstrRegs, downstrRegs)
  }
  if (length(result) == 0) {
    message('[', Sys.time(), '] 0 regions is extracted! Returning empty ',
            'GRanges')
    return(result)
  }
  
  # split by gene
  geneToTxMap <- mcolsResult[,.(transcript_id, gene_id, gene_name)]
  geneToTxMap <- geneToTxMap[!duplicated(geneToTxMap)]
  setkey(geneToTxMap, 'transcript_id')
  mcols(result) <- geneToTxMap[names(result)][,.(gene_id, gene_name)]
  result$rCode <- rCode
  result <- split(result, result$gene_id)
  
  # reduce and disjoin, if needed
  if (reduceDisjoin) {
    result <- unlist(disjoin(reduce(result)))
    geneNameToIDmap <- geneToTxMap[,.(gene_id, gene_name)]
    geneNameToIDmap <- geneNameToIDmap[!duplicated(geneNameToIDmap)]
    setkey(geneNameToIDmap, 'gene_id')
    mcols(result) <- data.table(gene_id = names(result),
                                gene_name = geneNameToIDmap[names(result)]$gene_name,
                                'rCode' = rCode)
    result <- split(result, result$gene_id)
  }
  
  # remove chr from chromosome names for simplicity of handling
  seqlevelsStyle(result) <- "NCBI"
  
  result
}

#' extractRegionsByCode_gr
#' @description Extracts coordinates of desired type of genomic regions from 
#' Granges object derived from custom bed file.
#' @author Maria Litovchenko
#' @param rCode region code, can be any string. ss and promoter are a special 
#' codes.
#' @param grObj GRanges object 
#' @param upstr upstream length to extract. For splice sites it will be length
#' of intronic region upstream of acceptor site. For promoters it will be 
#' length upstream TSS.
#' @param downstr dowstream length to extract. For splice sites it would be 
#' length of intronic region downstream of donor site. For promoters it will be 
#' length downstream TSS.
#' @param reduceDisjoin boolean, indicates, whatever regions needed to be 
#' reduced and disjoined on by gene level
#' @return list of genomic ranges. One item = one genes. Extra columns: 
#' gene_id, gene_name, rCode
extractRegionsByCode_gr <- function(rCode, grObj, upstr = 0, downstr = 0, 
                                    reduceDisjoin = T) {
  upstr <- as.integer(upstr)
  downstr <- as.integer(downstr)
  message('[', Sys.time(), '] Extracting ', rCode, ' upstream length: ', 
          upstr, ', downstream length: ', downstr)
  msg <-  paste0('[', Sys.time(), '] Special code (', rCode ,') is used on ',
                 'bed data. REPLACE_HERE If you want standard behavior ',
                 '(i.e. just match rCode from bed to gr_code from ',
                 'the analysis inventory table0, change gr_code in the ',
                 'inventory table from ss to something else and change rCode ',
                 'in bed respectively.')
  
  if (rCode %in% c('ss')) {# splice sites, %in% to work with NA
    rplMsg <- paste('Will extract splice sites as intronic regions in',
                    'between regions in submitted bed file.')
    message(gsub('REPLACE_HERE', rplMsg, msg))
    result <- extract_ss_Granges(grObj, upstr, downstr)
  }
  
  if (rCode %in% c('promoter')) {# promoters
    rplMsg <- paste('Will extract promoter as region around start site of ',
                    'the first "exon" belonging to a gene.')
    message(gsub('REPLACE_HERE', rplMsg, msg))
    result <- extract_promoters_Granges(grObj, upstr, downstr)
  }
  
  if (!rCode %in% c('ss', 'promoter')) {
    if (!is.na(rCode)) {
      result <- grObj[grObj$rCode == rCode]
    } else {
      message('[', Sys.time(), '] NA as rCode is detected. All regions will ',
              'pass')
      result <- copy(grObj)
    }
    
    # add upstream and downstream, if requested (except ss and promoter)
    if (upstr != 0 | downstr != 0) {
      upstrRegs <- NULL
      downstrRegs <- NULL
      if (upstr != 0) {
        upstrRegs <- shift_upstream(result, upstr) %>% anchor_5p() %>% 
          mutate(width = upstr)
      }
      if (downstr != 0) {
        downstrRegs <- shift_downstream(result, downstr) %>% anchor_3p() %>% 
          mutate(width = downstr)
      }
      result <- c(upstrRegs, downstrRegs)
    }
  }
  
  # reduce and disjoin, if needed
  if (reduceDisjoin) {
    geneNameToIDmap <- mcols(result)[, c('gene_id', 'gene_name')]
    geneNameToIDmap <- as.data.table(geneNameToIDmap)
    geneNameToIDmap <- geneNameToIDmap[!duplicated(geneNameToIDmap)]
    setkey(geneNameToIDmap, 'gene_id')
    
    result <- unlist(disjoin(reduce(split(result, result$gene_id))))
    mcols(result) <- data.table(gene_id = names(result),
                                gene_name = geneNameToIDmap[names(result)]$gene_name,
                                'rCode' = rCode)
  }
  
  # split by gene
  result <- split(result, result$gene_id)
  # remove chr from chromosome names for simplicity of handling
  seqlevelsStyle(result) <- "NCBI"
  result
}

# FUNCTIONS: removing regions to exclude  -------------------------------------
#' reduceGRkeepMcols 
#' @description Performs reduce & disjoin of list of genomic regions on by gene
#' bases, keeping mcols information
#' @author Maria Litovchenko
#' @param grlist List of GRanges objects, each containing gene_id, gene_name 
#' and rCode as columns in mcols
#' @return GRanges object with reduce & disjoined regions. rCodes of 
#' overlapping original regions are put separated by comma.
reduceGRkeepMcols <- function(grlist) {
  GRs <- unlist(GRangesList(lapply(grlist, unlist)))
  names(GRs) <- NULL
  
  # reduce & disjoin on gene basis
  result <- disjoin(reduce(split(GRs, GRs$gene_id)))
  result <- unlist(result)
  # add gene ID
  result$gene_id_r <- names(result)
  # get overlap with initial regions to retrieve mcols
  ovrl <- as.data.table(findOverlaps(result, GRs, type = "any", 
                                     select = 'all',ignore.strand = F))
  colnames(ovrl) <- c('resultIdx', 'targetIdx')
  # add mcols to result. Both all result and all GRs are in ovrl, because 
  # we just did disjoin reduce.
  result <- cbind(as.data.table(result)[ovrl$resultIdx],
                  as.data.table(GRs)[ovrl$targetIdx][,.(gene_id, gene_name,
                                                        rCode)])
  # against overlapping genes
  result <- result[gene_id_r == gene_id]
  result[, gene_id_r := NULL] 
  result[, width := NULL]
  # collapse, so that we have a unique genomic region - gene pair.
  result <- result[,.(gene_name = unique(gene_name),
                      rCode = paste(unique(rCode), collapse = ', ')),
                   by = .(seqnames, start, end, strand, gene_id)]
  result <- makeGRangesFromDataFrame(result, keep.extra.columns = T)
  result
}

#' excludeFromTarget
#' @description Excludes regions from targetGR if they overlap with targetGR 
#' in case of black listed regions (type = 'black') or only includes regions 
#' overlapping between targetGR and excludeGR in case of white listed regions
#' (type = 'white')
#' @author Maria Litovchenko
#' @param targetGR GRanges List or GRanges object, target regions
#' @param excludeGR GRanges List or GRanges object, regions to exclude or over-
#' lap with 
#' @param type character, white for whitelisted and black for blacklisted
#' @param ignore.strand boolean, indicates whatever strand information in 
#' targetGR and excludeGR should be taken into account. In case a strand is not
#' set for a region in excludeGR, it will be matched against both positive and
#' negative strands.
#' @return GRanges object
excludeFromTarget <- function(targetGR, excludeGR, type, ignore.strand = T) {
  if (!type %in% c('white', 'black')) {
    stop('[', Sys.time(), '] excludeFromTarget: type should be black or white')
  }
  
  # unlist, if needed
  if (class(targetGR) == 'CompressedGRangesList') {
    targetGR <- unlist(targetGR)
  }
  if (class(excludeGR) == 'CompressedGRangesList') {
    excludeGR <- unlist(excludeGR)
  }
  seqlevelsStyle(excludeGR) <- intersect(c('UCSC', 'NCBI'), 
                                         seqlevelsStyle(targetGR))
   
  # if one of targetGR or excludeGR do not have strand assigned, then result
  # of setdiff or intersect will also have no strand. In order to preserve the
  # strand information we need to 1) in case of ignoring strand completely
  # we need to duplicate regions to exclude and give them 3 types of strand
  # (+, - and *) strands to cater for all possible strands in 2) in case of not
  # ignoring strand in regions to exclude we only need to duplicate regions 
  # to exclude which have * strand and assign them + or - strand, because 
  # essentially they can be on any.
  if (ignore.strand) {
    excludeGR_plus <- copy(excludeGR)
    strand(excludeGR_plus) <- '+'
    excludeGR_minus <- copy(excludeGR)
    strand(excludeGR_minus) <- '-'
    excludeGR_none <- copy(excludeGR)
    strand(excludeGR_none) <- '*'
    excludeGR <- c(excludeGR_plus, excludeGR_minus, excludeGR_none)
  } else {
    excludeGR_none <- excludeGR[strand(excludeGR) == '*']
    excludeGR_none_plus <- copy(excludeGR_none)
    strand(excludeGR_none_plus) <- '+'
    excludeGR_none_minus <- copy(excludeGR_none)
    strand(excludeGR_none_minus) <- '-'
    excludeGR <- c(excludeGR[strand(excludeGR) != '*'], excludeGR_none_plus,
                   excludeGR_none_minus)
  }
  
  # here ignore.strand is F because we gave all possible strand combinations
  # of strand in excludeGR and we want to preserve strand in targetGR
  result <- switch(type, 
                   'black' = setdiff(targetGR, excludeGR, ignore.strand = F), 
                   'white' = intersect(targetGR, excludeGR, ignore.strand = F))
  if (length(result) == 0) {
    return(result)
  }

  # return mcols to result
  ovrl <- findOverlaps(result, targetGR, type = 'within', ignore.strand = F, 
                       select = 'all')
  ovrl <- as.data.table(ovrl)
  colnames(ovrl) <- c('result_idx', 'target_idx')
  if (length(setdiff(1:length(result), ovrl$result_idx)) != 0) {
    # experience showed that if one of the result region's ends overlaps with
    # the start/end of targetGR, and it is fully inside, findOverlaps with
    # type = 'within' won't recognize it! So, let's give it a chance without
    # type = 'within'.
    rogueRegs <- result[setdiff(1:length(result), ovrl$result_idx)]
    ovrl_patch <- findOverlaps(rogueRegs, targetGR, ignore.strand = F,  
                               select = 'all')
    ovrl_patch <- as.data.table(ovrl_patch)
    colnames(ovrl_patch) <- c('result_idx', 'target_idx')
    # well, now, if something isn't found, than it's definetely a major bug
    if (length(setdiff(1:length(rogueRegs), ovrl_patch$result_idx)) != 0) {
      stop('[', Sys.time(), '] MAJOR ERROR: excludeFromTarget created ',
           'genomics regions which are not in targets')
    }
  }
  # in case regions with strand '*' are submitted to findOverlaps function, it
  # will overlap them with all regions, even if ignore.strand = F. Therefore,
  # we need to cake care of the strand
  ovrl[, result_strand := as.character(strand(result))[result_idx]]
  ovrl[, target_strand := as.character(strand(targetGR))[target_idx]]
  ovrl <- ovrl[result_strand == target_strand]
  result <- result[ovrl$result_idx]
  mcols(result) <- mcols(targetGR)[ovrl$target_idx, ]
  names(result) <- NULL
  
  # print some information
  lenBefore <- sum(width(targetGR))
  lenAfter <- sum(width(result))
  message('[', Sys.time(), '] ',
          ifelse(type == 'black', 'Excluded ', 'Overlapped '),
          paste(unique(excludeGR$rCode), collapse = ', '), ' regions ',
          ifelse(type == 'black', 'from ', 'with '),
          paste(unique(targetGR$rCode), collapse = '+'), ' region set. ',
          'Removed ', lenBefore - lenAfter, ' from ', lenAfter, ' bases (',
          round(100 * (lenBefore - lenAfter) / lenBefore), '%).')
  
  result
}

#' getCleanRegions
#' @description Returns GRanges containing target genomic regions without 
#' regions we want to remove from them. Essentially a wrapper around 
#' excludeFromTarget. Does not handle black-/white- listed regions
#' @author Maria Litovchenko 
#' @param inventoryDT an inventory data table which determines target regions 
#' and regions to remove for one gr_id. Essential columns: gr_code, gr_file, 
#' gr_upstr, gr_downstr, gr_excl_code, gr_excl_file, gr_excl_upstr, 
#' gr_excl_downstr.
#' @param allAvailableRegs named list of lists all available genomic regions. 
#' First list - file, nested list - genomic region.
#' @param ignore.strand boolean, indicates whatever strand information in 
#' targetGR and excludeGR should be taken into account. In case a strand is not
#' set for a region in excludeGR, it will be matched against both positive and
#' negative strands.
#' @return GenomicRanges of target regions, cleaned from regions we want to 
#' remove.
getCleanRegions <- function(inventoryDT, allAvailableRegs, ignore.strand = T) {
  if (length(unique(inventoryDT$gr_id)) != 1) {
    stop('[', Sys.time(), '] Only 1 gr_id is allowed into function ',
         'getCleanRegions')
  }
  
  # target genomic regions
  targetsDT <- inventoryDT[,.(gr_id, gr_file, gr_code, gr_upstr, gr_downstr)]
  targetsDT <- targetsDT[!duplicated(targetsDT)]
  targetsDT[, id := apply(targetsDT[,.(gr_code, gr_upstr, gr_downstr)], 1, 
                          paste, collapse = '-')]
  targetsDT[, id := gsub(' ', '', id)]
  targetsGR <- apply(targetsDT, 1, 
                     function(x) allAvailableRegs[[x['gr_file']]][[x['id']]])
  if (length(targetsGR) > 1) {
    targetsGR <- reduceGRkeepMcols(targetsGR)
  } else {
    targetsGR <- unlist(GRangesList(lapply(targetsGR, unlist)))
  }
  targetsGR$rCode_print <- copy(targetsGR$rCode)
  targetsGR$rCode <- unique(inventoryDT$gr_id)
  
  # genomic regions to remove
  toRemoveDT <- inventoryDT[,.(gr_excl_id, gr_excl_file, gr_excl_code,
                               gr_excl_upstr, gr_excl_downstr)]
  toRemoveDT <- toRemoveDT[!duplicated(toRemoveDT)]
  toRemoveDT <- toRemoveDT[!is.na(gr_excl_code)]
  toRemoveGR <- NULL
  if (nrow(toRemoveDT) != 0) {
    toRemoveDT[, id := apply(toRemoveDT[,.(gr_excl_code, gr_excl_upstr, 
                                           gr_excl_downstr)], 1,
                             paste, collapse = '-')]
    toRemoveDT[, id := gsub(' ', '', id)]
    toRemoveDT[, print_id := apply(toRemoveDT[,.(gr_excl_id, gr_excl_code,
                                                 gr_excl_upstr, 
                                                 gr_excl_downstr)], 1,
                                  function(x) paste0(x['gr_excl_id'], '(', 
                                                     x['gr_excl_code'], ',u:',
                                                     x['gr_excl_upstr'], ',d:',
                                                     x['gr_excl_downstr'], 
                                                     ')'))]
    toRemoveDT[, print_id := gsub(' ', '', print_id)]
    
    toRemoveGR <- apply(toRemoveDT, 1, 
                        function(x) allAvailableRegs[[x['gr_excl_file']]][[x['id']]])
    toRemoveGR <- toRemoveGR[sapply(toRemoveGR, function(x) length(x) != 0)]
    toRemoveGR <- unlist(GRangesList(lapply(toRemoveGR, unlist)))
    # we don't care about mcols for the regions we're going to remove.
    toRemoveGR <- disjoin(reduce(toRemoveGR))
    toRemoveGR$rCode <- paste(toRemoveDT$print_id, collapse = ', ')
  }
  
  if (!is.null(toRemoveGR)) {
    result <- excludeFromTarget(targetsGR, toRemoveGR, 'black', ignore.strand)
  } else {
    result <- targetsGR
  }
  result$rCode <- result$rCode_print
  result$rCode_print <- NULL
  
  result
}

#' filterBWregions
#' @description 
#' @author Maria Litovchenko
#' @param bwGRs GRanges object representation of black and white regions
#' @param bwScoreCol string, name of the column containing scores on which 
#' regions should be filtred.
#' @param bwScoreMin numeric, minimum value of a score
#' @param bwScoreMax numeric, maximum value of a score
#' @return filtered Granges
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

# FUNCTIONS : identification & processing of overlapping regions --------------
#' calcPercOverlap
#' @description Calculates percentage of bases of a region from targetGR 
#'              overlapping with any region in backgroundGR. In case targetGR
#'              and backgroundGR are the same, overlap of a region to itself
#'              will not be counted.
#' @author Maria Litovchenko
#' @param targetGR GRanges object, target regions.
#' @param backgroundGR GRanges object, background regions.
#' @param ignore.strand boolean, whatever strand should be ignored
#' @param select when select is "all" (the default), then all target to 
#'               background pairs will be returned. max_target - for target  
#'               genomic regions pair with the maximum percentage of 
#'               intersection with background is reported, min_target - same,
#'               but minimum. max_background and min_background - same, but
#'               for background.
#' @return data table with columns targetIdx, bgIdx, intersectionWidth,
#'         intersectionPercTarget, intersectionPercBg. targetIdx - index in 
#'         targetGR, bgIdx - index in backgroundGR, intersectionWidth - number
#'         of basepairs in intersection, intersectionPercTarget - percentage
#'         of target region which overlaps with the background region, 
#'         intersectionPercBg - percentage of background region which overlaps
#'         with the target region. Values in intersectionPerc* range from 0 to
#'         100.
calcPercOverlap <- function(targetGR, backgroundGR, remove.self = T, 
                            ignore.strand = F,
                            select = 'all') {
  if (!select %in% c('all', 'max_target', 'min_target', 'max_background', 
                     'min_background')) {
    stop('[', Sys.time(), '] Error in calcPercOverlap: select can only be ',
         'one of all, max_target, min_target, max_background, min_background.')
  }
  
  # Step 1: find overlapping pairs of regions
  ovrl <- findOverlaps(targetGR, backgroundGR, ignore.strand = ignore.strand)
  ovrl <- as.data.table(ovrl)
  colnames(ovrl) <- c('targetIdx', 'backgroundIdx')
  
  if (!identical(targetGR, backgroundGR) & remove.self) {
    message('[', Sys.time(), '] remove.self was set to True, although ',
            'targetGR and backgroundGR are not the same. Will ignore it.')
  }
  
  if (identical(targetGR, backgroundGR) & remove.self) {
    # this line covers the case target and background are actually the 
    # same genomic regions and we want to compute overlap with the other  
    # regions except self
    ovrl <- ovrl[targetIdx != backgroundIdx]
    # also, in case target and background are actually the same there could be 
    # duplicated pairs of overlaps, i.e. 1 - 2 and 2 - 1 
    ovrl[, id := apply(ovrl, 1, function(x) paste(sort(x), collapse = '_'))]
    ovrl <- ovrl[!duplicated(id)]
    ovrl[, id := NULL]
  }
  
  # Step 2: for each pair, get region in intersect
  # ... first, we need to create GRangesList each element of which holds 
  # overlapping region pair
  ovrl[, pairIdx := apply(ovrl, 1, paste, collapse = '--')]
  ovrlPairs_targets <- targetGR[ovrl$targetIdx]
  mcols(ovrlPairs_targets) <- cbind(mcols(ovrlPairs_targets), 
                                    pairIdx = ovrl$pairIdx)
  ovrlPairs_bg <- backgroundGR[ovrl$backgroundIdx]
  mcols(ovrlPairs_bg) <- cbind(mcols(ovrlPairs_bg), pairIdx = ovrl$pairIdx)
  ovrlPairs <- c(ovrlPairs_targets, ovrlPairs_bg)
  rm(ovrlPairs_targets, ovrlPairs_bg)
  ovrlPairs <- split(ovrlPairs, ovrlPairs$pairIdx)
  # now, we can get region of intersection. As intersect/pintersect functions
  # do not operate on lists, we have to be creative with reduce. Reduce will
  # produce 3 regions. Then sorted, the middle one will contain overlap region.
  # If there is only 1 region after disjoin - regions overlap fully, if there
  # are 2 regions - one of the pair is inside the other one and the width of
  # the smaller one needs to be reported.
  disjoinedGR <- sort(disjoin(ovrlPairs, ignore.strand = ignore.strand))
  disjoinedGRwidth <- width(disjoinedGR)
  disjoinedGRlen <- sapply(disjoinedGRwidth, length)
  
  intersectLen <- lapply(disjoinedGRwidth, function(x) x[2])
  intersectLen[disjoinedGRlen != 3] <- lapply(sort(width(ovrlPairs[disjoinedGRlen != 3])),
                                              function(x) x[1])
  intersectLen <- unlist(intersectLen)
  
  result <- data.table(targetIdx = sapply(names(intersectLen), 
                                          function(x) gsub('--.*', '', x)))
  result[, targetIdx := as.integer(targetIdx)]
  result[, bgIdx := sapply(names(intersectLen), 
                           function(x) gsub('.*--', '', x))]
  result[, bgIdx := as.integer(bgIdx)]
  result[, intersectionWidth := intersectLen]
  # compute percentage of intersection (from target)
  result[, intersectionPercTarget := sapply(width(ovrlPairs), 
                                            function(x) x[1])]
  result[, intersectionPercTarget := 100 * intersectionWidth / 
           intersectionPercTarget]
  # compute percentage of intersection (from background)
  result[, intersectionPercBg := sapply(width(ovrlPairs), function(x) x[2])]
  result[, intersectionPercBg := 100 * intersectionWidth / 
           intersectionPercBg]
  
  if (select == 'max_target') {
    result <- result[order(-intersectionPercTarget)][,.SD[1], by = targetIdx]
  }
  if (select == 'min_target') {
    result <- result[order(intersectionPercTarget)][,.SD[1], by = targetIdx]
  }
  if (select == 'max_background') {
    result <- result[order(-intersectionPercBg)][,.SD[1], by = bgIdx]
  }
  if (select == 'min_background') {
    result <- result[order(intersectionPercBg)][,.SD[1], by = bgIdx]
  }
  result <- result[order(targetIdx, bgIdx)]
  result
}

#' calcPercOverlap_byGeneID 
#' @description 
#' @author 
#' @author Maria Litovchenko
#' @param targetGR GRanges object, target regions.
#' @param backgroundGR GRanges object, background regions.
#' @param ignore.strand boolean, whatever strand should be ignored
#' @param select when select is "all" (the default), then all target to 
#'               background pairs will be returned. max_target - for target  
#'               genomic regions pair with the maximum percentage of 
#'               intersection with background is reported, min_target - same,
#'               but minimum. max_background and min_background - same, but
#'               for background.
#' @return data table with columns target_gene_id, target_gene_name, 
#'         bg_gene_id, bg_gene_name, intersectionWidth, target_len, bg_len, 
#'         intersectionPercTarget, intersectionPercBg. target_gene_id - gene id
#'         of a target gene, target_gene_name - gene name of a target gene. 
#'         bg_gene_id, bg_gene_name - same, but for background genes. 
#'         target_len - sum of lengths of of all regions for taget gene. bg_len
#'         the same, but for all regions for background gene. intersectionWidth 
#'         - number of basepairs in intersection, intersectionPercTarget - 
#'         percentage of target region which overlaps with the background 
#'         region, intersectionPercBg - percentage of background region which 
#'         overlaps with the target region. Values in intersectionPerc* range 
#'         from 0 to 100.
calcPercOverlap_byGeneID <- function(targetGR, backgroundGR, ignore.strand = F, 
                                     remove.self = T, select = 'all') {
  # Step 1: find overlaps between individual regions
  result <- calcPercOverlap(targetGR, backgroundGR, remove.self, ignore.strand, 
                            select = 'all')
  
  # Step 2: since per gene we can have several individual regions, we will have
  # to recompute percentage of overlapping bases for background and target 
  # taking into account all regions which belong to a gene
  result <- result[,.(targetIdx, bgIdx, intersectionWidth)]
  # add target gene name and id
  result <- cbind(result, as.data.table(mcols(targetGR)[result$targetIdx, 
                                                        c('gene_id', 
                                                          'gene_name')]))
  setnames(result, c('gene_id', 'gene_name'), 
           c('target_gene_id', 'target_gene_name'))
  # same for background
  result <- cbind(result, as.data.table(mcols(backgroundGR)[result$bgIdx, 
                                                            c('gene_id',
                                                              'gene_name')]))
  setnames(result, c('gene_id', 'gene_name'), c('bg_gene_id', 'bg_gene_name'))
  # sum lengths of all intersecting regions for each target - background pair
  result <- result[,.(target_gene_id, target_gene_name, bg_gene_id, 
                      bg_gene_name, intersectionWidth)]
  result <- result[, intersectionWidth := sum(intersectionWidth), 
                   by = .(target_gene_id, target_gene_name, 
                          bg_gene_id, bg_gene_name)]
  result <- unique(result)
  # compute total length of target regions per gene & add it
  regsTotalLen <- data.table(gene_id = targetGR$gene_id, 
                             gene_name = targetGR$gene_name, 
                             target_len = width(targetGR))
  regsTotalLen <- unique(regsTotalLen[, target_len := sum(target_len), 
                                      by = .(gene_id, gene_name)])
  setnames(regsTotalLen, c('gene_id', 'gene_name'), 
           c('target_gene_id', 'target_gene_name'))
  result <- merge(result, regsTotalLen,  all.x = T,
                  by = c('target_gene_id', 'target_gene_name'))
  # same for background
  regsTotalLen <- data.table(gene_id = backgroundGR$gene_id, 
                             gene_name = backgroundGR$gene_name,
                             bg_len = width(backgroundGR))
  regsTotalLen <- unique(regsTotalLen[, bg_len := sum(bg_len), 
                                      by = .(gene_id, gene_name)])
  setnames(regsTotalLen, c('gene_id', 'gene_name'), 
           c('bg_gene_id', 'bg_gene_name'))
  result <- merge(result, regsTotalLen,  all.x = T,
                  by = c('bg_gene_id', 'bg_gene_name'))
  
  result[, intersectionPercTarget := 100 * intersectionWidth / target_len]
  result[, intersectionPercBg := 100 * intersectionWidth / bg_len]
  
  if (select == 'max_target') {
    result <- result[order(-intersectionPercTarget)][,.SD[1],
                                                     by = .(target_gene_id, 
                                                            target_gene_name)]
  }
  if (select == 'min_target') {
    result <- result[order(intersectionPercTarget)][,.SD[1],
                                                    by = .(target_gene_id, 
                                                           target_gene_name)]
  }
  if (select == 'max_background') {
    result <- result[order(-intersectionPercBg)][,.SD[1],
                                                 by = .(bg_gene_id, 
                                                        bg_gene_name)]
  }
  if (select == 'min_background') {
    result <- result[order(intersectionPercBg)][,.SD[1],
                                                by = .(bg_gene_id, 
                                                       bg_gene_name)]
  }
  setcolorder(result, c('target_gene_id', 'target_gene_name', 'bg_gene_id', 
                        'bg_gene_name', 'intersectionWidth', 'target_len', 
                        'bg_len', 'intersectionPercTarget',
                        'intersectionPercBg'))
  result 
}

#' cluster_genes_by_location
#' @description Clusters genes by physical overlap based on matrix of 
#'              intersection percentages.
#' @author Maria Litovchenko
#' @param gene_intersect_perc_matr data frame with genes as row names and 
#'        column names values in cells = intersection percentages. Values 
#'        should range from 0 to 100.
#' @param threshold_cut cutoff on distance between genes
#' @return data table with columns gene_id, group_id, n_genes, min_ovrl_perc -
#'         where group_id is group id, n_genes - number of genes in a group,
#'         min_ovrl_perc - minimal overlap (in %, from 0 to 100) between genes
#'         in a group.
cluster_genes_by_location <- function(gene_intersect_perc_matr, threshold_cut) {
  # convert intersecion percentage into distance
  gene_dist <- abs(gene_intersect_perc_matr - 100)
  # if NA is percentage of intersection, will set 100 to it
  gene_dist[is.na(gene_dist)] <- 100
  
  # group genes based on that distance. As cutoff use the min of thresholds
  gene_tree <- hclust(as.dist(as.matrix(gene_dist)))
  result <- cutree(gene_tree, h = threshold_cut)
  result <- data.table(gene_id = names(result), group_id = result)
  result <- result[, n_genes := .N, by = group_id]
  result <- result[n_genes > 1]
  result <- result[order(group_id, gene_id)]
  
  result
}

#' merge_mcols
#' @description Merges meta information of genomic ranges based on group_id.
#' @author Maria Litovchenko
#' @param inGR input genomic ranges with any columns + group_id. Merge of mcols
#'        will be performed by group_id.
#' @return data table with the same columns + group_id.
merge_mcols <- function(inGR) {
  result <- as.data.table(mcols(inGR))
  result <- split(result, by = 'group_id')
  result <- lapply(result, 
                   function(x) apply(x, 2, 
                                     function(y) paste0(unique(y),
                                                        collapse = '__and__')))
  result <- lapply(result, function(x) as.data.table(t(x)))
  result <- do.call(rbind, result)
  result[, group_id := as.character(group_id)]
  setkey(result, 'group_id')
  result
}

#' repeat_groupid
#' @description Repeats names of GRanges list number of times corresponding to
#'              number of ranges in each element of the list.
#' @author Maria Litovchenko
#' @param gr_list named GRanges list
#' @return vector of strings
repeat_groupid <- function(gr_list) {
  n_regs_per_item <- sapply(gr_list, length)
  names(n_regs_per_item) <- as.character(names(n_regs_per_item))
  result <- unlist(lapply(names(gr_list),
                          function(x) rep(x, n_regs_per_item[x])))
  result
}

#' preserve_strand_info
#' @description 
#' @author Maria Litovchenko
#' @param inGR
#' @return 
preserve_strand_info <- function(inGR) {
  result <- lapply(strand(inGR), function(x) unique(as.character(x)))
  result <- sapply(result, 
                   function(x) if (length(x) == 1) {return(x)} 
                   else {return('*')})
  result
}

#' self_intersect_list
#' @description Returns intersections of all elements in a list, per element.
#'              (Intersection is also done per element)
#' @author Maria Litovchenko
#' @param inGR input GRanges with gene_id, gene_name, rCode, group_id as 
#'             additional info.
#' @param ignore.strand boolean, whatever or not strand should be ignored
#' @return GRanges
self_intersect_list <- function(inGR, ignore.strand = T) {
  # we'll find intersection though disjoin + findOverlaps
  disjoined <- split(inGR, inGR$group_id)
  disjoined <- disjoin(disjoined, ignore.strand = ignore.strand)
  # add group_id to disjoined 
  disjoin_groupid <- repeat_groupid(disjoined)
  disjoined <- unlist(disjoined)
  disjoined$group_id <- disjoin_groupid
  
  # overlap disjoined with original regions and select ones which overlap all
  # genes
  ovrl <- findOverlaps(disjoined, inGR, ignore.strand = ignore.strand)
  ovrl <- as.data.table(ovrl)
  colnames(ovrl) <- c('disjoinedIdx', 'inputIdx')
  ovrl[, disjoinedGroup := disjoined[disjoinedIdx]$group_id]
  ovrl[, originalGroup := inGR[inputIdx]$group_id]
  ovrl <- ovrl[disjoinedGroup == originalGroup]
  ovrl <- ovrl[,.(disjoinedIdx, inputIdx, originalGroup)]
  setnames(ovrl, 'originalGroup', 'group_id')
  ovrl[, original_gene_id := inGR[inputIdx]$gene_id]
  ovrl <- merge(ovrl, as.data.table(mcols(inGR))[,.(length(unique(gene_id))),
                                                 by = group_id],
                all.x = T, by = 'group_id')
  ovrl[, n_genes := .(length(unique(original_gene_id))), by = disjoinedIdx]
  ovrl <- ovrl[n_genes == V1]
  
  disjoined <- disjoined[unique(ovrl$disjoinedIdx)]
  disjoined_mcols <- merge_mcols(inGR)
  mcols(disjoined) <- disjoined_mcols[disjoined$group_id]
  
  disjoined
}

#' union_regions_by_genes
#' @description Returns union of all elements in a list, per element.
#' @author Maria Litovchenko
#' @param inGR input GRanges with gene_id, gene_name, rCode, group_id as 
#'             additional info.
#' @param ignore.strand boolean, whatever or not strand should be ignored
#' @return GRanges
union_regions_by_genes <- function(inGR, ignore.strand = T) {
  result <- inGR[order(inGR$group_id)]
  
  # first of all, merge information in mcols
  mcols_upd <- merge_mcols(result)
  setkey(mcols_upd, 'group_id')
  
  # perform the actual merge of genomic ranges
  result <- split(result, result$group_id)
  # ... 1. if strands are all the same (i.e. all + or -), keep strand
  result_strand <- preserve_strand_info(result)
  names(result_strand) <- as.character(names(result_strand))
  # reduce
  result <- reduce(result, ignore.strand = ignore.strand)
  # add group_id
  reduced_groupid <- repeat_groupid(result)
  result <- unlist(result)
  result$group_id <- reduced_groupid
  # add strand
  strand(result) <- result_strand[as.character(result$group_id)]
  # add mcols
  mcols(result) <- mcols_upd[result$group_id]
  
  result
}

#' process_overlapping_gregions
#' @description Merges or intersects overlapping genomic regions of interest,
#'              i.e. promoters, in PER GENE ID manner.
#' @author Maria Litovchenko
#' @param inGR GRanges object holding genomic regions of interest
#' @param ignore.strand boolean, indicates, whatever or not strand should be
#'        ignored in overlapping & merging & intersecting operations
#' @param unionTreshold double, minimal % of overlap between regions 
#'        corresponding to separate genes in order for the regions to be merged
#'        (merged = reduced, no bp left behind)
#' @param intersectTreshold double, minimal % of overlap between regions 
#'        corresponding to separate genes in order for the regions to be 
#'        overlapped (base pairs, specific to each gene which do not overlap, 
#'        are discarded.)
#' @return GRanges object.
#' @note If unionTreshold < intersectTreshold then genes with similarity
#'       value (% overlap) in range unionTreshold - intersectTreshold will be
#'       merged, and genes with similarity value in range >= intersectTreshold
#'       will be intersected and the other way around.
process_overlapping_gregions <- function(inGR, ignore.strand = T, 
                                         unionTreshold = NA, 
                                         intersectTreshold = NA) {
  if (is.na(unionTreshold) & is.na(intersectTreshold)) {
    message('[', Sys.time(), '] Both unionTreshold and intersectTreshold are ',
            'NULL, return original GRanges.')
    result <- copy()
    return(inGR)
  }
  if (is.na(unionTreshold)) {
    unionTreshold <- 100
  }
  if (is.na(intersectTreshold)) {
    intersectTreshold <- 100
  }
  if (unionTreshold == intersectTreshold) {
    stop('[', Sys.time(), '] Merge and intersect thresholds can not be the ',
         'same value.')
  }
  # sort out thresholds
  tresholds <- sort(c('merge' = unionTreshold, 
                      'intersect' = intersectTreshold))
  
  # calculate percentage of all bases corresponding to a gene in overlap with
  # the other gene. If gene pair is not present in the table, it means that
  # they do not overlap.
  perc_ovrl <- calcPercOverlap_byGeneID(inGR, inGR, remove.self = F,
                                        ignore.strand = ignore.strand, 
                                        select = 'all')
  # compute minimum of percentage overlap between "target" and "background" 
  perc_ovrl[, minIntersect := apply(perc_ovrl[,.(intersectionPercTarget,
                                                 intersectionPercBg)], 1, min)]
  if (max(perc_ovrl[target_gene_id != bg_gene_id]$minIntersect,
          na.rm = T) < tresholds[1]) {
    return(inGR)
  }
  
  # use minIntersect as distance between genes and cluster genes based on it
  gene_dist <- as.data.frame(dcast(perc_ovrl, target_gene_id  ~ bg_gene_id,
                                   value.var = "minIntersect"))
  rownames(gene_dist) <- gene_dist$target_gene_id
  gene_dist$target_gene_id <- NULL
  gene_groups <- cluster_genes_by_location(gene_dist,
                                           threshold_cut = tresholds[1])
  # add minimal percentage of overlap per gene group
  gene_groups[, min_ovrl_perc := -1]
  for (groupID in unique(gene_groups$group_id)) {
    group_perc_ovrl <- perc_ovrl[target_gene_id %in% 
                                   gene_groups[group_id == groupID]$gene_id &
                                   bg_gene_id %in% 
                                   gene_groups[group_id == groupID]$gene_id & 
                                   target_gene_id != bg_gene_id]
    group_perc_ovrl <- min(group_perc_ovrl$minIntersect)
    gene_groups[group_id == groupID]$min_ovrl_perc <- group_perc_ovrl
  }
  # add what an action to perform on that gene group
  gene_groups[, to_do := names(tresholds)[1]]
  gene_groups[min_ovrl_perc >= tresholds[2]]$to_do <- names(tresholds)[2]
  
  # do union
  united_gr <- GRanges()
  if (nrow(gene_groups[to_do == 'merge']) != 0) {
    gr_to_unite <- inGR[inGR$gene_id %in% 
                          gene_groups[to_do == 'merge']$gene_id]
    # ... add index in inGR to ensure that entries are not mixed in
    gr_to_unite$idx <- 1:length(gr_to_unite)
    # ... add info about groups
    mcols(gr_to_unite) <- merge(as.data.table(mcols(gr_to_unite)), 
                                gene_groups[,.(gene_id, group_id)], 
                                by = 'gene_id', all.x = T)[order(idx)]
    gr_to_unite$idx <- NULL
    united_gr <- union_regions_by_genes(gr_to_unite)
  }
  # do reduction 
  intersected_gr <- GRanges()
  if (nrow(gene_groups[to_do == 'intersect']) != 0) {
    gr_to_intersect <- inGR[inGR$gene_id %in% 
                              gene_groups[to_do == 'intersect']$gene_id]
    # ... add index in inGR to ensure that entries are not mixed in
    gr_to_intersect$idx <- 1:length(gr_to_intersect)
    # ... add info about groups
    mcols(gr_to_intersect) <- merge(as.data.table(mcols(gr_to_intersect)),
                                    gene_groups[,.(gene_id, group_id)], 
                                    by = 'gene_id', all.x = T)[order(idx)]
    gr_to_intersect$idx <- NULL
    intersected_gr <- self_intersect_list(gr_to_intersect)
  }
  
  result <- c(inGR[!inGR$gene_id %in% gene_groups$gene_id], intersected_gr, 
              united_gr)
  result$group_id <- NULL
  result <- sort(result)
  result
}

# FUNCTIONS: informing user ---------------------------------------------------
#' summurizeGenomicRegions
#' @description Puts basic statistics about genomic ranges to data table 
#' @author Maria Litovchenko
#' @gr GRanges object to sum up, columns rCode and gene_id are required
#' @return data table with columns gr_id, n_regions, n_genes, min_reg_len, 
#' mean_reg_len, max_reg_len, min_gene_len, mean_gene_len, max_gene_len,
#' total_len, perc_genome
summurizeGenomicRegions <- function(gr, gr_codename = NULL) {
  gr_lens <- width(gr)
  gene_lens <- split(gr, gr$gene_id)
  gene_lens <- sum(width(gene_lens))
  result <- data.table(gr_id = gr_codename,
                       rCodes = paste(unique(gr$rCode), collapse = ','),
                       n_regions = length(gr), 
                       n_genes = length(unique(gr$gene_id)), 
                       min_reg_len = min(gr_lens), 
                       mean_reg_len = round(mean(gr_lens)),
                       sd_reg_len = sd(gr_lens),
                       max_reg_len = max(gr_lens), 
                       min_gene_len = min(gene_lens),
                       mean_gene_len = round(mean(gene_lens)),
                       sd_gene_len = sd(gene_lens),
                       max_gene_len = max(gene_lens), 
                       total_len = sum(gr_lens))
  result[, perc_genome := round(100 * total_len / (3 * 10^9), 2)]
  result
}
# FUNCTIONS: writing ---------------------------------------------------------
#' gRangesToBed12
#' @description Converts GRanges object to bed12-formatted data table
#' @author Maria Litovchenko
#' @param gr GRanges object with additional columns gene_id, gene_name, rCode,
#' gr_id. gr_id should have one and only one value.
#' @param genomeversion genome version of regions' genome
#' @return data table with columns chr, start, end, name, score, strand,
#' thickStart (same value as start), thickEnd (same value as end), itemRgb 
#' (always 0), blockCount, blockSizes, blockStarts. name will be constructed as
#' genomeversion-gr_id-gene_id-gene_name. score will be 1000, unless rCode is
#' CDS, 3primeUTR, 5primeUTR.
#' @note https://genome.ucsc.edu/FAQ/FAQformat.html#format1
gRangesToBed12 <- function(gr, genomeversion) {
  grList <- unlist(GRangesList(gr))
  if (length(unique(grList$gr_id)) != 1) {
    stop('[', Sys.time(), '] Can not handle > 1 gr_id')
  }
  
  # get data table with scores
  scoresDT <- as.data.table(mcols(grList)[, c('gene_id', 'gene_name', 
                                              'rCode', 'gr_id')])
  scoresDT <- scoresDT[!duplicated(scoresDT)]
  scoresDT[, score := ifelse(rCode %in% c('3primeUTR', '5primeUTR', 'CDS'), 
                             0, 1000)]
  scoresDT[, rCode := NULL]
  
  # re-split by gene
  grList <- split(grList, grList$gene_id)
  
  # construct result
  result <- data.table(chr = unlist(unique(seqnames(grList))), 
                       start = min(start(grList)), end = max(end(grList)),
                       strand = unlist(unique(strand(grList))),
                       gene_id = names(grList))
  result[, blocksStart := paste(start(grList) - min(start(grList)), 
                                collapse = ',')]
  result[, blocksStart := paste0(blocksStart, ',')]
  result[, blockSizes := paste(end(grList) - start(grList) + 1, 
                               collapse = ',')]
  result[, blockSizes := paste0(blockSizes, ',')]
  result[, nBlocks := sapply(start(grList) - min(start(grList)), length)]
  result <- merge(result, scoresDT, by = 'gene_id')
  result[, gene_id := apply(result[,.(gr_id, gene_id, gene_name)], 1,
                            paste, collapse = '--')]
  result[, gene_id := paste0(genomeversion, '--', gene_id)]
  result[, itemRgb := 0]
  result <- result[,.(chr, start, end, gene_id, score, strand, start, 
                      end, itemRgb, nBlocks, blockSizes, blocksStart)]
  colnames(result) <- c('chr', 'start', 'end', 'name', 'score', 'strand',
                        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                        'blockSizes', 'blockStarts')
  
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'create_input_genomic_regions_files.R')

analysisHelp <- paste('Path to inventory table containing details of the',
                      'future analysis to be conducted. Minimal columns:',
                      'tumor_subtype,', 'software,', 'gr_id,', 'gr_code,', 
                      'gr_file,', 'gr_upstr,', 'gr_downstr,', 'gr_genome,', 
                      'gr_excl_id,', 'gr_excl_code,', 'gr_excl_file,',
                      'gr_excl_upstr,', 'gr_excl_downstr,', 'gr_excl_genome,',
                      'blacklisted_codes.')
parser$add_argument("-a", "--inventory_analysis", required = T, 
                    type = 'character', help = analysisHelp)

blackListHelp <- paste('Path to patientsInv table containing details of the',
                       'black&white lists to use. Minimal columns:',
                       'list_name,', 'file_path,', 'file_genome,', 
                       'file_type,', 'score_column,', 'min_value,',
                       'max_value.')
parser$add_argument("-b", "--inventory_blacklisted", required = F, 
                    type = 'character', default = NULL, help = blackListHelp)

targetGenomePathHelp <- paste('Path to the fasta file, genome version of',
                              'which should be same as ',
                              '--target_genome_version.')
parser$add_argument("-f", "--target_genome_path", required = T,
                    default = '', type = 'character',
                    help = targetGenomePathHelp)

targetGenomeChrHelp <- paste('Path to the tab-separated file containing ', 
                             'chromosomal lengths of the target genome. ',
                             'Must have 3 columns: chr, start(1) and length ',
                             'of the chromosome. No header.')
parser$add_argument("-cl", "--target_genome_chr_len", required = F,
                    default = '', type = 'character', 
                    help = targetGenomeChrHelp)

targetGenomeHelp <- paste("Genome version, i.e. hg19, to which input variants",
                          "files for software should be brought.",
                          "Default: hg19.")
parser$add_argument("-g", "--target_genome_version", required = F,
                    default = 'hg19', type = 'character',
                    help = targetGenomeHelp)

chainHelp <- paste('Path to chain file in case genome version of mutations is',
                   'not the same as --target_genome_version')
parser$add_argument("-l", "--chain", required = F, default = NULL,
                    type = 'character', help = chainHelp)

strandHelp <- paste('Boolean, indicating, whatever or not strand information ',
                    'should be taken into account while removing regions to ',
                    'exclude from the target regions. Default: T')
parser$add_argument("-s", "--ignore_strand", required = F, default = 'T',
                    choices = c('T', 'F'), type = 'character', 
                    help = strandHelp)

minRegLenHelp <- paste('Minimal region length. Regions with length less than',
                       'that will be removed.')
parser$add_argument("-m", "--min_reg_len", required = T, default = 5,
                    type = 'integer', help = minRegLenHelp)

coresHelp <- 'How many cores the script should use. Default: 1.'
parser$add_argument("-n", "--cores", required = F, type = 'integer', 
                    default = 1, help = coresHelp)

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output folder")

args <- parser$parse_args()
args$ignore_strand <- as.logical(args$ignore_strand)
check_input_arguments(args, outputType = 'folder')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

dir.create(args$output, recursive = T)

# Test arguments --------------------------------------------------------------
# args <- list(inventory_analysis = '../data/inventory/inventory_analysis_tcga.csv',
#              inventory_blacklisted = '../data/inventory/inventory_blacklist_tcga.csv', 
#              target_genome_version = 'hg19',
#              target_genome_path = '../data/assets/reference_genome/Homo_sapiens_assembly19.fasta',
#              target_genome_chr_len = '../data/assets/reference_genome/Homo_sapiens_assembly19.chrom.sizes.bed',
#              chain = '../data/assets/reference_genome/hg38ToHg19.over.chain',
#              ignore_strand = T, min_reg_len = 5, cores = 2, output = 'test/')

# READ inventories ------------------------------------------------------------
analysisInv <- readAnalysisInventory(args$inventory_analysis, args$cores)
# check that if gr_genome in analysisInv is not the same as target one
# chain file is submitted
if (any(!unique(c(analysisInv$gr_genome, analysisInv$gr_excl_genome)) %in% 
        args$target_genome_version) & is.null(args$chain)) {
  stop('[', Sys.time(), '] genome version of some genomic regions is not ',
       'the same as --target_genome_version, but no chain file is provided')
}
message('[', Sys.time(), '] Read --inventory_analysis: ', 
        args$inventory_analysis)

if (!is.null(args$inventory_blacklisted)) {
  bwInv <- readBlacklistInventory(args$inventory_blacklisted, args$cores)
  if (is.null(bwInv)) {
    args$inventory_blacklisted <- NULL
    message('[', Sys.time(), '] Black & white lists inventory file ',
            '(--inventory_blacklisted) was empty. Proceeding without it.')
  } else {
    message('[', Sys.time(), '] Read --inventory_blacklisted: ', 
            args$inventory_blacklisted)
  }
}

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

# read in chain file
chain <- NULL
if (!is.null(args$chain)) {
  chain <- import.chain(args$chain)
}

# determine, if chromosomal names in fasta file have 'chr' in front of them or
# not. This will help us to output variants with chromosomal names in the
# proper format matching reference genome
outChrStyle <- getSeqlevelsStyle(args$target_genome_path)

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Determine, if liftover is needed at the reading of files --------------------
# Basically, lift over to target version will be done just after file(s) will
# be loaded, unless all submitted files are on non-target genome version. It 
# may not be the most optimal strategy, but it simplifies the logic a lot.
filesGenomeVersion <- rbind(setNames(analysisInv[, .(gr_file, gr_genome)], 
                                     c('filePath', 'genomeVersion')),
                            setNames(analysisInv[, .(gr_excl_file, 
                                                     gr_excl_genome)], 
                                     c('filePath', 'genomeVersion')))
if (!is.null(bwInv)) {
  filesGenomeVersion <- rbind(filesGenomeVersion,
                              setNames(bwInv[,.(file_path, file_genome)],
                                       c('filePath', 'genomeVersion')))
}
filesGenomeVersion <- unique(filesGenomeVersion)
filesGenomeVersion <- filesGenomeVersion[!is.na(filePath)]
if (all(filesGenomeVersion$genomeVersion == args$target_genome_version)) {
  liftOverIsNeeded <- 'not'
} else {
  if (all(filesGenomeVersion$genomeVersion != args$target_genome_version)) {
    liftOverIsNeeded <- 'end'
  } else {
    liftOverIsNeeded <- 'start'
  }
}

# EXTRACT target and exclusion genomic regions (& liftover, if needed) --------
# get unique combos of files and gr_code so that if in one analysis a region,
# i.e. CDS is target, and in the other it is excluded we wouldn't extract it 
# twice.
grIDtoFileDT <- rbind(setNames(analysisInv[,.(gr_code, gr_file, gr_upstr, 
                                              gr_downstr, gr_genome)],
                               c('gr_code', 'filePath', 'upLen', 'downLen', 
                                 'gr_genome')),
                      setNames(analysisInv[,.(gr_excl_code, gr_excl_file, 
                                              gr_excl_upstr, gr_excl_downstr,
                                              gr_excl_genome)],
                               c('gr_code', 'filePath', 'upLen', 'downLen', 
                                 'gr_genome')))
grIDtoFileDT <- unique(grIDtoFileDT)
grIDtoFileDT <- grIDtoFileDT[!is.na(filePath)]
# add file extension = file type
grIDtoFileDT[, file_type := gsub('.*[.]', '', gsub('.gz$', '', filePath))]
# add unique id for regions
grIDtoFileDT[, id := apply(grIDtoFileDT[,.(gr_code, upLen, downLen)], 1, 
                           paste, collapse = '-')]
grIDtoFileDT[, id := gsub(' ', '', id)]

# split grIDtoFileDT by file path so that we would read each file just once.
grIDtoFileDT <- split(grIDtoFileDT, grIDtoFileDT$filePath)
# read files and extract regions sequentially so that the process is not too
# hard on memory
allGenomicRegions <- list()
for (fileIdx in 1:length(grIDtoFileDT)) {
  # check, that file needs liftover now
  fileNeedsLiftOver <- unique(grIDtoFileDT[[fileIdx]]$gr_genome)
  fileNeedsLiftOver <- fileNeedsLiftOver != args$target_genome_version
  fileNeedsLiftOver <- fileNeedsLiftOver & liftOverIsNeeded == 'start'
  # read file
  genomeAnno <- switch(grIDtoFileDT[[fileIdx]]$file_type[1],
                       'gtf' = gtfToTxDB(gtfPath = grIDtoFileDT[[fileIdx]]$filePath[1],
                                         chrStyle = outChrStyle, 
                                         acceptedChrCodes = acceptedChrNames, 
                                         doLiftOver = fileNeedsLiftOver,
                                         chainObj = chain),
                       'bed' = bedToGR(grIDtoFileDT[[fileIdx]]$filePath[1],
                                       chrStyle = outChrStyle, 
                                       acceptedChrCodes = acceptedChrNames,
                                       doLiftOver = fileNeedsLiftOver,
                                       chainObj = chain))
  # extract regions
  regs <- list()
  for (regIdx in 1:nrow(grIDtoFileDT[[fileIdx]])) {
    regInfo <- unlist(grIDtoFileDT[[fileIdx]][regIdx, ])
    
    if (grIDtoFileDT[[fileIdx]]$file_type[1] == 'gtf') {
      regs[[regIdx]] <- extractRegionsByCode_txdb(rCode = regInfo['gr_code'], 
                                                  txdbObj = genomeAnno$gtfTxdb,
                                                  geneBiotypeDT = genomeAnno$geneMap,
                                                  upstr = regInfo['upLen'],
                                                  downstr = regInfo['downLen'])
    } else {
      regs[[regIdx]] <- extractRegionsByCode_gr(rCode = regInfo['gr_code'], 
                                                grObj = genomeAnno,
                                                upstr = regInfo['upLen'],
                                                downstr = regInfo['downLen'])
      
    }
    # add source column - where the region came from 
    if (length(regs[[regIdx]]) != 0) {
      regs[[regIdx]] <- unlist(regs[[regIdx]])
      regs[[regIdx]]$source <- basename(unique(grIDtoFileDT[[fileIdx]]$filePath))
      regs[[regIdx]] <- split(regs[[regIdx]], regs[[regIdx]]$gene_id)
    }
  }
  names(regs) <- grIDtoFileDT[[fileIdx]]$id
  allGenomicRegions[[fileIdx]] <- regs
  
  rm(genomeAnno)
  rm(regs)
  suppressMessages(suppressWarnings(gc()))
}
names(allGenomicRegions) <- names(grIDtoFileDT)

# SUBSTRACT exclusion regions from target regions -----------------------------
analysisInv <- split(analysisInv, analysisInv$gr_id)
cleanTargets <- lapply(analysisInv, getCleanRegions, allGenomicRegions,
                       args$ignore_strand)
# add gr_id column cleanTargets
for (i in 1:length(cleanTargets)) {
  cleanTargets[[i]]$gr_id <- names(cleanTargets)[i]
}

rm(allGenomicRegions)
suppressMessages(suppressWarnings(gc()))

# Liftover (if all files were on other than target genome) --------------------
# ... if it was not the case, liftover was performed at the reading of all 
# files.
if (liftOverIsNeeded == 'end') {
  cleanTargets <- lapply(cleanTargets,
                         function(x) unlist(liftOverGenomicRegion(x, chain)))
}

# REMOVE black-/include only white- listed regions from/into target regions----
# read all black- and white- listed files
if (!is.null(bwInv)) {
  setkey(bwInv, 'list_name')
  bwGR <- lapply(bwInv$list_name, 
                 function(x) readInAndFilterBWregions(bwInv[x]$file_path, 
                                                      chrStyle = outChrStyle,
                                                      bwScoreCol = bwInv[x]$score_column, 
                                                      bwScoreMin = bwInv[x]$min_value,
                                                      bwScoreMax = bwInv[x]$max_value))
  names(bwGR) <- bwInv$list_name
  for (bwCode in names(bwGR)) {
    bwGR[[bwCode]]$rCode <- bwCode
  }
  
  # perform removal of black- regions and filter in white ones
  for (grID in names(analysisInv)) {
    removeBW <- unlist(analysisInv[[grID]]$blacklisted_codes[1])
    removeBW <- intersect(removeBW, names(bwGR))
    if (any(!is.na(removeBW))) {
      # do it with loop, because where could be both black- and white- listed 
      # regions
      for (bwCode in removeBW) {
        cleanTargets[[grID]] <- excludeFromTarget(targetGR = cleanTargets[[grID]],  
                                                  excludeGR = bwGR[[bwCode]],
                                                  type = bwInv[bwCode]$file_type, 
                                                  ignore.strand = args$ignore_strand)
      }
    }
  }
  
  rm(bwGR)
  suppressMessages(suppressWarnings(gc()))
}
# Perform UNION or INTERSECTION of overlapping regions of different genes -----
doUnionIntersect <- sapply(analysisInv, 
                           function(x) any(c(!is.na(x$union_percentage), 
                                             c(!is.na(x$intersect_percentage)))))
doUnionIntersect <- names(doUnionIntersect[doUnionIntersect == T])

if (length(doUnionIntersect) != 0) {
  message('[', Sys.time(), '] union_percentage and/or intersect_percentage ',
          'was given in analysis inventory table for gr_id(s): ', 
          paste(doUnionIntersect, collapse = ', '), '. Will perform union or ',
          'intersection for noncoding regions. It will not be performed ',
          'for coding regions due to presence of internal structure (codons) ',
          'inside coding regions.')
  for (gr_id_i in doUnionIntersect) {
    message('[', Sys.time(), '] Started performing union or intersection of ',
            'overlapping regions of different genes on ', gr_id_i)
    cleanTargets[[gr_id_i]] <- process_overlapping_gregions(inGR = cleanTargets[[gr_id_i]], 
                                                            ignore.strand = args$ignore_strand,
                                                            unionTreshold = unique(analysisInv[[gr_id_i]]$union_percentage),
                                                            intersectTreshold = unique(analysisInv[[gr_id_i]]$intersect_percentage))
    
    message('[', Sys.time(), '] Finished performing union or intersection of ',
            'overlapping regions of different genes on ', gr_id_i)
  }
}

suppressMessages(suppressWarnings(gc()))

# REMOVE regions outside of chromosomal lengths -------------------------------
if (!is.null(args$target_genome_chr_len)) {
  message('[', Sys.time(), '] --target_genome_chr_len file is given. Will ',
          'remove any regions surpassing chromosomal length')
  
  chrLensGR <- fread(args$target_genome_chr_len, stringsAsFactors = F, 
                     header = F,  col.names = c('chr', 'start', 'end'))
  chrLensGR <- makeGRangesFromDataFrame(chrLensGR)
  
  # perform removal of black- regions and filter in white ones
  for (grID in names(analysisInv)) { 
    cleanTargets[[grID]] <- excludeFromTarget(cleanTargets[[grID]], chrLensGR,
                                              'white', args$ignore_strand)
  }
  
  message('[', Sys.time(), '] Removed any regions surpassing chromosomal ',
          'length')
}

# REMOVE too short regions ----------------------------------------------------
if (args$min_reg_len > 1) {
  for (grID in names(cleanTargets)) {
    nbefore <- length(cleanTargets[[grID]])
    longEnough <- width(cleanTargets[[grID]]) >= args$min_reg_len
    cleanTargets[[grID]] <- cleanTargets[[grID]][longEnough]
    nafter <- length(cleanTargets[[grID]])
    message('[', Sys.time(), '] Removed ', nbefore - nafter, '(',
            round(100 * (nbefore - nafter) / nbefore, 2), '%) regions from ', 
            grID, ' due length < ', args$min_reg_len)
    
    suppressMessages(suppressWarnings(gc()))
  }
}

# Print summary of extracted regions to screen --------------------------------
cleanTargsSum <- lapply(names(cleanTargets), 
                        function(x) summurizeGenomicRegions(cleanTargets[[x]],
                                                            x))
cleanTargsSum <- do.call(rbind, cleanTargsSum)

# to print
summaryToPrint <- copy(cleanTargsSum)
summaryToPrint[, `N.regs(k)` := n_regions / 1000]
summaryToPrint[, `N.genes(k)` := n_genes / 1000]
summaryToPrint[, `Mean.per.gene,kb` := mean_gene_len / 1000]
summaryToPrint[, `Max.per.gene,kb` := max_gene_len / 1000]
summaryToPrint[, `Total.len,Mb` := total_len / 10^6]
setnames(summaryToPrint, 
         c('min_reg_len', 'mean_reg_len', 'max_reg_len', 'min_gene_len',
           'perc_genome'),
         c('Min.reg,bp', 'Mean.reg,bp', 'Max.reg,bp', 'Min.per.gene,bp', 
           '% genome'))
summaryToPrint <- summaryToPrint[, c("gr_id", 'rCodes', "N.regs(k)", 
                                     "N.genes(k)", "Min.reg,bp", "Mean.reg,bp",
                                     "Max.reg,bp", "Min.per.gene,bp", 
                                     "Mean.per.gene,kb", "Max.per.gene,kb",
                                     "Total.len,Mb", "% genome")]
summaryToPrint <- summaryToPrint[order(-`% genome`)]
message('[', Sys.time(), '] Overview of the extracted genomic regions: ')
message(paste0(capture.output(knitr::kable(summaryToPrint,
                                           format = "markdown")),
               collapse = '\n'))

# [SAVE] BED12 file with all genomic regions to scan --------------------------
# these files (one per tumor type) will be very useful further down the 
# pipeline for overlapping mutations and regions.
bed12Regs <- lapply(cleanTargets, gRangesToBed12, args$target_genome_version)

# output to file(s)
outBed12inv <- lapply(analysisInv, function(x) x[,.(tumor_subtype, gr_id)])
outBed12inv <- lapply(outBed12inv, function(x) x[!duplicated(x)])
outBed12inv <- do.call(rbind, outBed12inv)
setkey(outBed12inv, 'tumor_subtype')

for (tumSubt in unique(outBed12inv$tumor_subtype)) {
  outfile <- paste0(args$output, '/inputGR-', tumSubt, '-', 
                    args$target_genome_version, '.bed')
  if (file.exists(outfile)) {
    file.remove(outfile)
  }
  for (gr_id_idx in outBed12inv[tumSubt]$gr_id) {
    regsForOut <- copy(bed12Regs[[gr_id_idx]])
    regsForOut[, chr := factor(chr, orderChromosomes(chr))]
    regsForOut <- regsForOut[order(chr, start)]
    write.table(regsForOut, outfile, col.names = F, row.names = F, quote = F,
                sep = '\t', append = T)
    rm(regsForOut)
  }
}

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')