#!/usr/bin/env Rscript
# FILE: calculate_mutation_rates.R --------------------------------------------
#
# DESCRIPTION: Rscript to calculate local [hyper]mutation rate to use it 
#              further in filtering
#
# USAGE: Rscript --vanilla calculate_mutation_rates.R [options]
#
# OPTIONS: Run in terminal: 
#          Rscript --vanilla calculate_mutation_rates.R -h 
#          to see available options and help message. Or just scroll down to
#          parser definition.
#
# REQUIREMENTS: argparse, data.table, GenomicRanges, rtracklayer
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  22.01.2022
# REVISION: 22.09.2023

box::use(./custom_functions[...])

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges)))
suppressWarnings(suppressPackageStartupMessages(library(maftools)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(VariantAnnotation)))
options(scipen = 999)

# FUNCTIONS: reading ---------------------------------------------------------
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

# FUNCTIONS: common -----------------------------------------------------------
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

#' getUniqQuantileBreaks
#' @description Creates data table containing unique breaks and their labels
#' based on data vector
#' @param x numeric vector
#' @param nQuants number of quantiles
#' @return data table with columns topBound (use it as breaks) and quant (use 
#' it as labels)
getUniqQuantileBreaks <- function(x, nQuants = 100) {
  # add quantile of local mutation rate as well
  quantCuts <- data.table(quant = seq(0, 1, length.out = nQuants), 
                          topBound = quantile(x, seq(0, 1,
                                                     length.out = nQuants)))
  quantCuts <- quantCuts[,.(quant = max(quant)), by = topBound]
  quantCuts[topBound == max(topBound)]$topBound <- Inf
  quantCuts[, quant := round(100 * quant, 2)]
  quantCuts
}

# FUNCTIONS: matching variants to regions -------------------------------------
#' matchVariantsToRegions
#' @description Matches mutations to genomic regions of interest. It is 
#'              possible to submit mutations either as maf file via mafPath 
#'              argument or as data table via varsDT argument. Same goes 
#'              towards genomic regions, which can be submitted either though
#'              file (bed12Path) or genomic regions (GRs)
#' @author Maria Litovchenko
#' @param mafPath path to MAF file with mutations
#' @param varsDT data table with essential columns chr, start, end
#' @param bed12Path path to BED12 file with genomic regions
#' @param GRs genomic ranges object with gr_name column. Other annotation 
#'            columns may be present too.
#' @param nCores number of cores to use
#' @return data table with columns gr_name + all annotation columns from 
#'         genomic regions file + all columns from variants data table
matchVariantsToRegions <- function(mafPath = NULL, varsDT = NULL,
                                   bed12Path = NULL, GRs = NULL,
                                   nCores = 1) {
  if ((is.null(bed12Path) & is.null(GRs)) | 
      (!is.null(bed12Path) & !is.null(GRs))) {
    stop('[', Sys.time(), '] matchVariantsToRegions: please submit either ',
         'bed12Path or GRs')
  }
  if ((is.null(mafPath) & is.null(varsDT)) | 
      (!is.null(mafPath) & !is.null(varsDT))) {
    stop('[', Sys.time(), '] matchVariantsToRegions: please submit either ',
         'mafPath or varsDT')
  }
  
  # read genomic regions from BED12
  if (!is.null(bed12Path)) { 
    GRs <- import(bed12Path)
    GRs <- unlist(GRangesList(blocks(GRs)))
    GR_parsed <- parseBED12regName(unique(names(GR)))
    setkey(GR_parsed, name)
    mcols(GR) <- GR_parsed[names(GR)]
    rm(GR_parsed)
    message('[', Sys.time(), '] Read ', bed12Path)
    message('[', Sys.time(), '] rtracklayer import function while reading ',
            'bed12 assumes that it is 0-based. Therefore it adds 1 to all ',
            'coordinates. We will correct it.')
    start(GRs) <- start(GRs) - 1
    end(GRs) <- end(GRs) - 1
  }
  
  # read variants from MAF
  if (!is.null(mafPath)) {
    message('[', Sys.time(), '] Started reading ', mafPath)
    colsToRead <- c('Tumor_Sample_Barcode', 'key', 'Chromosome',
                    'Start_Position', 'End_Position', 'Gene.refGene', 
                    'Variant_Classification', 'Variant_Type', 'mut_len', 
                    'struct_type', 'patient_tumor_subtype')
    newColNames <- c('participant_id', 'key', 'chr', 'start', 'end', 
                     'gene_name_var', 'var_class', 'var_type', 'mut_len',
                     'struct_type', 'patient_tumor_subtype')
    varsDT <- fread(mafPath, header = T, stringsAsFactors = F, 
                    nThread = nCores, select = colsToRead)
    setnames(varsDT, colsToRead, newColNames, skip_absent = T) 
    message('[', Sys.time(), '] Finished reading ', mafPath)
  }
  
  if (!'gr_name' %in% colnames(mcols(GRs))) {
    stop('[', Sys.time(), '] Genomic regions should have gr_name column')
  }
  
  # find overlaps
  # ... step 1: find overlaps for variants with length = 1 (SNPs and one bp ins
  #             or dels). These variants have to be whole inside the regions. 
  #             This will also exclude match variants to 2 neighboring 
  #             regions, i.e. CDS and ss. In MAF format, end = start + 1, for
  #             us, to simplify, let's make end = start
  varsDT_1bp <- varsDT[end - start <= 1]
  result_1bp <- data.table()
  if (nrow(varsDT_1bp) > 0) {
    varsDT_1bp[, end := start]
    varsGR_1bp <- makeGRangesFromDataFrame(varsDT_1bp)
    seqlevelsStyle(varsGR_1bp) <- suppressWarnings(seqlevelsStyle(GRs)[1])
    ovrl_1bp <- findOverlaps(varsGR_1bp, GRs, type = "within")
    ovrl_1bp <- as.data.table(ovrl_1bp)
    colnames(ovrl_1bp) <- c('varIdx', 'grIdx')
    # inform, if variants were not mapped to any region
    notMappedVars <- setdiff(1:length(varsGR_1bp), unique(ovrl_1bp$varIdx))
    if (length(notMappedVars) > 0) {
      notMappedChr <- table(seqnames(varsGR_1bp[notMappedVars]))
      notMappedChr <- as.data.table(notMappedChr)
      colnames(notMappedChr) <- c('chr', 
                                  'N unmapped to genomic regions variants')
      message('[', Sys.time(), '] WARNING: ', length(notMappedVars), '(',
              round(100 * length(notMappedVars)/length(varsGR_1bp), 2), '%)',
              'SNV variants were not mapped to genomic tiles. Distribution ',
              'over chromosomes: ')
      message(paste0(capture.output(knitr::kable(notMappedChr, 
                                                 format = "markdown")),
                     collapse = '\n'))
      message('[', Sys.time(), '] Possible reasons: variants lifted over ',
              'outside chromosomal length or genomic regions which do not ',
              'cover whole genome')
    }
    rm(varsGR_1bp)
    result_1bp <- as.data.table(mcols(GRs[ovrl_1bp$grIdx]))
    result_1bp <- cbind(result_1bp, 
                        varsDT_1bp[ovrl_1bp$varIdx, 
                                   !colnames(varsDT_1bp) %in% c('chr', 'start',
                                                                'end'),  
                                   with = F])
  }
  
  # ... step 2: find overlaps for variants with length > 1 (Indels and MNPs). 
  #             For these variants we will not require full overlap with a 
  #             region. A variant can overlap multiple regions even. The one(s)
  #             with the maximum percentage overlap will be selected. In MAF 
  #             format, end = start + length(ref) + 1, so we will take that + 1
  #             away.
  varsDT_multBp <- varsDT[end - start > 1]
  result_MiltiBp <- data.table()
  if (nrow(varsDT_multBp) > 0) {
    varsGR_multBp <- makeGRangesFromDataFrame(varsDT_multBp)
    end(varsGR_multBp) <- end(varsGR_multBp) - 1
    seqlevelsStyle(varsGR_multBp) <- seqlevelsStyle(GRs)[1]
    ovrl_multBp <- findOverlaps(varsGR_multBp, GRs, type = "any", 
                                select = 'all')
    ovrl_multBp <- as.data.table(ovrl_multBp)
    colnames(ovrl_multBp) <- c('varIdx', 'grIdx')
    # inform, if variants were not mapped to any region
    notMappedVars <- setdiff(1:length(varsGR_multBp), 
                             unique(ovrl_multBp$varIdx))
    if (length(notMappedVars) >0) {
      notMappedChr <- table(seqnames(varsGR_multBp[notMappedVars]))
      notMappedChr <- as.data.table(notMappedChr)
      colnames(notMappedChr) <- c('chr', 
                                  'N unmapped to genomic regions variants')
      message('[', Sys.time(), '] WARNING: ', length(notMappedVars),  '(',
              round(100 * length(notMappedVars)/length(varsGR_multBp), 2), 
              '%) MNV variants were not mapped to genomic tiles. ',
              'Distribution over chromosomes: ')
      message(paste0(capture.output(knitr::kable(notMappedChr, 
                                                 format = "markdown")),
                     collapse = '\n'))
      message('[', Sys.time(), '] Possible reasons: variants lifted over ',
              'outside chromosomal length or genomic regions which do not ',
              'cover whole genome')
    }
    rm(varsGR_multBp)
    result_MiltiBp <- as.data.table(mcols(GRs[ovrl_multBp$grIdx]))
    result_MiltiBp <- cbind(result_MiltiBp, 
                            varsDT_multBp[ovrl_multBp$varIdx, 
                                          !colnames(varsDT_multBp) %in% 
                                            c('chr', 'start', 'end'), 
                                          with = F]) 
  }

  result <- rbind(result_1bp, result_MiltiBp)
  rm(result_1bp, result_MiltiBp)
  result <- result[!duplicated(result)]
  result
}

#' checkGeneMatchBetweenVarAnnoAndGR
#' @description Checks match between gene name of genomic region and gene name
#'              which was assigned to mutation by variant annotation tool
#' @author Maria Litovchenko
#' @param mutsToGRmap data table with each row corresponding to one mutation. 
#'                    Columns gene_name, gene_name_var, gr_id and key are 
#'                    required
#' @param synonymsDT data table of gene name synonyms. Should have columns idx
#'                   and gene_name. Synonymous gene names will have the same 
#'                   idx.
#' @return mutsToGRmap with added column matchedGeneName
checkGeneMatchBetweenVarAnnoAndGR <- function(mutsToGRmap, synonymsDT = NULL) {
  result <- copy(mutsToGRmap)
  
  # sometimes in MAF files if a mutations maps to 2 genomic regions 
  # simultaneously (this happens, for example, then 2 genes CDS overlap)
  # then both genes will be written in gene_name_var, separated by ;. Therefore
  # it is sensible to grep gene_name in gene_name_var (and other way around) to
  # confirm that genes are the same in gene_name and gene_name_var. However!
  # do not do grep as a word because there are antisense RNA which can be 
  # written as ITFG1-AS1.
  result[, gene_name_toGrep := paste0(gene_name, ';|;', gene_name)]
  result[, gene_name_var_toGrep := paste0(gene_name_var, ';|;', gene_name_var)]
  namesMatchByGrep <- apply(result, 1, 
                            function(x) grepl(x['gene_name_toGrep'], 
                                              x['gene_name_var']) |
                              grepl(x['gene_name_var_toGrep'],
                                    x['gene_name']))
  namesMatched <- (result$gene_name == result$gene_name_var) | namesMatchByGrep
  result[, matched := namesMatched]
  message('[', Sys.time(), '] Found exact match between genomic region gene ',
          'name and gene name assigned to mutation by variant annotation ',
          'tool for ', 
          round(100 * sum(result$matched) / nrow(result), 2), '% (',
                sum(result$matched), ' out of ', nrow(result), ') mutations.')
  
  # Let's check gene names synonyms, maybe genomic regions and variant 
  # annotations are actually meant the same gene, but known under different 
  # names.
  namesNotMatchedSyns <- getGeneSymbolSynonyms(result[matched == F]$gene_name, 
                                               synonymsDT)
  setnames(namesNotMatchedSyns, 'gene_name_orig', 'gene_name')
  # '^', gene_name_syn, '$|'  - is for exact match
  namesNotMatchedSyns[, gene_name_syn := paste0('^', gene_name_syn, '$|', 
                                                gene_name_syn, ';|;', 
                                                gene_name_syn)]
  namesNotMatchedSyns[, gene_name_syn := paste0(gene_name_syn, collapse = '|'),
                      by = gene_name]
  namesNotMatchedSyns <- unique(namesNotMatchedSyns)
  result <- merge(result, namesNotMatchedSyns, by = 'gene_name', 
                       all.x = T)
  result[!is.na(gene_name_syn) & 
           matched == F]$matched <- apply(result[!is.na(gene_name_syn) & 
                                                   matched == F], 1, 
                                          function(x) grepl(x['gene_name_syn'], 
                                                            x['gene_name_var']) |
                                            grepl(x['gene_name_var_toGrep'],
                                                  x['gene_name_syn']))
  message('[', Sys.time(), '] Found exact or synonymous match between ',
          'genomic region gene name and gene name assigned to mutation by ',
          'variant annotation tool for ', 
          round(100 * sum(result$matched) / nrow(result), 2), '% (',
          sum(result$matched), ' out of ', nrow(result), ') ',
          'mutations.')
  
  result <- result[, setdiff(colnames(result),
                             c('gene_name_syn', 'gene_name_toGrep', 
                               'gene_name_var_toGrep')), with = F]
  
  # print table summary of what was not annotated
  noMatchSumTab <- result[matched == F][, length(unique(key)),
                                        by = .(gr_id, var_class)]
  if (nrow(noMatchSumTab) > 0) {
    colnames(noMatchSumTab) <- c('gr_id', 
                                 'genome region type from variant annotator',
                                 'N. mutations')
    noMatchSumTab <- noMatchSumTab[order(-`N. mutations`)]
    noMatchSumTab[, gr_id := factor(gr_id, unique(gr_id))]
    noMatchSumTab <- noMatchSumTab[order(gr_id, -`N. mutations`)]
    message('[', Sys.time(), '] Overview of mutations with unmatched gene ',
            'from genomic region and annotation:')
    message(paste0(capture.output(knitr::kable(noMatchSumTab,
                                               format = "markdown")),
                   collapse = '\n'))
    message('[', Sys.time(), '] Absence of match between gene names of ',
            'genomic regions and genome region type from variant annotator ',
            'is OK for noncoding regions. It can be explained by custom ',
            'definitions of genomic regions submitted by user, overlap of ', 
            'genomic regions which belong to different genes (in such case ',
            'gene name will have format Gene1__and__Gene2) or by ',
            'preference of variant annotator of one region type over the ',
            'other. However!\n DO CHECK MISMATCH FOR CODING VARIANTS BECAUSE ',
            'MAJORITY OF CODING DRIVER CALLERS DO TAKE INFORMATION ABOUT ',
            'VARIANT TYPE, i.e. missence INTO ACCOUNT.')
    percCDSunmatched <- round(100 * nrow(result[matched == F &
                                                  gr_id == 'CDS']) / 
                                nrow(result), 2)
    message('[', Sys.time(), '] Did not find a match between genomic region ',
            'and genome region type from variant annotator for ', 
            nrow(result[matched == F & gr_id == 'CDS']), '(', percCDSunmatched,
            '%) coding mutations. Will re-annotate them with ',
            'VariantAnnotation package')
    if (percCDSunmatched > 1.00) {
      warning('[', Sys.time(), '] Please check your genome definition and 
            annotation files!')
    }
  }
  
  setnames(result, 'matched', 'matchedGeneName')
  result
}

#' checkCodVarClassBetweenVarAnnoAndGR
#' @description Checks that mutations assigned to CDS have accepted variant 
#'              class, i.e. Frame_Shift_Del and not Splice_Site. Also, if 
#'              gene names differ, will indicate that variant needs to be 
#'              re-annotated.
#' @author Maria Litovchenko
#' @param mutsToGRmap data table with each row corresponding to one mutation. 
#'                    Columns gr_id, var_class and matchedGeneName are 
#'                    required
#' @param cdsAcceptedClass variant annotation which are expected to be found in 
#'                      coding regions of the genome
#' @return mutsToGRmap with added column matchedVarClassCod
checkCodVarClassBetweenVarAnnoAndGR <- function(mutsToGRmap,
                                                cdsAcceptedClass = c("Frame_Shift_Del", 
                                                                     "Frame_Shift_Ins",
                                                                     "In_Frame_Del", 
                                                                     "In_Frame_Ins", 
                                                                     "Missense_Mutation", 
                                                                     "Nonsense_Mutation", 
                                                                     "Silent", 
                                                                     "Translation_Start_Site", 
                                                                     "Nonstop_Mutation", 
                                                                     "De_novo_Start_InFrame", 
                                                                     "De_novo_Start_OutOfFrame",
                                                                     'Unknown')) {
  result <- copy(mutsToGRmap)
  
  result[, matchedVarClassCod := T]
  # Variants in CDS, but variant class is not something expected from coding 
  # variant
  nbefore <- nrow(result)
  result[gr_id == 'CDS' & 
           !var_class %in% cdsAcceptedClass]$matchedVarClassCod <- F
  result[gr_id == 'CDS' & matchedGeneName == F]$matchedVarClassCod <- F
  nafter <- nrow(result[matchedVarClassCod == T])
  if (nafter != nbefore) {
    message('[', Sys.time(), '] Found ', nbefore - nafter, '(',
            round(100 * (nbefore - nafter)/nbefore , 2), '% from total) ',
            'coding mutations with the variant class unexpected for coding ',
            'variant. They may arise due to liftover and/or by preference ',
            'of variant annotator of one region type over the other. They ',
            'will be re-annotated.')
  }
  result
}

#' checkNCVarClassBetweenVarAnnoAndGR
#' @description Checks that mutations assigned to noncoding regions have 
#'              accepted variant class, i.e. Splice_Site and not 
#'              Missense_Mutation.
#' @author Maria Litovchenko
#' @param mutsToGRmap data table with each row corresponding to one mutation. 
#'                    Columns gr_id, var_class and matchedGeneName are 
#'                    required
#' @param ncAcceptedClass variant annotation which are expected to be found in 
#'                        noncoding regions of the genome
#' @return mutsToGRmap with added column matchedVarClassNC                    
#' @note In contrast to CDS, if gene names differ, will not report it.
checkNCVarClassBetweenVarAnnoAndGR <- function(mutsToGRmap,
                                               ncAcceptedClass = c("3'UTR",
                                                                   "3'Flank", 
                                                                   "5'UTR", 
                                                                   "5'Flank",
                                                                   "IGR",
                                                                   "Intron",
                                                                   "RNA", 
                                                                   "Targeted_Region", 
                                                                   "Splice_Site",
                                                                   'Unknown')) {
  result <- copy(mutsToGRmap)
  
  result[, matchedVarClassNC := T]
  nbefore <- nrow(result)
  result[gr_id != 'CDS' & 
           !var_class %in% ncAcceptedClass]$matchedVarClassNC <- F
  nafter <- nrow(result[matchedVarClassNC == T])
  if (nafter != nbefore) {
    message('[', Sys.time(), '] Found ', nbefore - nafter, '(',
            round(100 * (nbefore - nafter)/nbefore , 2), '% from total) ',
            'non-coding mutations with the variant class unexpected for ',
            'non-coding variant. Most often it happens to indels, as they ',
            'may be quite long and therefore span several neighboring ',
            'genomic entities, i.e. splice sites and CDS. SNPs assigned to ',
            'non-coding regions which are annotated as coding are most ',
            'likely the result of lost transcrips during liftover or not ',
            'excluding CDS from the set of noncoding regions. Variant class ',
            'of non-coding mutations with the variant class unexpected ',
            'for non-coding variant will be assigned the one of genomic ',
            'region.')
    message('[', Sys.time(), '] Overview of mutations distribution by ',
            'strural variant type:')
    msgTab <- result[matchedVarClassNC == F][,.N, by = .(var_type, gr_id)]
    message(paste0(capture.output(knitr::kable(msgTab, format = "markdown")),
                   collapse = '\n'))
  }
  
  result
}

# FUNCTIONS: coding variants re-annotation ------------------------------------
#' parseKeys
#' @description Parses vector of variant keys into data table with columns chr,
#'              start, end, REF, ALT
#' @author Maria Litovchenko
#' @param keysVector string vector of keys
#' @return data table with columns chr, start, end, REF, ALT
parseKeys <- function(keysVector) {
  result <- lapply(keysVector, function(x) unlist(strsplit(x, ':')))
  result <- do.call(rbind, result)
  colnames(result) <- c('chr', 'start', 'REF', 'ALT')
  result <- as.data.table(result)
  result[, start := as.numeric(start)]
  result[, end := start + nchar(REF)]
  result[, key := keysVector]
  setcolorder(result, c('key', 'chr', 'start', 'end', 'REF', 'ALT'))
  result
}

#' convertIndelsToVarAnnoFormat
#' @description Converts data table containing information about indels from 
#'              maf format to format used in VariantAnnotation package
#' @author Maria Litovchenko
#' @param indelsDT data table with columns chr, start, end, REF, ALT
#' @param bsGenome reference genome, instance of BSgenome
#' @return data table with modified start, end, REF and ALT column
convertIndelsToVarAnnoFormat <- function(indelsDT, bsGenome) {
  result <- copy(indelsDT)
  
  # move the variant 1 base pair left as we will use that base pair as new ALT
  # allele in case of 
  result$start <- result$start - 1
  result$end <- result$start + nchar(result$REF)
  # get REF from reference genome
  seqnames(bsGenome) <- gsub('chr', '', seqnames(bsGenome))
  if (grepl('^chr', result$chr[1])) {
    seqnames(bsGenome) <- paste0('chr', seqnames(bsGenome))
  }
  result[, REFupd := as.character(getSeq(bsGenome, result$chr, 
                                         start = result$start,
                                         end = result$start))]
  result[, REF := apply(result, 1, function(x) paste0(x['REFupd'], x['REF']))]
  result[, ALT := apply(result, 1, function(x) paste0(x['REFupd'], x['ALT']))]
  result[, ALT := gsub('[-]', '', ALT)]
  result[, REF := gsub('[-]', '', REF)]
  result[, REFupd := NULL]
  
  result
}

#' annotateCodVarsWithVariantAnnotation
#' @description Annotates variant with use of VariantAnnotation package.
#' @author Maria Litovchenko
#' @param varsDT data table with information about mutations. Columns chr, 
#'        start, key, REF and ALT are needed.
#' @param specieCode string denoting specie. I.e. Hsapiens for human and 
#'        Mmusculus for mouse
#' @param targetGenomeVersion string, target genome version code, i.e. hg19
#' @return data table with columns key, var_class_reanno, gene_name_reanno. 
#'         key = mutation ID (key), var_class_reanno - new annotation provided
#'         by VariantAnnotation package. gene_name_reanno - gene name.
annotateCodVarsWithVariantAnnotation <- function(varsDT, 
                                                 specieCode = 'Hsapiens', 
                                                 targetGenomeVersion = 'hg19'){
  # load reference genome and annotation from R packages
  bsGenomeLib <- paste0('BSgenome.', specieCode, '.UCSC.', targetGenomeVersion)
  txdbLib <- paste0('TxDb.', specieCode, '.UCSC.', targetGenomeVersion, 
                    '.knownGene')
  orgLib <- paste0('org.', substr(specieCode, start = 1, stop = 2), '.eg.db')
  message('[', Sys.time(), '] Loading: ', bsGenomeLib, ', ', orgLib, ' and ', 
          txdbLib)
  suppressPackageStartupMessages(library(bsGenomeLib, character.only = T))
  suppressPackageStartupMessages(library(txdbLib, character.only = T))
  suppressPackageStartupMessages(library(orgLib, character.only = T))
  txdb <- eval(parse(text = paste0(txdbLib, '::', txdbLib)))  
  bsGenome <- eval(parse(text = paste0(bsGenomeLib, '::', bsGenomeLib)))
  orgDb <- eval(parse(text = paste0(orgLib, '::org.', 
                                    substr(specieCode, start = 1, stop = 2),
                                    '.egSYMBOL')))
  message('[', Sys.time(), '] Finished loading: ', bsGenomeLib, ', ', orgLib,
          ' and ', txdbLib)
  
  varsToAnno <- unique(varsDT[, c('key', 'chr', 'start', 'REF', 'ALT'), 
                              with = F])
  varsToAnno[, end := start]
  # take care of indels - they need to be re-formatted
  indelsIdx <- which(!(toupper(varsToAnno$REF) %in% c('A', 'T', 'G', 'C') &
                         toupper(varsToAnno$ALT) %in% c('A', 'T', 'G', 'C')))
  if (length(indelsIdx) != 0) {
    indelsDT <- varsToAnno[indelsIdx]
    indelsDT <- convertIndelsToVarAnnoFormat(indelsDT, bsGenome)
    varsToAnno <- rbind(varsToAnno[setdiff(1:nrow(varsToAnno), indelsIdx)],
                        indelsDT)
  }
  varsToAnno <- makeGRangesFromDataFrame(varsToAnno, keep.extra.columns = T)
  varsToAnno$ALT <- DNAStringSet(varsToAnno$ALT)
  seqlevelsStyle(txdb) <- seqlevelsStyle(varsToAnno)
  seqlevelsStyle(bsGenome) <- seqlevelsStyle(varsToAnno)
  # perform annotation. lapply makes it slower, but otherwise predictCoding may
  # not annotate certain variants, but it does annotate them then the variants
  # are submitted individually. Therefore, we will first annotate without 
  # lapply what VariantAnnotation can, and then we will top up with lapply
  message('[', Sys.time(), '] Started re-annotation of coding variants')
  result <- predictCoding(varsToAnno, txdb, seqSource = bsGenome, 
                          varAllele = varsToAnno$ALT)
  notAnnotatedIdx <- which(!varsToAnno$key %in% result$key)
  if (sum(notAnnotatedIdx) > 0) {
    result_topup <- GRanges()
    for (notAnnotatedIdx_i in notAnnotatedIdx) {
      message('[', Sys.time(), '] ... ', 
              round(100 * which(notAnnotatedIdx == notAnnotatedIdx_i) /
                      length(notAnnotatedIdx), 2), '%')
      # tryCatch here because if a variant overlaps significantly noncoding 
      # regions, i.e. intron, predictCoding will give an erro.
      tryCatch({toAdd <- predictCoding(varsToAnno[notAnnotatedIdx_i], txdb,
                                       seqSource = bsGenome, 
                                       varAllele = varsToAnno[notAnnotatedIdx_i]$ALT)},
               error = function(e) {toAdd <<- GRanges()})
      result_topup <- c(result_topup, toAdd)
      toAdd <- NULL
    } 
    result <- c(result, result_topup)
  }
  message('[', Sys.time(), '] Finished re-annotation of coding variants')
  
  # sometimes GENEID is set to NA for one line of annotation for a variant 
  # despite being not NA in the other lines. Let's remove it.
  result <- result[!is.na(result$GENEID)]
  
  # deal with not annotated keys
  result <- as.data.table(result)
  if (nrow(result) == 0) {
    keysNotAnnotated <- varsToAnno$key
    result <- data.table(keyStr = character(), GENEID = character(), 
                         CONSEQUENCE = character())
    setnames(result, 'keyStr', 'key')
  } else {
    result <- result[, c('key', 'GENEID', 'CONSEQUENCE')]
    keysNotAnnotated <- setdiff(varsToAnno$key, result$key)
  }
  if (length(keysNotAnnotated) > 0) {
    message('[', Sys.time(), '] Could not re-annotated ', 
            length(keysNotAnnotated), ' (', 
            round(100 * length(keysNotAnnotated) / length(varsToAnno), 2),
            '%) variants with VariantAnnotation package. They will be set to ',
            'Unknown.')
    keysNotAnnotated <- data.table('keyI' = keysNotAnnotated, 
                                   GENEID = '123456789', 
                                   CONSEQUENCE = 'Unknown')
    setnames(keysNotAnnotated, 'keyI', 'key')
    result <- rbind(result, keysNotAnnotated)
  } 
  # convert gene id to gene name
  result[, GENEID := as.character(GENEID)]
  result[, gene_name_reanno := unlist(mget(GENEID, orgDb, 
                                           ifnotfound = NA))[GENEID]]
  result <- result[,.(key, gene_name_reanno, CONSEQUENCE)]
  setnames(result, 'CONSEQUENCE', 'var_class_reanno')
  
  # because gene can have multiple transcripts, we will count annotation by
  # gene and select the most frequent one
  result <- result[,.N, by = .(key, gene_name_reanno, var_class_reanno)]
  result <- result[order(-N)]
  result <- result[,.SD[1], by = .(key, gene_name_reanno, var_class_reanno)]
  result[, N := NULL]
  result
}

#' convertVarAnnoCodesToUserAnnoCodes
#' @description Converts protein coding annotations produced by 
#'              VariantAnnotation to the ones which were produced by software 
#'              used by a user
#' @author Maria Litovchenko
#' @param varsDTreanno data table with columns var_type and var_class_reanno
#' @param codesDT data table with columns: variantAnnotation_anno, var_type, 
#'                var_class which denotes conversion from variantAnnotation,
#'                terms of protein coding variant annotations (i.e. synonymous)
#'                to the ones given by tool used by a user to annotate their 
#'                variants (i.e. Silent in case of annovar). var_type should 
#'                be one of SNP, DEL, INS, MNP.
#' @param notAnnotatedVarCode String which should be used if variantAnnotation
#'                            package fails to re-annotate a variant, i.e. 
#'                            Unknown.
#' @return varsDTreanno data table with updated column var_class_reanno
convertVarAnnoCodesToUserAnnoCodes <- function(varsDTreanno, codesDT, 
                                               notAnnotatedVarCode = 'Unknown') {
  codesDTcopy <- copy(codesDT)
  setnames(codesDTcopy, c('variantAnnotation_anno', 'var_class'),
           c('var_class_reanno', 'var_class_reanno_upd'))
  
  result <- merge(varsDTreanno, codesDTcopy, all.x = T,
                  by = c('var_class_reanno', 'var_type'))
  result[is.na(var_class_reanno_upd)]$var_class_reanno_upd <- notAnnotatedVarCode
  result[, var_class_reanno := NULL]
  setnames(result, 'var_class_reanno_upd', 'var_class_reanno')
  result
}

#' reannotateCodVars
#' @description Re-annotates coding variants with use of VariantAnnotation 
#'              package and updates var_class field in the submitted table, and
#'              unifies it between VariantAnnotation terms and user submitted
#'              terms, if codesDT is given
#' @author Maria Litovchenko
#' @param varsDT data table with information about mutations. Columns chr, 
#'        start, key, REF, ALT, gene_name, var_class and var_type are needed.
#' @param specieCode string denoting specie. I.e. Hsapiens for human and 
#'        Mmusculus for mouse
#' @param targetGenomeVersion string, target genome version code, i.e. hg19
#' @param codesDT data table with columns: variantAnnotation_anno, var_type, 
#'                var_class which denotes conversion from variantAnnotation,
#'                terms of protein coding variant annotations (i.e. synonymous)
#'                to the ones given by tool used by a user to annotate their 
#'                variants (i.e. Silent in case of annovar). var_type should 
#'                be one of SNP, DEL, INS, MNP.
#' @param notAnnotatedVarCode String which should be used if variantAnnotation
#'                            package fails to re-annotate a variant, i.e. 
#'                            Unknown.
#' @return varsDT data table with updated column var_class. In case 
#'         re-annotation was not successful, notAnnotatedVarCode will be used
#'         to fill in var_class.
reannotateCodVars <- function(varsDT, specieCode = 'Hsapiens', 
                              targetGenomeVersion = 'hg19', 
                              codesDT = NULL, 
                              notAnnotatedVarCode = 'Unknown') {
  message('[', Sys.time(), '] Will re-annotate ', nrow(varsDT), ' coding ',
          'variants')
  result <- annotateCodVarsWithVariantAnnotation(parseKeys(varsDT$key),
                                                 specieCode, 
                                                 targetGenomeVersion)
  result <- merge(result, varsDT, by = 'key', allow.cartesian = T)
  # for variants for which gene_name_reanno and gene_name differ, we will
  # set var_class_reanno to annotation_failed and gene_name_reanno to NA
  geneNotMatchIdx <- which(result$gene_name != result$gene_name_reanno |
                             is.na(result$gene_name_reanno))
  result[geneNotMatchIdx]$gene_name_reanno <- NA
  result[geneNotMatchIdx]$var_class_reanno <- 'annotation_failed'
  result <- unique(result)
  # in case the same mutation was annotated to several genes, take the
  # one which matches gene_name from genomic regions
  result[, N := .N, by = .(key, gene_name)]
  result <- result[(N > 1 & gene_name_reanno == gene_name) | N == 1]
  result[, N := NULL]
  
  if (!is.null(codesDT)) {
    result <- convertVarAnnoCodesToUserAnnoCodes(result, codesDT, 
                                                 notAnnotatedVarCode)
  } else {
    message('[', Sys.time(), '] --varanno_conversion_table is not provided, ',
            'will not be able to unify VariantAnnotation codes with codes ',
            'from MAF. That can cause different names for the same mutation ',
            'type, i.e. synonymous and Silent.')
  }
  result[, var_class := var_class_reanno]
  result <- result[, colnames(varsDT), with = F]
  result
}

# FUNCTIONS: mutational rate --------------------------------------------------
#' calcMutRate
#' @description 
#' @author 
#' @param mutsToGRmap
#' @param allParticipantsID 
#' @param GRs
#' @param average_across_patients
#' @return 
calcMutRate <- function(mutsToGRmap, allParticipantsID, GRs,
                        average_across_patients = T) {
  # calculate length of each genomic region
  grsLen <- GRangesList(split(GRs, GRs$gr_name))
  grsLen <- sum(width(grsLen))
  grsLen <- data.table(gr_name = names(grsLen), gr_len = grsLen)
  
  # in order to take into account regions which are not mutated in certain 
  # participants, create a data table with all possible combinations of 
  # participants IDs and genomic regions IDs
  result <- expand.grid(unique(allParticipantsID), unique(GRs$gr_name))
  result <- as.data.table(result)
  colnames(result) <- c('participant_id', 'gr_name')
  
  result <- merge(result, mutsToGRmap[,.N, by = .(participant_id, gr_name)], 
                  by = c('participant_id', 'gr_name'), all.x = T)
  setnames(result, 'N', 'nVars')
  # if a participant didn't have any variants in a genomic bin
  result[is.na(nVars)]$nVars <- 0
  result <- merge(result, grsLen, by = 'gr_name', all.x = T)
  result[, mutRate := nVars/gr_len]
  
  # calculate, to which quantile of mutation rate each region belongs. Do not 
  # take into account regions with 0 mutation rate
  quantBreaks <- result[mutRate != 0]
  quantBreaks <- getUniqQuantileBreaks(quantBreaks$mutRate)
  result[, mutRateQuant := cut(mutRate, 
                               breaks = c(-1, quantBreaks$topBound),
                               labels = quantBreaks$quant)]
  result[, mutRateQuant := as.numeric(as.character(mutRateQuant))]
  result[mutRate == 0]$mutRateQuant <- NA
  
  if (average_across_patients) {
    meanResult <- result[,.(mean(mutRate)), by = gr_name]
    setnames(meanResult, 'V1', 'meanMutRate')
    # calculate, to which quantile of mutation rate each region belongs. Do not 
    # take into account regions with 0 mutation rate
    quantBreaks <- meanResult[meanMutRate != 0]
    quantBreaks <- getUniqQuantileBreaks(quantBreaks$meanMutRate)
    meanResult[, meanMutRateQuant := cut(meanMutRate, 
                                         breaks = c(-1, quantBreaks$topBound),
                                         labels = quantBreaks$quant)]
    meanResult[, meanMutRateQuant := as.character(meanMutRateQuant)]
    meanResult[, meanMutRateQuant := as.numeric(meanMutRateQuant)]
    meanResult[meanMutRate == 0]$meanMutRateQuant <- NA
    result <- merge(result, meanResult, by = 'gr_name')
  }
  
  result
}

#' mutRate_wrapper
#' @description A wrapper function around matchVariantsToRegions and 
#'              calcMutRate functions.
#' @author Maria Litovchenko
#' @param varsDT data table with essential columns chr, start, end. If
#'        checkAnnoMatch or calc_synonymous is requested, var_class is also
#'        needed
#' @param GRs genomic ranges object with gr_name column. Other annotation 
#'            columns may be present too.
#' @param checkAnnoMatch logical, indicates, if check that regions type (i.e. 
#'                      coding or noncoding) corresponds to variant annotation,
#'                      (i.e. if variant is annotated as 3'UTR it should not be
#'                       mapped to CDS) should be performed
#' @param synonymsDT data table with information about gene name synonyms.
#'                   Should have columns idx and gene_name. Synonymous gene
#'                   names will have the same idx. Used to check gene name 
#'                   match between genomic region and mutation. Optional.
#' @param calc_synonymous logical, indicates, if synonymous mutation rate 
#'                        should be computed
#' @param cdsAcceptedClass variant annotation which are expected to be found in 
#'                      coding regions of the genome
#' @param ncAcceptedClass variant annotation which are expected to be found in 
#'                      noncoding regions of the genome
#' @param specieCode string denoting specie. I.e. Hsapiens for human and 
#'        Mmusculus for mouse
#' @param targetGenomeVersion string, target genome version code, i.e. hg19
#' @return named list with values varsToRegsMap - result of 
#'         matchVariantsToRegions function, mutRatesDT and synMutRatesDT -
#'         optional.
mutRate_wrapper <- function(varsDT, GRs, checkAnnoMatch = F, 
                            synonymsDT = NULL, 
                            calc_synonymous = F, cdsAcceptedClass = NULL, 
                            ncAcceptedClass = NULL, specieCode = NULL,
                            targetGenomeVersion = NULL, 
                            codesConversionDT = NULL,
                            annoFailedCode = 'Unknown') {
  if (checkAnnoMatch & (is.null(cdsAcceptedClass) | is.null(ncAcceptedClass))) {
    stop('[', Sys.time(), '] checkAnnoMatch is requested, but ',
         'cdsAcceptedClass or ncAcceptedClass are not submitted.')
  }
  
  message('[', Sys.time(), '] Started mapping variants to genomic regions')
  varsToRegsMap <- matchVariantsToRegions(varsDT = varsDT, GRs = GRs)
  if (checkAnnoMatch) {
    varsToRegsMap <- checkGeneMatchBetweenVarAnnoAndGR(varsToRegsMap,  
                                                       synonymsDT)
    
    # check variant class (i.e. missence, nonsence, frameshift, etc) for 
    # coding mutations
    varsToRegsMap <- checkCodVarClassBetweenVarAnnoAndGR(varsToRegsMap, 
                                                         cdsAcceptedClass)
    varsToRegsMap[, reannotated := F]
    if (nrow(varsToRegsMap[matchedVarClassCod == F]) > 0) {
      varsToReAnnot <- varsToRegsMap[matchedVarClassCod == F]
      varsReAnnoted <- reannotateCodVars(varsToReAnnot, specieCode, 
                                         targetGenomeVersion, 
                                         codesConversionDT, annoFailedCode)
      varsReAnnoted[, reannotated := T]
      varsToRegsMap <- rbind(varsToRegsMap[matchedVarClassCod == T],
                             varsReAnnoted)
    }
    varsToRegsMap[, matchedVarClassCod := NULL]
    
    # check variant class for noncoding mutations
    varsToRegsMap <- checkNCVarClassBetweenVarAnnoAndGR(varsToRegsMap,
                                                        ncAcceptedClass)
    varsToRegsMap[matchedVarClassNC == F]$reannotated <- T
    varsToRegsMap[matchedVarClassNC == F]$var_class <- annoFailedCode
    varsToRegsMap[, matchedVarClassNC := NULL]
  } 
  setcolorder(varsToRegsMap, 
              intersect(c('gr_id', 'gr_name', 'gene_id', 'gene_name', 
                          'participant_id', 'key', 'var_type', 'var_class'), 
                        colnames(varsToRegsMap)))
  message('[', Sys.time(), '] Finished mapping variants to genomic regions')
  
  message('[', Sys.time(), '] Started calculating averaged across ',
          'participants mutation rate')
  mutRatesDT <- calcMutRate(mutsToGRmap = varsToRegsMap, 
                            allParticipantsID = unique(varsDT$participant_id),
                            GRs = GRs, average_across_patients = T)
  message('[', Sys.time(), '] Finished calculating averaged across ',
          'participants mutation rate')
  
  # add info about quantiles of mut.rates to variants to GR map
  varsToRegsMap <- merge(varsToRegsMap, 
                         mutRatesDT[,.(gr_name, participant_id, mutRate,
                                       mutRateQuant, meanMutRate, 
                                       meanMutRateQuant)],
                         by = c('participant_id', 'gr_name'), all.x = T)
  
  if (calc_synonymous) {
    message('[', Sys.time(), '] Started calculating averaged across ',
            'participants synonymous mutation rate')
    synMutRatesDT <- calcMutRate(mutsToGRmap = varsToRegsMap[var_class == 'Silent'], 
                                 allParticipantsID = unique(varsDT$participant_id),
                                 GRs = GRs, average_across_patients = T)
    setnames(synMutRatesDT, 
             c('mutRate', 'mutRateQuant', 'meanMutRate', 'meanMutRateQuant'),
             c('synMutRate', 'synMutRateQuant', 'synMeanMutRate', 
               'synMeanMutRateQuant'))
    message('[', Sys.time(), '] Finished calculating averaged across ',
            'participants synonymous mutation rate')
    varsToRegsMap <- merge(varsToRegsMap, 
                           synMutRatesDT[,.(gr_name, participant_id, 
                                            synMutRate, synMutRateQuant,
                                            synMeanMutRate,
                                            synMeanMutRateQuant)],
                           by = c('participant_id', 'gr_name'), all.x = T)
  }
  result <- list('varsToRegsMap' = varsToRegsMap, 'mutRatesDT' = mutRatesDT)
  if (calc_synonymous) {
    result[['synMutRatesDT']] <- synMutRatesDT
  }
  
  result
}

#' matchLocalMutRateToRegions
#' @description Matches local mutation rate (computed via genomic bins) with
#'              genomic regions of interest
#' @author Maria Litovchenko
#' @param GRs GenomicRanges with genomic regions, columns name, gene_name,
#'            gr_id, gene_id, target_genome_version are required. Can be given
#'            if bed12Path is not available. 
#' @param genomicTilesGR GenomicRanges with genomic coordinates of tiles across
#'        the genome
#' @param genomeTilesMutDT data table with local mutation rates. Columns chr, 
#'        start, end, meanMutRate are required 
#' @return GRs updated with 2 columns meanMutRate_local and  
#'         meanMutRateQuant_local
matchLocalMutRateToRegions <- function(GRs, genomicTilesGR, genomeTilesMutDT) {
  seqlevelsStyle(genomicTilesGR) <- seqlevelsStyle(GRs)
  
  # add mean local mutation rate to genomic regions. If a region overlaps 
  # several tiles from genomicTilesGR, the weighted mean value will be taken
  ovrl <- as.data.table(findOverlaps(GRs, genomicTilesGR))
  colnames(ovrl) <- c('grIdx', 'mutRateIdx')
  # compute % of overlap
  overlappingRegs <- pintersect(GRs[ovrl$grIdx], 
                                genomicTilesGR[ovrl$mutRateIdx])
  percentOvrl <- width(overlappingRegs) / width(GRs[ovrl$grIdx])
  ovrl[, percOvrl := percentOvrl]
  rm(overlappingRegs)
  rm(percentOvrl)
  # add gr_name
  ovrl[, gr_name := GRs[ovrl$grIdx]$gr_name]
  ovrl[, tiles_name := genomicTilesGR[ovrl$mutRateIdx]$gr_name]
  # add local mutation rate & calculate weighted mean
  setkey(genomeTilesMutDT, 'gr_name')
  ovrl[, meanMutRate_local := genomeTilesMutDT[ovrl$tiles_name]$meanMutRate]
  ovrl[, meanMutRate_local := sum(meanMutRate_local * percOvrl), by = grIdx]
  ovrl <- unique(ovrl[,.(gr_name, grIdx, meanMutRate_local)])
  # now, each individual genomic region has mutation rate assigned. But the
  # genomic region we score may consist out of several individual regions, 
  # i.e. CDS out of several exons. We want one estimation of local mutation
  # rate per one genomic region we scanned, so let's take a mean of the
  # individual ones
  ovrl <- ovrl[,.(mean(meanMutRate_local)), by = gr_name]
  setnames(ovrl, 'V1', 'meanMutRate_local')
  
  setkey(ovrl, 'gr_name')
  result <- copy(GRs)
  result$meanMutRate_local <- ovrl[result$gr_name]$meanMutRate_local
  # calculate, to which quantile of meanMutRate each region belongs. Do not 
  # take into account regions with 0 meanMutRate
  quantBreaks <- genomeTilesMutDT[meanMutRate != 0]$meanMutRate
  quantBreaks <- getUniqQuantileBreaks(quantBreaks)
  result$meanMutRateQuant_local <- cut(result$meanMutRate_local, 
                                       breaks = c(-1, quantBreaks$topBound),
                                       labels = quantBreaks$quant)
  result$meanMutRateQuant_local <- as.character(result$meanMutRateQuant_local)
  result$meanMutRateQuant_local <- as.numeric(result$meanMutRateQuant_local)
  result[which(result$meanMutRate_local == 0)]$meanMutRateQuant_local <- NA
  
  result
}

# FUNCTIONS: enrichment of certain structural types ---------------------------
#' checkVarTypeEnrichment
#' @description Checks if a certain type of mutations (i.e. SNP, MNP, INDEL of
#' 1bp length, INDEL of 2bp lenght, etc) is enriched in a gene in comparison to
#' a background via binomial test
#' @author Maria Litovchenko
#' @param varsToGRmap data table with columns key, gr_id, gene_id, gene_name, 
#' participant_id, var_type, mut_len
#' @return data table with columns gr_id, var_cat, gene_id, gene_name, Nvars,
#'         N_gene, bgPerc, binomP, binomPadj
checkVarTypeEnrichment <- function(varsToGRmap) {
  varsDT <- varsToGRmap[,.(key, gr_id, gene_id, gene_name, participant_id, 
                           var_type, mut_len)]
  
  # assign every variant its structural type category
  varsDT[, var_cat := var_type]
  varsDT[var_type %in% c('INS', 'DEL')]$var_cat <- 'INDEL, '
  varsDT[var_cat == 'INDEL, ' & mut_len == 1]$var_cat <- 'INDEL, 1bp'
  varsDT[var_cat == 'INDEL, ' & 
           mut_len > 1 & mut_len < 6]$var_cat <- 'INDEL, 2-5bp'
  varsDT[var_cat == 'INDEL, ' & 
           mut_len > 5 & mut_len < 10]$var_cat <- 'INDEL, 6-9bp'
  varsDT[var_cat == 'INDEL, ' & mut_len > 9]$var_cat <- 'INDEL, 10+bp'
  varsDT[, var_cat := factor(var_cat, 
                             levels = c('SNP', 'MNP', 'INDEL, 1bp', 
                                        'INDEL, 2-5bp', 'INDEL, 6-9bp',
                                        'INDEL, 10+bp'))]
  
  # 1. get background distribution of indels and snps across all genes
  bgDist <- varsDT[,.N, by = .(gr_id, var_cat)] 
  setnames(bgDist, 'N', 'Nvars')
  bgDist <- merge(bgDist, varsDT[,.(length(unique(gene_id))), by = gr_id],
                  by = 'gr_id', all = T)
  setnames(bgDist, 'V1', 'Ngenes')
  bgDist[, bgPerc := Nvars / sum(Nvars), by = .(gr_id)]
  
  # 2. compute enrichment of a variant category per gene in comparison to 
  # background
  geneDist <- varsDT[,.N, by = .(gr_id, gene_id, gene_name, var_cat)]
  setnames(geneDist, 'N', 'Nvars')
  geneDist[, N_gene := sum(Nvars), by = .(gr_id, gene_id, gene_name)]
  geneDist <- merge(geneDist, bgDist[,.(gr_id, var_cat, bgPerc)],
                    by = c('gr_id', 'var_cat'), all = T)
  binomP <- apply(geneDist[,.(Nvars, N_gene, bgPerc)], 1, 
                  function(x) binom.test(x[1], x[2], x[3], 
                                         alternative = 'greater')$p.value)
  geneDist[, binomP := binomP]
  geneDist[, binomPadj := p.adjust(binomP, method = 'BH'), by = gr_id]
  geneDist
}

# Parse input arguments -------------------------------------------------------
parser <- ArgumentParser(prog = '4_calculate_mutation_rates.R')

variantsHelp <- paste('Path to maf file containing genomic variants',
                      '(mutations)')
parser$add_argument("-v", "--variants", required = T, type = 'character',
                    help = variantsHelp)

grHelp <- paste('Part to file containing genomic ranges of region of ',
                'interest, bed12 format')
parser$add_argument("-g", "--genomic_regions", required = T, 
                    type = 'character', help = grHelp)

synTabHelp <- paste('Path to table containing information about gene names ',
                    '(symbols) and their synonyms. Required columns: idx and ',
                    'gene_name.')
parser$add_argument("-m", "--gene_name_synonyms", required = F, default = NULL,  
                    type = 'character', help = synTabHelp)

binHelp <- paste('Length of bin (in bp) in which local mutation rates should',
                 'be computed.')
parser$add_argument("-b", "--bin_len", required = F, default = NULL,
                    type = 'integer', help = binHelp)

subtypeHelp <- paste('Tumor subtype for which mutation rates are computed.',
                     'If given, a column tumor_subtype will be added to ',
                     'all the output.')
parser$add_argument("-t", "--tumor_subtype", required = T, default = NULL,
                    type = 'character', help = subtypeHelp)

synHelp <- paste('Tumor subtype for which mutation rates are computed. If ',
                 'given, a column tumor_subtype will be added to output.')
parser$add_argument("-s", "--calc_synonymous", required = F, default = T,
                    type = 'logical', help = synHelp)

targetGenomeChrHelp <- paste('Path to the tab-separated file containing ', 
                             'chromosomal lengths of the target genome. ',
                             'Must have 2 columns: chr and length. No header.')
parser$add_argument("-l", "--target_genome_chr_len", required = T, 
                    default = '', type = 'character', 
                    help = targetGenomeChrHelp)

cdsClassHelp <- paste('Variant_Classification-s from MAF format which are ',
                      'acceptable as annotation of coding variants. ',
                      'Suggested values: Frame_Shift_Del, Frame_Shift_Ins',
                      "In_Frame_Del, In_Frame_Ins, Missense_Mutation, ",
                      "Nonsense_Mutation, Silent, Translation_Start_Site, ",
                      "Nonstop_Mutation, De_novo_Start_InFrame,",
                      "De_novo_Start_OutOfFrame, Unknown. It is higly ",
                      'reccomended to include Unknown to handle MNPs')
parser$add_argument("-c", "--cdsAcceptedClass", required = T, nargs = '+',
                    type = 'character', help = cdsClassHelp)

ncClassHelp <- paste('Variant_Classification-s from MAF format which are ',
                     'acceptable as annotation of noncoding variants. ',
                     "Suggested values: 3primeUTR, 3primeFlank, 5primeUTR, ",
                     "5primeFlank, IGR, Intron, RNA, Targeted_Region, ", 
                     "Splice_Site, Unknown")
parser$add_argument("-n", "--ncAcceptedClass", required = T, nargs = '+',
                    type = 'character', help = ncClassHelp)

targetVersionHelp <- paste0('Target genome version')
parser$add_argument("-gv", "--target_genome_version", required = T, nargs = 1, 
                    type = 'character', help = targetVersionHelp)

conversionTabHelp <- paste0('Path to table with columns: ',
                            'variantAnnotation_anno, var_type, var_class ',
                            'which denotes conversion from variantAnnotation ',
                            'terms of protein coding variant annotations (i.',
                            'e. synonymous) to the ones given by tool used ',
                            'by a user to annotate their variants (i.e. ',
                            'Silent in case of annovar). var_type should be ',
                            'one of SNP, DEL, INS, MNP.')
parser$add_argument("-a", "--varanno_conversion_table", required = F,
                    nargs = 1, type = 'character', help = conversionTabHelp)

notAnnotatedVarCodeHelp <- paste0('String which should be used if ',
                                  'variantAnnotation package fails to ',
                                  're-annotate a variant, i.e. Unknown.')
parser$add_argument("-f", "--annotation_failed_code", required = T, nargs = 1,
                    type = 'character', help = notAnnotatedVarCodeHelp)

outputHelp <- paste("Prefix to use in the output files.",
                    "All files will start with this.")
parser$add_argument("-p", "--output", required = T, type = 'character',
                    help = outputHelp)
args <- parser$parse_args()

# Test arguments --------------------------------------------------------------
args <- list('variants' = '../work/46/4f35fcb3f86675007be24a347c6f5b/inputs/inputMutations-LUSC-hg19.maf',
             'genomic_regions' = '../work/8c/3bb3b24a9af0a805dc0cc73ca1a6de/inputs/inputGR-LUSC-hg19.bed',
             'tumor_subtype' = 'LUSC', 'target_genome_version' = 'hg19',
             'target_genome_chr_len' = '../data/assets/reference_genome/Homo_sapiens_assembly19.chrom.sizes',
             'gene_name_synonyms' = '../data/assets/gene_names_synonyms/hgnc_complete_set_2022-07-01_proc.csv',
             'varanno_conversion_table' = '../data/assets/variantAnnotation_to_annovar_conversion.txt',
             'bin_len' = 5*(10^4), 'calc_synonymous' = T,
             'cdsAcceptedClass' = list("Frame_Shift_Del", "Frame_Shift_Ins",
                                       "In_Frame_Del", "In_Frame_Ins", 
                                       "Missense_Mutation", 
                                       "Nonsense_Mutation", "Silent", 
                                       "Translation_Start_Site", 
                                       "Nonstop_Mutation", 
                                       "De_novo_Start_InFrame", 
                                       "De_novo_Start_OutOfFrame", 
                                       "Unknown"),
             'ncAcceptedClass' = list("3'UTR", "3'Flank", "5'UTR", 
                                      "5'Flank", "IGR", "Intron", "RNA", 
                                      "Targeted_Region", "Splice_Site", 
                                      'Unknown'), 
             'annotation_failed_code' = 'Unknown',
             'output' = 'mutRate-MET_PANCAN-hg19-')
 
args$ncAcceptedClass <- gsub('prime', "'", args$ncAcceptedClass)

check_input_arguments(args, outputType = 'file')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# READ in mutation, genome region and chr lengths files -----------------------
message('[', Sys.time(), '] Started reading input mutation file')
message('[', Sys.time(), '] Started reading ', args$variants)
colsToKeep <- c('Tumor_Sample_Barcode', 'key', 'Chromosome',
                'Start_Position', 'End_Position', 'Gene.refGene', 
                'Variant_Classification', 'Variant_Type', 'mut_len', 
                'struct_type', 'patient_tumor_subtype')
updColNames <- c('participant_id', 'key', 'chr', 'start', 'end', 
                 'gene_name_var', 'var_class', 'var_type', 'mut_len',
                 'struct_type', 'patient_tumor_subtype')
allVars <- fread(args$variants, header = T, stringsAsFactors = F, 
                 select = colsToKeep)
setnames(allVars, colsToKeep, updColNames, skip_absent = T)
message('[', Sys.time(), '] Finished reading ', args$variants)

message('[', Sys.time(), '] Started reading input genomic ranges file')
# GR = regions of interest
GR <- import(args$genomic_regions)
GR <- unlist(blocks(GR))
GR_parsed <- parseBED12regName(unique(names(GR)))
setnames(GR_parsed, 'name', 'gr_name')
setkey(GR_parsed, gr_name)
mcols(GR) <- GR_parsed[names(GR)]
rm(GR_parsed)
message('[', Sys.time(), '] Finished reading input genomic ranges file')
message('[', Sys.time(), '] rtracklayer import function while reading ',
        'bed12 assumes that it is 0-based. Therefore it adds 1 to all ',
        'coordinates. We will correct it.')
start(GR) <- start(GR) - 1
end(GR) <- end(GR) - 1

message('[', Sys.time(), '] Started reading chromosomal length file')
chrLensDT <- fread(args$target_genome_chr_len, header = F, 
                   stringsAsFactors = F)
chrLensVect <- chrLensDT$V2
names(chrLensVect) <- chrLensDT$V1
message('[', Sys.time(), '] Finished reading chromosomal length file')

# READ in table with gene name synonyms & variant codes conversion ------------
symbolSynsDT <- NULL
if (!is.null(args$gene_name_synonyms)) {
  message('[', Sys.time(), '] Started reading gene name synonyms table')
  symbolSynsDT <- fread(args$gene_name_synonyms, header = T, 
                        stringsAsFactors = F, select = c('idx', 'gene_name'))
  message('[', Sys.time(), '] Finished reading gene name synonyms table')
}

codesConvertDT <- NULL
if (!is.null(args$varanno_conversion_table)) {
  message('[', Sys.time(), '] Started reading variant annotation conversion ',
          'table')
  codesConvertDT <- fread(args$varanno_conversion_table, header = T, 
                          stringsAsFactors = F)
  message('[', Sys.time(), '] Finished reading variant annotation conversion ',
          'table')
}

# Calculate number of variants per participant ---------------------------------
nMutsPerPart <- allVars[,.N, by = participant_id][order(N)]
nMutsPerPart[, tumor_subtype := args$tumor_subtype]

# Bin genome on tiles ---------------------------------------------------------
if (!is.null(args$bin_len)) {
  message('[', Sys.time(), '] Started binning genome on bins sized ', 
          args$bin_len, 'bp.')
  # bin genome 
  genomeTilesGR <- tileGenome(chrLensVect, tilewidth = args$bin_len, 
                              cut.last.tile.in.chrom = T)
  mcols(genomeTilesGR)$gr_name <- apply(as.data.table(genomeTilesGR), 1,
                                        function(x) paste(x['seqnames'],
                                                          x['start'], x['end'],
                                                          sep = '--'))
  mcols(genomeTilesGR)$gr_name <- gsub(' ', '', mcols(genomeTilesGR)$gr_name)
  message('[', Sys.time(), '] Finished binning genome on bins sized ', 
          args$bin_len, 'bp.')
}

# Match mutations to genomic regions, calculate mut.rates ---------------------
if (!is.null(args$bin_len)) {# for tiles over genome
  message('[', Sys.time(), '] Working with genomic regions created as bins')
  mutRates_GenomeTiles <- mutRate_wrapper(varsDT = allVars, 
                                          GRs = genomeTilesGR, 
                                          checkAnnoMatch = F, 
                                          calc_synonymous = T)
}
if (!is.null(args$genomic_regions)) {# for custom genomic regions
  message('[', Sys.time(), '] Working with genomic regions from bed12 file')
  mutRates_GR <- mutRate_wrapper(varsDT = allVars, GRs = GR, 
                                 checkAnnoMatch = T, 
                                 synonymsDT = symbolSynsDT,
                                 calc_synonymous = T,
                                 cdsAcceptedClass = args$cdsAcceptedClass, 
                                 ncAcceptedClass = args$ncAcceptedClass, 
                                 specieCode = 'Hsapiens',
                                 targetGenomeVersion = args$target_genome_version,
                                 codesConversionDT = codesConvertDT, 
                                 annoFailedCode = args$annotation_failed_code)
}

# Matched local mut.rates(tiled genome) to custom regions ---------------------
if (!is.null(args$genomic_regions) & !is.null(args$bin_len)) {
  GR <- matchLocalMutRateToRegions(GR, genomeTilesGR, 
                                   unique(mutRates_GenomeTiles$mutRatesDT[,.(gr_name,
                                                                             meanMutRate)]))
  message('[', Sys.time(), '] Because local mutation rates are available, ',
          'will also add them to the variats to genomic regions map in ',
          'columns meanMutRate_local and meanMutRateQuant_local')
  GRdt <- as.data.table(mcols(GR))
  GRdt <- GRdt[,.(gr_name, meanMutRate_local, meanMutRateQuant_local)]
  GRdt <- unique(GRdt)
  mutRates_GR$varsToRegsMap <- merge(mutRates_GR$varsToRegsMap, GRdt,
                                     by = 'gr_name', all.x = T)
}

# Check if there is a significant enrichment of variant types (i.e.INDELs)-----
if (!is.null(args$genomic_regions)) {
  message('[', Sys.time(), '] Started calculations of variant types ',
          'enrichment across genes')
  varCatEnrich <- checkVarTypeEnrichment(mutRates_GR$varsToRegsMap)
  message('[', Sys.time(), '] Finished calculations of variant types ',
          'enrichment across genes')
}

# [OUTPUT] to files -----------------------------------------------------------
write.table(nMutsPerPart, quote = F, sep = '\t', row.names = F, col.names = T,
            paste0(args$output, '-nMutsPerParticipant.csv'), append = F)

if (!is.null(args$bin_len)) {
  write.table(cbind('tumor_subtype' = args$tumor_subtype, 
                    mutRates_GenomeTiles$mutRatesDT), 
              sep = '\t', col.names = T, row.names = F, append = F, quote = F,
              file = paste0(args$output, '-meanMutRatePerBin.csv'))
  
  
  write.table(cbind('tumor_subtype' = args$tumor_subtype, 
                    mutRates_GenomeTiles$varsToRegsMap), 
              sep = '\t', col.names = T, row.names = F, append = F, quote = F,
              file = paste0(args$output, '-mutMapToGenomeBins.csv'))
}

if (!is.null(args$genomic_regions)) {
  write.table(cbind('tumor_subtype' = args$tumor_subtype, 
                    mutRates_GR$mutRatesDT), 
              sep = '\t', row.names = F, col.names = T, append = F, quote = F,
              file = paste0(args$output, '-meanMutRatePerGR.csv'))
  write.table(cbind('tumor_subtype' = args$tumor_subtype, 
                    mutRates_GR$varsToRegsMap), 
              sep = '\t', row.names = F, col.names = T, append = F, quote = F,
              file = paste0(args$output, '-mutMapToGR.csv'))
  write.table(cbind('tumor_subtype' = args$tumor_subtype, varCatEnrich), 
              sep = '\t', row.names = F, col.names = T, append = F, quote = F,
              file = paste0(args$output, '-varCatEnrich.csv'))
}

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')