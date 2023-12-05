#!/usr/bin/env Rscript
# FILE: combine_p_values_and_annotate_with_metadata.R -------------------------
#
# DESCRIPTION: A script to 
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
# CREATED:  10.12.2020
# REVISION: 16.11.2023

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
suppressWarnings(suppressPackageStartupMessages(library(EmpiricalBrownsMethod)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(poolr)))

# Functions : Reading results of de-novo cancer driver genes calling  ---------
#' read_DIGdriver_results
#' @description Reads raw results files from DIGdriver
#' @author Maria Litovchenko
#' @param filePath path to the file
#' @return data.table with columns 
read_DIGdriver_results <- function(filePath) {
  colsToRead <- c('ELT', 'PVAL_MUT_BURDEN', 'PVAL_SNV_BURDEN',
                  'PVAL_INDEL_BURDEN')
  result <- fread(filePath, header = T, stringsAsFactors = F, 
                  select = colsToRead)
  if (ncol(result) < 2) {
    stop('[', Sys.time(), '] DigDriver output lacks all 3 columns: ',
         'PVAL_SNV_BURDEN, PVAL_INDEL_BURDEN, PVAL_MUT_BURDEN') 
  }
  result <- result[, intersect(colnames(result), colsToRead)[1:2], with = F]
  setnames(result, colnames(result), c('binID', 'raw_p'), skip_absent = T)
  result <- cbind(result, parseBED12regName(result$binID))
  result <- result[, !colnames(result) %in% c('binID', 'name', 
                                              'target_genome_version'),
                   with = F]
  setnames(result, 'gr_id', 'grID')
  
  result
}

#' read_dNdScv_results
#' @description Reads raw results files from dNdScv
#' @author Maria Litovchenko
#' @param filePath path to the file
#' @return data.table with columns gene_name and raw_p.
read_dNdScv_results <- function(filePath) {
  colsToRead <- c('gene_name', 'pallsubs_cv', 'pglobal_cv')
  result <- fread(filePath, header = T, stringsAsFactors = F,
                  select = colsToRead)
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
  setnames(result, 'pglobal_cv', 'raw_p')

  result
}

#' read_MutPanning_results
#' @description Reads raw results files from MutPanning
#' @author Maria Litovchenko
#' @param filePath path to the file
#' @return data.table with columns gene_name and raw_p.
read_MutPanning_results <- function(filePath) {
  result <- fread(filePath, header = T, stringsAsFactors = F,
                  select = c('Name', 'Significance'))
  setnames(result, c('Name', 'Significance'), c('gene_name', 'raw_p'),
           skip_absent = T)
  result
}

#' read_NBR_results
#' @description Reads raw results files from NBR
#' @author Maria Litovchenko
#' @param filePath path to the file
#' @return data.table with columns gene_id, raw_p.
read_NBR_results <- function(filePath) {
  colsToRead <- c('region', 'pval_subs_CV', 'pval_indels_CV', 'pval_both_CV')
  result <- fread(filePath, header = T, stringsAsFactors = F, 
                  select = colsToRead)
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
  result
}

#' read_OncodriveFML_results
#' @description Reads raw results files from OncodriveFML
#' @author Maria Litovchenko
#' @param filePath path to the file
#' @return data.table with columns gene_id, gene_name, raw_p.
read_OncodriveFML_results <- function(filePath) {
  result <- fread(filePath, header = T, stringsAsFactors = F,
                  select = c('GENE_ID', 'SYMBOL', 'P_VALUE'))
  setnames(result, c('GENE_ID', 'SYMBOL', 'P_VALUE'), 
           c('gene_id', 'gene_name', 'raw_p'), skip_absent = T)
  result
}

#' readSoftwareResults
#' @description Reads in result files from dndscv, driverpower, mutpanning, 
#' nbr, oncodrivefml into uniformal data table.
#' @author Maria Litovchenko
#' @param softName name of the software filePath corresponds to. One of dndscv,
#' driverpower, mutpanning, nbr, oncodrivefml
#' @param acceptedSoftware vector of accepted software
#' @param filePath path to result file
#' @param infoStr named vector with at least items cancer_subtype, gr_id
#' @return data table with columns software, tumor_subtype, gr_id
readSoftwareResults <- function(softName, acceptedSoftware, filePath, infoStr){
  if (!file.exists(filePath)) {
    message('[', Sys.time(), '] File not found: ', filePath, '. Return empty ',
            'data table.')
    return(data.table())
  }
  if (!softName %in% acceptedSoftware) {
    stop('[', Sys.time(), '] readSoftwareResults: softName should be one of ',
         paste0(acceptedSoftware, collapse = ','), '. Current value: ',
         softName)
  }
  
  if (!all(c('cancer_subtype', 'gr_id') %in% names(infoStr))) {
    stop('[', Sys.time(), '] readSoftwareResults: infoStr vector should have ',
         'items under names cancer_subtype, gr_id')
  }
  infoStr <- infoStr[c('cancer_subtype', 'gr_id')]
  names(infoStr)[1] <- 'tumor_subtype'
  infoStr <- as.data.table(t(infoStr))
  
  message('[', Sys.time(), '] Reading: ', filePath)
  result <- switch(softName,
                   "digdriver"    = read_DIGdriver_results(filePath), 
                   "dndscv"       = read_dNdScv_results(filePath),
                   "mutpanning"   = read_MutPanning_results(filePath),
                   "nbr"          = read_NBR_results(filePath),
                   "oncodrivefml" = read_OncodriveFML_results(filePath))
  result <- cbind(result, infoStr)
  
  if (softName == 'digdriver') {
    result <- result[grID == gr_id]
    result[, grID := NULL]
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
  
  if (nrow(unique(result)) != nrow(result)) {
    duplicatedDT <- unique(result[duplicated(result)])
    message('[', Sys.time(), '] Found duplicated entries for in ', filePath, 
            ' for ', length(duplicatedDT$gene_name), ' genes: ',
            paste0(duplicatedDT$gene_name, collapse = ', '), 
            '. Will removed them.')
    result <- result[!duplicated(result)]
  }
  
  result
}

# Functions : gene ID and gene name unification across software ---------------
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

# Functions : Brown method of combining p-values ------------------------------
#' transformData
#' Derived from EmpiricalBrownsMethod package, modified to handle NA values
#' @description Transforms (aka normalizes) data_vector with mean and standard
#' deviation. 
#' @param data_vector numerical vector
#' @return numerical vector
transformData <- function(data_vector) {
  dataVectorNoNA <- data_vector[!is.na(data_vector)]
  
  dvm <- mean(dataVectorNoNA)
  dvsd <- var(dataVectorNoNA) 
  dvsd <- dvsd * (length(dataVectorNoNA) - 1) / length(dataVectorNoNA)
  dvsd <- sqrt(dvsd)
  s <- (dataVectorNoNA - dvm) / dvsd
  distr = ecdf(dataVectorNoNA)
  sapply(data_vector, function(a) -2*log(distr(a)))
}

#' combinePValues
#' Derived from EmpiricalBrownsMethod package, modified to handle NA
#' @param covar_matrix A m x m numpy array of covariances between transformed
#'        data vectors 
#' @param p_values vector of m p-values to combine.
#' @param extra_info == True: also returns the p-value from Fisher's method,
#'        the scale factor c, and the new degrees of freedom from Brown's 
#'        Method
#' @return a combined P-value.
combinePValues <- function(covar_matrix, p_values, extra_info = FALSE){
  N = ncol(covar_matrix) # number of samples
  df_fisher = 2.0*N
  Expected  = 2.0*N
  cov_sum <- (2*sum(covar_matrix[lower.tri(covar_matrix, diag=FALSE)]))
  Var = 4.0*N+cov_sum
  c = Var/(2.0*Expected)
  df_brown = (2.0*Expected^2)/Var
  if (df_brown > df_fisher) {
    df_brown = df_fisher
    c = 1.0
  }
  x = 2.0*sum( -log(p_values[!is.na(p_values)]) )
  
  p_brown = pchisq(df=df_brown, q=x/c, lower.tail=FALSE)
  p_fisher = pchisq(df=df_fisher, q=x, lower.tail=FALSE)
  
  if (extra_info) {
    return(list(P_test=p_brown, P_Fisher=p_fisher, Scale_Factor_C=c, DF=df_brown))
  }
  else {
    return(p_brown)
  }
}

#' empiricalBrownsMethod_V
#' Vectorised version of empiricalBrownsMethod from EmpiricalBrownsMethod 
#' package which is derived from EmpiricalBrownsMethod package.
#' @param covar_matrix covariance matrix of p-values. Compute as following:
#'  cov(apply(p_val_matr, 1, transformData)), where p_val_matr is a value
#'  matrix. Genes in columns, methods in rows.
#' @param p_values_list a vector of p values from which brown should be 
#'        computed
#' @param extra_info whatever or not fisher, c and such should be returned
#' @return list
empiricalBrownsMethod_V <- function(covar_matrix, p_values_list, 
                                    extra_info = F) {
  result <- lapply(p_values_list, 
                   function(x) combinePValues(covar_matrix, x, extra_info))
  result <- lapply(result, unlist)
  result <- do.call(rbind, result)
  return(result)
}

# Functions : Combining p-values ----------------------------------------------
#' meltToWide
#' @description Melts result data table into wide format where each software
#' has its own column with raw p values.
#' @author Maria Litovchenko
#' @param resDT data table, columns software and raw_p are essential.
#' @return data table with all the column in resDT which are not software or 
#' raw_p + one column per software with raw p-vlaues.
meltToWide <- function(resDT) {
  frml <- paste(paste(setdiff(colnames(resDT), c('software', 'raw_p')),
                      collapse = '+'), "~ software")
  result <- dcast(resDT, as.formula(frml), value.var = 'raw_p', )
  result <- as.data.table(result)
  setnames(result, sort(unique(resDT$software)),
           paste0(sort(unique(resDT$software)), '.raw_p'))
  result
}

#' harmonic_mean
#' @description Computes harmonic mean
#' @author Maria Litovchenko 
#' @param x numerical vector
#' @return numeric value
harmonic_mean <- function(x) {
  x <- x[!is.na(x)]
  result <- length(x) / sum(1/x)
  result
}

#' createUID
#' @description Creates gene unique identifier based on data table with 
#'              tumor_subtype, gr_id, gene_name and gene_id as columns.
#' @author Maria Litovchenko
#' @param geneInfoDT data table with columns tumor_subtype, gr_id, gene_name, 
#'        gene_id
#' @param separator string to use as a divider for future unique IDs
#' @return vector of unique IDs made by pasting tumor_subtype, gr_id, 
#'         gene_name and gene_id together.
createUID <- function(geneInfoDT, separator = '--') {
  result <- apply(geneInfoDT[,.(tumor_subtype, gr_id, gene_name, 
                                gene_id)], 1, 
                  function(x) gsub(' ', '', paste(x, collapse = separator)))
  result
}

#' parseUID
#' @description Parses unique gene identifier to data table
#' @author Maria Litovchenko
#' @param UIDs vector/list of unique gene identifiers
#' @param separator string to use as a divider for future unique IDs
#' @return data table with columns tumor_subtype, gr_id, gene_id, gene_name
parseUID <- function(UIDs, separator = '--') {
  result <- lapply(UIDs, strsplit, split = separator)
  result <- lapply(result, unlist)
  result <- as.data.table(do.call(rbind, result))
  colnames(result) <- c('tumor_subtype', 'gr_id', 'gene_name', 'gene_id')
  result
}

#' combineRawPvalues
#' @description Combines raw p values from various software into one single
#' p-value.
#' @author Maria Litovchenko
#' @param resDTwide data table with raw results of various software runs. 
#' Columns gene_id as well as columns ending with .raw_p holding p values of
#' individual softwares are required. 
#' @return resDTwide with added columns brown.comb_p, fisher.comb_p, 
#' stouffer.comb_p  and harmonic.comb_p
combineRawPvalues <- function(resDTwide) {
  pValMatr <- resDTwide[,(grepl('.raw_p$', colnames(resDTwide))), with = F]
  
  # sometimes there are situations then there is just one non-NA value for 
  # one of the driver calling methods. In such case covariance matrix can't be
  # calculated and that entry should be removed.
  if(sum(complete.cases(pValMatr)) == 1) {
    nNonNa <- apply(pValMatr, 2, function(x) sum(!is.na(x)))
    toRemove <- which(nNonNa == min(nNonNa))
    message('[', Sys.time(), '] Removed ', colnames(pValMatr)[toRemove], ' ',
            'from matrix of raw p values for ', unique(resDTwide$gr_id), 
            ' in ', unique(resDTwide$tumor_subtype), ' because there was ',
            'just one non-NA entry for that software.')
    pValMatr <- pValMatr[, setdiff(1:ncol(pValMatr), toRemove), with = F]
  }
  
  if (any(apply(pValMatr, 1, function(x) any(x[!is.na(x)] == 0)))) {
    stop('[', Sys.time(), '] Found 0 in p-value matrix. 0 can not be handled ',
         'by harmonic mean calculations. Please replace 0 with suitable value')
  }
  
  # brown method function understands only matrices
  pValMatr <- as.matrix(pValMatr)
  rownames(pValMatr) <- createUID(resDTwide)
  pValMatr <- t(pValMatr)
  # Brown method of combining p-values
  brown <- empiricalBrownsMethod_V(cov(apply(pValMatr, 1, transformData),
                                       use = "pairwise.complete.obs"), 
                                   lapply(1:ncol(pValMatr), 
                                          function(y) pValMatr[, y]), 
                                   extra_info = T)
  brown <- as.data.table(brown)[,.(P_test, P_Fisher)]
  setnames(brown, c('P_test', 'P_Fisher'), c('brown.comb_p', 'fisher.comb_p'))
  # return back all the meta info
  brownMetaInfo <- parseUID(colnames(pValMatr))
  brown <- cbind(brown, brownMetaInfo)
  
  # Stouffer method of combining p-values
  stoufferP <- apply(pValMatr, 2, function(y) stouffer(y[!is.na(y)])$p)
  stoufferP <- cbind(parseUID(colnames(pValMatr)),
                     'stouffer.comb_p' = stoufferP)
  
  # Harmonic method of combining p-values
  harmonic <- apply(pValMatr, 2, harmonic_mean)
  harmonic <- cbind(parseUID(colnames(pValMatr)), 'harmonic.comb_p' = harmonic)
  
  result <- merge(resDTwide, brown, by = c('tumor_subtype', 'gr_id',
                                           'gene_name', 'gene_id'), all.x = T)
  result <- merge(result, stoufferP, by = c('tumor_subtype', 'gr_id',
                                            'gene_name', 'gene_id'), all.x = T)
  result <- merge(result, harmonic, by = c('tumor_subtype', 'gr_id',
                                           'gene_name', 'gene_id'), all.x = T)
  result
}

# Functions : Annotations with expression -------------------------------------
#' annotateWithExpression
#' @description Annotates table containing genes with expression status.
#' @author Maria Litovchenko
#' @param dt data table with essential columns gene_id, gene_name and 
#' tumor_subtype. Data should be only for 1 tumor subtype.
#' @param expression_InvPath path to the file with expression inventory. The
#' file must contain columns tumor_subtype, expr_db and min_expr.
#' @param expression_filePath path to processed table with expression data for
#' one expression data base. Columns gene_id, gene_name and a column 
#' containing expression data for the tumor_subtype (from dt) should be
#' present.
#' @return data table dt with added column expr_in_{expr_db} indicating state
#' of expression (T or F).
annotateWithExpression <- function(dt, expression_InvPath, 
                                   expression_filePath) {
  dt_tumorSubtype <- unique(dt$tumor_subtype)
  if (length(dt_tumorSubtype) > 1) {
    stop('[', Sys.time(), '] > 1 unique values of tumor_subtype column in ',
         'table to be annotated.')
  }
  
  inventory_expression <- fread(expression_InvPath, header = T, 
                                stringsAsFactors = F)
  message('[', Sys.time(), '] Started annotation with ', 
          unique(inventory_expression$expr_db))
  if (!dt_tumorSubtype %in% inventory_expression$tumor_subtype) {
    message('[', Sys.time(), '] No expression data is available for ',
            paste0(dt_tumorSubtype, collapse = ', '), ', skipping.')
  }
  min_expr <- inventory_expression[tumor_subtype == dt_tumorSubtype]$min_expr
  min_expr <- as.numeric(min_expr)
  
  message('[', Sys.time(), '] Started reading expression table ',
          expression_filePath)
  exprDT <- fread(expression_filePath, header = T, stringsAsFactors = F,
                  select = c('gene_id', 'gene_name', dt_tumorSubtype))
  message('[', Sys.time(), '] Finished reading expression table ',
          expression_filePath)
  
  result <- mergeBasedOnAtLeastOneKey(dt, exprDT, 
                                      keyCols = c('gene_id', 'gene_name'))
  setnames(result, dt_tumorSubtype, 'expr_in')
  result[, expr_in := ifelse(expr_in > min_expr, T, F)]
  setnames(result, 'expr_in', 
           paste0('expr_in_', unique(inventory_expression$expr_db)))
  result
}
# Functions : Misc ------------------------------------------------------------
#' mergeBasedOnAtLeastOneKey
#' @description Performs merge of two data tables based on at least one 
#'              overlapping key column
#' @param baseDT base data table. All rows and columns on this data table will
#'        be preserved
#' @param toAddDT data table from which columns will be added to baseDT. 
#' @param keyCols names of key columns
#' @return data table baseDT with added columns from toAddDT matched based on
#'         at least one column from keyCols
mergeBasedOnAtLeastOneKey <- function(baseDT, toAddDT, keyCols) {
  result <- merge(baseDT, toAddDT, by = keyCols)
  notKeys <- setdiff(colnames(toAddDT), keyCols)
  for (oneKeyCol in keyCols) {
    notMatched <- merge(baseDT, unique(cbind(result[, keyCols, with = F],
                                             matched = T)), 
                        by = keyCols, all.x = T)
    notMatched <- notMatched[is.na(matched)]
    notMatched[, matched := NULL]
    result <- rbind(result, merge(notMatched, 
                                  toAddDT[, c(oneKeyCol, notKeys), with = F], 
                                  by = oneKeyCol))
  }
  # add lines for which match wasn't found at all
  notMatched <- merge(baseDT, unique(cbind(result[, keyCols, with = F], 
                                           matched = T)), 
                      by = keyCols, all.x = T)
  notMatched <- notMatched[is.na(matched)]
  notMatched[, matched := NULL]
  result <- rbind(result, notMatched, fill = T)
  # return back the original order
  result <- merge(baseDT, result, 
                  by = intersect(colnames(baseDT), colnames(result)))
  result
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'combine_p_values_and_annotate_with_metadata.R')

subtypeHelp <- paste('A cancer subtype to select from patientsInv table. Only',
                     'mutations from patients with that cancer type will be',
                     'selected. In case an analysis of several cancer types',
                     'needed to be performed please run this script ',
                     'separetedly for each cancer type.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    default = NULL, help = subtypeHelp)

parser$add_argument("-g", "--gr_id", required = T, type = 'character',
                    default = NULL, help = 'IDs of a genomic region')

softHelp <- paste('A list of software names run on current cancer subtype',
                  'genomic region pair.')
parser$add_argument("-s", "--software", required = T, type = 'character',
                    nargs = '+', default = NULL, help = softHelp)

resHelp <- paste('A list of file paths to results of software runs on current',
                 'cancer subtype genomic region pair. The order of paths',
                 'should match the order of software names in --software',
                 'argument')
parser$add_argument("-d", "--run_results", required = T, type = 'character',
                    nargs = '+', default = NULL, help = resHelp)

mutRhelp <- paste('A path to file containing mutations rates for the',
                  'selected --cancer_subtype. Usually this file is created',
                  'by 4_calculate_mutation_rates.R script and has a prefix',
                  'meanMutRatePerGR-')
parser$add_argument("-r", "--mut_rate", required = T, type = 'character',
                    default = NULL, help = mutRhelp)

synTabHelp <- paste('Path to table containing information about gene names',
                    '(symbols) and their synonyms. Required columns: idx and',
                    'gene_name.')
parser$add_argument("-m", "--gene_name_synonyms", required = F, default = NULL,
                    type = 'character', help = synTabHelp)

knownHelp <- paste('Path to table containing information about known cancer',
                   'associated genes. Must have columns gene_id and gene_name')
parser$add_argument("-k", '--known_cancer_genes', required = F, default = NULL,
                    type = 'character', help = knownHelp)

olfactHelp <- paste('Path to table containing information about known',
                    'olfactory genes. Must have columns gene_id and gene_name')
parser$add_argument("-of", '--olfactory_genes', required = F, default = NULL,
                    type = 'character', help = olfactHelp)

capHelp <- paste('Cap on raw p-values: minimum p-value values below which',
                 'would be set to that value. I.e. if rawP_cap is 1e-8, all',
                 'p values below 1e-8 (i.e.  1e-16,  1e-15, 0) will be set',
                 'to 1e-8.')
parser$add_argument("-t", '--rawP_cap', required = T, default = NULL,
                    type = 'character', help = capHelp)

expIhelp <- paste('A list of paths to inventories outlining expression data',
                  'to be used in driver annotation. Must have columns',
                  'expr_db, tumor_subtype and min_expr.')
parser$add_argument('-ei', '--expression_inventories', required = F, 
                    nargs = '+', default = NULL, type = 'character',
                    help = expIhelp)

expHelp <- paste('A list of file paths to expression tables. The order of',
                 'paths should match the order of inventories in',
                 '--expression_inventories argument')
parser$add_argument('-e', '--expression', required = F, default = NULL,
                    nargs = '+', type = 'character', help = expHelp)

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output file")

args <- parser$parse_args()
# check_input_arguments_postproc(args, outputType = 'file')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
names(args[['run_results']]) <- unlist(args$software)
args$rawP_cap <- as.numeric(args$rawP_cap)
printArgs(args)

# Test inputs -----------------------------------------------------------------
# args <- list(cancer_subtype = 'LUAD', gr_id = '5primeUTR',
#              software = list('digdriver', 'dndscv', 'mutpanning', 'nbr',
#                              'oncodrivefml'),
#              run_results = list('../TEST/results/digdriver/digdriverResults-LUAD-5primeUTR-hg19.csv', 
#                                 '../TEST/results/dndscv/dndscvResults-LUAD-5primeUTR-hg19.csv', 
#                                 '../TEST/results/mutpanning/mutpanningResults-LUAD-5primeUTR-hg19.csv', 
#                                 '../TEST/results/nbr/nbrResults-LUAD-5primeUTR-hg19.csv', 
#                                 '../TEST/results/oncodrivefml/oncodrivefmlResults-LUAD-5primeUTR-hg19.csv'),
#              mut_rate = '../TEST/results/mut_rates/meanMutRatePerGR-LUAD--hg19.csv',
#              gene_name_synonyms = '../data/assets/hgnc_complete_set_processed.csv',
#              known_cancer_genes = '../data/assets/cgc_knownCancerGenes.csv',
#              olfactory_genes = '../data/assets/olfactory_barnes_2020.csv',
#             rawP_cap = 10^(-8),
#              expression_inventories = list('../data/inventory/inventory_expression_gtex.csv',
#                                            '../data/inventory/inventory_expression_tcga.csv'),
#              expression = list('../data/assets/GTEx_expression.csv',
#                                '../data/assets/TCGA_expression.csv'),
#              output = '../TEST/results/tables/combined_p_values/combinedP_LUAD-5primeUTR_hg19.csv')

# Read in raw results of software runs ----------------------------------------
message('[', Sys.time(), '] Started reading results of de-novo driver ',
        'discovery')
rawPs <- lapply(names(args$run_results), 
                function(x) readSoftwareResults(x, names(SOFTWARE_GR_CODES),
                                                args$run_results[[x]],
                                                unlist(args)))
rawPs <- do.call(rbind, rawPs)
message('[', Sys.time(), '] Finished reading results of de-novo driver ',
        'discovery')

# Read in mutation rate table -------------------------------------------------
message('[', Sys.time(), '] Started reading mutation rate table')
mutRate <- fread(args$mut_rate, header = T, stringsAsFactors = F)
message('[', Sys.time(), '] Finished reading mutation rate table')

# Read in gene name synonyms --------------------------------------------------
GENE_NAME_SYNS <- NULL
if (!is.null(args$gene_name_synonyms)) {
  message('[', Sys.time(), '] Started reading gene names synonyms table')
  GENE_NAME_SYNS <- fread(args$gene_name_synonyms, header = T, 
                          stringsAsFactors = F)
  message('[', Sys.time(), '] Finished reading gene names synonyms table')
}

# Unify gene_id and gene_name across various software -------------------------
# create gene_id to gene_name map from bed12 and from rawPs
GENE_ID_TO_NAME <- rawPs[,.(gene_id, gene_name)]
GENE_ID_TO_NAME <- GENE_ID_TO_NAME[complete.cases(GENE_ID_TO_NAME)]
GENE_ID_TO_NAME <- rbind(GENE_ID_TO_NAME, mutRate[,.(gene_name, gene_id)])
GENE_ID_TO_NAME <- unique(GENE_ID_TO_NAME) 
rawPs <- fillInGeneIDs(DT = rawPs, id_To_Name_Map = GENE_ID_TO_NAME, 
                       name_syns = GENE_NAME_SYNS)
rawPs <- fillInGeneName(rawPs, GENE_ID_TO_NAME)

rm(GENE_ID_TO_NAME, GENE_NAME_SYNS)
gc()

nbefore <- nrow(rawPs)
nafter <- nrow(rawPs[!is.na(gene_id)])
if (nafter != nbefore) {
  nRemovedGenes <- length(unique(rawPs[is.na(gene_id)]$gene_name))
  nRemovedGenesBySoft <- rawPs[is.na(gene_id)][,.(length(unique(gene_name))), 
                                               by = software]
  setnames(nRemovedGenesBySoft, 'V1', 'N genes')
  message('[', Sys.time(), '] Removed ', nbefore - nafter, ' entries out of ',
          nbefore, ' (', 100 * round((nbefore - nafter) / nbefore, 4), '%) ',
          'corresponding to ', nRemovedGenes, ' (',
          100 * round(nRemovedGenes / length(unique(rawPs$gene_name)), 4),
          '%) genes from the raw results table due to the absence of gene_id.',
          ' This can happen due to removal of genes from the set of genomic ',
          'regions of interest due to overlap with blacklisted regions. Yet, ',
          'MutPanning regions are not controlled by the pipeline, and ',
          'therefore those regions could be retained:')
  message(paste0(capture.output(knitr::kable(nRemovedGenesBySoft, 
                                             format = "markdown")), 
                 collapse = '\n'))
  rawPs <- rawPs[!is.na(gene_id)]
}

nbefore <- nrow(rawPs)
nafter <- nrow(rawPs[!is.na(gene_name)])
if (nafter != nbefore) {
  nRemovedGenes <- length(unique(rawPs[is.na(gene_name)]$gene_id))
  nRemovedGenesBySoft <- rawPs[is.na(gene_name)][,.(length(unique(gene_id))), 
                                                 by = software]
  setnames(nRemovedGenesBySoft, 'V1', 'N genes')
  message('[', Sys.time(), '] Removed ', nbefore - nafter, ' entries out of ',
          nbefore, ' (', 100 * round((nbefore - nafter) / nbefore, 4), '%) ',
          'corresponding to ', nRemovedGenes, ' (',
          100 * round(nRemovedGenes / length(unique(rawPs$gene_id)), 4),
          '%) genes from the raw results table due to the absence of ',
          'gene_name. This can happen due to removal of genes from the set of ',
          'genomic regions of interest due to overlap with blacklisted ',
          'regions. Yet, MutPanning regions are not controlled by the ',
          'pipeline, and therefore those regions could be retained:')
  message(paste0(capture.output(knitr::kable(nRemovedGenesBySoft, 
                                             format = "markdown")), 
                 collapse = '\n'))
  rawPs <- rawPs[!is.na(gene_name)]
}

# Cap raw p-values ------------------------------------------------------------
# As stated in PCAWG paper methods: raw P values smaller than 10-16 were 
# trimmed to that value before proceeding.
if (!is.null(args$rawP_cap)) {
  message('[', Sys.time(), '] Capping raw p values < ', args$rawP_cap, ' to ',
          args$rawP_cap)
  rawPs[raw_p < args$rawP_cap]$raw_p <- args$rawP_cap
}

# Combine p-values: Fisher, Brown, Stouffer, harmonic -------------------------
combinedPs <- meltToWide(rawPs)
if (length(unique(rawPs[!is.na(raw_p)]$software)) == 1) {
  message('[', Sys.time(), '] Results from only 1 software (',
          unique(rawPs[!is.na(raw_p)]$software), ') are present. Can not ',
          'perform p-values combination.')
} else {
  # remove entries where only 1 gene was found significant, because combined 
  # p-values can't be calculated on 1x1 matrices
  if (nrow(rawPs) == 1) {
    stop('[', Sys.time(), '] Detected gene only 1 gene. Can not perform ',
         'p-values combination.')
  }
  message('[', Sys.time(), '] Started combining p-values')
  combinedPs <- as.data.table(combineRawPvalues(combinedPs))
  message('[', Sys.time(), '] Finished combining p-values')
}

# Read & annotate with known cancer genes -------------------------------------
if (!is.null(args$known_cancer_genes)) {
  message('[', Sys.time(), '] Started reading known cancer genes table')
  known_cancer <- fread(args$known_cancer_genes, header = T, 
                        stringsAsFactors = F)
  known_cancer <- known_cancer[!duplicated(known_cancer[,.(gene_id,
                                                           gene_name)])]
  message('[', Sys.time(), '] Finished reading known cancer genes table')
  message('[', Sys.time(), '] Started annotation with known cancer genes ',
          'table')
  combinedPs <- mergeBasedOnAtLeastOneKey(combinedPs, known_cancer, 
                                          keyCols = c('gene_id', 'gene_name'))
  combinedPs[is.na(is_known_cancer)]$is_known_cancer <- F
  message('[', Sys.time(), '] Finished annotation with known cancer genes ',
          'table')
} else {
  message('[', Sys.time(), '] Known cancer genes table is not given. Will ',
          'set value of column is_known_cancer to F.')
  combinedPs[, is_known_cancer := F]
}

# Annotate with olfactory gene status -----------------------------------------
if (!is.null(args$olfactory_genes)) {
  message('[', Sys.time(), '] Started reading olfactory genes table')
  olfactory <- fread(args$olfactory_genes, header = T, sep = '\t', 
                     stringsAsFactors = F, select = c('gene_id', 'gene_name'))
  olfactory <- olfactory[!duplicated(olfactory[,.(gene_id, gene_name)])]
  olfactory[, is_olfactory := T]
  message('[', Sys.time(), '] Finished reading olfactory genes table')
  combinedPs <- mergeBasedOnAtLeastOneKey(combinedPs, olfactory, 
                                          keyCols = c('gene_id', 'gene_name'))
  combinedPs[is.na(is_olfactory)]$is_olfactory <- F
}

# Annotate with length & quantile of length, N mutations, etc -----------------
combinedPs <- merge(combinedPs, mutRate[, setdiff(colnames(mutRate), 
                                                  c('gr_name', 'gene_name')), 
                                        with = F], allow.cartesian = T,
                    all.x = T, by = c('gr_id', 'gene_id'))

mismatched <- combinedPs[is.na(nMuts)]
mismatched <- mismatched[apply(mismatched[, grep('.raw_p$', 
                                                 colnames(mismatched)), 
                                          with = F], 1, 
                               function(x) any(x[!is.na(x)] < 0.05))]
mismatched <- mismatched[, c('tumor_subtype', 'gr_id', 'gene_id', 'gene_name',
                             grep('raw_p$', colnames(mismatched), value = T)), 
                         with = F]
mismatched <- melt(mismatched, id.vars = c('tumor_subtype', 'gr_id',
                                           'gene_id', 'gene_name'))
setnames(mismatched, 'variable', 'software')
mismatched <- mismatched[!is.na(value)]
mismatched[, software := gsub('.raw_p$', '', software)]
mismatched <- unique(mismatched[,.(tumor_subtype, gr_id, gene_name, software)])

if (nrow(mismatched) > 0) {
  message('[', Sys.time(), '] Found ', nrow(mismatched), '(',
          round(100 * nrow(mismatched) / 
                  nrow(unique(rawPs[,.(tumor_subtype, gr_id, gene_name, 
                                       software)])), 4), '%) ', 
          'genomic region - software combinations where raw p value is ',
          'significant, but no mutations were assigned to that region. This ',
          'can be cause by slight misalignment of annotations between ',
          'softwares and genomic ranges. dNdScv and DIGdriver are prone to ',
          'that due to RefCDS creation and inability to remove positions ',
          'partically affected by low mappability, etc. MutPanning is prone ',
          'to that due do genome annotation being fixed inside the software. ',
          'Overview of distribution across software: ',
          paste0(capture.output(knitr::kable(mismatched[,.N, by = software], 
                                             format = "markdown")), 
                 collapse = '\n'), 
          '\n\nOverview of distribution across genomic regions: ',
          paste0(capture.output(knitr::kable(mismatched[,.N, 
                                                        by = gr_id][order(-N)], 
                                             format = "markdown")), 
                 collapse = '\n'))
}

rm(mutRate)
gc()
message('[', Sys.time(), '] Annotated with mutation rate statistics')

# Annotate with expression status (True/False/NA) -----------------------------
if (!is.null(args$expression)) {
  for (eIdx in 1:length(args$expression)) {
    combinedPs <- annotateWithExpression(combinedPs, 
                                         args$expression_inventories[[eIdx]],
                                         args$expression[[eIdx]])
  }
}

# Save raw results as table ---------------------------------------------------
write.table(combinedPs, args$output, append = F, quote = F,  sep = '\t', 
            row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote raw drivers table to ', args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')