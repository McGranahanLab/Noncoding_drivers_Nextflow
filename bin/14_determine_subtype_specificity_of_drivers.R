#!/usr/bin/env Rscript
# FILE: determine_subtype_specificity_of_drivers.R ----------------------------
#
# DESCRIPTION: Determines tumor subtype specificity (specific, preferential or
# not enough evidence) of cancer driver genomic elements detected in the tumor
# subtypes of uniform histological origin. This method is implemented as binary
# method of assigning tumor subtype specificity does not work in case cohort 
# sizes are very different.
#
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla determine_subtype_specificity_of_drivers.R -h
#          to see the full list of options and their descriptions.
# EXAMPLE: 
# REQUIREMENTS: 
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  17.12.2020
# REVISION: 27.12.2023

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
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))

# Functions : comparison of dNdScv/NBR selection rates post-hoc (LRT test) ----
#' compare_dNdScv_selection
#' @description Compares strength of selection for a gene in two different 
#' tumor subtypes or two genes in the same tumor subtype. For coding regions
#' only.
#' @author Maria Litovchenko & Ariana Huebner
#' @param dt data table with columns n_mis, n_non, n_spl, exp_mis, exp_non, 
#' exp_spl, wmis_global, wtru_global and 2 rows describing the same gene in two
#' different tumor subtypes or two different genes in the same tumor subtypes
#' @return data.table with columns p.value, r_mis, r_tru
compare_dNdScv_selection <- function(dt) {
  if (any(is.na(dt$wmis_global))) {
    return(data.table(p.value = NA, r_mis = NA, r_tru = NA))
  }
  
  # MLE dN/dS ratios using the uniform model under H0 and H1
  wmis_mle0 <- apply(dt[,.(exp_mis, wmis_global)], 1, 
                     function(x) as.numeric(x[1])*as.numeric(x[2]))
  wmis_mle0 = sum(dt$n_mis) / sum(wmis_mle0)
  
  wmis_mle1 <- apply(dt[,.(exp_mis, wmis_global)], 1,
                     function(x) as.numeric(x[1])*as.numeric(x[2]))
  wmis_mle1 = dt$n_mis / wmis_mle1
  
  wtru_mle0 <- apply(dt[,.(exp_non, exp_spl, wtru_global)], 1,
                     function(x) as.numeric(x[1:2])*as.numeric(x[3]))
  wtru_mle0 <- sum(rowSums(dt[,.(n_non, n_spl)])) / sum(wtru_mle0)
  
  wtru_mle1 <- apply(dt[,.(exp_non, exp_spl, wtru_global)], 1,
                     function(x) sum(as.numeric(x[1:2])*as.numeric(x[3])))
  wtru_mle1 <- apply(dt[,.(n_non, n_spl)], 1, sum) / wtru_mle1
  
  # Observed and predicted counts under H0 and H1
  obs <- c(dt$n_mis[1], sum(dt[1,.(n_non, n_spl)]), 
           dt$n_mis[2], sum(dt[2,.(n_non, n_spl)]))
  exp0 <- c(dt$exp_mis[1]*dt$wmis_global[1]*wmis_mle0, 
            sum(dt[1,.(exp_non, exp_spl)])*dt$wtru_global[1]*wtru_mle0, 
            dt$exp_mis[2]*dt$wmis_global[2]*wmis_mle0, 
            sum(dt[2,.(exp_non, exp_spl)])*dt$wtru_global[2]*wtru_mle0)
  # Note that exp1 == obs (we only have this line here for confirmation purposes)
  exp1 <- c(dt$exp_mis[1]*dt$wmis_global[1]*wmis_mle1[1], 
            sum(dt[1,.(exp_non, exp_spl)])*dt$wtru_global[1]*wtru_mle1[1],
            dt$exp_mis[2]*dt$wmis_global[2]*wmis_mle1[2], 
            sum(dt[2,.(exp_non, exp_spl)])*dt$wtru_global[2]*wtru_mle1[2])
  
  ll0 <- c(sum(dpois(x = obs[c(1,3)], lambda = exp0[c(1,3)], log = T)),
           sum(dpois(x = obs[c(2,4)], lambda = exp0[c(2,4)], log = T)))
  ll1 <- c(sum(dpois(x = obs[c(1,3)], lambda = exp1[c(1,3)], log = T)), 
           sum(dpois(x = obs[c(2,4)], lambda = exp1[c(2,4)], log = T)))
  
  # One-sided p-values
  pvals = (1 - pchisq(2 * (ll1 - ll0), df = 1))
  if (any(is.na(pvals))) {
    result <- data.table(p.value = NA, r_mis = NA, r_tru = NA)
  } else {
    if (wmis_mle1[1] < wmis_mle1[2]) { pvals[1] = 1 } 
    else { pvals[1] = pvals[1]/2 }
    if (wtru_mle1[1] < wtru_mle1[2]) { pvals[2] = 1 } 
    else { pvals[2] = pvals[2]/2 }
  }
  
  result <- data.table(p.value = 1 - pchisq(-2 * sum(log(pvals)), df = 4),
                       r_mis = wmis_mle1[1] / wmis_mle1[2], 
                       r_tru = wtru_mle1[1] / wtru_mle1[2])
  result
}

#' compare_nbr_selection
#' @description Compares strength of selection for a gene in two different 
#' tumor subtypes or two genes in the same tumor subtype. For noncoding regions
#' only.
#' @author Maria Litovchenko & Ariana Huebner
#' @param dt data table with columns n_sub, exp_mis, exp_sub, wsub_global and 2
#' rows describing the same gene in twodifferent tumor subtypes or two 
#' different genes in the same tumor subtypes
#' @return data.table with columns p.value, r_sub
compare_nbr_selection <- function(dt) {
  if (any(is.na(dt$wsub_global))) {
    return(data.table(p.value = NA, r_sub = NA))
  }
  
  # MLE dN/dS ratios using the uniform model under H0 and H1
  wsub_mle0 <- apply(dt[,.(exp_subs, wsub_global)], 1, 
                     function(x) as.numeric(x[1])*as.numeric(x[2]))
  wsub_mle0 = sum(dt$n_sub) / sum(wsub_mle0)
  
  wsub_mle1 <- apply(dt[,.(exp_subs, wsub_global)], 1,
                     function(x) as.numeric(x[1])*as.numeric(x[2]))
  wsub_mle1 = dt$n_sub / wsub_mle1
  
  # Observed and predicted counts under H0 and H1
  obs <- c(dt$n_sub[1], dt$n_sub[2])
  exp0 <- c(dt$exp_subs[1]*dt$wsub_global[1]*wsub_mle0,  
            dt$exp_subs[2]*dt$wsub_global[2]*wsub_mle0)
  # Note that exp1 == obs (we only have this line here for confirmation purposes)
  exp1 <- c(dt$exp_subs[1]*dt$wsub_global[1]*wsub_mle1[1],  
            dt$exp_subs[2]*dt$wsub_global[2]*wsub_mle1[2])
  
  ll0 <- c(sum(dpois(x = obs[c(1)], lambda = exp0[c(1)], log = T)),
           sum(dpois(x = obs[c(2)], lambda = exp0[c(2)], log = T)))
  ll1 <- c(sum(dpois(x = obs[c(1)], lambda = exp1[c(1)], log = T)), 
           sum(dpois(x = obs[c(2)], lambda = exp1[c(2)], log = T)))
  
  # One-sided p-values
  pvals = (1 - pchisq(2 * (ll1-ll0), df = 1))
  if (any(is.na(pvals))) {
    result <- data.table(p.value = NA, r_sub = NA)
  } else {
    if (wsub_mle1[1] < wsub_mle1[2]) { pvals[1] = 1 } 
    else { pvals[1] = pvals[1]/2 } 
    
    result <- data.table(p.value = 1 - pchisq(-2 * sum(log(pvals)), df = 4),
                         r_sub = wsub_mle1[1] / wsub_mle1[2])
    
  }
  result
}

#' compare_selection
#' @description Compares strength of selection for a gene in two different 
#' tumor subtypes or two genes in the same tumor subtype
#' @author Maria Litovchenko & Ariana Huebner
#' @param rawResDT data table with columns tumor_subtype, gr_id, gene_id, 
#' n_mis/n_sub, n_non, n_spl, exp_mis/exp_sub, exp_non, exp_spl, 
#' wmis_global, wsub_global, wtru_global and 2 rows describing the same
#' gene in two different tumor subtypes or two different genes in the same 
#' tumor subtypes
#' @return data.table with columns tumor_subtype_1/tumor_subtype_2 or 
#' gene_id_1/gene_id_2, gr_id, p.value, r_mis, r_tru
#' @note We can implement a simple LRT model based on the uniform dNdS model.
#' This is different from the Fisher test in that it uses synonymous mutations
#' (i.e. dN/dS ratios) instead of comparing the contribution of nonsyn muts of 
#' a gene relative to other genes. Being a uniform model it assumes no 
#' considerable changes in the mutation rate variation or coverage across genes
#' in both datasets. But takes into account signature and rate variation 
#' between two datasets.
#' H0: wmis1==wmis2 & wtru1==wtru2
#' H1: wmis1!=wmis2 & wtru1!=wtru2
compare_selection <- function(rawResDT) {
  if (nrow(rawResDT) != 2) {
    stop('[', Sys.time(), '] compare_dNdScv_selection: should have 2 rows')
  }
  dndscvCols <- c('n_mis', 'n_non', 'n_spl', 'exp_mis','exp_non', 'exp_spl', 
                  'wmis_global', 'wtru_global')
  nbrCols <- c('n_subs', 'exp_subs','exp_indels', 'wsub_global')
  if (length(setdiff(dndscvCols, colnames(rawResDT))) > 0 & 
      length(setdiff(nbrCols, colnames(rawResDT))) > 0 ) {
    stop('[', Sys.time(), '] compare_selection: needs following columns: ',
         paste0(c('n_mis', 'n_non', 'n_spl', 'exp_mis','exp_non', 'exp_spl',
                  'wmis_global', 'wtru_global'), collapse = ', '), ' OR ',
         paste0(c('n_subs', 'exp_subs','exp_indels', 'wsub_global'),
                collapse = ', '))
  }
  
  if (all(!complete.cases(rawResDT[, nbrCols, with = F]))) {
    result <- compare_dNdScv_selection(rawResDT)
  } else {
    result <- compare_nbr_selection(rawResDT)
  }
  
  # add info about tumor subtype, gene_id, etc
  if (length(unique(rawResDT$tumor_subtype)) != 1) {
    result <- cbind(tumor_subtype_1 = rawResDT$tumor_subtype[1], 
                    tumor_subtype_2 = rawResDT$tumor_subtype[2], 
                    gene_id = unique(rawResDT$gene_id),
                    gr_id = unique(rawResDT$gr_id), result)
  } else {
    if (length(unique(rawResDT$gene_id)) != 1) {
      result <- cbind(gene_id_1 = rawResDT$gene_id[1], 
                      gene_id_2 = rawResDT$gene_id[2], 
                      tumor_subtype = unique(rawResDT$tumor_subtype),
                      gr_id = unique(rawResDT$gr_id), result)
    } else {
      stop('[', Sys.time(), '] a data table of selection values of two ',
           'different genes in two different tissues is submitted. ',
           'Such genes should not be compared.')
    }
  }
  setcolorder(result,
              intersect(colnames(result),
                        c('tumor_subtype', 'tumor_subtype_1', 
                          'tumor_subtype_2', 'gene_id', 'gene_id_1', 
                          'gene_id_2', 'gr_id', 'p.value', 'r_mis', 'r_tru')))
  result
}

#' assign_specificity
#' @description Assigns tumor subtype specificity for the pair wise comparison 
#' of 1 gene across 2 tumor subtypes
#' @param pwCompDT data table with 2 rows and essential columns: 
#' tumor_subtype_1, is_driver_ts_1, p.value
#' @param pCutOff numeric, cut off on p-value
#' @return  save data table with added columns tumor_subtype_spec and 
#' specificity_mode. tumor_subtype_spec lists tumor subtype to which driver 
#' has specificity, and specificity_mode tells the mode of specificity, one of
#' 'not enough evidence', 'common', 'preferential', 'specific'. 
#' tumor_subtype_spec will be empty ('') if specificity_mode is 
#' 'not enough evidence' or 'common'.
assign_specificity <- function(pwCompDT, pCutOff = 0.05) {
  result <- copy(pwCompDT)
  result[, tumor_subtype_spec := '']
  result[, specificity_mode := '']
  
  if (sum(pwCompDT$p.value < pCutOff, na.rm = T) == 0) {
    # no need to check is_driver_ts_2 because it's the same gene in 2 tumor
    # subtypes
    res_spec <- switch(sum(pwCompDT$is_driver_ts_1), 'not enough evidence',
                       'common')
    if (!is.null(res_spec)) {
      result[, specificity_mode := res_spec]
    }
  }
  if (sum(pwCompDT$p.value < pCutOff, na.rm = T) == 1) {
    res_spec <- switch(sum(pwCompDT$is_driver_ts_1), 'specific', 
                       'preferential')
    if (!is.null(res_spec)) {
      if (res_spec == 'specific') {
        result[, tumor_subtype_spec := tumor_subtype_1[is_driver_ts_1]]
        result[, specificity_mode := 'specific']
      } else {
        result[, tumor_subtype_spec := tumor_subtype_1[p.value < pCutOff]]
        result[, specificity_mode := 'preferential']
      }
    }
  }
  result
}


# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'determine_subtype_specificity_of_drivers.R')

selectHelp <- paste('Path(s) to files containing calculated selection rates',
                     'for all cancer driver genomic elements detected in',
                     'tumor subtypes of single histological origin')
parser$add_argument("-s", "--selection_rates", required = T, nargs = '+', 
                    type = 'character', help = selectHelp)

pHelp <- paste('P value below which test is considered significant')
parser$add_argument("-p", "--p_val_max", required = T, default = 0.05,
                    type = 'double', help = pHelp)

outputHelp <- paste('Path to the output file')
parser$add_argument("-o", "--output", required = T, 
                    type = 'character', help = outputHelp)

args <- parser$parse_args()
# check_input_arguments_postproc(args, outputType = 'file')
names(args$run_result) <- unlist(args$gr_id)

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test inputs -----------------------------------------------------------------
# args <- list(selection_rates = c('selectionRates-Adenocarcinoma--hg19.csv',
#                                  'selectionRates-Squamous_cell--hg19.csv'),
#              p_val_max = 0.05, output = '')

# Read selection rates --------------------------------------------------------
selection_rates <- lapply(args$selection_rates, fread, header = T,
                          stringsAsFactors = F)
selection_rates <- do.call(rbind, selection_rates)

# during production of selectionRates-.* files driver genomic elements found
# in several histologically uniform subtypes will have their selection rates
# computed multiple times. Therefore, in order to perform tests on selection
# proper amount of times, it is better to separate here information about
# in which tumor subtypes genomic element was found a driver, and its 
# selection strenth in different subtypes
colsToKeep <- c('tumor_subtype', 'gr_id', 'gene_id', 'gene_name', 'tier',
                'is_known_cancer', 'known_in_tumor_subtype')
driverInfo <- selection_rates[, intersect(colsToKeep, 
                                          colnames(selection_rates)), with = F]
driverInfo <- unique(driverInfo)
selection_rates <- selection_rates[, setdiff(colnames(selection_rates),
                                             c('tumor_subtype', 'tier',
                                               'is_known_cancer', 
                                               'known_in_tumor_subtype')),
                                   with = F]
selection_rates <- unique(selection_rates)
setnames(selection_rates, 'scored_in_tumor_subtype', 'tumor_subtype')

# All pairwise combinations of tumor subtypes for for each driver element ----
subtype_combs <- selection_rates[,.(gr_id, gene_id, gene_name, tumor_subtype)]
subtype_combs <- unique(subtype_combs)
subtype_combs <- split(subtype_combs, by = c('gr_id', 'gene_id', 'gene_name'),
                       drop = T)
subtype_combs <- lapply(subtype_combs, 
                        function(x) expand.grid(gr_id = unique(x$gr_id),
                                                gene_id = unique(x$gene_id),
                                                gene_name = unique(x$gene_name),
                                                tumor_subtype_1 = unique(x$tumor_subtype),
                                                tumor_subtype_2 = unique(x$tumor_subtype)))
subtype_combs <- as.data.table(do.call(rbind, subtype_combs))
subtype_combs <- subtype_combs[tumor_subtype_1 != tumor_subtype_2]
subtype_combs[, id := paste(sort(unique(c(tumor_subtype_1, tumor_subtype_2))),
                            collapse = '-'), by = .(gr_id, gene_id, gene_name)]
subtype_combs <- subtype_combs[order(gr_id, gene_id, gene_name, 
                                     tumor_subtype_1, tumor_subtype_2)]
subtype_combs <- subtype_combs[!duplicated(subtype_combs[,.(gr_id, gene_id, 
                                                            gene_name, id)])]
subtype_combs[, id := NULL] 

# Compare selection strenth ---------------------------------------------------
subtype_combs <- split(subtype_combs, by = c('gr_id', 'gene_id', 'gene_name'),
                       drop = T)
selection_rates <- lapply(subtype_combs, 
                          function(x) merge(x[,.(gr_id, gene_id, gene_name)],
                                            selection_rates[tumor_subtype %in%
                                                              unlist(x)], 
                                            by = c('gr_id', 'gene_id',
                                                   'gene_name')))
# the test is one sides, that's why both direction should be scorred
select_comp <- lapply(selection_rates, 
                      function(x) rbind(compare_selection(x[order(tumor_subtype)]),
                                        compare_selection(x[order(-tumor_subtype)])))
select_comp <- as.data.table(do.call(rbind.fill, select_comp))

# Add back information about drivers ------------------------------------------
# add information, if gene is driver in one/two tumor subtypes
setnames(select_comp, 'tumor_subtype_1', 'tumor_subtype')
select_comp <- merge(select_comp, 
                     cbind(driverInfo[,.(tumor_subtype, gr_id, gene_id,
                                         gene_name)], is_driver_ts_1 = T),
                     all.x = T, by = c('tumor_subtype', 'gr_id', 'gene_id'))
setnames(select_comp, 'tumor_subtype', 'tumor_subtype_1')
select_comp[is.na(is_driver_ts_1)]$is_driver_ts_1 <- F
setnames(select_comp, 'tumor_subtype_2', 'tumor_subtype')
select_comp <- merge(select_comp, 
                     cbind(driverInfo[,.(tumor_subtype, gr_id, gene_id)], 
                           is_driver_ts_2 = T),  all.x = T,
                     by = c('tumor_subtype', 'gr_id', 'gene_id'))
setnames(select_comp, 'tumor_subtype', 'tumor_subtype_2')
select_comp[is.na(is_driver_ts_2)]$is_driver_ts_2 <- F
setcolorder(select_comp,
            c('tumor_subtype_1', 'tumor_subtype_2', 'gene_id', 'gr_id'))

# assign specificity based on pair wise comparisons 
select_comp[, tumor_subtype_pair := apply(select_comp[,.(tumor_subtype_1,
                                                         tumor_subtype_2)], 1, 
                                          function(x) paste(sort(x), 
                                                            collapse = '--'))]
select_comp <- split(select_comp, drop = T,
                     by = c('gr_id', 'gene_id', 'gene_name',
                            'tumor_subtype_pair'))
select_comp <- lapply(select_comp, assign_specificity, args$p_val_max)
select_comp <- do.call(rbind, select_comp)

# Write to output -------------------------------------------------------------
write.table(select_comp, args$output, append = F, quote = F, sep = '\t',
            row.names = F, col.names = T)
message('[', Sys.time(), '] Wrote tumor subtype specificity of detected in ',
        'histologically uniform cancer subtypes driver genomic elements to ',
        args$output)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')