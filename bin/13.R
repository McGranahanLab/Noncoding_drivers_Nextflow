# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))

# Functions : reading from dNdScv & NBR ---------------------------------------
#' readExpObsNmuts
#' @description Reads in result files from dndscv & nbr into uniformal data 
#' table while extracting columns containing information about expected and 
#' observed number of mutations. 
#' @author Maria Litovchenko
#' @param filePath path to result file
#' @param softName name of the software filePath corresponds to. dndscv or nbr
#' @param infoStr named vector with at least items cancer_subtype, gr_id
#' @return data table with columns software, tumor_subtype, gr_id +
#' n_syn, n_mis, n_non, n_spl, n_ind, mis_mle, tru_mle, mis_low, tru_low, 
#' mis_high, tru_high, ind_mle, ind_low, ind_high
readExpObsNmuts <- function(filePath, softName, infoStr) {
  if (!file.exists(filePath)) {
    message('[', Sys.time(), '] File not found: ', filePath, '. Return empty ',
            'data table.')
    return(data.table())
  }
  if (!softName %in% c('dndscv', 'nbr')) {
    stop('[', Sys.time(), '] readExpectedNmuts: softName should be one of ',
         'dndscv or nbr.', 'Current value: ', softName)
  }
  
  if (!all(c('cancer_subtype', 'gr_id') %in% names(infoStr))) {
    stop('[', Sys.time(), '] readExpectedNmuts: infoStr vector should have ',
         'items under names cancer_subtype, gr_id')
  }
  infoStr <- infoStr[c('cancer_subtype', 'gr_id')]
  infoStr <- as.data.table(t(infoStr))
  setnames(infoStr, 'cancer_subtype', 'tumor_subtype')
  
  if (softName == 'dndscv') {
    result <- fread(filePath, header = T, stringsAsFactors = F)
    result <- result[,.(gene_name, n_syn, n_mis, n_non, n_spl, n_ind, mis_mle,
                        tru_mle, mis_low, tru_low, mis_high, tru_high, ind_mle,
                        ind_low, ind_high, exp_syn, exp_mis, exp_non, exp_spl)]
    setnames(result, 'gene_name', 'gene_id')
    result <- cbind(result, infoStr)
  }
  
  if (softName == 'nbr') {
    result <- fread(filePath, header = T, stringsAsFactors = F,
                    select = c('region', 'obs_subs', 'obs_indels',
                               'obsexp_subs_mle', 'obsexp_subs_low',
                               'obsexp_subs_high', 'obsexp_indels_mle', 
                               'obsexp_indels_low', 'obsexp_indels_high',
                               'exp_subs', 'exp_indels'))
    setnames(result, 
             c('region', 'obs_subs', 'obs_indels', 'obsexp_subs_mle', 
               'obsexp_subs_low', 'obsexp_subs_high', 'obsexp_indels_mle', 
               'obsexp_indels_low', 'obsexp_indels_high'), 
             c('gene_id', 'n_subs', 'n_ind', 'subs_mle', 'subs_low', 
               'subs_high', 'ind_mle', 'ind_low', 'ind_high'), skip_absent = T)
    result <- cbind(result, infoStr)
  }
  
  result[, software := softName]
  
  result
}

# Functions : reading CHASM+ --------------------------------------------------
#' readCHASMplus
#' @description Reads results of CHASMplus into data table
#' @author Maria Litovchenko
#' @param chasmResPath path to results of CHASM+ 
#' @return data table with columns chr, pos, ref, alt, chasmPval, chasmScore
readCHASMplus <- function(chasmResPath) {
  message('[', Sys.time(), '] Started reading & processing ', chasmResPath)
  # read and parse results of chasmplus
  result <- fread(chasmResPath, header = T, stringsAsFactors = F, 
                  select = c('Chrom', 'Pos', 'Ref Base', 'Alt Base',
                             'P-value', 'Score'),
                  col.names = c('chr', 'pos', 'ref', 'alt', 'chasmPval',
                                'chasmScore'))
  result <- result[!is.na(chasmScore)]
  message('[', Sys.time(), '] Finished reading & processing ', chasmResPath)
  result
}

# Functions : % of driver mutations estimations -------------------------------
#' estimatePercOfDriverMuts
#' @description Estimates number and percentage of driver mutations per gene
#' @author Maria Litovchenko
#' @param inDT data table, result of run of either dNdScv or NBR. Should have
#' columns gene_id, tumor_subtype, gr_id, software, obsexp_subs_mle, 
#' obsexp_subs_low, obsexp_subs_high, obsexp_indels_mle, obsexp_indels_low, 
#' obsexp_indels_high for result of NBR and mis_mle, tru_mle, mis_low, tru_low,
#' mis_high, tru_high for result of dNdScv. 
#' @return data table with columns gene_id, tumor_subtype, gr_id, software,
#' var_type (indicates type of mutations, i.e. mis - missense), dN/dS ratios 
#' for each gene and mutation type (high - top level of 95th percentile 
#' confidence interval, low - bottom level of 95th percentile confidence 
#' interval, and mle - maximum likelihood estimation), observed_n - observed 
#' number of mutations, prob_is_driver_low - bottom level of 95th percentile 
#' confidence interval of estimated percentage of driver mutations, 
#' prob_is_driver_mle - maximum likelihood of estimated percentage of driver 
#' mutations, prob_is_driver_high - top level of 95th percentile confidence 
#' interval of estimated percentage of driver mutations
estimatePercOfDriverMuts <- function(inDT) {
  stopMsg <- paste('[', Sys.time(), '] estimatePercOfDriverMuts: please ',
                   'submit results of ONE software run. Software should be ',
                   'either dndscv or nbr.')
  if (length(unique(inDT$software)) != 1 |
      !all(inDT$software %in% c('dndscv', 'nbr'))) {
    stop('[', Sys.time(), '] estimatePercOfDriverMuts: please  submit ',
         'results of ONE software run. Software should be either dndscv or ',
         'nbr.')
  } 
  
  # columns identifying a gene driver
  idCols <- c('gene_id', 'tumor_subtype', 'gr_id', 'software')
  # columns with MLE
  mleCols <- '_mle$|_low$|_high$'
  # columns with observed number of mutations
  nObsCols <- '^n_'
  
  # retrieve MLE
  # columns with MLE of dNdS/dNdS-like (from NBR) ratio 
  colsOI <- paste0(paste0('^', idCols, '$|', collapse = ''), mleCols)
  mleDT <- inDT[, grepl(colsOI, colnames(inDT)), with = F]
  mleDT <- as.data.table(melt(mleDT, id.vars = idCols))
  # tru stands for nonsense in the original dNdScv results. Checked with 
  # function geneci
  setnames(mleDT, 'variable', 'measurement')
  mleDT[, measurement := gsub('tru', 'non', measurement)]
  # var_type will indicate a type of mutation
  mleDT[, var_type := gsub(mleCols, '',  measurement)]
  mleDT[, measurement := gsub('.*_', '', measurement)]
  # convert back to wide, it prevents a lot of bugs
  mleDT <- dcast(mleDT, gene_id + tumor_subtype + gr_id + 
                   var_type + software ~ measurement, value.var = 'value')
  
  # retrieve number of observed mutations
  nObsVars <- inDT[, grepl(paste0(paste0('^', idCols, '$|', collapse = ''), 
                                  nObsCols, collapse = ''), colnames(inDT)), 
                   with = F]
  nObsVars <- nObsVars[, !apply(nObsVars, 2, function(x) all(is.na(x))),
                       with = F]
  nObsVars <- as.data.table(melt(nObsVars, id.vars = idCols))
  setnames(nObsVars, 'value', 'observed_n')
  nObsVars[, var_type := gsub('.*_', '', variable)]
  nObsVars[, var_type := gsub('indused', 'ind', var_type)]
  nObsVars[, variable := NULL] 
  
  # merge them together and calculate % of driver mutations and number of 
  # driver mutations. For the origin of (value - 1) / value see dNdS paper,
  # section 'Estimating the number of substitutions fixed by positive selection
  # from dN/dS'
  result <- merge(mleDT, nObsVars, by = c(idCols, 'var_type'))
  result[, prob_is_driver_low := (low - 1) / low]
  result[prob_is_driver_low < 0 | 
           is.nan(prob_is_driver_low)]$prob_is_driver_low <- NA
  result[, prob_is_driver_mle := (mle - 1) / mle]
  result[prob_is_driver_mle < 0 | 
           is.nan(prob_is_driver_mle)]$prob_is_driver_mle <- NA
  result[, prob_is_driver_high := (high - 1) / high]
  result[prob_is_driver_high < 0 | 
           is.nan(prob_is_driver_high)]$prob_is_driver_high <- NA
  
  result[var_type == 'ind']$var_type <- 'indels'
  
  result
}

# Functions : number of tumors with a driver mutation(s) ----------------------
#' generateMutCombsWithMinMaxNtums
#' @description Generates mutation combinations which will correspond to 
#' minimum number of patients with driver mutation and to the maximum number of
#' patients with driver mutations.
#' @author Maria Litovchenko
#' @param data.table with columns key, n_tumors, where key is the unique ID of
#' a mutation and n_tumors is number of tumors in which it is found.
#' @param n_drivers integer, number of driver mutations.
#' @return data table with 2 columns: 1st column gives indices of mutations 
#' which in combination give the least number of patients with driver mutations
#' 2nd column the biggest number of patients with driver mutations
generateMutCombsWithMinMaxNtums <- function(tums_per_mut, n_drivers) {
  result <- copy(tums_per_mut)
  result[, idx_orig := 1:nrow(result)]
  result <- result[order(n_tumors)]
  result <- cbind(head(result, n_drivers)$idx_orig, 
                  tail(result, n_drivers)$idx_orig)
  result
}

#' generateAllMutCombs
#' @description Generates all possible combinations of mutations 
#' @param data.table with columns key, n_tumors, where key is the unique ID of
#' a mutation and n_tumors is number of tumors in which it is found.
#' @param n_drivers integer, number of driver mutations.
#' @return data table with number of columns equal 
generateAllMutCombs <- function(tums_per_mut, n_drivers) {
  result <- tryCatch({
    combn(1:nrow(tums_per_mut), n_drivers)
  }, error = function(cond) {
    message('[', Sys.time(), '] Can not compute all combinations of ',
            n_drivers, ' from ', nrow(tums_per_mut), '. Will return 2 ',
            'combinations: with maximum and minimum number of patients with ',
            'driver mutations.')
    generateMutCombsWithMinMaxNtums(tums_per_mut, n_drivers)
  })
  result
}

#' estimateNtumorsWithDriverMut
#' @description Estimates number of patients with driver mutations (using 
#' analytical solution!)
#' @author Maria Litovchenko
#' @param muts_in_gene data table with columns tumor_subtype, gr_id, gene_id,
#' gene_name, var_type, key, prob_is_driver_mle, participant_id, where key is
#' a unique mutations identifier, prob_is_driver_mle is essentially percentage
#' of driver mutations in gene estimated by dNdScv/NBR. The data table should
#' contains mutations of only one structural type, i.e. mis/non/sub/indels.
#' @return data table with columns tumor_subtype, gr_id, gene_id, gene_name and
#' n_tums_w_driver_mut - numerical, estimation of number of patients with
#' driver mutations.
#' @note description of analytical solution:
#' Let's say we have n mutations in a gene and we also know that k of them are 
#' drivers. However, we do not know which of mutations are drivers. Also, one
#' mutation of the set is present in 2 patients. Therefore, number of patients
#' with driver mutations can range from k to k+1. But this is just the range,
#' we don't know yet how many times it is possible that number of patients will
#' be k and how many k+1. So,
#' 1. number of mutations' combinations which will produce k patients with 
#' driver mutations:
#'    \binom{n}{k}
#' 2. number of mutations' combinations which will produce k+1 patients with 
#' driver mutations:
#'    \binom{n-1}{k-1} (because we essentially need to select k mutations from
#'    n - 1 set, as we already fixed mutation present in 2 patients to be in
#'    our set)
#' 3. Then % of all possible mutations' combinations which will produce k+1
#' patient is:
#'    \binom{n}{k}/\binom{n-1}{k-1} = k/n
#' 4. Therefore, mean number of patients with driver mutation can be computed
#' as:
#'    (k/n)*(k+1) + (1 - k/n)*k
#' This is an analytical solution for the case then there is just 1 mutation
#' present in > 1 patient is in the set. If there are more mutations like this,
#' then it all gets more complicated, especially if one mutation is present in,
#' for example 3 patients and the other one in 2 patients. Therefore, we 
#' address such cases by solving them programmatically. 
estimateNtumorsWithDriverMut <- function(muts_in_gene) {
  # number of estimated driver mutations
  n_driver_muts <- unique(muts_in_gene$prob_is_driver_mle)
  n_driver_muts <- round(n_driver_muts * length(unique(muts_in_gene$key)))
  if (n_driver_muts > 0) {
    # number of tumors in which each mutation is present
    n_tumors_per_mut <- unique(muts_in_gene[,.(key, participant_id)])
    n_tumors_per_mut <- data.table(table(n_tumors_per_mut$key))
    setnames(n_tumors_per_mut, c('V1', 'N'), c('key', 'n_tumors'))
    
    if (max(n_tumors_per_mut$n_tumors) != 1) {
      if (nrow(n_tumors_per_mut[n_tumors > 1]) == 1) {
        # analytical solution to preserve memory and CPU
        # percentage of mutations combinations which will have mutation present
        # in > 1 patient. This result is easy to prove by dividing C(n-1)(k-1)
        # by C(n)(k)
        perc_mut_more1 <- n_driver_muts/nrow(n_tumors_per_mut)
        # minimal and maximum n patients:
        mut_combs <- generateMutCombsWithMinMaxNtums(n_tumors_per_mut, 
                                                     n_driver_muts)
        n_tumors_per_mut_v <- n_tumors_per_mut$n_tumors
        result <- apply(mut_combs, 2, function(x) sum(n_tumors_per_mut_v[x]))
        rm(mut_combs)
        result <- perc_mut_more1*result[1] + (1 - perc_mut_more1)*result[2]
      } else {
        # all possible combinations of n_driver_muts number of mutations from 
        # all mutation in the gene
        mut_combs <- generateAllMutCombs(n_tumors_per_mut, n_driver_muts)
        if (!is.vector(mut_combs)) {
          n_tumors_per_mut_v <- n_tumors_per_mut$n_tumors
          # all possible numbers of patients with driver mutations
          result <- apply(mut_combs, 2, function(x) sum(n_tumors_per_mut_v[x]))
          rm(mut_combs)
          result <- median(result)
        } else {
          result <- n_driver_muts
        }
      }
    } else {
      result <- n_driver_muts
    }
  } else {
    result <- 0
  }
  result <- cbind(unique(muts_in_gene[,.(tumor_subtype, gr_id, gene_id, 
                                         gene_name, var_type)]), 
                  'n_tums_w_driver_mut' = result)
  result
}

# Functions : comparison of dNdScv/NBR selection rates post-hoc (LRT test) ----
#' calculateGlobalMisTru
#' @description computes global dN/dS ratios from all genes except gene of 
#' interest (to normalise the differences for the genes being tested in 
#' function compare_dNdScv_selection)
#' @author Maria Litovchenko & Ariana Huebner
#' @param raw_results_DT
#' @param goi gene_id of interest
#' @param goi_gr_id genomic region of regions of interest
#' @return data table with columns tumor_subtype, software, gr_id, gene_id, 
#' n_mis/n_sub, n_non, n_spl/n_ind, exp_mis/exp_sub, exp_non, 
#' exp_spl/exp_indels, wmis_global/wsub_global, wtru_global/wind_global
calculateGlobalMisTru <- function(raw_results_DT, goi, goi_gr_id) {
  result <- data.table(gene_id = goi, gr_id = goi_gr_id) 
  if ('exp_non' %in% colnames(raw_results_DT)) { #dNdScv
    wmis <- sum(raw_results_DT[gene_id != goi]$n_mis)
    wmis <- wmis / sum(raw_results_DT[gene_id != goi]$exp_mis)
    wtru <- sum(rowSums(raw_results_DT[gene_id != goi][,.(n_non, n_spl)]))
    wtru <- wtru / sum(rowSums(raw_results_DT[gene_id != goi][,.(exp_non, 
                                                                 exp_spl)]))
    result <- cbind(result, wmis_global = wmis, wtru_global = wtru)
  } else { # NBR
    wsub <- sum(raw_results_DT[gene_id != goi]$n_subs)
    wsub <- wsub / sum(raw_results_DT[gene_id != goi]$exp_subs)
    wind <- sum(raw_results_DT[gene_id != goi]$n_indels)
    wind <- wind / sum(raw_results_DT[gene_id != goi]$exp_indels)
    result <- cbind(result, wsub_global = wsub, wind_global = wind)
  }
  colsToKeep <- intersect(c('tumor_subtype', 'software', 'gene_id',
                            'n_mis', 'n_non', 'n_spl', 'n_subs', 'n_ind',
                            'exp_mis', 'exp_non', 'exp_spl', 'exp_subs', 
                            'exp_indels'), colnames(raw_results_DT))
  result <- merge(raw_results_DT[, colsToKeep, with = F], result,
                  by = 'gene_id')
  result
}

#' calculateGlobalMisTruList
#' @description Applies function calculateGlobalMisTru to list of data tables
#' containing raw results of driver discovery by dNdScv/NBR on one genomic 
#' region and one tumor subtype
#' @author Maria Litovchenko & Ariana Huebner
#' @param goi gene_id of interest
#' @param goi_gr_id genomic region of regions of interest
#' @param raw_results_list 
#' @return data table with columns tumor_subtype, software, gr_id, gene_id, 
#' n_mis, exp_mis, n_non, exp_non, n_spl, exp_spl, wmis_global, wtru_global
calculateGlobalMisTruList <- function(goi, goi_gr_id, raw_results_list) {
  raw_results_goi <- lapply(raw_results_list, 
                            function(x) x[gr_id == goi_gr_id])
  raw_results_goi <- raw_results_goi[sapply(raw_results_goi, 
                                            function(x) nrow(x) != 0)]
  raw_results_goi <- raw_results_goi[sapply(raw_results_goi, 
                                            function(x) goi %in% x$gene_id)]
  result <- lapply(raw_results_goi, calculateGlobalMisTru, goi, goi_gr_id)
  result <- as.data.table(do.call(rbind.fill, result))
  result
}

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
#' software, n_mis/n_sub, n_non, n_spl, exp_mis/exp_sub, exp_non, exp_spl, 
#' wmis_global, wsub_global, wtru_global and 2 rows describing the same
#' gene in two different tumor subtypes or two different genes in the same 
#' tumor subtypes
#' @return data.table with columns tumor_subtype_1/tumor_subtype_2 or 
#' gene_id_1/gene_id_2, software, gr_id, p.value, r_mis, r_tru
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
  codingCols <- c('n_mis', 'n_non', 'n_spl', 'exp_mis','exp_non', 'exp_spl', 
                  'wmis_global', 'wtru_global')
  noncodingCols <- c('n_subs', 'exp_subs','exp_indels', 'wsub_global')
  if (length(setdiff(codingCols, colnames(rawResDT))) > 0 & 
      length(setdiff(noncodingCols, colnames(rawResDT))) > 0 ) {
    stop('[', Sys.time(), '] compare_dNdScv_selection: needs following ',
         'columns: ', 
         paste0(c('n_mis', 'n_non', 'n_spl', 'exp_mis','exp_non', 
                  'exp_spl', 'wmis_global', 'wtru_global'), collapse = ', '),
         ' OR ',
         paste0(c('n_subs', 'exp_subs','exp_indels', 'wsub_global'),
                collapse = ', '))
  }
  
  if (unique(rawResDT$software) == 'dndscv') {
    result <- compare_dNdScv_selection(rawResDT)
  } else {
    if (unique(rawResDT$software) == 'nbr') {
      result <- compare_nbr_selection(rawResDT)
    } else {
      stop('[', Sys.time(), '] Can not handle data from software: ',
           unique(rawResDT$software))
    }
  }
  
  # add info about tumor subtype, gene_id, etc
  if (length(unique(rawResDT$tumor_subtype)) != 1) {
    result <- cbind(tumor_subtype_1 = rawResDT$tumor_subtype[1], 
                    tumor_subtype_2 = rawResDT$tumor_subtype[2], 
                    gene_id = unique(rawResDT$gene_id),
                    gr_id = unique(rawResDT$gr_id),
                    software = unique(rawResDT$software),
                    result)
  } else {
    if (length(unique(rawResDT$gene_id)) != 1) {
      result <- cbind(gene_id_1 = rawResDT$gene_id[1], 
                      gene_id_2 = rawResDT$gene_id[2], 
                      tumor_subtype = unique(rawResDT$tumor_subtype),
                      gr_id = unique(rawResDT$gr_id),
                      software = unique(rawResDT$software),
                      result)
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
                          'gene_id_2', 'gr_id', 'software','p.value', 'r_mis', 
                          'r_tru')))
  
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

# Functions : information messages --------------------------------------------
#' display_driverMut_update_msg
#' @description Prints info message about how many of driver mutations were 
#'              assigned to a certain confidence level
#' @author Maria Litovchenko
#' @param patientsInvDT data table, inventory of patients with participant_id
#'        column
#' @param driveMutsDT data table with driver mutations with column 
#'        confidenceLvl
#' @param confLvl string, confidence level
#' @return void
display_driverMut_update_msg <- function(patientsInvDT, driveMutsDT, confLvl) {
  # all available patients
  total_nParticip <- length(unique(patientsInvDT$participant_id))
  # number of patients for which a driver mutation at confLvl level of 
  # confidence is identified
  nExplained_lvl <- driveMutsDT[confidenceLvl == confLvl]
  nExplained_lvl <- length(unique(nExplained_lvl$participant_id))
  # total number of patients for which driver mutation was found
  nExplained_total <- length(unique(driveMutsDT[confidenceLvl != 'passenger']$participant_id))
  
  message('[', Sys.time(), '] Found driver mutations at ', confLvl, ' level ',
          'of confidence in ', nExplained_lvl, ' participants. Total current ',
          'number of participants with identified driver mutation: ',
          nExplained_total, '(', 
          100 * round(nExplained_total/total_nParticip, 4), '%)')
}

# Test arguments --------------------------------------------------------------
args <- list(inventory_analysis = 'data/inventory/inventory_analysis.csv', 
             inventory_patients = 'data/inventory/inventory_patients.csv',
             cancer_subtype = 'Panlung', gr_id = 'CDS', software = 'dndscv',
             run_result = 'test_fixed_dndscv_Panlung_CDS.csv', 
             chasmplus = 'chasmplus-results-PANCAN-CDS-hg19.csv', 
             drivers = 'completed_runs/2023-12-14/results/tables/drivers/drivers-Panlung--hg19.csv',
             inferred_biotype = 'completed_runs/2023-12-14/results/',
             muts_to_gr = 'completed_runs/2023-12-14/results/mut_rates/mutMapToGR-Panlung--hg19.csv',
             synAcceptedClass = 'Silent', 
             chasm_padj = 0.05, chasm_score_min = 0.5, max_fp = 0.05)
# known_driver_mutations = ''

if (!args$software %in% c('dndscv', 'nbr')) {
  stop('[', Sys.time(), '] Can not perform inference of number of driver ',
       'mutations on ', args$software, '. Please use results of dNdScv',
       '(coding regions) or NBR intead')
}

# Read in patients inventory --------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients, 1)
patientsInv <- patientsInv[tumor_subtype %in% args$cancer_subtype]
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)

# Read in analysis inventory --------------------------------------------------
analysisInv <- readAnalysisInventory(args$inventory_analysis)
message('[', Sys.time(), '] Read --inventory_analysis: ', 
        args$inventory_analysis)
analysisInv <- analysisInv[tumor_subtype == args$cancer_subtype]
# assign coding genomic regions - anything which contains CDS as gr_code
coding_gr_id <- unique(analysisInv[gr_code == 'CDS']$gr_id)
if (length(coding_gr_id) != 0) {
  message('[', Sys.time(), '] Following genomic regions: ', 
          paste0(coding_gr_id, collapse = ', '), ', will be considered as ',
          'coding.')
}
rm(analysisInv)
gc()

# Read driver genes -----------------------------------------------------------
if (!is.null(args$drivers)) {
  message('[', Sys.time(), '] Started reading ', args$drivers)
  drivers <- fread(args$drivers, header = T, stringsAsFactors = F, 
                   select = c('gr_id', 'gene_id', 'gene_name', 'FILTER', 
                              'tier', 'is_known_cancer', 
                              'known_in_tumor_subtype', 
                              'known_cancer_biotype'))
  drivers <- unique(drivers)
  drivers <- drivers[FILTER == 'PASS' & !is.na(tier)]
  
  message('[', Sys.time(), '] Finished reading ', args$drivers)
  if (nrow(drivers) == 0) {
    stop('[', Sys.time(), '] no significant (FILTER is PASS and tier is not ',
         'NA) driver genes is found in ', args$drivers, ' table.')
  }
  
  
  drivers <- drivers[gr_id %in% args$gr_id]
  
  "SOFT STOOOOOP"
}

# Read in observed and estimated number of driver mutations - dNdScv, NBR -----
mleDriveMuts <- readExpObsNmuts(args$run_result, args$software, unlist(args))
message('[', Sys.time(), '] Read ', args$run_result)
# restrict to drivers
mleDriveMuts <- mleDriveMuts[gene_id %in% unique(drivers$gene_id)]

# Estimate percentage and number of driver mutations --------------------------
mleDriveMuts <- estimatePercOfDriverMuts(mleDriveMuts)
message('[', Sys.time(), '] Finished estimating number of driver mutations')
mleDriveMuts[, tumor_subtype := NULL]

# Read mutation map to genomic regions & restrict to driver genes -------------
varsToGRmap <- fread(args$muts_to_gr, header = T, stringsAsFactors = F, 
                     select = c('participant_id', 'gr_id', 'gene_id', 
                                'gene_name', 'key', 'key_orig', 'var_class',
                                'struct_type'))
message('[', Sys.time(), '] Read ', args$muts_to_gr)

# restrict to driver genes
varsToGRmap <- merge(varsToGRmap, drivers[,.(gr_id, gene_id, gene_name)],
                     by = c('gr_id', 'gene_id', 'gene_name'))
message('[', Sys.time(), '] Mutations were restricted to mutations ',
        'which belong to identified driver genomic regions.')

# remove silent mutations, if requested
if (!is.null(args$synAcceptedClass)) {
  nbefore <- nrow(varsToGRmap)
  varsToGRmap <- varsToGRmap[!var_class %in% args$synAcceptedClass]
  nafter <- nrow(varsToGRmap)
  message('[', Sys.time(), '] Removed ', nbefore - nafter, ' mutations out ',
          'of ', nbefore, '(', 100*round((nbefore - nafter)/nbefore, 4), '%) ',
          'because they fell into ', 
          paste(args$synAcceptedClass, collapse = ', '), ' class(es).')
}

# sync var_type of varsToGRmap to the future mleDriveMuts derived from NBR/dNsScv
varsToGRmap[, var_type := character()]
varsToGRmap[struct_type %in% c('INS', 'DEL')]$var_type <- 'indels'
# ... cds 
varsToGRmap[gr_id %in% coding_gr_id & struct_type == 'SNP']$var_type <- 'mis'
varsToGRmap[var_class == 'Nonsense_Mutation']$var_type <- 'non'
# ... noncoding 
varsToGRmap[!gr_id %in% coding_gr_id & struct_type == 'SNP']$var_type <- 'subs'
# ... for MNPs - set up var_type to mnp, it won't be used for drivers as NBR/
# dNdScv does not score them
varsToGRmap[struct_type == 'MNP' & var_class == 'Unknown']$var_type <- 'mnp'

# remove MNPs from consideration
varsToGRmap <- varsToGRmap[!var_type %in% 'mnp']
message('[', Sys.time(), '] Removed MNPs from considereation')

# Assign a probability of being a driver to every mutation --------------------
mleDriveMuts <- merge(varsToGRmap, mleDriveMuts, all.x = T, 
                 by = c('gr_id', 'gene_id', 'var_type'))
rm(varsToGRmap)
gc()
message('[', Sys.time(), '] Assigned a probability of being a driver ',
        'mutation to mutations in driver genes.')

# deal with NAs
mleDriveMuts[is.na(prob_is_driver_low)]$prob_is_driver_low <- 0
mleDriveMuts[is.na(prob_is_driver_mle)]$prob_is_driver_mle <- 0
mleDriveMuts[is.na(prob_is_driver_high)]$prob_is_driver_high <- 0

# Read results of CHASMplus & add scores to mutations -------------------------
if (!is.null(args$chasmplus)) {
  message('[', Sys.time(), '] Started reading results of CHASM+ analysis')
  
  chasmDT <- readCHASMplus(args$chasmplus)
  chasmDT[, chasmPadj := p.adjust(chasmPval, method = 'BH')]
  chasmDT[, key := apply(chasmDT[,.(chr, pos, ref, alt)], 1, paste0, 
                          collapse = ':')]
  chasmDT[, key := gsub(' ', '', key)]
  chasmDT <- chasmDT[,.(key, chasmScore, chasmPadj)]
  chasmDT[, chasmScore := as.numeric(chasmScore)]
  chasmDT[, chasmPadj := as.numeric(chasmPadj)]
  
  # assign CHASMplus scores to mutations
  mleDriveMuts <- merge(mleDriveMuts, chasmDT, by = 'key', all.x = T)
  # deal with NAs
  mleDriveMuts[is.na(chasmScore)]$chasmScore <- 0
  mleDriveMuts[is.na(chasmPadj)]$chasmPadj <- 1
  rm(chasmDT)
  gc()
  message('[', Sys.time(), '] Finished reading results of CHASM+ analysis')
} else {
  mleDriveMuts[, chasmScore := 0]
  mleDriveMuts[, chasmPadj := 1]
}

# Read known driver mutations & annotate mutations with them ------------------
if (!is.null(args$known_driver_mutations)) {
  message('[', Sys.time(), '] Started annotation with known driver mutations')
  if (args$known_mut_drivers_ovrl_mode == 'location') {
    mleDriveMuts[, key_loc := gsub(':[ATGC-].*', '', key)]
    mleDriveMuts <- merge(mleDriveMuts, 
                         unique(knownDriverMuts[,.(gene_name, key_loc,
                                                   is_known_driver_mut)]),
                         all.x = T, by = c('gene_name', 'key_loc'))
  } else {
    mleDriveMuts <- merge(mleDriveMuts, knownDriverMuts, all.x = T,
                         by = c('gene_name', 'key'))
  }
  mleDriveMuts[var_type == 'indels']$is_known_driver_mut <- F
  mleDriveMuts[is.na(is_known_driver_mut)]$is_known_driver_mut <- F
  message('[', Sys.time(), '] Finished annotation with known driver mutations')
  gc()
} else {
  mleDriveMuts[, is_known_driver_mut := F]
}

# Initialize driver mutations data table --------------------------------------
driverMutations <- data.table(mkey = character(), gr_id = character(),
                              gene_id = character(), gene_name = character(),
                              var_type = character(), 
                              participant_id = character(), 
                              key_orig = character(), var_class = character(),
                              struct_type = character(),
                              software = character(), high = numeric(), 
                              low = numeric(), mle = numeric(), 
                              observed_n = numeric(),
                              prob_is_driver_low = numeric(),  
                              prob_is_driver_mle = numeric(), 
                              prob_is_driver_high = numeric(), 
                              chasmScore = numeric(), chasmScore = numeric(),
                              is_known_driver_mut = logical(), 
                              confidenceLvl = character())
setnames(driverMutations, 'mkey', 'key')

# Identify MISSENSE DRIVER mutation based on dNdScv & CHASM & known muts ------
# Select missence mutations. We do selection based on var_type = 'mis' and not
# by var_class, because var_class can be Missense_Mutation, Nonstop_Mutation,
# Translation_Start_Site or Unknown for var_type = 'mis'. Yet, the same dNdScv
# probability as for Missense_Mutation will be assigned to for example 
# Nonstop_Mutation. It could even be scored by CHASM+, which normally works on
#  missence mutations only!
missen_muts <- mleDriveMuts[gr_id %in% coding_gr_id & var_type == 'mis']

# compute total number of missense mutations per driver gene & preserve mle
# probabilities of being found as driver.
missen_muts[, missen_total := length(unique(key)),
            by = .(gr_id, gene_id, gene_name)]

# calculate number of missence driver mutations identified by CHASM+. CHASM+
# scores only missence mutations.
missen_muts[, inCHASM := F]
missen_muts[chasmScore >= args$chasm_score_min & 
              chasmPadj <= args$chasm_padj]$inCHASM <- T
missen_muts[, chasm_n_driver := length(unique(key[inCHASM])),
            by = .(gr_id, gene_id, gene_name)]
missen_muts[, inCHASM := NULL]

# calculate number of known driver mutations, NOT identified by CHASM+ in the
# genes
missen_muts[, known_n_driver := 0]
if (!is.null(args$known_driver_mutations)) {
  missen_muts[, known_not_chasm := F]
  missen_muts[is_known_driver_mut %in% T & 
                (chasmScore < args$chasm_score_min |
                   chasmPadj > args$chasm_padj)]$known_not_chasm <- T
  missen_muts[, known_n_driver := length(unique(key[known_not_chasm])),
              by = .(gr_id, gene_id, gene_name)] 
  missen_muts[, known_not_chasm := NULL]
}

# calculate estimates of number of missense DRIVER mutations ACROSS ALL 
# PATIENTS for each gene (since it's from probabilities - it is from dNdScv)
missen_muts[, dnds_n_driver_low := prob_is_driver_low * missen_total]
missen_muts[, dnds_n_driver_mle := prob_is_driver_mle * missen_total]
missen_muts[, dnds_n_driver_high := prob_is_driver_high * missen_total]

# Now, ACROSS ALL PATIENTS we have 2 numbers per gene: number of drivers from 
# dNdScv and sum of number of drivers from CHASM+ and known driver mutation 
# (not identified as significant by CHASM+). We want to produce conservative
# estimation. Therefore, for genes for which number of driver mutations from 
# dNdScv is bigger than sum of number of driver mutations from CHASM+ and known
# driver mutations, we will take as drivers mutations mutations identified by 
# CHASM+ & known driver mutations.
missen_muts[, status := '']
missen_muts[dnds_n_driver_mle > 0 & chasm_n_driver + known_n_driver != 0 &
              chasm_n_driver + known_n_driver <= dnds_n_driver_mle &
              status == '']$status <- 'chasm_less_dnds'
if (nrow(missen_muts[status == 'chasm_less_dnds']) > 0) {
  chasm_less_dnds <- missen_muts[status == 'chasm_less_dnds']
  chasm_less_dnds[, confidenceLvl := 'passenger']
  chasm_less_dnds[(chasmScore >= args$chasm_score_min &
                     chasmPadj <= args$chasm_padj) |
                    is_known_driver_mut %in% T]$confidenceLvl <- 'dNdScv&CHASM+&known'
  chasm_less_dnds[(chasmScore >= args$chasm_score_min &
                     chasmPadj <= args$chasm_padj) |
                    is_known_driver_mut %in% F]$confidenceLvl <- 'dNdScv&CHASM+'
  # remove extra columns
  chasm_less_dnds <- chasm_less_dnds[, c(colnames(mleDriveMuts),
                                         'confidenceLvl'), with = F]
  
  driverMutations <- rbind(driverMutations, chasm_less_dnds, fill = T)

  # info message
  driverMutations[, participant_id := as.character(participant_id)]
  display_driverMut_update_msg(patientsInv, driverMutations, 
                               confLvl = 'dNdScv&CHASM+&known')
  display_driverMut_update_msg(patientsInv, driverMutations,
                               confLvl = 'dNdScv&CHASM+')
}

# If it's the other way around (CHASM & known driver mutations > dNdScv), then
# we can order known driver mutations and CHASM+ results by significance and 
# select out of them top N, there N is the number from dNdScv.
missen_muts[dnds_n_driver_mle > 0 & 
              chasm_n_driver + known_n_driver > dnds_n_driver_mle &
              status == '']$status <- 'chasm_more_dnds'
if (nrow(missen_muts[status == 'chasm_more_dnds']) > 0) {
  dnds_less_chasm <- missen_muts[status == 'chasm_more_dnds']
  dnds_less_chasm[, dnds_n_driver_mle := round(dnds_n_driver_mle)]

  # sort mutation by probability of being drivers
  dnds_less_chasm[, confidenceLvl := 'passenger']
  dnds_less_chasm[is_known_driver_mut %in% T]$confidenceLvl <- 'known'
  dnds_less_chasm[chasmScore >= args$chasm_score_min &
                     chasmPadj <= args$chasm_padj]$confidenceLvl <- 'CHASM+'
  dnds_less_chasm[, confidenceLvl := factor(confidenceLvl, 
                                            levels = c('passenger', 'CHASM+',
                                                       'known'))]
  dnds_less_chasm <- dnds_less_chasm[order(-is_known_driver_mut, -chasmScore, 
                                           -chasmPadj, -confidenceLvl)]
  dnds_less_chasm <- split(dnds_less_chasm, drop = T, 
                           by = c('gr_id', 'gene_id', 'gene_name'))
  for(dlc_idx in 1:length(dnds_less_chasm)) {
    n_driver_muts_dls <- unique(dnds_less_chasm[[dlc_idx]]$dnds_n_driver_mle)
    dnds_less_chasm[[dlc_idx]][1:n_driver_muts_dls]$confidenceLvl <- 'driver'
  }
  dnds_less_chasm <- do.call(rbind, dnds_less_chasm)
  dnds_less_chasm[confidenceLvl != 'driver']$confidenceLvl <- 'passenger'
  dnds_less_chasm[confidenceLvl == 'driver']$confidenceLvl <- 'dNdScv&CHASM+'
  
  dnds_less_chasm <- dnds_less_chasm[, c(colnames(mleDriveMuts),
                                         'confidenceLvl'), with = F]
  driverMutations <- rbind(driverMutations, dnds_less_chasm, fill = T)
  # info message
  driverMutations[, participant_id := as.character(participant_id)]
  display_driverMut_update_msg(patientsInv, driverMutations,
                               confLvl = 'dNdScv&CHASM+')
}

# If for some reason number of sum of driver mutations identified by CHASM+ and
# known driver mutations is 0, we can still pinpoint driver mutation if their
# percentage is high enough
missen_muts[, perc_diff := 1 - round(dnds_n_driver_mle)/missen_total]
missen_muts[round(dnds_n_driver_mle) > 0 & chasm_n_driver == 0 & 
              known_n_driver == 0 & perc_diff <= args$max_fp & 
              status == '']$status <- 'high_perc'
missen_muts[, perc_diff := NULL]
if (nrow(missen_muts[status == 'high_perc']) > 0) {
  dnds_only <- missen_muts[status == 'high_perc']
  dnds_only[, confidenceLvl := 'dNdScv']
  dnds_only <- dnds_only[, c(colnames(mleDriveMuts), 'confidenceLvl'), 
                         with = F]
  driverMutations <- rbind(driverMutations, dnds_less_chasm, fill = T)
  
  # info message
  driverMutations[, participant_id := as.character(participant_id)]
  display_driverMut_update_msg(patientsInv, driverMutations,
                               confLvl = 'dNdScv')
}

# If gene was not identified as driver by dNdScv, it is likely that the 
# expected number of driver mutations for it will be 0. However, gene can have
# known driver mutations.
missen_muts[dnds_n_driver_mle == 0 & known_n_driver != 0 &
              status == '']$status <- 'known mut., gene is not driver'
if (nrow(missen_muts[status == 'known mut., gene is not driver']) > 0) {
  known_only <- missen_muts[status == 'known mut., gene is not driver']
  known_only[, confidenceLvl := 'known mut., gene is not driver']
  known_only <- known_only[, c(colnames(mleDriveMuts), 'confidenceLvl'), 
                           with = F]
  driverMutations <- rbind(driverMutations, known_only, fill = T)
  # info message
  driverMutations[, participant_id := as.character(participant_id)]
  display_driverMut_update_msg(patientsInv, driverMutations,
                               confLvl = 'known mut., gene is not driver')
}


# If gene was not identified as driver by dNdScv, it is likely that the
# expected number of driver mutations for it will be 0. However, CHASM+ can 
# identify potential driver mutations. 
missen_muts[dnds_n_driver_mle == 0 & chasm_n_driver != 0  &
              status == '']$status <- 'CHASM+'
if (nrow(missen_muts[status == 'CHASM+']) > 0) {
  chasm_only <- missen_muts[status == 'CHASM+']
  chasm_only[, confidenceLvl := 'dNdScv']
  chasm_only <- chasm_only[, c(colnames(mleDriveMuts), 'confidenceLvl'), 
                           with = F]
  driverMutations <- rbind(driverMutations, chasm_only, fill = T)
  # info message
  driverMutations[, participant_id := as.character(participant_id)]
  display_driverMut_update_msg(patientsInv, driverMutations,
                               confLvl = 'CHASM+')
}

missen_muts[status == '']

