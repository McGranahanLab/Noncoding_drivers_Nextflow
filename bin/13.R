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

# Functions : identification of driver mutations ------------------------------
#' findDriverMutations
#' @description Distinguishes between driver and passenger SNV and small indel
#' mutations 
#' @author Maria Litovchenko
#' @param mutDT data table with mutations of one type (i.e. missense, indel,
#' nonsense or noncoding.). Essential columns: participant_id, key, gr_id,
#' gene_id, gene_name, prob_is_driver_low, prob_is_driver_mle, 
#' prob_is_driver_high, chasmScore, chasmPadj, is_known_driver_mut
#' @param patientsInvDT data table with information about patients. Essential
#' columns: participant_id.
#' @param varType string, type of variants under consideration
#' @param chasm_score_min minimal score for CHASM+ result to be significant
#' @param chasm_padj maximum adjusted p-value of CHASM+ result  to be 
#' significant
#' @return mutDT with added column confidenceLvl
findDriverMutations <- function(mutDT, patientsInvDT, varType, 
                                chasm_score_min, chasm_padj) {
  mutInfoDT <- copy(mutDT)
  mutInfoDT[, status := '']
  result <- data.table()
  
  if (nrow(mutDT) == 0) {
    return(result)
  }
  
  # compute total number of driver mutations per driver gene & preserve mle
  # probabilities of being found as driver.
  mutInfoDT[, n_total := length(unique(key)), by =.(gr_id, gene_id, gene_name)]
  
  # calculate estimates of number of DRIVER mutations ACROSS ALL PATIENTS for
  # each gene (from dNdScv or NBR)
  mutInfoDT[, n_driverMut_low  := prob_is_driver_low  * n_total]
  mutInfoDT[, n_driverMut_mle  := prob_is_driver_mle  * n_total]
  mutInfoDT[, n_driverMut_high := prob_is_driver_high * n_total]
  
  # calculate number of driver mutations identified by CHASM+. CHASM+ scores
  # only missence mutations. Sometimes nonsense mutations will also be scored
  # due to misalignment of annotation between CHASM+ and input mutations. 
  # However, this number will usually be 0 for nonsense, indel or noncoding 
  # mutations
  mutInfoDT[, inCHASM := F]
  mutInfoDT[chasmScore >= chasm_score_min & 
              chasmPadj <= chasm_padj]$inCHASM <- T
  mutInfoDT[, chasm_n_driver := length(unique(key[inCHASM])),
            by = .(gr_id, gene_id, gene_name)]
  mutInfoDT[, inCHASM := NULL]
  
  # calculate number of known driver mutations, NOT identified by CHASM+ in the
  # genes
  mutInfoDT[, known_not_chasm := F]
  mutInfoDT[is_known_driver_mut %in% T & 
              (chasmScore < chasm_score_min |
                 chasmPadj > chasm_padj)]$known_not_chasm <- T
  mutInfoDT[, known_n_driver := length(unique(key[known_not_chasm])),
            by = .(gr_id, gene_id, gene_name)] 
  mutInfoDT[, known_not_chasm := NULL]
  
  # Now, ACROSS ALL PATIENTS we have 2 numbers per gene: number of drivers from 
  # dNdScv and sum of number of drivers from CHASM+ and known driver mutations
  # (not identified as significant by CHASM+). We want to produce a
  # conservative estimation. Therefore, for genes for which number of driver
  # mutations from dNdScv is bigger than sum of number of driver mutations from
  # CHASM+ and known driver mutations, we will take as drivers mutations 
  # mutations identified by CHASM+ & known driver mutations. This will only be
  # executed for missense mutations, because CHASM+ will not score other types
  # of mutations.
  mutInfoDT[, status := '']
  mutInfoDT[n_driverMut_mle > 0 & chasm_n_driver + known_n_driver != 0 &
              chasm_n_driver + known_n_driver <= n_driverMut_mle &
              status == '']$status <- 'chasmLess'
  if (nrow(mutInfoDT[status == 'chasmLess']) > 0) {
    chasmLess <- mutInfoDT[status == 'chasmLess']
    chasmLess[, confidenceLvl := 'passenger']
    cLvl <- paste0(varType, ',', 'dNdScv&CHASM+&known')
    chasmLess[(chasmScore >= args$chasm_score_min &
                 chasmPadj <= args$chasm_padj) |
                is_known_driver_mut %in% T]$confidenceLvl <- cLvl
    cLvl <- paste0(varType, ',', 'dNdScv&CHASM+')
    chasmLess[(chasmScore >= args$chasm_score_min &
                 chasmPadj <= args$chasm_padj) |
                is_known_driver_mut %in% F]$confidenceLvl <- cLvl
    
    result <- rbind(result, chasmLess, fill = T)
    # info message
    display_driverMut_update_msg(patientsInvDT, result, 
                                 confLvl = paste0(varType, ',', 
                                                  'dNdScv&CHASM+&known'))
    display_driverMut_update_msg(patientsInvDT, result,
                                 confLvl = paste0(varType, ',dNdScv&CHASM+'))
  }
  
  # If it's the other way around (CHASM & known driver mutations > dNdScv), then
  # we can order known driver mutations and CHASM+ results by significance and 
  # select out of them top N, there N is the number from dNdScv.
  mutInfoDT[n_driverMut_mle > 0 & 
              chasm_n_driver + known_n_driver > n_driverMut_mle &
              status == '']$status <- 'chasmMore'
  if (nrow(mutInfoDT[status == 'chasmMore']) > 0) {
    chasmMore <- mutInfoDT[status == 'chasmMore']
    chasmMore[, n_driverMut_mle := round(n_driverMut_mle)]
    
    # sort mutation by probability of being drivers
    chasmMore[, confidenceLvl := 'passenger']
    chasmMore[is_known_driver_mut %in% T]$confidenceLvl <- 'known'
    chasmMore[chasmScore >= args$chasm_score_min &
                chasmPadj <= args$chasm_padj]$confidenceLvl <- 'CHASM+'
    chasmMore[, confidenceLvl := factor(confidenceLvl, 
                                        levels = c('passenger', 'CHASM+', 
                                                   'known'))]
    chasmMore <- chasmMore[order(-is_known_driver_mut, -chasmScore, 
                                 -chasmPadj, -confidenceLvl)]
    chasmMore <- split(chasmMore, drop = T, 
                       by = c('gr_id', 'gene_id', 'gene_name'))
    for(dlc_idx in 1:length(chasmMore)) {
      n_driver_muts_dls <- unique(chasmMore[[dlc_idx]]$n_driverMut_mle)
      chasmMore[[dlc_idx]][1:n_driver_muts_dls]$confidenceLvl <- 'driver'
    }
    chasmMore <- do.call(rbind, chasmMore)
    chasmMore[confidenceLvl != 'driver']$confidenceLvl <- 'passenger'
    cLvl <- paste0(varType, ',', 'dNdScv&CHASM+')
    chasmMore[confidenceLvl == 'driver']$confidenceLvl <- cLvl
    
    result <- rbind(result, chasmMore, fill = T)
    # info message
    display_driverMut_update_msg(patientsInvDT, result,
                                 confLvl = 'dNdScv&CHASM+')
  }
  
  # If for some reason number of sum of driver mutations identified by CHASM+ 
  # andknown driver mutations is 0, we can still pinpoint driver mutation if 
  # their percentage is high enough. This section will be executed for nonsense
  # indels and noncoding drivers
  mutInfoDT[, perc_diff := 1 - round(n_driverMut_mle)/n_total]
  mutInfoDT[round(n_driverMut_mle) > 0 & chasm_n_driver == 0 & 
              known_n_driver == 0 & perc_diff <= args$max_fp & 
              status == '']$status <- 'high_perc'
  mutInfoDT[, perc_diff := NULL]
  if (nrow(mutInfoDT[status == 'high_perc']) > 0) {
    # although the name of object here is "dnds_only", in case of noncoding
    # if should be "nbr_only"
    dnds_only <- mutInfoDT[status == 'high_perc']
    cLvl <- ifelse(varType == 'noncoding', 'NBR', 'dNdScv')
    dnds_only[, confidenceLvl := cLvl]
    
    result <- rbind(result, dnds_only, fill = T)
    # info message
    display_driverMut_update_msg(patientsInvDT, result, confLvl = cLvl)
  }
  
  # If gene was not identified as driver by dNdScv/NBR (but it was identified
  # as such by the other software, otherwise it wouldn't be here), it is likely
  # that theexpected number of driver mutations for it will be 0. However, gene
  #  can have known driver mutations.
  mutInfoDT[n_driverMut_mle == 0 & known_n_driver != 0 &
              status == '']$status <- 'known mut.'
  if (nrow(mutInfoDT[status == 'known mut.']) > 0) {
    known_only <- mutInfoDT[status == 'known mut.']
    known_only[, confidenceLvl := 'known mut.']
    
    result <- rbind(result, known_only, fill = T)
    # info message
    display_driverMut_update_msg(patientsInvDT, result, confLvl = 'known mut.')
  }
  
  # If gene was not identified as driver by dNdScv, it is likely that the
  # expected number of driver mutations for it will be 0. However, CHASM+ can 
  # identify potential driver mutations. 
  mutInfoDT[n_driverMut_mle == 0 & chasm_n_driver != 0  &
              status == '']$status <- 'CHASM+'
  if (nrow(mutInfoDT[status == 'CHASM+']) > 0) {
    chasm_only <- mutInfoDT[status == 'CHASM+']
    chasm_only[, confidenceLvl := 'CHASM+']
    
    result <- rbind(result, chasm_only, fill = T)
    # info message
    display_driverMut_update_msg(patientsInvDT, result, confLvl = 'CHASM+')
  }
  
  # If at this point status is still empty, it means that we can only 
  # reliably estimate number of missence drivers in the gene
  if (nrow(mutInfoDT[status == '']) != 0) {
    n_only <- mutInfoDT[status == '']
    n_only[, confidenceLvl := 'only n. driver mutations'] 
    
    result <- rbind(result, n_only, fill = T)
    display_driverMut_update_msg(patientsInvDT, result,
                                 confLvl = 'only n. driver mutations')
  }
  
  # simple check that all mutations were accounted for
  result <- result[, status := NULL]
  result <- unique(result)
  result[, participant_id := as.character(participant_id)]
  resultCheck <- mutDT[,.(key, gr_id, gene_id, var_type, gene_name,
                          participant_id)]
  resultCheck[, participant_id := as.character(participant_id)]
  resultCheck <- merge(resultCheck, result, all.x = T,
                       by = c('key', 'gr_id', 'gene_id', 'var_type', 
                              'gene_name', 'participant_id'))
  if (any(is.na(resultCheck$confidenceLvl))) {
    stop('[', Sys.time(), '] Error in function : not all mutations are accounted for.')
  }
  
  result
}

# Functions : information messages --------------------------------------------
#' display_driverMut_update_msg
#' @description Prints info message about how many of driver mutations were 
#'              assigned to a certain confidence level
#' @author Maria Litovchenko
#' @param patientsInvDT data table, inventory of patients with participant_id
#' column
#' @param driveMutsDT data table with driver mutations with column 
#' confidenceLvl
#' @param confLvl string, confidence level
#' @return void
display_driverMut_update_msg <- function(patientsInvDT, driveMutsDT, confLvl) {
  patientsInvDT[, participant_id := as.character(participant_id)]
  driveMutsDT[, participant_id := as.character(participant_id)]
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
             muts_to_gr = 'completed_runs/2023-12-14/results/mut_rates/mutMapToGR-Panlung--hg19.csv',
             inferred_biotype = 'completed_runs/2023-12-14/results/',
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

# Identify DRIVER mutation based on dNdScv, NBR, CHASM & known muts -----------
# Select missence mutations. We do selection based on var_type = 'mis' and not
# by var_class, because var_class can be Missense_Mutation, Nonstop_Mutation,
# Translation_Start_Site or Unknown for var_type = 'mis'. Yet, the same dNdScv
# probability as for Missense_Mutation will be assigned to for example 
# Nonstop_Mutation. It could even be scored by CHASM+, which normally works on
#  missence mutations only!
missen_muts <- mleDriveMuts[gr_id %in% coding_gr_id & var_type == 'mis']
missenDrivers <- findDriverMutations(missen_muts, patientsInv, 
                                     varType = 'missense',
                                     chasm_score_min = args$chasm_score_min,
                                     chasm_padj = args$chasm_padj)
nonsense_muts <- mleDriveMuts[gr_id %in% coding_gr_id & var_type == 'non']
nonsenDrivers <- findDriverMutations(nonsense_muts, patientsInv, 
                                     varType = 'nonsense',
                                     chasm_score_min = args$chasm_score_min,
                                     chasm_padj = args$chasm_padj)
indels_muts <- mleDriveMuts[gr_id %in% coding_gr_id & var_type == 'indels']
indeDrivers <- findDriverMutations(indels_muts, patientsInv, 
                                   varType = 'indels',
                                   chasm_score_min = args$chasm_score_min,
                                   chasm_padj = args$chasm_padj)
noncoding_muts <- mleDriveMuts[!gr_id %in% coding_gr_id]
noncodDrivers <- findDriverMutations(noncoding_muts, patientsInv, 
                                     varType = 'noncoding',
                                     chasm_score_min = args$chasm_score_min,
                                     chasm_padj = args$chasm_padj)
