args <- list(tieredPs = '../TEST/results/tables/rawDrivers_LUAD--_hg19.csv',
             codingGRid = list('CDS'), min_n_soft_cod = 3,
             min_n_soft_noncod = 2, min_n_muts = 3, min_n_patients = 3,
             max_local_mut_rate_q = 0.99, use_synonumous_mut_rate = T,
             max_gr_mut_rate_q = 0.99, max_gr_len_q = 0.99)

# Read in tiered p-values -----------------------------------------------------
tieredPs <- fread(args$tieredPs, header = T, stringsAsFactors = F)

# Initialize FILTER -----------------------------------------------------------
tieredPs[, FILTER := 'PASS']

# Mark genomic elements not scanned with enough of software -----------------
# get columns with raw P values
rawPcols <- grep('[.]raw_p$', colnames(tieredPs), value = T)
tieredPs[, nSoft := apply(tieredPs[, rawPcols, with = F], 1, 
                          function(x) sum(!is.na(x)))]
tieredPs[, minNsoft := ifelse(gr_id %in% args$codingGRid, args$min_n_soft_cod, 
                              args$min_n_soft_noncod)]
pass <- tieredPs$nSoft >= tieredPs$minNsoft
tieredPs[!pass]$FILTER <- paste0(tieredPs[!pass]$FILTER,
                                 '; low N. software')
tieredPs[, nSoft := NULL]
tieredPs[, minNsoft := NULL]

# Mark genomic elements with not enough of N mutations or patients ----------
pass <- tieredPs$nMuts_total >= args$min_n_muts & !is.na(tieredPs$nMuts_total)
tieredPs[!pass]$FILTER <- paste0(tieredPs[!pass]$FILTER, '; low N. mutations')

pass <- tieredPs$nParts_total >= args$min_n_patients & 
  !is.na(tieredPs$nParts_total)
tieredPs[!pass]$FILTER <- paste0(tieredPs[!pass]$FILTER,
                                 '; low N. participants') 

# Mark genes with high local mutational rate  -------------------------------
if ('meanMutRateQuant_local' %in% colnames(tieredPs)) {
  if (args$max_local_mut_rate_q < 1) {
    pass <- tieredPs$meanMutRateQuant < 100*args$max_local_mut_rate_q &
      !is.na(tieredPs$meanMutRateQuant_local)
    tieredPs[!pass]$FILTER <- paste0(tieredPs[!pass]$FILTER,
                                     '; high local mut.rate')
  }
} else {
  if (!'meanMutRateQuant_local' %in% colnames(tieredPs)) {
    message('[', Sys.time(), '] Filtering based on mean mutation rate ',
            'quantile was not performed, because column ',
            'meanMutRateQuant_local was not found.')
  }
}

# Mark by genomic regions specific mutational rate/synonymous mut.rate ------
colToUse <- 'meanMutRateQuant'
filterName <-  '; high genomic region mut.rate'
if (args$use_synonumous_mut_rate) {
  colToUse <- 'synMeanMutRateQuant'
  filterName <-  '; high syn gr mut.rate'
}

if (colToUse %in% colnames(tieredPs)) {
  if (args$max_gr_mut_rate_q < 1) {
    grMutRate <- unlist(tieredPs[, colToUse, with = F])
    pass <- grMutRate < 100*args$max_gr_mut_rate_q & !is.na(grMutRate)
    tieredPs[!pass]$FILTER <- paste0(tieredPs[!pass]$FILTER, filterName)
  }
} else {
  if (!colToUse %in% colnames(tieredPs)) {
    msg <- 'Filtering based on mean genome region '
    if (args$use_synonumous_mut_rate) {
      msg <- paste(msg, 'SYNONYMOUS')
    }
    msg <- paste(msg, 'mutation rate quantile was not performed, because ',
                 'column ', colToUse, ' was not found.')
    message('[', Sys.time(), '] ', msg)
  }
}

# Mark olfactory genes --------------------------------------------------------
if ('is_olfactory' %in% colnames(tieredPs)) {
  pass <- !tieredPs$is_olfactory
  tieredPs[!pass]$FILTER <- paste0(tieredPs[!pass]$FILTER, '; olfactory') 
}

# Mark by quantile of length --------------------------------------------------
if (args$max_gr_len_q < 1) {
  pass <- tieredPs$gr_lenQuant < 100*args$max_gr_len_q
  tieredPs[!pass]$FILTER <- paste0(tieredPs[!pass]$FILTER, '; long gene') 
}

# Mark out gene based on being not expressed ----------------------------------
exprCols <- grep('expr_in', colnames(tieredPs), value = T)
if (length(exprCols) > 0) {
  for (eCol in exprCols) {
    eVals <- unlist(tieredPs[, eCol, with = F])
    pass <- eVals | is.na(eVals)
    filterName <- paste0('; not expr. in ', gsub('expr_in_', '', eCol))
    tieredPs[!pass]$FILTER <- paste0(tieredPs[!pass]$FILTER, filterName) 
  }
} 

# Remove PASS at the start of FILTER for entries which didn't pass ------------
tieredPs[, FILTER := gsub('^PASS; ', '', FILTER)]

# [SAVE] Tiered results as a table --------------------------------------------
colOrderToPrint <- c('tumor_subtype', 'gr_id', 'gene_id', 'gene_name', 
                     'is_known_cancer', 'FILTER', 'participant_tumor_subtype', 
                     'nMuts', 'nMutsUniq', 'nParts', 'nMuts_total', 
                     'nMutsUniq_total', 'nParts_total', 
                     'meanMutRate', 'meanMutRateQuant',
                     sort(grep('.raw_p$', colnames(tieredPs), 
                               value = T)),
                     sort(grep('.comb_p$', colnames(tieredPs), 
                               value = T)),
                     sort(grep('.bh_p$', colnames(tieredPs), 
                               value = T)),
                     sort(grep('.tier$', colnames(tieredPs), value = T)))
setcolorder(tieredPs, colOrderToPrint)
write.table(tieredPs, args$output, append = F, quote = F,  sep = '\t', 
            row.names = F, col.names = T)

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')