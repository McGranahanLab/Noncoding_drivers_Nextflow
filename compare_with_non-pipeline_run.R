#!/usr/bin/env Rscript
# FILE: compare_with_non-pipeline_run.R ---------------------------------------
#
# DESCRIPTION:
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
# CREATED:  07.12.2020
# REVISION: 07.12.2023

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

# Inputs ----------------------------------------------------------------------
currentRunPath <- 'completed_runs/06_12_2023/results/tables/drivers/drivers-Panlung--hg19.csv'
tumor_subtype_current <- 'Panlung'
oldRunPath <- '../Noncoding_drivers/results/2023-03-21_synRemoved/tables/tiered_drivers_CGConly.csv'
tumor_subtype_old <- 'PANCAN'

keyCols <- c('gr_id', 'gene_id', 'gene_name')
colsToCheck <- c('FILTER' = 'character', 
                 
                 'nMuts_total' = 'numeric', 'nParts_total' = 'numeric',
                 'digdriver.raw_p' = 'numeric', 'dndscv.raw_p' = 'numeric', 
                 'mutpanning.raw_p' = 'numeric', 'nbr.raw_p' = 'numeric', 
                 'oncodrivefml.raw_p' = 'numeric',
                 
                 'digdriver.bh_p' = 'numeric', 'dndscv.bh_p'= 'numeric', 
                 'mutpanning.bh_p' = 'numeric', 'nbr.bh_p' = 'numeric', 
                 'oncodrivefml.bh_p' = 'numeric',
                 
                 'brown.comb_p' = 'numeric', 'brown.bh_p' = 'numeric')

# Read tables -----------------------------------------------------------------
current <- fread(currentRunPath, header = T, stringsAsFactors = F)
current <- current[tumor_subtype == tumor_subtype_current]
old <- fread(oldRunPath, header = T, stringsAsFactors = F)
old <- old[tumor_subtype == tumor_subtype_old]
# fix some region names
fixNeeded <- grepl('__and__', old$gene_name)
fixedNames <- unique(old$gene_name[fixNeeded])
names(fixedNames) <- fixedNames
fixedNames <- sapply(fixedNames, 
                     function(x) paste0(sort(strsplit(x, '__and__')[[1]]),
                                        collapse = '__and__'))
old[fixNeeded]$gene_name <- fixedNames[old[fixNeeded]$gene_name]
# fix some region ids
fixNeeded <- grepl('__and__', old$gene_id)
fixedIDs <- unique(old$gene_id[fixNeeded])
names(fixedIDs) <- fixedIDs
fixedIDs <- sapply(fixedIDs, 
                   function(x) paste0(sort(strsplit(x, '__and__')[[1]]),
                                      collapse = '__and__'))
old[fixNeeded]$gene_id <- fixedIDs[old[fixNeeded]$gene_id]
  
# Find drivers not coherent drivers -------------------------------------------
current_drivers <- current[FILTER == 'PASS' & !is.na(tier)]
current_drivers <- unique(current_drivers[, keyCols, with = F])
old_drivers <- old[FILTER == 'PASS' & !is.na(brown.tier)]
old_drivers <- unique(old_drivers[, keyCols, with = F])

notCoherent <- merge(cbind(old_drivers, status = 'old'),
                     cbind(current_drivers, status = 'current'), by = keyCols,
                     all = T)
notCoherent[, status := 'ok']
notCoherent[is.na(status.x)]$status <- 'current_only'
notCoherent[is.na(status.y)]$status <- 'old_only'
notCoherent <- notCoherent[, c(keyCols, 'status'), with = F]
notCoherent <- notCoherent[order(gr_id, gene_name)]

if (all(notCoherent$status == 'ok')) {
  stop('All drivers match')
} else {
  message('[', Sys.time(), '] Incoherent drivers: ')
  print(table(notCoherent[status != 'ok']$gr_id, 
              notCoherent[status != 'ok']$status))
}

# Merge results of old and current runs ---------------------------------------
merged <- merge(unique(old[, c(keyCols, names(colsToCheck)), with = F]),
                unique(current[, c(keyCols, names(colsToCheck)), with = F]), 
                by = keyCols, all = T)
colnames(merged) <- gsub('[.]x$', '.old', colnames(merged))
colnames(merged) <- gsub('[.]y$', '.current', colnames(merged))
merged <- merge(merged, notCoherent, by = keyCols, all = T)
merged[is.na(status)]$status <- 'not driver'

# Plot correlations between p-values ------------------------------------------
corrVals <- c()
for (varToCheck in names(colsToCheck[colsToCheck == 'numeric'])) {
  corrDT <- merged[, c(paste0(varToCheck, '.current'),
                       paste0(varToCheck, '.old')), with = F]
  corrDT <- corrDT[complete.cases(corrDT)]
  corrVals <- c(corrVals, cor(corrDT[, 1], corrDT[, 2]))
  names(corrVals)[length(corrVals)] <- varToCheck
  
  compPlot <- ggplot(merged,
                     aes(x = get(paste0(varToCheck, '.current')),
                         y = get(paste0(varToCheck, '.old')), 
                         color = status)) + 
    geom_point() + coord_equal() + xlab(paste0(varToCheck, ', current')) +
    ylab(paste0(varToCheck, ', old')) + geom_abline(slope = 1)
  print(compPlot)
}
message('[', Sys.time(), '] Look at warnings! They are helpful to see if ',
        'there is an influx of NA')
warnings()
rm(varToCheck)

# Plot correlations between brown p-values for each genomic region ------------
varToCheck <- 'brown.comb_p'
corrVals_grID <- c()
for (grID in sort(unique(merged$gr_id))) {
  corrDT <- merged[gr_id == grID][, c(paste0(varToCheck, '.current'),
                                      paste0(varToCheck, '.old')), with = F]
  corrDT <- corrDT[complete.cases(corrDT)]
  corrVals_grID <- c(corrVals_grID, cor(corrDT[, 1], corrDT[, 2]))
  names(corrVals_grID)[length(corrVals_grID)] <- grID
  
  compPlot <- ggplot(merged[gr_id == grID],
                     aes(x = get(paste0(varToCheck, '.current')),
                         y = get(paste0(varToCheck, '.old')), 
                         color = status)) + 
    geom_point() + coord_equal() + xlab(paste0(varToCheck, ', current')) +
    ylab(paste0(varToCheck, ', old')) + geom_abline(slope = 1) +
    ggtitle(grID)
  print(compPlot)
}
message('[', Sys.time(), '] Look at warnings! They are helpful to see if ',
        'there is an influx of NA')
warnings()

# Find where previously brown.comb_p was not computed, but now it is ----------
noBrownInOld <- merged[is.na(brown.comb_p.old) & !is.na(brown.comb_p.current)]
noBrownInOld <- merge(noBrownInOld[, keyCols, with = F], old, by = keyCols)

