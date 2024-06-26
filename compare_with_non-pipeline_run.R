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
oldRunPath <- '../Noncoding_drivers/results/2023-03-21_synRemoved/tables/tiered_drivers_CGConly_USE_ME.csv'
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
compPlotList <- list()
for (varToCheck in names(colsToCheck[colsToCheck == 'numeric'])) {
  corrDT <- merged[, c(paste0(varToCheck, '.current'),
                       paste0(varToCheck, '.old')), with = F]
  corrDT <- corrDT[complete.cases(corrDT)]
  corrVals <- c(corrVals, cor(corrDT[, 1], corrDT[, 2]))
  names(corrVals)[length(corrVals)] <- varToCheck
  
  plotData <- merged[, c('gr_id', paste0(varToCheck, '.current'), 
                        paste0(varToCheck, '.old')), with = F]
  plotData <- plotData[!is.na(unlist(plotData[, 2])) | 
                         !is.na(unlist(plotData[, 3]))]
  compPlot <- ggplot(plotData,
                     aes(x = get(paste0(varToCheck, '.current')),
                         y = get(paste0(varToCheck, '.old')), 
                         color = gr_id)) + 
    geom_point() + coord_equal() + xlab(paste0(varToCheck, ', current')) +
    ylab(paste0(varToCheck, ', old')) + geom_abline(slope = 1) + 
    theme_classic() + xlab('pipeline') + ylab('previously') +
    theme(legend.position = 'bottom', legend.direction = 'horizontal') +
    guides(color = guide_legend(title = "Genomic region")) +
    ggtitle(paste('Version comparison for', varToCheck,
                  '. R = ', round(corrVals[length(corrVals)], 4)))
  compPlotList[[length(compPlotList) + 1]] <- compPlot
  rm(plotData)
  gc()
}
message('[', Sys.time(), '] Look at warnings! They are helpful to see if ',
        'there is an influx of NA')
warnings()
rm(varToCheck)

# Find where raw p-values are computed in one version, but not in the other----
noRawPold <- merged[(is.na(digdriver.raw_p.old) & !is.na(digdriver.raw_p.current)) |
                      (is.na(dndscv.raw_p.old) & !is.na(dndscv.raw_p.current)) | 
                      (is.na(mutpanning.raw_p.old) & !is.na(mutpanning.raw_p.current)) |
                      (is.na(nbr.raw_p.old) & !is.na(nbr.raw_p.current)) |
                      (is.na(oncodrivefml.raw_p.old) & !is.na(oncodrivefml.raw_p.current))]
message('[', Sys.time(), '] Number of entries, where previously raw p-value ',
        'was not computed: ', dim(noRawPold))
noRawPcurrent <- merged[(!is.na(digdriver.raw_p.old) & is.na(digdriver.raw_p.current)) |
                          (!is.na(dndscv.raw_p.old) & is.na(dndscv.raw_p.current)) | 
                          (!is.na(mutpanning.raw_p.old) & is.na(mutpanning.raw_p.current)) |
                          (!is.na(nbr.raw_p.old) & is.na(nbr.raw_p.current)) |
                          (!is.na(oncodrivefml.raw_p.old) & is.na(oncodrivefml.raw_p.current))]
message('[', Sys.time(), '] Number of entries, where currently raw p-value ',
        'is not computed: ', dim(noRawPcurrent))

# Find where previously brown.comb_p was not computed, but now it is ----------
noBrownInOld <- merged[is.na(brown.comb_p.old) & !is.na(brown.comb_p.current)]
noBrownInOld <- merge(noBrownInOld[, keyCols, with = F], old, by = keyCols)

noBrownInNew <- merged[!is.na(brown.comb_p.old) & is.na(brown.comb_p.current)]
noBrownInNew <- merge(noBrownInNew[, keyCols, with = F], old, by = keyCols)
