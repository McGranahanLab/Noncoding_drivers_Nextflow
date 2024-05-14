library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)

source('bin/custom_functions.R')

# Functions -------------------------------------------------------------------
readAndSumUpDriverMutations <- function(filePath, tumor_subtype, 
                                       allDetectedDrivers) {
  colsToKeep <- c("is_driver", "gr_id", "gene_id", "gene_name", "var_type",
                  "key", "participant_id", "n_total", "n_driverMut_low", 
                  "n_driverMut_mle", "n_driverMut_high",
                  "n_tums_w_driver_mut")
  driverMuts <- fread(filePath, header = T, select = colsToKeep,
                      stringsAsFactors = F)
  
  driverMuts <- driverMuts[var_type %in% NON_CNV_VAR_TYPES]
  driverMuts[, n_participants := length(unique(participant_id)),
             by = .(gr_id, gene_id, gene_name)]
  driverMuts[, n_muts := length(unique(key)), 
             by = .(gr_id, gene_id, gene_name)]
  driverMuts[, participant_id := NULL]
  driverMuts[, key := NULL]
  driverMuts <- unique(driverMuts)
  
  driverMuts <- merge(driverMuts, allDetectedDrivers, all.y = T,
                      by = c("gr_id", "gene_id", "gene_name"))
  
  # sum across var_type
  driverMuts[, n_total := sum(n_total), by = .(gr_id, gene_id, gene_name)]
  driverMuts[, perc_driverMut_low := 100 * sum(n_driverMut_low, na.rm = T)/n_total,
             by = .(gr_id, gene_id, gene_name)]
  driverMuts[, perc_driverMut_mle := 100 * sum(n_driverMut_mle, na.rm = T)/n_total,
             by = .(gr_id, gene_id, gene_name)]
  driverMuts[, perc_driverMut_high := 100 * sum(n_driverMut_high, na.rm = T)/n_total,
             by = .(gr_id, gene_id, gene_name)]
  driverMuts[, n_tums_w_driver_mut := sum(n_tums_w_driver_mut),
             by = .(gr_id, gene_id, gene_name)]
  
  driverMuts[, tumor_subtype := tumor_subtype]
  driverMuts <- driverMuts[,.(tumor_subtype, is_driver, gr_id, gene_id,
                              gene_name, is_known_cancer, perc_driverMut_low, 
                              perc_driverMut_mle, perc_driverMut_high, 
                              n_participants, n_muts, n_tums_w_driver_mut)]
  driverMuts <- unique(driverMuts)
  driverMuts[is.na(driverMuts)] <- 0
  driverMuts[, is_driver := as.logical(is_driver)]
  driverMuts
}

# Test input arguments --------------------------------------------------------
args <- list(inventory_patients = 'data/inventory/inventory_patients.csv',
             cancer_subtype_1 = 'Adenocarcinoma',
             drivers_1 = "completed_runs/2023_12_25/results/tables/drivers/drivers-Adenocarcinoma--hg19.csv",
             driver_mutations_1 = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Adenocarcinoma--hg19.csv",
             cancer_subtype_2 = 'Squamous_cell',
             driver_mutations_2 = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Squamous_cell--hg19.csv",
             drivers_2 = "completed_runs/2023_12_25/results/tables/drivers/drivers-Squamous_cell--hg19.csv",
             subtype_specificity_file = 'completed_runs/2023_12_25/results/tables/subtype_specificity/subtypeSpecificity---hg19.csv',
             min_n_patients = 3)

NON_CNV_VAR_TYPES = c("subs", "mis", "non", "indels")
MODE_UPD <- "found as driver in a tumor subtype,\nbut is preferential/specific to the other subtype"

# Read in patients inventory --------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients, 1)
patientsInv <- patientsInv[tumor_subtype %in% c(args$cancer_subtype_1, 
                                                args$cancer_subtype_2)]
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)

# Calculate cohort sizes ------------------------------------------------------
cohort_sizes <- patientsInv[,.(length(unique(participant_id))), 
                            by = tumor_subtype]
setnames(cohort_sizes, 'V1', 'cohort_size')

# Read results of driver discovery ----------------------------------------------
drivers_1 <- fread(args$drivers_1, header = T, stringsAsFactors = F, 
                   select = c("gr_id", "gene_id", "gene_name", "FILTER", 
                              "tier", "is_known_cancer"))
drivers_1 <- drivers_1[FILTER == "PASS" & !is.na(tier)]
message('[', Sys.time(), '] Read --drivers_1: ',  args$drivers_1)

drivers_2 <- fread(args$drivers_2, header = T, stringsAsFactors = F, 
                   select = c("gr_id", "gene_id", "gene_name", "FILTER", 
                              "tier", "is_known_cancer"))
drivers_2 <- drivers_2[FILTER == "PASS" & !is.na(tier)]
message('[', Sys.time(), '] Read --drivers_2: ',  args$drivers_2)

all_drivers <- rbind(drivers_1, drivers_2)
all_drivers[, FILTER := NULL]
all_drivers[, tier := NULL]
all_drivers <- unique(all_drivers)

# Read results of subtype specificity analysis --------------------------------
subtypeSpec <- fread(args$subtype_specificity_file, header = T, 
                     stringsAsFactors = F, 
                     select = c('tumor_subtype_1', 'tumor_subtype_2', 'gr_id',
                                'gene_id', 'gene_name', 'tumor_subtype_spec',
                                'specificity_mode', 'is_driver_ts_1',
                                'is_driver_ts_2'))
subtypeSpec <- subtypeSpec[tumor_subtype_1 %in% args[c('cancer_subtype_1', 
                                                       'cancer_subtype_2')] &
                             tumor_subtype_2 %in% args[c('cancer_subtype_1', 
                                                         'cancer_subtype_2')]]
subtypeSpec <- subtypeSpec[is_driver_ts_1 == T | is_driver_ts_2 == T]

subtypeSpec <- subtypeSpec[,.(gr_id, gene_id, gene_name, tumor_subtype_spec,
                              specificity_mode)]
subtypeSpec <- unique(subtypeSpec)

# Read driver mutations -------------------------------------------------------
excessMuts <- rbind(readAndSumUpDriverMutations(args$driver_mutations_1,
                                                args$cancer_subtype_1,
                                                all_drivers),
                    readAndSumUpDriverMutations(args$driver_mutations_2,
                                                args$cancer_subtype_2,
                                                all_drivers))
subtypeSpec <- merge(excessMuts, subtypeSpec, all = T,
                     by = c('gr_id', 'gene_id', 'gene_name'))
rm(excessMuts)
subtypeSpec <- subtypeSpec[n_participants >= args$min_n_patients]

# Prepare % excess of driver mutations for plotting ---------------------------
subtypeSpec[tumor_subtype_spec != tumor_subtype & 
              tumor_subtype_spec %in% 
              c(args$cancer_subtype_1, 
                args$cancer_subtype_2)]$specificity_mode <- MODE_UPD
subtypeSpec[tumor_subtype_spec != tumor_subtype & 
              tumor_subtype_spec %in% 
              c(args$cancer_subtype_1, 
                args$cancer_subtype_2)]$tumor_subtype_spec <- ""
subtypeSpec[, specificity := paste0(tumor_subtype_spec, '-', specificity_mode)]
subtypeSpec[, specificity := gsub('^-', '', specificity)]

# sort so that tumor subtype with more specific/preferential drivers would 
# appear first
subtypeOrder <- subtypeSpec[tumor_subtype_spec %in%
                              c(args$cancer_subtype_1, args$cancer_subtype_2)]

if (nrow(subtypeOrder) > 0) {
  subtypeOrder <- subtypeOrder[,.N, by = tumor_subtype_spec]
  subtypeOrder <- subtypeOrder[order(-N)]
  subtypeOrder <- subtypeOrder$tumor_subtype_spec
} else {
  subtypeOrder <- sort(subtypeOrder, descreasing = T)
}
specificity_lvls <- c(paste0(subtypeOrder[2], "-specific"),
                      paste0(subtypeOrder[2], "-preferential"),
                      "common", "not enough evidence",
                      paste0(subtypeOrder[1], "-preferential"),
                      paste0(subtypeOrder[1], "-specific"),
                      MODE_UPD)
subtypeSpec[, specificity := factor(specificity, specificity_lvls, 
                                    ordered = T)]
empty_specificity <- which(subtypeSpec$tumor_subtype_spec == "")
subtypeSpec[empty_specificity]$tumor_subtype_spec <- subtypeSpec[empty_specificity]$specificity_mode
subtypeSpec[, tumor_subtype_spec := factor(tumor_subtype_spec, 
                                           c(subtypeOrder[2], "common", 
                                             "not enough evidence",
                                             subtypeOrder[1], MODE_UPD))]
subtypeSpec[, specificity_mode := NULL]

subtypeSpec <- subtypeSpec[order(-perc_driverMut_mle)]
subtypeSpec[, gr_id := factor(gr_id, rev(unique(gr_id)))]
subtypeSpec <- subtypeSpec[order(gr_id, tumor_subtype_spec, 
                                 perc_driverMut_mle)]
subtypeSpec[, gene_plots_id := paste(gene_name, gr_id)]
subtypeSpec[, gene_plots_id := factor(gene_plots_id, unique(gene_plots_id))]

# Prepare % of patients with driver mutations for plotting --------------------
subtypeSpec <- merge(subtypeSpec, cohort_sizes, by = "tumor_subtype")
subtypeSpec[, percTumsWdriverMut := 100 * n_tums_w_driver_mut/cohort_size]

percTumsWdriverMutBreaks <- c(-0.01, 0.5, 5, 10)
percTumsWdriverMutLabels <- c('<1%', '5%', '10%')
if (max(subtypeSpec$percTumsWdriverMut) > 0.1) {
  maxPerc <- max(subtypeSpec$percTumsWdriverMut)
  percTumsWdriverMutBreaks <- c(percTumsWdriverMutBreaks,  
                                seq(20, maxPerc, by = 10))
  percTumsWdriverMutLabels <- c(percTumsWdriverMutLabels,
                                paste0(seq(20, maxPerc, by = 10), '%'))
}

subtypeSpec[, percTumsWdriverMut_bin := cut(percTumsWdriverMut, 
                                            percTumsWdriverMutBreaks, 
                                            percTumsWdriverMutLabels)]
maxLvl <- percTumsWdriverMutLabels[length(percTumsWdriverMutLabels)]
subtypeSpec[is.na(percTumsWdriverMut_bin)]$percTumsWdriverMut_bin <- maxLvl

percTumsWdriverMutColor <- colorRampPalette(percTumorsWithDriverMutsPalette)
percTumsWdriverMutColor <- percTumsWdriverMutColor(length(percTumsWdriverMutLabels))
names(percTumsWdriverMutColor) <- percTumsWdriverMutLabels

# X axis labels and color palletes --------------------------------------------
xAxisLabs <- format_xLabels(subtypeSpec, 'gene_plots_id', 'gene_name',
                            'is_known_cancer', highlightColor = 'red')

names(subtypeSpecificityPalette) <- gsub('leftSubtype', subtypeOrder[1], 
                                         names(subtypeSpecificityPalette))
names(subtypeSpecificityPalette) <- gsub('rightSubtype', subtypeOrder[2], 
                                         names(subtypeSpecificityPalette))
names(subtypeSpecificityPalette) <- gsub('MODE_UPD', MODE_UPD, 
                                         names(subtypeSpecificityPalette))

# Left plot -------------------------------------------------------------------
LEFT_PLOTData <- subtypeSpec[tumor_subtype == subtypeOrder[1]]
LEFT_PLOT <- ggplot()

# first of all, if there are several genomic regons, create shadowing
if (length(unique(LEFT_PLOTData$gr_id)) > 1) {
  # compute coordinates for future lines showing genome region 
  grShadeCoords <- get_gr_change_coords(plotDT = LEFT_PLOTData, 
                                        xLabsDT = xAxisLabs$dt,
                                        xOrderCol = "gene_plots_id")
  grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                  nrow(grShadeCoords$gr_coords),
                                                  by = 2)]
  LEFT_PLOT <- LEFT_PLOT + 
    geom_rect(data = grShadeCoordsCut, fill = DIVIDER_COLOR, inherit.aes = F,
              aes(ymin = start - 0.5, ymax = end + 0.5, xmin = -Inf, 
                  xmax = Inf), alpha = 0.5)
}

LEFT_PLOT <- LEFT_PLOT +
  geom_bar(data = LEFT_PLOTData,
           mapping = aes(y = gene_plots_id, x = -perc_driverMut_mle,
                         fill = specificity, color = is_driver,
                         alpha = is_driver), 
           stat = 'identity', linewidth = 0.25) +
  geom_point(data = LEFT_PLOTData,
             mapping = aes(y = gene_plots_id, x = -perc_driverMut_mle),
             color = 'black', size = 1) +
  geom_errorbar(data = LEFT_PLOTData,
                mapping = aes(y = gene_plots_id, xmin = -perc_driverMut_low,
                              xmax = -perc_driverMut_high),
                color = 'black', linewidth = 0.25, width = 0.25) +
  geom_bar(data = LEFT_PLOTData,
           mapping = aes(y = gene_plots_id, x = 15, 
                         fill = percTumsWdriverMut_bin), 
           stat = 'identity', linewidth = 0.25) +
  geom_text(data = LEFT_PLOTData,
            mapping = aes(y = gene_plots_id, x = 8, 
                          label = round(percTumsWdriverMut)),
            size = 2) +
  scale_fill_manual(values = c(percTumsWdriverMutColor,
                               subtypeSpecificityPalette), drop = F) +
  scale_color_manual(values = c('FALSE' = '#CCCCCC', 'TRUE' = 'black')) +
  scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 1)) +
  scale_x_continuous(breaks = seq(0, -100, by = -25), 
                     labels = seq(0, 100, by = 25),
                     expand = c(0,0)) +
  scale_y_discrete(drop = F, position = 'right') +
  xlab('% excess of mutations') + ylab('') + 
  labs(color = 'is driver\n in tumor subtype', 
       alpha = 'is driver\n in tumor subtype') +
  customGgplot2Theme + ggtitle(subtypeOrder[1]) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        axis.text.y = element_blank(), axis.line.y = element_blank(),
        plot.title = element_text(size = 8))

# Right plot ------------------------------------------------------------------------
RIGHT_PLOTData <- subtypeSpec[tumor_subtype == subtypeOrder[2]]
RIGHT_PLOT <- ggplot()

# first of all, if there are several genomic regons, create shadowing
if (length(unique(RIGHT_PLOTData$gr_id)) > 1) {
  # compute coordinates for future lines showing genome region 
  grShadeCoords <- get_gr_change_coords(plotDT = RIGHT_PLOTData, 
                                        xLabsDT = xAxisLabs$dt,
                                        xOrderCol = "gene_plots_id")
  grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                  nrow(grShadeCoords$gr_coords),
                                                  by = 2)]
  RIGHT_PLOT <- RIGHT_PLOT + 
    geom_rect(data = grShadeCoordsCut, fill = DIVIDER_COLOR, inherit.aes = F,
              aes(ymin = start - 0.5, ymax = end + 0.5, xmin = -Inf, 
                  xmax = Inf), alpha = 0.5) +
    geom_segment(data = grShadeCoords$gr_coords, inherit.aes = F,
                 aes(y = start, yend = end, 
                     x = 100 * grShadeCoords$lab_coords$y,
                     xend = 100 * grShadeCoords$lab_coords$y)) +
    geom_text(mapping = aes(y = x, x = 101 * y, label = gr, angle = 0, 
                            hjust = -0.01), inherit.aes = F,
              data = grShadeCoords$lab_coords, size = 2)
}

RIGHT_PLOT <- RIGHT_PLOT +
  geom_bar(data = RIGHT_PLOTData,
           mapping = aes(y = gene_plots_id, x = perc_driverMut_mle,
                         fill = specificity, color = is_driver,
                         alpha = is_driver), 
           stat = 'identity', linewidth = 0.25) +
  geom_point(data = RIGHT_PLOTData,
             mapping = aes(y = gene_plots_id, x = perc_driverMut_mle),
             color = 'black', size = 1) +
  geom_errorbar(data = RIGHT_PLOTData,
                mapping = aes(y = gene_plots_id, xmin = perc_driverMut_low,
                              xmax = perc_driverMut_high),
                color = 'black', linewidth = 0.25, width = 0.25) +
  geom_bar(data = RIGHT_PLOTData,
           mapping = aes(y = gene_plots_id, x = -15, 
                         fill = percTumsWdriverMut_bin), 
           stat = 'identity', linewidth = 0.25) +
  geom_text(data = RIGHT_PLOTData,
            mapping = aes(y = gene_plots_id, x = -8, 
                          label = round(percTumsWdriverMut)),
            size = 2) +
  scale_fill_manual(values = c(percTumsWdriverMutColor,
                               subtypeSpecificityPalette), drop = F) +
  scale_color_manual(values = c('FALSE' = '#CCCCCC', 'TRUE' = 'black')) +
  scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 1)) +
  scale_x_continuous(breaks = seq(0, 100, by = 25), 
                     labels = seq(0, 100, by = 25),
                     expand = c(0,0)) +
  scale_y_discrete(labels = xAxisLabs$label, drop = F) +
  xlab('% excess of mutations') + ylab('') + 
  labs(color = 'is driver\n in tumor subtype', 
       alpha = 'is driver\n in tumor subtype') +
  customGgplot2Theme + ggtitle(subtypeOrder[2]) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        axis.line.y = element_blank(),
        axis.text.y = element_text(color = xAxisLabs$color),
        plot.title = element_text(size = 8))

# Arrange plots ---------------------------------------------------------------
# width 3 heigth 5
ggarrange(LEFT_PLOT + theme(legend.position = 'none'),
          RIGHT_PLOT + theme(legend.position = 'none'), 
          nrow = 1, align = 'h', widths = c(3, 4))

plot(extract_legend(LEFT_PLOT))