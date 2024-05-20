#!/usr/bin/env Rscript
# FILE: figures_create_subtype_specificity_plot.R -----------------------------
#
# DESCRIPTION: An R script which creates a subtype specificity plot showing 
# drivers' preference to occur in one tumor subtype vs the other one.
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla figures_create_subtype_specificity_plot.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: R v4.1.0, data.table, ggplot2, ggpubr, plyr
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  01.02.2023
# REVISION: 20.05.2024

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
Sys.sleep(sample(1:15, 1))
source(paste0(srcDir, '/custom_functions_for_figures.R'))

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ggpubr)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))

# Functions -------------------------------------------------------------------
#' readAndSumUpDriverMutations
#' @description Reads file with detected driver mutations and summarises it 
#' (i.e. MLE % of driver mutations in driver genomic element)
#' @author Maria Litovchenko
#' @param filePath string, path to detected driver mutations file, columns
#' is_driver, gr_id, gene_id, gene_name, var_type, key, participant_id, 
#' n_total, n_driverMut_low, n_driverMut_mle, n_driverMut_high, 
#' n_tums_w_driver_mut are needed.
#' @param tumor_subtype string, tumor subtype in which driver mutations were 
#' detected
#' @param allDetectedDrivers data table with drivers which were detected in at
#' least one tumor subtype participating in analysis. Columns gr_id, gene_id,
#' gene_name are needed. No duplicated rows.
#' @return data table with columns tumor_subtype, is_driver, gr_id, gene_id,
#' gene_name, is_known_cancer, perc_driverMut_low,perc_driverMut_mle, 
#' perc_driverMut_high, n_participants, n_muts, n_tums_w_driver_mut
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

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'figures_create_subtype_specificity_plot.R')

patientsHelp <- paste('Path to table listing information about patients,',
                      'their cancer types and mutation files. Minimal',
                      'columns: participant_id, tumor_subtype,',
                      'participant_tumor_subtype, somatic_path,',
                      'somatic_genome, cohort_name.')
parser$add_argument("-p", "--inventory_patients", required = T, 
                    type = 'character', help = patientsHelp)

subtype1Help <- paste('A name of the first tumor subtype in the comparison',
                      'pair.')
parser$add_argument("-c1", "--cancer_subtype_1", required = T, 
                    type = 'character', help = subtype1Help)

subtype1DriversHelp <- paste('Path to file containing de-novo detected',
                              'cancer drivers in the first cancer subtype (',
                             '--cancer_subtype_1)')
parser$add_argument("-d1", "--drivers_1", required = T, type = 'character',
                    default = NULL, help = subtype1DriversHelp)

mutations1Help <- paste('Path to detected driver mutations file in the first',
                        'cancer subtype (--cancer_subtype_1). Columns',
                        'is_driver, gr_id, gene_id, gene_name, var_type, key,',
                        'participant_id, n_total, n_driverMut_low,',
                        'n_driverMut_mle, n_driverMut_high,',
                        'n_tums_w_driver_mut are needed.')
parser$add_argument("-m1", "--driver_mutations_1", required = T, 
                    type = 'character', default = NULL, help = mutations1Help)

subtype2Help <- paste('A name of the second tumor subtype in the comparison',
                      'pair.')
parser$add_argument("-c2", "--cancer_subtype_2", required = T, 
                    type = 'character', help = subtype2Help)

subtype2DriversHelp <- paste('Path to file containing de-novo detected',
                             'cancer drivers in the first cancer subtype (',
                             '--cancer_subtype_2)')
parser$add_argument("-d2", "--drivers_2", required = T, type = 'character',
                    default = NULL, help = subtype2DriversHelp)

mutations2Help <- paste('Path to detected driver mutations file in the second',
                        'cancer subtype (--cancer_subtype_2). Columns',
                        'is_driver, gr_id, gene_id, gene_name, var_type, key,',
                        'participant_id, n_total, n_driverMut_low,',
                        'n_driverMut_mle, n_driverMut_high,',
                        'n_tums_w_driver_mut are needed.')
parser$add_argument("-m2", "--driver_mutations_2", required = T,
                    type = 'character', default = NULL, help = mutations2Help)

subtypeSpecificityHelp <- paste('Path to the file with subtype specificity of',
                                'detected drivers. Required columns:',
                                'tumor_subtype_1, tumor_subtype_2, gr_id,',
                                'gene_id, gene_name, tumor_subtype_spec,',
                                'specificity_mode, is_driver_ts_1,',
                                'is_driver_ts_2.')
parser$add_argument("-s", "--subtype_specificity_file", required = T,
                    type = 'character', default = NULL, 
                    help = subtypeSpecificityHelp)

minPatientsHelp <- paste('Minimal number of patients with a mutation in',
                         'genomic region.')
parser$add_argument("-mp", "--min_n_patients", required = F, type = 'integer',
                    default = 0, help = minPatientsHelp)

jsonHelp <- paste('Path to a JSON file containing visual parameters to be',
                  'used in the plot.')
parser$add_argument("-v", "--visuals_json", required = T, 
                    type = 'character', help = jsonHelp)

outputTypeHelp <- paste('Type of image to create: pdf or png')
parser$add_argument("-ot", "--output_type", required = F, 
                    type = 'character', default = 'pdf',  
                    choices = c('pdf', 'png'), help = outputTypeHelp)

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output file")

args <- parser$parse_args()

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

NON_CNV_VAR_TYPES = c("subs", "mis", "non", "indels")
MODE_UPD <- "found as driver in a tumor subtype,\nbut is preferential/specific to the other subtype"

# Test input arguments --------------------------------------------------------
# args <- list(inventory_patients = 'data/inventory/inventory_patients.csv',
#              cancer_subtype_1 = 'Adenocarcinoma',
#              drivers_1 = "completed_runs/2023_12_25/results/tables/drivers/drivers-Adenocarcinoma--hg19.csv",
#              driver_mutations_1 = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Adenocarcinoma--hg19.csv",
#              cancer_subtype_2 = 'Squamous_cell',
#              driver_mutations_2 = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Squamous_cell--hg19.csv",
#              drivers_2 = "completed_runs/2023_12_25/results/tables/drivers/drivers-Squamous_cell--hg19.csv",
#              subtype_specificity_file = 'completed_runs/2023_12_25/results/tables/subtype_specificity/subtypeSpecificity---hg19.csv',
#              min_n_patients = 3)

# Read visuals JSON -----------------------------------------------------------
ESSENTIAL_VISUAL_NAMES <- c('ggplot2_theme', 'color_divider', 
                            'colors_perc_tumor_with_driver_mutation',
                            'colors_subtype_specificity', 
                            'subtype_specificity_plot_width', 
                            'subtype_specificity_plot_heigth')

visualParams <- readJsonWithVisualParameters(args$visuals_json)
message('[', Sys.time(), '] Read --visuals_json: ', args$visuals_json)

notFoundVisuals <- setdiff(ESSENTIAL_VISUAL_NAMES, names(visualParams))
if (length(notFoundVisuals)) {
  stop('[', Sys.time(), '] Following visuals: ', 
       paste(notFoundVisuals, collapse = ', '), ' not found in JSON.')
}

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
message('[', Sys.time(), '] Read --subtype_specificity_file: ',  
        args$subtype_specificity_file)
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
message('[', Sys.time(), '] Read --driver_mutations_1: ',  
        args$driver_mutations_1)
message('[', Sys.time(), '] Read --driver_mutations_2: ',  
        args$driver_mutations_2)
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

visualParams$colors_perc_tumor_with_driver_mutation <- colorRampPalette(visualParams$colors_perc_tumor_with_driver_mutation)(length(percTumsWdriverMutLabels))
names(visualParams$colors_perc_tumor_with_driver_mutation) <- percTumsWdriverMutLabels

# X axis labels and color palettes --------------------------------------------
xAxisLabs <- formatAxisLabels(subtypeSpec, 'gene_plots_id', 'gene_name',
                              'is_known_cancer', highlightColor = 'red')

names(visualParams$colors_subtype_specificity) <- gsub('leftSubtype', 
                                                       subtypeOrder[1], 
                                                       names(visualParams$colors_subtype_specificity))
names(visualParams$colors_subtype_specificity) <- gsub('rightSubtype', 
                                                       subtypeOrder[2], 
                                                       names(visualParams$colors_subtype_specificity))
names(visualParams$colors_subtype_specificity) <- gsub('MODE_UPD', 
                                                       MODE_UPD, 
                                                       names(visualParams$colors_subtype_specificity))

# Left plot -------------------------------------------------------------------
LEFT_PLOTData <- subtypeSpec[tumor_subtype == subtypeOrder[1]]
LEFT_PLOT <- ggplot()

# first of all, if there are several genomic regons, create shadowing
if (length(unique(LEFT_PLOTData$gr_id)) > 1) {
  # compute coordinates for future lines showing genome region 
  grShadeCoords <- get_gr_change_coords(plotDT = LEFT_PLOTData, 
                                        axisOrderCol = "gene_plots_id",
                                        yCol = 'perc_driverMut_mle',
                                        labsDT = xAxisLabs$dt)
  grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                  nrow(grShadeCoords$gr_coords),
                                                  by = 2)]
  LEFT_PLOT <- LEFT_PLOT + 
    geom_rect(data = grShadeCoordsCut, fill = visualParams$color_divider,
              inherit.aes = F,
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
  scale_fill_manual(values = c(visualParams$colors_perc_tumor_with_driver_mutation,
                               visualParams$colors_subtype_specificity), 
                    drop = F) +
  scale_color_manual(values = c('FALSE' = '#CCCCCC', 'TRUE' = 'black')) +
  scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 1)) +
  scale_x_continuous(breaks = seq(0, -100, by = -25), 
                     labels = seq(0, 100, by = 25),
                     expand = c(0,0)) +
  scale_y_discrete(drop = F, position = 'right') +
  xlab('% excess of mutations') + ylab('') + 
  labs(color = 'is driver\n in tumor subtype', 
       alpha = 'is driver\n in tumor subtype') +
  visualParams$ggplot2_theme + ggtitle(subtypeOrder[1]) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        axis.text.y = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom", legend.direction = 'horizontal',
        plot.title = element_text(size = 8))

# Right plot ------------------------------------------------------------------------
RIGHT_PLOTData <- subtypeSpec[tumor_subtype == subtypeOrder[2]]
RIGHT_PLOT <- ggplot()

# first of all, if there are several genomic regons, create shadowing
if (length(unique(RIGHT_PLOTData$gr_id)) > 1) {
  # compute coordinates for future lines showing genome region 
  grShadeCoords <- get_gr_change_coords(plotDT = RIGHT_PLOTData, 
                                        axisOrderCol = "gene_plots_id",
                                        yCol = 'perc_driverMut_mle',
                                        labsDT = xAxisLabs$dt)
  grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                  nrow(grShadeCoords$gr_coords),
                                                  by = 2)]
  RIGHT_PLOT <- RIGHT_PLOT + 
    geom_rect(data = grShadeCoordsCut, fill = visualParams$color_divider, 
              inherit.aes = F,
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
  scale_fill_manual(values = c(visualParams$colors_perc_tumor_with_driver_mutation,
                               visualParams$colors_subtype_specificity),
                    drop = F) +
  scale_color_manual(values = c('FALSE' = '#CCCCCC', 'TRUE' = 'black')) +
  scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 1)) +
  scale_x_continuous(breaks = seq(0, 100, by = 25), 
                     labels = seq(0, 100, by = 25),
                     expand = c(0,0)) +
  scale_y_discrete(labels = xAxisLabs$label, drop = F) +
  xlab('% excess of mutations') + ylab('') + 
  labs(color = 'is driver\n in tumor subtype', 
       alpha = 'is driver\n in tumor subtype') +
  visualParams$ggplot2_theme + ggtitle(subtypeOrder[2]) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        axis.line.y = element_blank(),
        axis.text.y = element_text(color = xAxisLabs$color),
        legend.position = "bottom", legend.direction = 'horizontal',
        plot.title = element_text(size = 8))

# Arrange plots ---------------------------------------------------------------
# width 3 height 5
SUBTYPE_SPECIFICITY_PLOT <- ggarrange(ggarrange(LEFT_PLOT + 
                                        theme(legend.position = 'none'),
                                      RIGHT_PLOT + 
                                        theme(legend.position = 'none'), 
                                      nrow = 1, align = 'h', widths = c(3, 4)),
                                      extract_legend(LEFT_PLOT), ncol = 1,
                                      heights = c(1, 0.15))

# Output to file --------------------------------------------------------------
if (args$output_type == 'pdf') {
  pdf(args$output, width = visualParams$subtype_specificity_plot_width, 
      height = visualParams$subtype_specificity_plot_heigth)
} else {
  png(args$output, width = visualParams$subtype_specificity_plot_width, 
      height = visualParams$subtype_specificity_plot_heigth)
}
SUBTYPE_SPECIFICITY_PLOT
dev.off()

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')