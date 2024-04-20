library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)

source('bin/custom_functions.R')

# Functions : genomic regions coloring ----------------------------------------
#' format_xLabels
#' @description Formats labels for x Axis, i.e. color and actual string to use
#'              as label
#' @author Maria Litovchenko
#' @param plotDT data table with columns gr_id, FILTER, and columns with names
#'               encoded in xOrderCol and xLabelCol
#' @param xLabelCol column name containing label to be on X axis
#' @param xOrderCol column name containing order of X axis
#' @param annoCol column containing annotation, i.e. COSMIC 
#' @param defaultColor default label color
#' @param highlightColor highlight label color(s). If several distinct values
#'        (levels) exist in annoCol and the same number of colors is given,
#'        each lavel will have distinct color.
#' @return list with members dt, labels, color
format_xLabels <- function(plotDT, xOrderCol, xLabelCol, annoCol, 
                           defaultColor = 'black', highlightColor = 'red') {
  # x axis labels
  colsToGet <- intersect(c('gr_id', 'FILTER', xOrderCol, xLabelCol),
                         colnames(plotDT))
  xLabsDT <- plotDT[, colsToGet, with = F]
  xLabsDT <- unique(xLabsDT)
  xLabs <- unlist(xLabsDT[, xLabelCol, with = F])
  names(xLabs) <- unlist(xLabsDT[, xOrderCol, with = F])
  
  # x axis labels color
  xColor <- defaultColor
  if (!is.null(annoCol)) {
    annoDT <- plotDT[, c(xLabelCol, xOrderCol, annoCol), with = F]
    annoDT <- unique(annoDT)
    annoDT <- annoDT[complete.cases(annoDT)]
    setkeyv(annoDT, xOrderCol)
    annoDT <- annoDT[order(unlist(annoDT[, xOrderCol, with = F]))]
    if (length(unique(annoDT)) > length(highlightColor)) {
      message('[', Sys.time(), '] There is more unique non NA values in ',
              annoCol, ' than colors in highlightColor. Will use the first ',
              'color to highlight all.')
      xColor <- ifelse(unlist(annoDT[, annoCol, with = F]), defaultColor,
                       highlightColor[1])
    } else {
      xColor <- unlist(annoDT[, annoCol, with = F])
      annoVals_Lvls <- sort(unique(unlist(annoDT[, annoCol, with = F])))
      print(annoVals_Lvls)
      print(sort(names(highlightColor)))
      if (identical(sort(as.character(annoVals_Lvls)),
                    sort(names(highlightColor)))) {
        for (annoUniq_i in annoVals_Lvls) {
          xColor <- gsub(paste0('^', annoUniq_i, '$'),
                         highlightColor[annoUniq_i], xColor)
        }
      } else {
        for (idx in 1:length(annoVals_Lvls)) {
          xColor <- gsub(annoVals_Lvls[idx], highlightColor[idx], xColor)
        }
      }
    }
    names(xColor) <- unlist(annoDT[, xOrderCol, with = F])
    
    annoDT <- plotDT[, c(xLabelCol, xOrderCol, annoCol), with = F]
    annoDT <- unique(annoDT)
    annoDT <- annoDT[!complete.cases(annoDT)]
    toAdd <- rep(defaultColor, nrow(annoDT))
    names(toAdd) <-  unlist(annoDT[, xOrderCol, with = F])
    xColor <- c(xColor, toAdd)
  }
  xColor <- xColor[as.character(sort(unlist(xLabsDT[, xOrderCol, with = F])))]
  result <- list(dt = xLabsDT, labels = xLabs, color = xColor)
  result
}

#' get_gr_change_coords
#' @description Create coordinates for future shadowed regions showing 
#'              different genomic region
#' @author Maria Litovchenko
#' @param plotDT data table with column with the name submitted as xOrderCol 
#'               and column with the name submitted as yCol
#' @param xOrderCol column name containing order of X axis
#' @param yCol column name containing data for Y axis
#' @param xLabsDT data table, result of format_xLabels, field dt
#' @return list with members gr_coords and lab_coords
get_gr_change_coords <- function(plotDT, xOrderCol, yCol, xLabsDT) {
  # x coordinates for changes of genomic regions
  xLabsDT <- xLabsDT[order(unlist(xLabsDT[, xOrderCol, with = F]))]
  regChange <- sapply(2:nrow(xLabsDT),
                      function(x) xLabsDT$gr_id[x] != xLabsDT$gr_id[x - 1])
  regChange <- which(regChange)
  x_coords <- data.table(start = c(1, regChange + 1))
  x_coords[, end := c(start[-1] - 1, nrow(xLabsDT))]
  
  # y coordinates of labels for genomic regions
  yVals <- unlist(plotDT[, yCol, with = F])
  if (class(yVals) == 'factor') {
    yMax <- factor(max(levels(yVals)), levels = levels(yVals))
  } else {
    if (class(yVals) %in% c('integer', 'numberic', 'double')) {
      plotDT[, toSumUp := yVals]
      yMax <- 1.01 * max(plotDT[,.(sum(toSumUp)), by = eval(xOrderCol)]$V1)
      plotDT[, toSumUp := NULL]
    } else {
      yMax <- 1.01
    }
  } 
  
  grLabs <- sapply(2:nrow(xLabsDT), 
                   function(x) xLabsDT$gr_id[x - 1] != xLabsDT$gr_id[x])
  grLabs <- xLabsDT$gr_id[c(T, grLabs)]
  lab_coords <- data.table(x = rowMeans(x_coords), y = yMax, gr = grLabs)
  
  result <- list('gr_coords' = x_coords, 'lab_coords' = lab_coords)
  result
}

#' get_filtered_out_gene_coords
#' @description Create coordinates for future shadowed regions showing 
#'              genes drivers which were filtered out
#' @author Maria Litovchenko
#' @param xLabsDT data table, result of format_xLabels, field dt
#' @param xOrderCol column name containing order of X axis
#' @return data table with columns start and end
get_filtered_out_gene_coords <- function(xLabsDT, xOrderCol) {
  xLabsDT <- xLabsDT[order(unlist(xLabsDT[, xOrderCol, with = F]))]
  xLabsDT[, removed := FILTER == 'PASS']
  rmChange <- sapply(2:nrow(xLabsDT),
                     function(x) xLabsDT$removed[x] != xLabsDT$removed[x - 1])
  rmChange <- which(rmChange)
  result <- data.table(start = c(1, rmChange + 1))
  result[, end := c(start[-1] - 1, nrow(xLabsDT))]
  result <- result[apply(result, 1, 
                         function(x) all(!xLabsDT[seq(x['start'], 
                                                      x['end'])]$removed))]
  result
}

# Function: dealing with legends ----------------------------------------------
#' extract_legend 
#' @description Extracts legend from ggplot
#' @author The internet
#' @param my_ggp ggplot2 object
#' @return legend as ggplot2
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  step3
}

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

# Visuals setup ---------------------------------------------------------------
font_config <- list('family' = 'Helvetica', 'base' = 6, 'text' = 6, 
                    'axis_title' = 8, 'plot_title' = 10, 'legend_text' = 6, 
                    'legend_title' = 8)

customGgplot2Theme <- list(
  theme_classic(base_size = font_config$base, 
                base_family = font_config$family) +
    theme(axis.line = element_line(colour = 'black', linewidth = 0.25,
                                   linetype = 'solid'),
          axis.ticks = element_line(colour = 'black', linewidth = 0.25,
                                    linetype = 'solid'),
          axis.text = element_text(colour = 'black', 
                                   family = font_config$family,
                                   size = font_config$text),
          axis.title = element_text(colour = 'black', 
                                    family = font_config$family,
                                    size = font_config$axis_title),
          plot.title = element_text(colour = 'black', 
                                    family = font_config$family,
                                    size = font_config$plot_title),
          legend.key.size = unit(0.25, 'cm'),
          legend.text = element_text(colour = 'black', 
                                     family = font_config$family,
                                     size = font_config$legend_text), 
          legend.title = element_text(colour = 'black', 
                                      family = font_config$family,
                                      size = font_config$legend_title),
          legend.background = element_blank(), legend.direction = "vertical",
          legend.box.background = element_blank(),
          panel.grid.major = element_line(colour = "#CCCCCC", 
                                          linewidth = 0.25, linetype = 2),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm")))

DIVIDER_COLOR <-  '#eed9c4'

subtypeSpecificityPalette <- c('leftSubtype-specific' = '#CC3311', 
                               'leftSubtype-preferential' = '#f8b8a9', 
                               'common' = '#9970AB', 
                               'MODE_UPD' = '#AAAAAA',
                               'not enough evidence' = '#F0E6B2',
                               'rightSubtype-preferential' = '#a7dfff',
                               'rightSubtype-specific' = '#0077BB')
percTumorsWithDriverMutsPalette <- c('#FFFFFF', '#D9F0D3', '#ACD39E',
                                     '#5AAE61', '#1B7837')

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
                                           c(subtypeOrder[2], 
                                             "common", 
                                             "not enough evidence",
                                             subtypeOrder[1],
                                             MODE_UPD))]
subtypeSpec[, specificity_mode := NULL]

subtypeSpec <- subtypeSpec[order(-perc_driverMut_mle)]
subtypeSpec[, gr_id := factor(gr_id, rev(unique(gr_id)))]
subtypeSpec <- subtypeSpec[order(gr_id, tumor_subtype_spec, perc_driverMut_mle)]
subtypeSpec[, gene_plots_id := paste(gene_name, gr_id)]
subtypeSpec[, gene_plots_id := factor(gene_plots_id, unique(gene_plots_id))]

# Prepare % of patients with driver mutations for plotting --------------------
subtypeSpec <- merge(subtypeSpec, cohort_sizes, by = "tumor_subtype")
subtypeSpec[, perc_tums_w_driver_mut := 100 * n_tums_w_driver_mut / cohort_size]

perc_tums_w_driver_mut_breaks <- c(-0.01, 0.5, 5, 10)
perc_tums_w_driver_mut_labels <- c('<1%', '5%', '10%')
if (max(subtypeSpec$perc_tums_w_driver_mut) > 0.1) {
  maxPerc <- max(subtypeSpec$perc_tums_w_driver_mut)
  perc_tums_w_driver_mut_breaks <- c(perc_tums_w_driver_mut_breaks, 
                                     seq(20, maxPerc, by = 10))
  perc_tums_w_driver_mut_labels <- c(perc_tums_w_driver_mut_labels,
                                     paste0(seq(20, maxPerc, by = 10), '%'))
}

subtypeSpec[, perc_tums_w_driver_mut_bin := cut(perc_tums_w_driver_mut, 
                                                perc_tums_w_driver_mut_breaks, 
                                                perc_tums_w_driver_mut_labels)]
maxLvl <- perc_tums_w_driver_mut_labels[length(perc_tums_w_driver_mut_labels)]
subtypeSpec[is.na(perc_tums_w_driver_mut_bin)]$perc_tums_w_driver_mut_bin <- maxLvl

perc_tums_w_driver_mut_color <- colorRampPalette(percTumorsWithDriverMutsPalette)
perc_tums_w_driver_mut_color <- perc_tums_w_driver_mut_color(length(perc_tums_w_driver_mut_labels))
names(perc_tums_w_driver_mut_color) <- perc_tums_w_driver_mut_labels

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
leftPlotData <- subtypeSpec[tumor_subtype == subtypeOrder[1]]
leftPlot <- ggplot()

# first of all, if there are several genomic regons, create shadowing
if (length(unique(leftPlotData$gr_id)) > 1) {
  # compute coordinates for future lines showing genome region 
  grShadeCoords <- get_gr_change_coords(plotDT = leftPlotData, 
                                        xLabsDT = xAxisLabs$dt,
                                        xOrderCol = "gene_plots_id")
  grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                  nrow(grShadeCoords$gr_coords),
                                                  by = 2)]
  leftPlot <- leftPlot + 
    geom_rect(data = grShadeCoordsCut, fill = DIVIDER_COLOR, inherit.aes = F,
              aes(ymin = start - 0.5, ymax = end + 0.5, xmin = -Inf, 
                  xmax = Inf), alpha = 0.5)
}

leftPlot <- leftPlot +
  geom_bar(data = leftPlotData,
           mapping = aes(y = gene_plots_id, x = -perc_driverMut_mle,
                         fill = specificity, color = is_driver,
                         alpha = is_driver), 
           stat = 'identity', linewidth = 0.25) +
  geom_point(data = leftPlotData,
             mapping = aes(y = gene_plots_id, x = -perc_driverMut_mle),
             color = 'black', size = 1) +
  geom_errorbar(data = leftPlotData,
                mapping = aes(y = gene_plots_id, xmin = -perc_driverMut_low,
                              xmax = -perc_driverMut_high),
                color = 'black', linewidth = 0.25, width = 0.25) +
  geom_bar(data = leftPlotData,
           mapping = aes(y = gene_plots_id, x = 15, 
                         fill = perc_tums_w_driver_mut_bin), 
           stat = 'identity', linewidth = 0.25) +
  geom_text(data = leftPlotData,
            mapping = aes(y = gene_plots_id, x = 8, 
                          label = round(perc_tums_w_driver_mut)),
            size = 2) +
  scale_fill_manual(values = c(perc_tums_w_driver_mut_color,
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
rightPlotData <- subtypeSpec[tumor_subtype == subtypeOrder[2]]
rightPlot <- ggplot()

# first of all, if there are several genomic regons, create shadowing
if (length(unique(rightPlotData$gr_id)) > 1) {
  # compute coordinates for future lines showing genome region 
  grShadeCoords <- get_gr_change_coords(plotDT = rightPlotData, 
                                        xLabsDT = xAxisLabs$dt,
                                        xOrderCol = "gene_plots_id")
  grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                  nrow(grShadeCoords$gr_coords),
                                                  by = 2)]
  rightPlot <- rightPlot + 
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

rightPlot <- rightPlot +
  geom_bar(data = rightPlotData,
           mapping = aes(y = gene_plots_id, x = perc_driverMut_mle,
                         fill = specificity, color = is_driver,
                         alpha = is_driver), 
           stat = 'identity', linewidth = 0.25) +
  geom_point(data = rightPlotData,
             mapping = aes(y = gene_plots_id, x = perc_driverMut_mle),
             color = 'black', size = 1) +
  geom_errorbar(data = rightPlotData,
                mapping = aes(y = gene_plots_id, xmin = perc_driverMut_low,
                              xmax = perc_driverMut_high),
                color = 'black', linewidth = 0.25, width = 0.25) +
  geom_bar(data = rightPlotData,
           mapping = aes(y = gene_plots_id, x = -15, 
                         fill = perc_tums_w_driver_mut_bin), 
           stat = 'identity', linewidth = 0.25) +
  geom_text(data = rightPlotData,
            mapping = aes(y = gene_plots_id, x = -8, 
                          label = round(perc_tums_w_driver_mut)),
            size = 2) +
  scale_fill_manual(values = c(perc_tums_w_driver_mut_color,
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
ggarrange(leftPlot + theme(legend.position = 'none'),
          rightPlot + theme(legend.position = 'none'), 
          nrow = 1, align = 'h', widths = c(3, 4))

plot(extract_legend(leftPlot))