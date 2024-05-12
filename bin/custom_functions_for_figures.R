#!/usr/bin/env Rscript
# FILE: custom_functions_for_figures.R ----------------------------------------
#
# DESCRIPTION: Custom created R functions which are used in creation of plots
#             
# USAGE: In R script, insert:
#        source('custom_functions_for_figures.R')
#
# OPTIONS: None
#
# REQUIREMENTS: R v4.1.0
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  10.05.2024
# REVISION: 10.05.2024

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ggpubr)))
suppressWarnings(suppressPackageStartupMessages(library(plyr))) 

# Functions : delineating genomic regions with use of colored background ------
#' get_gr_change_coords
#' @description Create coordinates for future shadowed regions showing 
#'              different genomic region
#' @author Maria Litovchenko
#' @param plotDT data table with column with the name submitted as axisOrderCol 
#'               and column with the name submitted as yCol
#' @param axisOrderCol column name containing order of X axis
#' @param yCol column name containing data for Y axis
#' @param labsDT data table, result of formatAxisLabels, field dt
#' @return list with members gr_coords and lab_coords
get_gr_change_coords <- function(plotDT, axisOrderCol, yCol, labsDT) {
  # x coordinates for changes of genomic regions
  labsDT <- labsDT[order(unlist(labsDT[, axisOrderCol, with = F]))]
  regChange <- sapply(2:nrow(labsDT),
                      function(x) labsDT$gr_id[x] != labsDT$gr_id[x - 1])
  regChange <- which(regChange)
  x_coords <- data.table(start = c(1, regChange + 1))
  x_coords[, end := c(start[-1] - 1, nrow(labsDT))]
  
  # y coordinates of labels for genomic regions
  yVals <- unlist(plotDT[, yCol, with = F])
  if (class(yVals) == 'factor') {
    yMax <- factor(max(levels(yVals)), levels = levels(yVals))
  } else {
    if (class(yVals) %in% c('integer', 'numberic', 'double')) {
      plotDT[, toSumUp := yVals]
      yMax <- 1.01 * max(plotDT[,.(sum(toSumUp)), by = eval(axisOrderCol)]$V1)
      plotDT[, toSumUp := NULL]
    } else {
      yMax <- 1.01
    }
  } 
  
  grLabs <- sapply(2:nrow(labsDT), 
                   function(x) labsDT$gr_id[x - 1] != labsDT$gr_id[x])
  grLabs <- labsDT$gr_id[c(T, grLabs)]
  lab_coords <- data.table(x = rowMeans(x_coords), y = yMax, gr = grLabs)
  
  result <- list('gr_coords' = x_coords, 'lab_coords' = lab_coords)
  result
}

#' get_filtered_out_gene_coords
#' @description Create coordinates for future shadowed regions showing 
#'              genes drivers which were filtered out
#' @author Maria Litovchenko
#' @param labsDT data table, result of formatAxisLabels, field dt
#' @param axisOrderCol column name containing order of X axis
#' @return data table with columns start and end
get_filtered_out_gene_coords <- function(labsDT, axisOrderCol) {
  labsDT <- labsDT[order(unlist(labsDT[, axisOrderCol, with = F]))]
  labsDT[, removed := FILTER == 'PASS']
  rmChange <- sapply(2:nrow(labsDT),
                     function(x) labsDT$removed[x] != labsDT$removed[x - 1])
  rmChange <- which(rmChange)
  result <- data.table(start = c(1, rmChange + 1))
  result[, end := c(start[-1] - 1, nrow(labsDT))]
  result <- result[apply(result, 1, 
                         function(x) all(!labsDT[seq(x['start'], 
                                                     x['end'])]$removed))]
  result
}

#' create_plot_background
#' @description Creates "background" for the future plots of coding and 
#' noncoding driver genetic elements. It creates shades to distinguish between
#' different types of genetic elements (i.e. CDS and promoters) as well as 
#' creates a separate shade for the filtered out genetic elements.
#' @author Maria Litovchenko
#' @param plotDT data table with processed data, ready to be plotted (main 
#' plot)
#' @param axisLabs result of function formatAxisLabels
#' @param axisOrderCol column name containing order of X axis
#' @param yCol column name containing data for Y axis
#' @param divider_color hex color which should be used for separating different
#' genomic regions
#' @param keep_region_labels boolean, whatever or not labels of regions on the
#' top of the plot should be kept.
#' @param direction string, one of "vertical" and "horizontal"
#' @return ggplot2 object
create_plot_background <- function(plotDT, axisLabs, axisOrderCol, yCol,
                                   divider_color = '#eed9c4', 
                                   direction = 'vertical') {
  if (!direction %in% c('vertical', 'horizontal')) {
    stop('[', Sys.time(), '] direction should be vertical or horizontal')
  }
  
  # first of all, if there are several genomic regions, create shadowing
  result <- ggplot()
  if (length(unique(plotDT$gr_id)) > 1) {
    # compute coordinates for future lines showing genome region 
    grShadeCoords <- get_gr_change_coords(plotDT, axisOrderCol, yCol, 
                                          axisLabs$dt)
    nGRs <- nrow(grShadeCoords$gr_coords)
    grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, nGRs, by = 2)]
    if (direction == 'vertical') {
      result <- result + 
        geom_rect(data = grShadeCoordsCut, fill = divider_color, alpha = 0.5,
                  mapping = aes(xmin = start - 0.5, xmax = end + 0.5, 
                                ymin = -Inf, ymax = Inf), inherit.aes = F)
    } else {
      result <- result + 
        geom_rect(data = grShadeCoordsCut, fill = divider_color, alpha = 0.5,
                  mapping = aes(ymin = start - 0.5, ymax = end + 0.5, 
                                xmin = -Inf, xmax = Inf), inherit.aes = F)
    }
  }
  
  # and if there are some filtered out drivers needed to be shown - shadowed 
  # areas for them
  if (any(plotDT$FILTER != 'PASS')) {
    filteredOutShadeCoords <- get_filtered_out_gene_coords(axisLabs$dt, axisOrderCol)
    filteredOutShadeCoordsCut <- filteredOutShadeCoords[seq(1,
                                                            nrow(filteredOutShadeCoords),
                                                            by = 2)]
    if (direction == 'vertical') {
      result <- result + 
        geom_rect(data = filteredOutShadeCoordsCut, fill = 'red',
                  inherit.aes = F, 
                  mapping = aes(xmin = start - 0.5, xmax = end + 0.5,
                                ymin = -Inf, ymax = Inf), alpha = 0.15)
    } else {
      result <- result + 
        geom_rect(data = filteredOutShadeCoordsCut, fill = 'red', 
                  inherit.aes = F, 
                  mapping = aes(ymin = start - 0.5, ymax = end + 0.5, 
                                xmin = -Inf, xmax = Inf), alpha = 0.15)
    }
  }
  
  result
}

# Functions : genomic regions coloring ----------------------------------------
#' formatAxisLabels
#' @description Formats labels for one axis, i.e. color and actual string to 
#' use as label
#' @author Maria Litovchenko
#' @param plotDT data table with columns gr_id, FILTER, and columns with names
#' encoded in axisOrderCol and axisLabelCol
#' @param axisLabelCol column name containing label to be on X axis
#' @param axisOrderCol column name containing order of X axis
#' @param annoCol column containing annotation, i.e. COSMIC 
#' @param defaultColor default label color
#' @param highlightColor highlight label color(s). If several distinct values
#' (levels) exist in annoCol and the same number of colors is given, each level
#' will have distinct color.
#' @return list with members dt, labels, color
formatAxisLabels <- function(plotDT, axisOrderCol, axisLabelCol, annoCol, 
                             defaultColor = 'black', highlightColor = 'red') {
  colsToGet <- intersect(c('gr_id', 'FILTER', axisOrderCol, axisLabelCol),
                         colnames(plotDT))
  labsDT <- plotDT[, colsToGet, with = F]
  labsDT <- unique(labsDT)
  xLabs <- unlist(labsDT[, axisLabelCol, with = F])
  names(xLabs) <- unlist(labsDT[, axisOrderCol, with = F])
  
  # labels color
  xColor <- defaultColor
  if (!is.null(annoCol)) {
    annoDT <- plotDT[, c(axisLabelCol, axisOrderCol, annoCol), with = F]
    annoDT <- unique(annoDT)
    annoDT <- annoDT[complete.cases(annoDT)]
    setkeyv(annoDT, axisOrderCol)
    annoDT <- annoDT[order(unlist(annoDT[, axisOrderCol, with = F]))]
    if (length(unique(annoDT)) > length(highlightColor)) {
      message('[', Sys.time(), '] There is more unique non NA values in ',
              annoCol, ' than colors in highlightColor. Will use the first ',
              'color to highlight all.')
      xColor <- ifelse(unlist(annoDT[, annoCol, with = F]), defaultColor,
                       highlightColor[1])
    } else {
      xColor <- unlist(annoDT[, annoCol, with = F])
      annoVals_Lvls <- sort(unique(unlist(annoDT[, annoCol, with = F])))
      if (all(as.character(annoVals_Lvls) %in% names(highlightColor))) {
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
    names(xColor) <- unlist(annoDT[, axisOrderCol, with = F])
    
    annoDT <- plotDT[, c(axisLabelCol, axisOrderCol, annoCol), with = F]
    annoDT <- unique(annoDT)
    annoDT <- annoDT[!complete.cases(annoDT)]
    toAdd <- rep(defaultColor, nrow(annoDT))
    names(toAdd) <-  unlist(annoDT[, axisOrderCol, with = F])
    xColor <- c(xColor, toAdd)
  }
  xColor <- xColor[as.character(sort(unlist(labsDT[, axisOrderCol, 
                                                   with = F])))]
  result <- list(dt = labsDT, labels = xLabs, color = xColor)
  result
}


# Function: dealing with legends ----------------------------------------------
#' extract_legend 
#' @description Extracts legend from ggplot
#' @author The internet
#' @param ggplotObj ggplot2 object
#' @return legend as ggplot2
extract_legend <- function(ggplotObj) {
  step1 <- ggplot_gtable(ggplot_build(ggplotObj))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  step3
}

# Functions : overview of all detected driver ---------------------------------
#' get_bar_transparency
#' @description 
#' @author Maria Litovchenko
#' @param plotDT data table which holds data for plotting. Have to have columns
#'        tumor_subtype, tier, participant_tumor_subtype and a column with the 
#'        name identical to axisOrderCol
#' @param pancanCode string encoding tumor subtype code for pancancer analysis
#' @param axisOrderCol name of the column which contains IDs of genes to plot and
#'        by which they will be ordered on the plot
#' @note Driver status of the gene will control the opacity of the bars. 
#'       1) If gene is a driver in a particular tumor type(s), but not in 
#'       pancancer analysis, the bar for it (them) should be opaque and for the
#'       rest of the tumor types - transparent.
#'       2) If gene is a driver in pancancer analysis, all bars should be 
#'       opaque. However, in the plotDT table now we do not have a breakdown of 
#'       number of patients per tumor type in case gene is a driver in the 
#'       pancancer analysis. Let's add the breakdown by retrieving number of 
#'       patients from plotDT They will be lines which may not have significant
#'       p-values in case a gene was found significant only!) in pancancer 
#'       analysis.
get_bar_transparency <- function(plotDT, pancanCode, axisOrderCol) {
  result <- copy(plotDT)
  # transparent, if significant in individual tumor subtypes
  result[, transparent := ifelse(tier == 0, 1, 0)]
  # all tumor_subtypes opaque, if significant in pancan analysis
  signPancan <- result[tumor_subtype == pancanCode & tier == 1]
  signPancan <- unlist(signPancan[, axisOrderCol, with = F])
  signPancan <- unlist(result[, axisOrderCol, with = F]) %in% signPancan
  result[signPancan]$transparent <- 0
  # no need to keep pancan in the table anymore
  result <- split(result, by = axisOrderCol, drop = T)
  result <- lapply(result, 
                   function(x) x[tumor_subtype != pancanCode | 
                                   (tumor_subtype == pancanCode & 
                                      !participant_tumor_subtype %in% tumor_subtype)])
  result <- do.call(rbind, result)
  result[, tumor_subtype := NULL]
  setnames(result, 'participant_tumor_subtype', 'tumor_subtype')
  if (class(plotDT$tumor_subtype) == 'factor') {
    result[, tumor_subtype := factor(tumor_subtype, 
                                     levels = levels(plotDT$tumor_subtype))]
  }
  
  result
}

#' create_tumortype_colorpalette
#' @description Creates color palette for tumor subtypes
#' @author Maria Litovchenko
#' @param plotDT data table with data for plotting. Have to have columns 
#'        tumor_subtype
#' @param colorpalette character vector, color palette submitted by user
#' @return named character vector with colors for tumor types. Colors are 
#'.        repeated for the case then driver is significant in that tumor type
#'         (--0) and not (--1)
create_tumortype_colorpalette <- function(plotDT, colorpalette) {
  if (is.null(colorpalette)) {
    if (is.null(levels(plotDT$tumor_subtype))) {
      tumor_types_lvl <- as.character(unique(plotDT$tumor_subtype))
    } else {
      tumor_types_lvl <- levels(plotDT$tumor_subtype)
    }
    colorpalette <- palette(rainbow(length(tumor_types_lvl)))
    names(colorpalette) <- tumor_types_lvl
  }
  
  colorsToUse <- colorpalette[names(colorpalette) %in% plotDT$tumor_subtype]
  result <- rep(colorsToUse, 2)
  names(result) <- c(paste0(names(colorsToUse), '--0'), 
                     paste0(names(colorsToUse), '--1'))
  result
}

#' create_pvalue_colorpalette
#' @description Creates binned color palette for p-values
#' @author Maria Litovchenko
#' @param plotDT data table with columns colorByLog10
#' @param doLog10 whatever -log10 operation on colorBy column should be 
#'                performed
#' @param scaleLims limits of the scale
#' @param nColors number of colors to use in the color palette 
#' @param colorPalette character vector of colors
#' @return list with members colorPalette and colorBreaks
create_pvalue_colorpalette <- function(plotDT, doLog10, scaleLims, nColors,
                                       colorPalette) {
  # color breaks
  colorBreaks <- round(seq(0, max(plotDT$colorByLog10, na.rm = T), 
                           length.out = nColors + 1), 2)
  if (!is.null(scaleLims)) {
    colorBreaks <- round(seq(scaleLims[1], scaleLims[2],
                             length.out = nColors + 1), 2)
  }
  if (doLog10 == T) {
    colorBreaks <- round(colorBreaks)
  }
  # actual colors
  if (is.null(colorPalette)) {
    colorPalette <- colorRampPalette(c("grey0", "grey80"))(nColors)
    if (doLog10 == T) {
      colorPalette <- rev(colorPalette)
    }
  }
  result <- list('palette' = colorPalette, 'colorBreaks' = colorBreaks)
  result
}

#' barplotNpatient
#' @description Creates stacked barplot with number of patients per gene with 
#' histological subtypes on Y axis. In case several regions are given, they are
#' groupped and indicated with horizontal bar.
#' @author Maria Litovchenko
#' @param xLabelCol column name containing label to be on X axis
#' @param colorPalette named vector of color for different tumour type
#' @param xOrderCol column name containing order of X axis
#' @param annoCol column containing annotation, i.e. COSMIC 
#' @param annoColors named vector of colors (in hex) in which X axis labels 
#'                   will be colored
#' @param divider_color hex color which should be used for separating different
#'        genomic regions
#' @param ggplot2Theme ggplot2 custom theme
#' @return ggplot2: barplot with number of patients per gene with histological
#' subtypes on Y axis
barplotNpatient <- function(driversDT, xLabelCol, colorPalette = NULL,
                            pancanCode = 'PANCAN', xOrderCol = NULL, 
                            annoCol = NULL, annoColors = 'red',
                            divider_color = '#eed9c4',
                            ggplot2Theme = NULL) {
  if (is.null(xOrderCol)) {
    xOrderCol <- xLabelCol
  }
  colsToGet <- unique(c('tumor_subtype', 'gr_id', 'gene_id', 'gene_name',
                        'tier', 'nParts', 'participant_tumor_subtype', 
                        'FILTER', xOrderCol, xLabelCol, annoCol))
  dt <- driversDT[, colsToGet, with = F]
  # we will not display exact tier information
  dt[, tier := ifelse(is.na(tier), 0, 1)]
  
  # assign transparency values based on significance
  dt <- get_bar_transparency(dt, pancanCode, xOrderCol)
  
  # in order to put first tumor types for which gene was found a driver, and
  # then the other one for which not, we have to modify tumor_subtype by adding 
  # is_driver_in_tt to it.
  dt[, tumor_subtype_upd := paste0(tumor_subtype, '--', transparent)]
  tumor_subtype_upd_lvls <- sort(unique(dt$tumor_subtype))
  if (class(dt$tumor_subtype) == 'factor') {
    tumor_subtype_upd_lvls <- levels(dt$tumor_subtype)
  }
  tumor_subtype_upd_lvls <- c(paste0(rev(tumor_subtype_upd_lvls), '--1'), 
                              paste0(rev(tumor_subtype_upd_lvls), '--0'))
  dt[, tumor_subtype_upd := factor(tumor_subtype_upd, tumor_subtype_upd_lvls)]
  
  # tumour type tier palette
  tumorType_tier_palette <- create_tumortype_colorpalette(dt, colorPalette)
  # x axis labels & color
  xAxisLabs <- formatAxisLabels(dt, xOrderCol, xLabelCol, annoCol, 
                                highlightColor = annoColors)
  # transparency scales
  alphaVals <- c('yes' = 1, 'no' = 0.5)
  if (all(dt$transparent == 1)) {
    alphaVals <- c('no' = 0.5)
  }
  if (all(dt$transparent == 0)) {
    alphaVals <- c('yes' = 1)
  }
  
  # first of all, if there are several genomic regions, create shadowing
  result <- ggplot()
  Y_MAX <- NULL
  if (length(unique(dt$gr_id)) > 1) {
    # compute coordinates for future lines showing genome region 
    grShadeCoords <- get_gr_change_coords(plotDT = dt, labsDT = xAxisLabs$dt,
                                          axisOrderCol = xOrderCol, 
                                          yCol = 'nParts')
    grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                    nrow(grShadeCoords$gr_coords),
                                                    by = 2)]
    result <- result + 
      geom_rect(data = grShadeCoordsCut, fill = divider_color, inherit.aes = F,
                aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, 
                    ymax = Inf), alpha = 0.5) + 
      geom_segment(data = grShadeCoords$gr_coords, inherit.aes = F,
                   aes(x = start, xend = end, 
                       y =  grShadeCoords$lab_coords$y,
                       yend = grShadeCoords$lab_coords$y)) +
      geom_text(mapping = aes(x = x, y = 1.01 * y, label = gr, angle = 90, 
                              hjust = -0.01), inherit.aes = F,
                data = grShadeCoords$lab_coords, size = 2)
    Y_MAX <- 1.03 * unique(grShadeCoords$lab_coords$y)
  }
  # and if there are some filtered out drivers needed to be shown - shadowed 
  # areas for them
  if (any(dt$FILTER != "PASS")) {
    filteredOutShadeCoords <- get_filtered_out_gene_coords(xAxisLabs$dt, xOrderCol)
    filteredOutShadeCoordsCut <- filteredOutShadeCoords[seq(1, nrow(filteredOutShadeCoords), 
                                                            by = 2)]
    result <- result + 
      geom_rect(data = filteredOutShadeCoordsCut, fill = 'red', inherit.aes = F, 
                aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, 
                    ymax = Inf), alpha = 0.15)
  }
  
  result <- result + 
    geom_bar(data = dt, stat = 'identity',
             mapping = aes_string(x = xOrderCol, y = 'nParts', 
                                  alpha = 'transparent',
                                  fill = 'tumor_subtype_upd')) +
    xlab('genomic element') + ylab('N. participants') + 
    scale_fill_manual(values = tumorType_tier_palette, 
                      labels = gsub('--0', '', 
                                    grep('--0', names(tumorType_tier_palette),
                                         value = T)),
                      breaks = grep('--0', names(tumorType_tier_palette),
                                    value = T)) +
    scale_alpha_continuous(breaks = alphaVals, range = alphaVals, 
                           labels = names(alphaVals)) + 
    scale_x_discrete(labels = xAxisLabs$label, expand = c(0, 0)) +
    labs(fill = 'Histological\nsubtype', alpha = 'Found as driver')
  
  if (!is.null(ggplot2Theme)) {
    result <- result + ggplot2Theme
  }
  
  result <- result + 
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 6,
                                     color = xAxisLabs$color))
  if (!is.null(Y_MAX)) {
    result <- result + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX))
  } else {
    result <- result + scale_y_continuous(expand = c(0, 0))
  }
  
  result
}

#' tilePlotPvalues
#' @description Plots p-value matrix.
#' @author Maria Litovchenko
#' @param pValsDT data table
#' @param xLabelCol column name containing label to be on X axis
#' @param yLabelCol column name containing label to be on Y axis
#' @param colorBy name of column containing values to be used as color value of
#'                cells
#' @param xOrderCol column name containing order of X axis
#' @param yOrderCol column name containing order of Y axis
#' @param annoCol column containing annotation, i.e. COSMIC
#' @param annoColors named vector of colors (in hex) in which X axis labels 
#'                   will be colored
#' @param divider_color hex color which should be used for separating different
#'        genomic regions
#' @param doLog10 whatever -log10 operation on colorBy column should be 
#'                performed
#' @param scaleLims 
#' @param nColors number of colors to use in the color palette 
#' @param colorPalette character vector of colors
#' @param colorPaletteName name of color legend
#' @param ggplot2Theme ggplot2 theme
tilePlotPvalues <- function(pValsDT, xLabelCol, yLabelCol, colorBy,
                            xOrderCol = NULL, yOrderCol = NULL, 
                            annoCol = NULL, annoColors = 'red', 
                            doLog10 = T, scaleLims = NULL, 
                            nColors = 5, colorPalette = NULL, 
                            divider_color = '#eed9c4',
                            ggplot2Theme = NULL) {
  if (!is.null(colorPalette) & nColors != length(colorPalette)) {
    stop('[', Sys.time(), '] Length of nColors and colorPalette should match')
  }
  dt <- copy(pValsDT)
  if (is.null(xOrderCol)) {
    xOrderCol <- xLabelCol
  }
  if (is.null(yOrderCol)) {
    yOrderCol <- yLabelCol
  }
  dt[, colorByLog10 := unlist(dt[, colorBy, with = F])]
  if (doLog10 == T) {
    dt[, colorByLog10 := -log10(colorByLog10 + 10^(-16))]
  }
  
  # p value palette
  pval_palette <- create_pvalue_colorpalette(dt, doLog10, scaleLims, nColors,
                                             colorPalette)
  # x axis labels & color
  xAxisLabs <- formatAxisLabels(dt, xOrderCol, xLabelCol, annoCol, 
                                highlightColor = annoColors)
  
  # first of all, if there are several genomic regons, create shadowing
  result <- ggplot()
  if (length(unique(dt$gr_id)) > 1) {
    # compute coordinates for future lines showing genome region 
    grShadeCoords <- get_gr_change_coords(plotDT = dt, labsDT = xAxisLabs$dt, 
                                          axisOrderCol = xOrderCol,
                                          yCol = yOrderCol)
    grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                    nrow(grShadeCoords$gr_coords),
                                                    by = 2)]
    result <- result + 
      geom_rect(data = grShadeCoordsCut, fill = divider_color, inherit.aes = F,
                aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, 
                    ymax = Inf), alpha = 0.5)
  }
  # and if there are some filtered out drivers needed to be shown - shadowed 
  # areas for them
  if (any(dt$FILTER != "PASS")) {
    filteredOutShadeCoords <- get_filtered_out_gene_coords(xAxisLabs$dt, xOrderCol)
    filteredOutShadeCoordsCut <- filteredOutShadeCoords[seq(1, 
                                                            nrow(filteredOutShadeCoords),
                                                            by = 2)]
    result <- result + 
      geom_rect(data = filteredOutShadeCoordsCut, fill = 'red', 
                inherit.aes = F, alpha = 0.15,
                aes(xmin = start - 0.5, xmax = end + 0.5, 
                    ymin = -Inf, ymax = Inf))
  } 
  
  result <- result + 
    geom_tile(data = dt,
              mapping = aes_string(x = xOrderCol, y = yOrderCol,
                                   fill = 'colorByLog10')) + 
    scale_x_discrete(labels = xAxisLabs$label)
  if (class(unlist(dt[, c('colorByLog10'), with = F])) %in% 
      c('numeric', 'integer', 'double')) {
    result <- result + scale_fill_stepsn(na.value = '#FFFFFF',
                                         breaks = pval_palette$colorBreaks, 
                                         colours = pval_palette$palette)
  } else {
    result <- result + scale_fill_manual(na.value = '#FFFFFF', 
                                         values = pval_palette$palette)
  }
  
  if (!is.null(ggplot2Theme)) {
    result <- result + ggplot2Theme
  }
  result <- result + xlab('genomic element') + ylab('') +
    scale_y_discrete(expand = c(0, 0)) +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, size = 6, hjust = 1, 
                                     color = xAxisLabs$color))
  result
}

# Visuals setup ---------------------------------------------------------------
font_config <- list('family' = 'Helvetica', 'base' = 6, 'text' = 6, 
                    'axis_title' = 8, 'plot_title' = 10, 'legend_text' = 6, 
                    'legend_title' = 8)

customGgplot2Theme <- list(
  theme_classic(base_size = font_config$base, 
                base_family = font_config$family) +
    theme(axis.line = element_line(colour = 'black', linewidth = 0.5,
                                   linetype = 'solid'),
          axis.ticks = element_line(colour = 'black', linewidth = 0.5,
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
                                          linewidth = 0.5, linetype = 2),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm")))

DIVIDER_COLOR <-  '#eed9c4'

# color blind friendly
tumourTypeColorPalette <- c('Adenocarcinoma' = '#CC3311',
                            'Adenocarcinoma_met' = '#CC6677',
                            'Adenosquamous' = '#332288',
                            'Carcinoid' = '#882255',
                            'Large_cell' = '#DDCC77',
                            'Mesothelioma' = '#999933',
                            'Neuroendocrine_carcinoma' = '#EE8866',
                            'Small_cell' = '#44AA99', 
                            'Small_cell_met' = '#bbe4dd', 
                            'Squamous_cell' = '#0077BB',
                            'Squamous_cell_met' = '#88CCEE',
                            'Other' = '#000000',
                            'Other_met' = '#CCCCCC') 
