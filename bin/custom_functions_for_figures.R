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
suppressWarnings(suppressPackageStartupMessages(library(rjson)))

# Function : setting up visual appearance of plots ----------------------------
#' readJsonWithVisualParameters
#' @description Reads JSON file determining visual appearance (i.e. colors & 
#' ggplot2 theme) of plots.
#' @author Maria Litovchenko
#' @param jsonPath string, path to JSON
#' @return named list with names ggplot2_theme, color_divider, 
#' colors_tumor_type, colors_biotype, colors_subtype_specificity, 
#' colors_perc_tumors_with_driver_mut, colors_divergent_palette
readJsonWithVisualParameters <- function(jsonPath) {
  json <- fromJSON(file = jsonPath)
  
  # customGgplot2Theme
  ggplot2_theme <- list(
    theme_classic(base_size = json$font_size_base, 
                  base_family = json$font_family) +
      theme(plot.title = element_text(colour = "black", 
                                      size = json$font_size_plot_title),
            axis.title = element_text(colour = "black", 
                                      size = json$font_size_axis_title),
            axis.text = element_text(colour = "black", 
                                     size =  json$font_size_text),
            axis.line = element_line(colour = "black", linewidth = 0.5,
                                     linetype = "solid"),
            axis.ticks = element_line(colour = "black", linewidth = 0.5,
                                      linetype = "solid"),
            legend.title = element_text(colour = "black", 
                                        size = json$font_size_legend_title),
            legend.text = element_text(colour = "black", 
                                       size = json$font_size_legend_text),
            legend.background = element_blank(),
            legend.direction = "vertical",
            legend.box.background = element_blank(),
            legend.key.size = unit(0.25, "cm"),
            panel.grid.major = element_line(colour = "#CCCCCC", 
                                            linewidth = 0.5, linetype = 2),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "cm")))
  
  color_divider <- json$color_divider
  
  # tumourTypeColorPalette
  colors_tumor_type <- unlist(json$colors_tumor_type)
  
  # biotypePalette
  colors_biotype <- unlist(json$colors_biotype)
  
  # subtypeSpecificityPalette
  colors_subtype_specificity <- unlist(json$colors_subtype_specificity)
  
  # percTumorsWithDriverMutsPalette
  colors_perc_tumors_with_driver_mut <- unlist(json$colors_perc_tumor_with_driver_mutation)
  
  # gradientColors
  colors_divergent_palette <- unlist(json$colors_divergent_palette)
  
  result <- list('ggplot2_theme' = ggplot2_theme, 
                 'color_divider' = color_divider,
                 'colors_tumor_type' = colors_tumor_type, 
                 'colors_biotype' = colors_biotype, 
                 'colors_subtype_specificity' = colors_subtype_specificity, 
                 'colors_perc_tumors_with_driver_mut' = colors_perc_tumors_with_driver_mut,
                 'colors_divergent_palette' = colors_divergent_palette)
  result
}

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

# Functions : coloring of axis labels -----------------------------------------
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

#  Function: folding splice sites into coding regions  ------------------------
#' foldSplicSiteDriversIntoCodingDrivers
#' @description Replaces gr_id of splice site driver genomic elements with 
#' gr_id corresponding to coding driver genomic elements if a genomic element
#' was detected as a driver both in a splice site region and in a coding region
#' @author Maria Litovchenko
#' @param codingGrIds vector of characters. gr_ids considered as coding genomic
#' elements
#' @param ssGrIds vector of characters. gr_ids considered as splice site 
#' genomic elements
#' @param driversDT data table with driver genomic elements. Must have columns:
#' gene_id, gene_name, gr_id
#' @return data table, updated driversDT
foldSplicSiteDriversIntoCodingDrivers <- function(codingGrIds, ssGrIds, 
                                                  driversDT) {
  if (length(codingGrIds) == 0 | length(ssGrIds) == 0) {
    message('[', Sys.time(), '] Either no coding or no splice site genomic ',
            'regions were detected. Will not perform folding of splice sites ',
            'into corresponding coding driver genetic elements.')
    return(driversDT)
  }
  
  # GE = genomic element
  driverGEs <- unique(driversDT[,.(gene_id, gene_name, gr_id)])
  driverGEs[, nCodingGrId := sum(unique(gr_id) %in% codingGrIds),
            by = .(gene_id, gene_name)]
  driverGEs[, nSsGrId := sum(unique(gr_id) %in% ssGrIds), 
            by = .(gene_id, gene_name)]
  
  driversToFold <- driverGEs[nCodingGrId >= 1 & nSsGrId >= 1]
  if (nrow(driversToFold) > 0) {
    message('[', Sys.time(), '] Will fold splice sites drivers detected in ',
            paste(unique(driversToFold$gene_name), collapse = ', '),
            ' into corresponding coding drivers.')
    
    if (any(driversToFold$nCodingGrId) > 1) {
      driversToFold <- driversToFold[nCodingGrId > 1]
      stop('[', Sys.time(), '] Found genes for which > 1 gr_id can be ',
           'considered as coding:\n', 
           paste0(colnames(driversToFold), collapse = '\t'), '\n', 
           paste0(apply(driversToFold, 1, paste0, collapse = '\t'), 
                  collapse = '\n'), '. Can not process such  situation.')
    }
    
    driversToFold <- driversToFold[,.(gene_name, gene_id, gr_id)]
    driversToFold <- split(driversToFold, by = c('gene_name', 'gene_id'),
                           drop = T)
    for (geneToFold in driversToFold) {
      geneToFoldIdx <- which(driversDT$gene_id %in% geneToFold$gene_id & 
                               driversDT$gene_name %in% geneToFold$gene_name &
                               driversDT$gr_id %in% ssGrIds)
      geneToFoldCodingGr <- geneToFold[gr_id %in% codingGrIds]$gr_id
      driversDT[geneToFoldIdx]$gr_id <- geneToFoldCodingGr
    }
  }
  
  driversDT
}
