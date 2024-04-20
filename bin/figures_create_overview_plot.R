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

#' create_plot_background
#' @description Creates "background" for the future plots of coding and 
#'              noncoding driver genetic elements. It creates shades to 
#'              distinguish between different types of genetic elements (i.e.
#'              CDS and promoters) as well as creates a separate shade for the
#'              filtered out genetic elements.
#' @author Maria Litovchenko
#' @param plotDT data table with processed data, ready to be plotted (main 
#'        plot)
#' @param axisLabs result of function format_xLabels
#' @param xOrderCol column name containing order of X axis
#' @param yCol column name containing data for Y axis
#' @param divider_color hex color which should be used for separating different
#'        genomic regions
#' @param keep_region_labels boolean, whatever or not labels of regions on the
#'        top of the plot should be kept.
#' @param direction string, one of "vertical" and "horizontal"
#' @return ggplot2 object
create_plot_background <- function(plotDT, axisLabs, xOrderCol, yCol,
                                   divider_color = '#eed9c4', 
                                   direction = 'vertical') {
  if (!direction %in% c('vertical', 'horizontal')) {
    stop('[', Sys.time(), '] direction should be vertical or horizontal')
  }
  
  # first of all, if there are several genomic regions, create shadowing
  result <- ggplot()
  if (length(unique(plotDT$gr_id)) > 1) {
    # compute coordinates for future lines showing genome region 
    grShadeCoords <- get_gr_change_coords(plotDT, xOrderCol, yCol, axisLabs$dt)
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
    filteredOutShadeCoords <- get_filtered_out_gene_coords(axisLabs$dt, xOrderCol)
    filteredOutShadeCoordsCut <- filteredOutShadeCoords[seq(1,
                                                            nrow(filteredOutShadeCoords),
                                                            by = 2)]
    if (direction == 'vertical') {
      result <- result + 
        geom_rect(data = filteredOutShadeCoordsCut, fill = 'red', inherit.aes = F, 
                  mapping = aes(xmin = start - 0.5, xmax = end + 0.5,
                                ymin = -Inf, ymax = Inf), alpha = 0.15)
    } else {
      result <- result + 
        geom_rect(data = filteredOutShadeCoordsCut, fill = 'red', inherit.aes = F, 
                  mapping = aes(ymin = start - 0.5, ymax = end + 0.5, 
                                xmin = -Inf, xmax = Inf), alpha = 0.15)
    }
  }
  
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

# Functions : Overview of all detected driver ---------------------------------
#' get_bar_transparency
#' @description 
#' @author Maria Litovchenko
#' @param plotDT data table which holds data for plotting. Have to have columns
#'        tumor_subtype, tier, participant_tumor_subtype and a column with the 
#'        name identical to xOrderCol
#' @param pancanCode string encoding tumor subtype code for pancancer analysis
#' @param xOrderCol name of the column which contains IDs of genes to plot and
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
get_bar_transparency <- function(plotDT, pancanCode, xOrderCol) {
  result <- copy(plotDT)
  # transparent, if significant in individual tumor subtypes
  result[, transparent := ifelse(tier == 0, 1, 0)]
  # all tumor_subtypes opaque, if significant in pancan analysis
  signPancan <- result[tumor_subtype == pancanCode & tier == 1]
  signPancan <- unlist(signPancan[, xOrderCol, with = F])
  signPancan <- unlist(result[, xOrderCol, with = F]) %in% signPancan
  result[signPancan]$transparent <- 0
  # no need to keep pancan in the table anymore
  result <- split(result, by = xOrderCol, drop = T)
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
  xAxisLabs <- format_xLabels(dt, xOrderCol, xLabelCol, annoCol, 
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
    grShadeCoords <- get_gr_change_coords(plotDT = dt, xLabsDT = xAxisLabs$dt,
                                          xOrderCol = xOrderCol, 
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
  xAxisLabs <- format_xLabels(dt, xOrderCol, xLabelCol, annoCol, 
                              highlightColor = annoColors)
  
  # first of all, if there are several genomic regons, create shadowing
  result <- ggplot()
  if (length(unique(dt$gr_id)) > 1) {
    # compute coordinates for future lines showing genome region 
    grShadeCoords <- get_gr_change_coords(plotDT = dt, xLabsDT = xAxisLabs$dt, 
                                          xOrderCol = xOrderCol,
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
    geom_tile(data = dt, mapping = aes_string(x = xOrderCol, y = yOrderCol,
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

# Test input arguments -------------------------------------------------------------
args <- list(cancer_subtype = 'Panlung',
             drivers = "completed_runs/2023_12_25/results/tables/drivers/drivers-Panlung--hg19.csv",
             inventory_patients = 'data/inventory/inventory_patients.csv',
             excluded_patients = "", 
             allowed_filter_values = list("PASS", "INDEL, 2-5bp"),
             drivers_uniform_subtypes = list("completed_runs/2023_12_25/results/tables/drivers/drivers-Adenocarcinoma--hg19.csv",
                                             "completed_runs/2023_12_25/results/tables/drivers/drivers-Adenocarcinoma_met--hg19.csv",
                                             "completed_runs/2023_12_25/results/tables/drivers/drivers-Adenosquamous--hg19.csv",
                                             "completed_runs/2023_12_25/results/tables/drivers/drivers-Carcinoid--hg19.csv",
                                             "completed_runs/2023_12_25/results/tables/drivers/drivers-Large_cell--hg19.csv",
                                             "completed_runs/2023_12_25/results/tables/drivers/drivers-Mesothelioma--hg19.csv",
                                             "completed_runs/2023_12_25/results/tables/drivers/drivers-Neuroendocrine_carcinoma--hg19.csv",
                                             "completed_runs/2023_12_25/results/tables/drivers/drivers-Small_cell--hg19.csv",
                                             "completed_runs/2023_12_25/results/tables/drivers/drivers-Squamous_cell--hg19.csv",
                                             "completed_runs/2023_12_25/results/tables/drivers/drivers-Squamous_cell_met--hg19.csv"),
             extra_studies = list("data/assets/intogene_detectedCancerGenes.csv",
                                  "data/assets/mc3_detectedCancerGenes.csv",
                                  "data/assets/cgc_knownCancerGenes.csv"),
             extra_studies_names = list("intogen", "mc3", "CGC"),
             extra_studies_tumorsubtype = list("LNET,LUAD,LUSC,NSCLC,SCLC", 
                                               "LUAD,LUSC", "nsclc,sclc,lung"))

if (is.null(args$drivers) & is.null(args$drivers_uniform_subtypes)){
  stop()
}
if (!is.null(args$extra_studies) & is.null(args$extra_studies_names)) {
  stop()
}
if (length(args$extra_studies) != length(args$extra_studies_names)) {
  stop() 
}
if (!is.null(args$extra_studies_tumorsubtype)) {
  if (length(args$extra_studies) != length(args$extra_studies_tumorsubtype)) {
    stop()
  }
  args$extra_studies_tumorsubtype <- lapply(args$extra_studies_tumorsubtype,
                                            function(x) unlist(strsplit(x, 
                                                                        ',')))
  args$extra_studies_tumorsubtype <- lapply(args$extra_studies_tumorsubtype,
                                            function(x) gsub(' ', '', x))
  args$extra_studies_tumorsubtype <- lapply(args$extra_studies_tumorsubtype,
                                            paste0, collapse = '|')
}

# Read in patients inventory --------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients, 1)
patientsInv <- patientsInv[tumor_subtype %in% args$cancer_subtype]
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)

# Read unfiltered drivers -----------------------------------------------------
if (!is.null(args$drivers)) {
  drivers_unfiltered <- fread(args$drivers, header = T, stringsAsFactors = F)
  drivers_unfiltered[, gene_name := gsub('__and__', '&', gene_name)]
  drivers_unfiltered[, gene_id := gsub('__and__', '&', gene_id)]
  message('[', Sys.time(), '] Read --drivers: ', args$drivers)
}

# Read drivers detected in uniform subtypes -----------------------------------
if (!is.null(args$drivers_uniform_subtypes)) {
  drivers_uniform_subtypes <- lapply(args$drivers_uniform_subtypes, fread, 
                                     header = T, stringsAsFactors = F)
  message('[', Sys.time(), '] Read --drivers_uniform_subtypes: ',
          paste0(args$drivers_uniform_subtypes, collapse = ', '))
  drivers_uniform_subtypes <- do.call(rbind.fill, drivers_uniform_subtypes)
  drivers_uniform_subtypes <- as.data.table(drivers_uniform_subtypes)
  drivers_uniform_subtypes[, gene_name := gsub('__and__', '&', gene_name)]
  drivers_uniform_subtypes[, gene_id := gsub('__and__', '&', gene_id)]
  
  if (sum(drivers_uniform_subtypes$tumor_subtype != 
          drivers_uniform_subtypes$participant_tumor_subtype, 
          na.rm = T) > 1) {
    stop()
  }
  
  uniform_subtypes <- unique(drivers_uniform_subtypes$tumor_subtype)
  part_of_composite <- unique(drivers_unfiltered$participant_tumor_subtype)
  not_in_composite <- setdiff(uniform_subtypes, part_of_composite)
  if (length(not_in_composite) > 0) {
      stop()
  }
}

# Read extra databases --------------------------------------------------------
if (!is.null(args$extra_studies)) {
  extra_studies <- lapply(args$extra_studies, fread, header = T,
                          stringsAsFactors = F, 
                          select = c('gene_name', 'gene_id', 
                                     'known_in_tumor_subtype'))
  names(extra_studies) <- unlist(args$extra_studies_names)
  for (i in 1:length(extra_studies)) {
    extra_studies[[i]][, known_in_current_tumor_subtype := grepl(tolower(args$extra_studies_tumorsubtype[i]),
                                                                 tolower(known_in_tumor_subtype))]
  }
}

# Select genes - drivers in args$tumor_subtype or in uniform subtypes ---------
drivers <- drivers_unfiltered[!is.na(tier) &
                                FILTER %in% args$allowed_filter_values]
drivers <- drivers[,.(gr_id, gene_id, gene_name)]
drivers <- unique(drivers)

if (!is.null(args$drivers_uniform_subtypes)) {
  drivers_uniform_subtypes_ids <- drivers_uniform_subtypes[!is.na(tier) & 
                                                             FILTER %in%
                                                             args$allowed_filter_values]
  drivers_uniform_subtypes_ids <- drivers_uniform_subtypes_ids[,.(tumor_subtype,
                                                                  gr_id, 
                                                                  gene_id,
                                                                  gene_name)]
  drivers_uniform_subtypes_ids <- unique(drivers_uniform_subtypes_ids)
  drivers_uniform_subtypes <- merge(drivers_uniform_subtypes,
                                    drivers_uniform_subtypes_ids,
                                    by = c("tumor_subtype", "gr_id",
                                           "gene_id", "gene_name"))
  
  drivers_uniform_subtypes_ids <- drivers_uniform_subtypes_ids[,.(gr_id, 
                                                                  gene_id,
                                                                  gene_name)]
  drivers_uniform_subtypes_ids <- unique(drivers_uniform_subtypes_ids)
  drivers <- rbind(drivers, drivers_uniform_subtypes_ids)
  drivers <- unique(drivers)
}

drivers <- merge(drivers_unfiltered, drivers, 
                 by = c("gr_id", "gene_id", "gene_name"))
if (!is.null(args$drivers_uniform_subtypes)) {
  drivers <- as.data.table(rbind.fill(drivers, drivers_uniform_subtypes))
}

if (nrow(drivers) == 0) {
  stop()
}

# Calculate cohort sizes ------------------------------------------------------
cohort_sizes <- patientsInv[,.(length(unique(participant_id))), 
                            by = participant_tumor_subtype]
setnames(cohort_sizes, c("participant_tumor_subtype", 'V1'), 
         c("tumor_subtype", 'cohort_size'))
cohort_sizes <- cohort_sizes[order(cohort_size, decreasing = T)]

# Set up levels(order) for genomic regions and tumor types --------------------
# sort tumor types by number of patients
TUMOR_SUBTYPE_ORDER <- cohort_sizes$tumor_subtype

# and genomic regions by maximum number of mutated patients
GR_ORDER <- drivers[FILTER == "PASS"][order(-nParts_total)]
GR_ORDER <- GR_ORDER[,.SD[1], by = gr_id]$gr_id

# in case plotting of 2-5bp filtered out drivers is requested too
GR_FILTERED_OUT_ORDER <- NULL
if (nrow(drivers[FILTER != "PASS"]) > 0) {
  GR_FILTERED_OUT_ORDER <- drivers[FILTER != "PASS"][order(-nParts_total)]
  GR_FILTERED_OUT_ORDER <- GR_FILTERED_OUT_ORDER[,.SD[1], by = gr_id]$gr_id
}

# Create gene_plots_id to uniquely ID gene + gr combo & order on plots --------
GENE_ORDER <- drivers[FILTER == "PASS"]
GENE_ORDER <- GENE_ORDER[,.(gr_id, gene_id, gene_name, nParts_total)]
GENE_ORDER <- unique(GENE_ORDER)
GENE_ORDER[, gr_id := factor(gr_id, GR_ORDER)]
GENE_ORDER <- GENE_ORDER[order(gr_id, -nParts_total, -gene_name)]
GENE_ORDER <- GENE_ORDER[,.(gr_id, gene_id, gene_name)]
# to handle multiple tumour subtypes
GENE_ORDER <- unique(GENE_ORDER)

# in case plotting of 2-5bp filtered out drivers is requested too
GENE_FILTERED_OUT_ORDER <- NULL
if (nrow(drivers[FILTER != "PASS"]) > 0) {
  GENE_FILTERED_OUT_ORDER <- drivers[FILTER != "PASS"]
  GENE_FILTERED_OUT_ORDER <- GENE_FILTERED_OUT_ORDER[,.(gr_id, gene_id, 
                                                        gene_name,
                                                        nParts_total)]
  GENE_FILTERED_OUT_ORDER <- unique(GENE_FILTERED_OUT_ORDER)
  GENE_FILTERED_OUT_ORDER[, gr_id := factor(gr_id, GR_FILTERED_OUT_ORDER)]
  GENE_FILTERED_OUT_ORDER <- GENE_FILTERED_OUT_ORDER[order(gr_id,
                                                           -nParts_total,
                                                           -gene_name)]
  GENE_FILTERED_OUT_ORDER <- GENE_FILTERED_OUT_ORDER[,.(gr_id, gene_id, 
                                                        gene_name)]
  # to handle multiple tumour subtypes
  GENE_FILTERED_OUT_ORDER <- unique(GENE_FILTERED_OUT_ORDER)
}

GENE_ORDER <- rbind(GENE_ORDER, GENE_FILTERED_OUT_ORDER)
GENE_ORDER[, gene_plots_id := paste(gene_name, gr_id)] 
GENE_ORDER <- GENE_ORDER$gene_plots_id

drivers[, gene_plots_id := paste(gene_name, gr_id)] 
drivers[, gene_plots_id := factor(gene_plots_id, GENE_ORDER)]
drivers[, participant_tumor_subtype := factor(participant_tumor_subtype, 
                                              TUMOR_SUBTYPE_ORDER)]
drivers[, tumor_subtype := factor(tumor_subtype, 
                                  unique(c(args$cancer_subtype, 
                                           TUMOR_SUBTYPE_ORDER)))]

# Overlap data from extra databases with detected drivers ---------------------
drivers_restricted <- drivers[!is.na(tier) & 
                                FILTER %in% args$allowed_filter_values]
drivers_restricted <- drivers_restricted[,.(gr_id, gene_name, gene_id, 
                                            gene_plots_id, FILTER)]
drivers_restricted <- unique(drivers_restricted)
if (max(drivers_restricted[,.N, by = gene_plots_id]$N) > 1) {
  stop()
}

for (i in 1:length(extra_studies)) {
  keys_i <- intersect(c("gene_name", "gene_id"), colnames(extra_studies[[i]]))
  if (!'is_known_cancer' %in% colnames(extra_studies[[i]])) {
    extra_studies[[i]][, is_known_cancer := T]
  }
  
  extra_studies[[i]] <- merge(drivers_restricted, extra_studies[[i]], 
                              by = keys_i, all.x = T)
  extra_studies[[i]][is.na(is_known_cancer)]$is_known_cancer <- F
  extra_studies[[i]][is.na(known_in_current_tumor_subtype)]$known_in_current_tumor_subtype <- F
  extra_studies[[i]] <- extra_studies[[i]][,.(gr_id, gene_name, gene_id, 
                                              gene_plots_id, FILTER, 
                                              is_known_cancer,
                                              known_in_current_tumor_subtype)]
  extra_studies[[i]] <- unique(extra_studies[[i]])
  extra_studies[[i]] <- extra_studies[[i]][order(gr_id, gene_name, gene_id, 
                                                 gene_plots_id, FILTER,
                                                 -is_known_cancer,
                                                 -known_in_current_tumor_subtype)]
  extra_studies[[i]] <- extra_studies[[i]][,.SD[1], 
                                           by = .(gr_id, gene_name, gene_id,
                                                  gene_plots_id, FILTER)]
  extra_studies[[i]][, dummy := names(extra_studies)[i]]
}

# Barplot number of patients -------------------------------------------------
# divider color
barNpats <- barplotNpatient(driversDT = drivers, 
                            xLabelCol = 'gene_name',
                            xOrderCol = 'gene_plots_id', 
                            annoCol = 'is_known_cancer', 
                            colorPalette = tumourTypeColorPalette, 
                            pancanCode = args$cancer_subtype, 
                            divider_color = DIVIDER_COLOR,
                            ggplot2Theme = customGgplot2Theme)

barNpats <- barNpats + theme(legend.position = 'bottom')
barNpats <- barNpats + guides(fill = guide_legend(nrow = 2, byrow = T)) +
    theme(legend.box = "horizontal")

# Plot annotation with extra databases ----------------------------------------
if (!is.null(args$extra_studies)) {
  extraDbPlots <- lapply(1:length(extra_studies),
                         function(x) tilePlotPvalues(extra_studies[[x]],
                                                     xLabelCol = 'gene_name',
                                                     yLabelCol = 'dummy', 
                                                     colorBy = 'is_known_cancer', 
                                                     xOrderCol = 'gene_plots_id',
                                                     doLog10 = F,
                                                     annoCol = 'is_known_cancer',
                                                     nColors = 2, 
                                                     divider_color = DIVIDER_COLOR,
                                                     colorPalette = c('#FFFFFF', '#8e918f'),
                                                     ggplot2Theme = customGgplot2Theme) +
                           geom_text(data = extra_studies[[x]][known_in_current_tumor_subtype == T], 
                                     mapping = aes(x = gene_plots_id, y = dummy, 
                                                   label = '*'),
                                     vjust = 0.75) + 
                           theme(legend.position = 'none',
                                 axis.title = element_blank(), 
                                 axis.ticks.y = element_blank()))
}

# Plot tile plot of p-values --------------------------------------------------
pValMatrTiles <- tilePlotPvalues(pValsDT = drivers[!is.na(tier) & 
                                                     FILTER %in% args$allowed_filter_values],
                                 xLabelCol = 'gene_name', 
                                 yLabelCol = 'tumor_subtype', 
                                 colorBy = 'brown.bh_p', 
                                 xOrderCol = 'gene_plots_id',
                                 annoCol = 'is_known_cancer',
                                 doLog10 = T,
                                 ggplot2Theme = customGgplot2Theme)
pValMatrTiles <- pValMatrTiles + labs(fill = '-log10(Brown p.value)') +
  theme(legend.direction = 'horizontal')

# Assemble --------------------------------------------------------------------
extraTheme <- list(theme(legend.position = 'none',
                         axis.title.x = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.line.x = element_blank(),
                         panel.grid.major.x = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         panel.grid.minor.y = element_blank())) 
ovrvFigPlotList <- list(barNpats + extraTheme)
if (!is.null(args$extra_studies)) {
  extraDbPlots <- lapply(extraDbPlots, function(x) x + extraTheme)
  ovrvFigPlotList <- append(ovrvFigPlotList, extraDbPlots)
}
ovrvFigPlotList[[length(ovrvFigPlotList) + 1]] <- pValMatrTiles + 
  theme(legend.position = 'none')

# main plot
plotHeights <- c(1)
if (!is.null(args$extra_studies)) {
  plotHeights <- c(plotHeights, rep(0.025, length(extra_studies)))
}
plotHeights <- c(plotHeights, 0.32)
DRIVERS_OVERVIEW_PLOT <- ggarrange(plotlist = ovrvFigPlotList, align = 'v',
                                   ncol = 1, heights = plotHeights)
DRIVERS_OVERVIEW_LEGENDS <- ggarrange(extract_legend(barNpats), 
                                      extract_legend(pValMatrTiles),
                                      ncol = 1)
DRIVERS_OVERVIEW_PLOT <- ggarrange(DRIVERS_OVERVIEW_PLOT, 
                                   DRIVERS_OVERVIEW_LEGENDS, ncol = 1,
                                   heights = c(1, 0.15))

