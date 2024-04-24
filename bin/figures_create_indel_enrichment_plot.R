library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)

source('bin/custom_functions.R')


# Functions : genomic regions coloring & annotation ---------------------------
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

#' annotatePlotWithGRlines
#' @description Annotates a ggplot with a lines and text deliniating groups of
#'              the same genomic regions.
#' @author Maria Litovchenko
#' @param basePlot a ggplot2 object to which annotation should be added
#' @param plotDT data table with processed data, ready to be plotted (main 
#'        plot)
#' @param axisLabs result of function format_xLabels
#' @param xOrderCol column name containing order of X axis
#' @param yCol column name containing data for Y axis
#' @param direction string, one of "vertical" and "horizontal"
#' @param ... other parameters to pass to geom_text
#' @return ggplot2 object
annotatePlotWithGRlines <- function(basePlot, plotDT, axisLabs, xOrderCol,
                                    yCol, direction = 'vertical', ...) {
  result <- basePlot
  
  if (!direction %in% c('vertical', 'horizontal')) {
    stop('[', Sys.time(), '] direction should be vertical or horizontal')
  }
  
  if (length(unique(plotDT$gr_id)) > 1) {
    # compute coordinates for future lines showing genome region 
    grShadeCoords <- get_gr_change_coords(plotDT, xOrderCol, yCol, axisLabs$dt)
    nGRs <- nrow(grShadeCoords$gr_coords)
    grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, nGRs, by = 2)]
    
    if (direction == 'vertical') {
      if (is.factor(grShadeCoords$lab_coords$y)) {
        result <- basePlot + 
          geom_segment(data = grShadeCoords$gr_coords,  
                       mapping = aes(x = start, xend = end, 
                                     y = grShadeCoords$lab_coords$y,
                                     yend = grShadeCoords$lab_coords$y)) +
          geom_text(mapping = aes(x = x, y = y, label = gr),
                    data = grShadeCoords$lab_coords, ...) +
          scale_y_discrete(drop = F)
      } else {
        result <- basePlot + 
          geom_segment(data = grShadeCoords$gr_coords, inherit.aes = F,
                       mapping = aes(x = start, xend = end, 
                                     y = grShadeCoords$lab_coords$y,
                                     yend = grShadeCoords$lab_coords$y)) +
          geom_text(mapping = aes(x = x, y = y, label = gr),
                    inherit.aes = F, data = grShadeCoords$lab_coords, ...)
      }
      
    } else {
      if (is.factor(grShadeCoords$lab_coords$y)) {
        result <- basePlot + 
          geom_segment(data = grShadeCoords$gr_coords, 
                       mapping = aes(y = start, yend = end, 
                                     x = grShadeCoords$lab_coords$y,
                                     xend = grShadeCoords$lab_coords$y)) +
          geom_text(mapping = aes(x = y, y = x, label = gr),
                    data = grShadeCoords$lab_coords, ...) +
          scale_x_discrete(drop = F)
      } else {
        result <- basePlot + 
          geom_segment(data = grShadeCoords$gr_coords,
                       mapping = aes(y = start, yend = end, 
                                     x = grShadeCoords$lab_coords$y,
                                     xend = grShadeCoords$lab_coords$y)) +
          geom_text(mapping = aes(x = y, y = x, label = gr),
                    data = grShadeCoords$lab_coords, ...)
      }
    }
  }
  
  result
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
                                          linewidth = 0.5, linetype = 2),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm")))

DIVIDER_COLOR <-  '#eed9c4'

tealPalette <- c('#125A56','#238F9D', '#60BCE9', '#C6DBED')

# Test input arguments --------------------------------------------------------
args <- list(cancer_subtype = 'Panlung',
             drivers = "completed_runs/2023_12_25/results/tables/drivers/drivers-Panlung--hg19.csv",
             allowed_filter_values = list("PASS", "INDEL, 2-5bp"),
             variant_categories_counts = "completed_runs/2023_12_25/results/mut_rates/varCatEnrich-Panlung--hg19.csv",
             categoty_of_interest = "INDEL, 2-5bp", p_abj = 0.05)

VAR_CAT_LEVELS <- rev(c("SNP", "INDEL, 1bp", "INDEL, 2-5bp", "INDEL, 6-9bp",
                        "INDEL, 10+bp", "MNP"))

# Read unfiltered drivers -----------------------------------------------------
drivers <- fread(args$drivers, header = T, stringsAsFactors = F)
drivers[, gene_name := gsub('__and__', '&', gene_name)]
drivers[, gene_id := gsub('__and__', '&', gene_id)]
message('[', Sys.time(), '] Read --drivers: ', args$drivers)

drivers <- drivers[!is.na(tier) & FILTER %in% args$allowed_filter_values]
drivers <- drivers[,.(gr_id, gene_id, gene_name, FILTER, nParts_total,
                      is_known_cancer)]
drivers <- unique(drivers)

if (nrow(drivers) == 0) {
  message('[', Sys.time(), '] No driver genomic regions passed requested ',
          'filters. Plotting will not be performed.')
  stop_quietly()
}

# Set up levels(order) for genomic regions ------------------------------------
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

# Read variant category enrichment -------------------------------------------
varCatEnrich <- fread(args$variant_categories_counts, header = T, 
                      stringsAsFactors = F, 
                      select = c("gr_id", "gene_id", "gene_name", "var_cat",
                                 "Nvars", "N_gene", "binomPadj"))
message('[', Sys.time(), '] Read --variant_categories_counts: ', 
        args$variant_categories_counts)
varCatEnrich[, gene_name := gsub('__and__', '&', gene_name)]
varCatEnrich[, gene_id := gsub('__and__', '&', gene_id)]

varCatEnrich <- varCatEnrich[gr_id %in% unique(drivers$gr_id)]
varCatEnrich <- merge(varCatEnrich, drivers, all.y = T,
                      by = c("gr_id", "gene_id", "gene_name"))
varCatEnrich[, perc_cat := 100 * Nvars/N_gene]
varCatEnrich[, Nvars := NULL]
varCatEnrich[, N_gene := NULL]

# Set up order of variant categories ------------------------------------------
varCatEnrich[, var_cat := factor(var_cat, VAR_CAT_LEVELS)]

# Create plot background (rectangles to split genomic regions) ----------------
# color of plot labels
xAxisLabs <- format_xLabels(varCatEnrich, 'gene_plots_id', 'gene_name',
                            'is_known_cancer')
xAxisLabs$labels <- gsub("__and__", "&", xAxisLabs$labels)
# background plot to split different genomic regions from each other
plotBackground <- create_plot_background(varCatEnrich, xAxisLabs, 
                                         'gene_plots_id', 'perc_cat', 
                                         direction = 'vertical')

# Plot bar plot enrichment of variant categories  -----------------------------
categoriesPallete <- colorRampPalette(tealPalette)
categoriesPallete <- categoriesPallete(length(unique(varCatEnrich$var_cat)))

VAR_CT_ENRICH_BAR <- plotBackground + 
  geom_bar(data = varCatEnrich,
           mapping = aes(x = gene_plots_id, y = perc_cat, fill = var_cat),
           position = "stack", stat = "identity") +
  geom_text(data = varCatEnrich[binomPadj <= args$p_abj & 
                                  var_cat %in% args$categoty_of_interest],
            mapping = aes(x = gene_plots_id, y = 101, label = "*")) + 
  scale_fill_manual(values = categoriesPallete) +
  scale_x_discrete(drop = F, labels = xAxisLabs$labels) + 
  xlab('driver genomic element') + ylab('% of variants') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  customGgplot2Theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   family = customGgplot2Theme[[1]]$axis.text$family,
                                   size = customGgplot2Theme[[1]]$axis.text$size,
                                   color = xAxisLabs$color),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'bottom', 
        legend.direction = "horizontal") +
  guides(fill = guide_legend(title = "variant category", 
                             nrow = 1, byrow = T))

VAR_CT_ENRICH_BAR <- annotatePlotWithGRlines(basePlot = VAR_CT_ENRICH_BAR,
                                             plotDT = cbind(dummy = 100,
                                                            varCatEnrich),
                                             axisLabs = xAxisLabs, 
                                             xOrderCol = 'gene_plots_id', 
                                             yCol = 'dummy', 
                                             direction = 'vertical',  
                                             angle = 45, vjust = -1, 
                                             size = 2)

