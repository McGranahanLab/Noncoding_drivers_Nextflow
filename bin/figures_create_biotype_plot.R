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

# Function: misc --------------------------------------------------------------
simplifyBiotype <- function(biotype_str) {
  simplified <- NA
  if (grepl("TSG", biotype_str) & grepl("oncogene", biotype_str)) {
    simplified <- 'TSG & OG'
  } else {
    if (grepl("TSG", biotype_str)) {
      simplified <- "TSG"
    }
    if (grepl("oncogene", biotype_str)) {
      simplified <- "OG"
    }
  }
  simplified
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
biotypePalette <- c('OG' = '#762A83', 'weak OG' = '#E7D4E8', 
                    'weak TSG' = '#D9F0D3', 'TSG' = '#1B7837', 
                    'TSG & OG' = '#FFEE99')

tealPalette <- c('#125A56','#238F9D', '#60BCE9', '#C6DBED')
orangePalette <- c('#F0E6B2', '#FD9A44','#F57634', '#E94C1F','#D11807', 
                   '#A01813')
neutralColor <- '#EEEEEE'

# Test input arguments --------------------------------------------------------
args <- list(cancer_subtype = 'Panlung',
             drivers_biotyped = "completed_runs/2023_12_25/results/tables/drivers_biotyped/driversBiotyped-Panlung--hg19.csv",
             min_n_tumors = 10)
# add args for 0.5 and 0.33
# min_n_tumors is not needed? Because it was already done

# Read driver biotype ---------------------------------------------------------
biotypes <- fread(args$drivers_biotyped, header = T, stringsAsFactors = F)
biotypes[, gene_plots_id := paste(gene_name, gr_id)]
# add filtering
biotypes[, perc_group := 100 * perc_group]
biotypes[, simplified_known_biotype := sapply(known_cancer_biotype, 
                                              simplifyBiotype)]

# Set up levels(order) for genomic regions ------------------------------------
# and genomic regions by maximum number of mutated patients (by only Mut type
# because CN is not used in the sorting of the regions of the overview plot)
GR_ORDER <- biotypes[alt_type == "Mut"]
if (nrow(GR_ORDER) != 0) {
  GR_ORDER <- GR_ORDER[order(-n_total)]
} 
GR_ORDER <- unique(GR_ORDER$gr_id)
biotypes[, gr_id := factor(gr_id, levels = GR_ORDER)]

# Set up levels(order) for val_group ------------------------------------------
val_group_order <- data.table(val_group = unique(biotypes$val_group))
val_group_order[, lower_bound := gsub(";.*", "", val_group)]
val_group_order[, lower_bound := gsub("^\\[", "", lower_bound)]
val_group_order[, lower_bound := gsub("[+].*", "", lower_bound)]
val_group_order[, lower_bound := as.numeric(lower_bound)]
val_group_order <- val_group_order[order(lower_bound)]
biotypes <- merge(biotypes, val_group_order, by = 'val_group', all = T)
biotypes[, val_group := factor(val_group, levels = val_group_order$val_group)]

# Set up levels(order) for genes ----------------------------------------------
gene_order_muts <- biotypes[leading == 'Mut'][order(gr_id, -percInact)]
gene_order_muts <- unique(gene_order_muts$gene_plots_id)

gene_order_cn <- biotypes[leading == 'CN'][order(gr_id, percAmp)]
gene_order_cn <- unique(gene_order_cn$gene_plots_id)

biotypes[, gene_plots_id := factor(gene_plots_id, 
                                   c(gene_order_muts, gene_order_cn))]
biotypes <- biotypes[order(gr_id, gene_plots_id)]
biotypes[, gene_plots_id := factor(gene_plots_id, unique(gene_plots_id))]

# Create plot background (rectangles to split genomic regions) ----------------
biotyped_drivers <- unique(biotypes[,.(gene_plots_id, gr_id, gene_name,
                                       is_known_cancer)])
# color of plot labels
xAxisLabs <- format_xLabels(biotyped_drivers, 'gene_plots_id', 'gene_name',
                            'is_known_cancer')
xAxisLabs$labels <- gsub("__and__", "&", xAxisLabs$labels)
# background plot to split different genomic regions from each other
plotBackground <- create_plot_background(cbind(dummy = 1, biotyped_drivers),
                                         xAxisLabs, 'gene_plots_id', 
                                         'dummy', direction = 'vertical')

# Inferred biotype plot -------------------------------------------------------
INFERRED_TILE <- plotBackground + 
  geom_tile(data = unique(biotypes[,.(gene_plots_id, inferred_biotype)]),
            mapping = aes(x = gene_plots_id, y = 1, fill = inferred_biotype)) + 
  scale_x_discrete(drop = F, labels = xAxisLabs$label) + 
  scale_y_continuous(breaks = 1, labels = 'inferred biotype',
                     expand = c(0, 0)) +
  scale_fill_manual(values = biotypePalette, na.value = '#FFFFFF') + 
  customGgplot2Theme + labs(fill = 'biotype') +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'bottom', legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1, byrow = T))

INFERRED_TILE <- annotatePlotWithGRlines(basePlot = INFERRED_TILE,
                                         plotDT = cbind(dummy = 1.5,
                                                        unique(biotypes[,.(gene_plots_id, gr_id,
                                                                           inferred_biotype)])),
                                         axisLabs = xAxisLabs, 
                                         xOrderCol = 'gene_plots_id', 
                                         yCol = 'dummy', 
                                         direction = 'vertical', 
                                         angle = 90, vjust = -1, 
                                         size = 2)

# Known biotype plot ----------------------------------------------------------
KNOWN_TILE <- plotBackground + 
  geom_tile(data = unique(biotypes[,.(gene_plots_id, simplified_known_biotype)]),
            mapping = aes(x = gene_plots_id, y = 1, 
                          fill = simplified_known_biotype)) + 
  scale_x_discrete(drop = F, labels = xAxisLabs$label) + 
  scale_y_continuous(breaks = 1, labels = 'CGC biotype',
                     expand = c(0, 0)) +
  scale_fill_manual(values = biotypePalette, na.value = '#FFFFFF') + 
  customGgplot2Theme + labs(fill = 'biotype') +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'bottom', legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1, byrow = T))

# Barplot for mutations -------------------------------------------------------
mutsPallete <- colorRampPalette(c(neutralColor, rev(tealPalette)))
mutsPallete <- mutsPallete(length(unique(biotypes[alt_type == 'Mut']$val_group)))

MUTS_BAR <- plotBackground + 
  geom_bar(data = biotypes[alt_type == 'Mut'],
           mapping = aes(x = gene_plots_id, y = perc_group, fill = val_group),
           position = "stack", stat = "identity") +
  scale_fill_manual(values = mutsPallete) +
  scale_x_discrete(drop = F, labels = xAxisLabs$labels) + 
  ylab('% of tumors with mut.') +
  scale_y_continuous() + 
  customGgplot2Theme + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'bottom', legend.direction = "horizontal") +
  guides(fill = guide_legend(title = "% mutated", nrow = 1, byrow = T)) +
  geom_hline(yintercept = c(33, 50), linetype = 2, linewidth = 0.5, 
             color = '#999999')

# Barplot for CNA -------------------------------------------------------------
delGroups <- unique(biotypes[alt_type == 'CN' & lower_bound < 1]$val_group)
ampGroups <- unique(biotypes[alt_type == 'CN' & lower_bound >= 1]$val_group)
neutralGroups <- setdiff(unique(biotypes[alt_type == 'CN']$val_group),
                         c(delGroups, ampGroups))
cnPallete <- c()
if (length(delGroups) > 0) {
  cnPallete <- c(cnPallete, colorRampPalette(tealPalette)(length(delGroups)))
}
if (length(neutralGroups) > 0) {
  cnPallete <- c(cnPallete, rep(neutralColor, length(neutralGroups)))
}
if (length(ampGroups) > 0) {
  cnPallete <- c(cnPallete, colorRampPalette(orangePalette)(length(ampGroups)))
}

CNA_BAR <- plotBackground + 
  geom_bar(data = biotypes[alt_type == 'CN'], 
           mapping = aes(x = gene_plots_id, y = perc_group, fill = val_group),
           position = "stack", stat = "identity") +
  scale_fill_manual(values = cnPallete) +
  xlab('genetic element') +
  scale_x_discrete(drop = F, labels = xAxisLabs$labels) + 
  ylab('% of tumors with CN') +
  scale_y_continuous(expand = c(0, 0)) + 
  customGgplot2Theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   family = customGgplot2Theme[[1]]$axis.text$family,
                                   size = customGgplot2Theme[[1]]$axis.text$size,
                                   color = xAxisLabs$color),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'bottom', 
        legend.direction = "horizontal") +
  guides(fill = guide_legend(title = "N.copies rel.ploidy", 
                             nrow = 1, byrow = T)) +
  geom_hline(yintercept = c(33, 50), linetype = 2, linewidth = 0.5, 
             color = '#999999')

# Aggange plots ---------------------------------------------------------------
LEGEND_MUTS_BAR <- extract_legend(MUTS_BAR)
LEGEND_CN_BAR <- extract_legend(CNA_BAR)
LEGEND_INFERRED_BIOTYPE <- extract_legend(INFERRED_TILE) 
LEGENDS <- ggarrange(LEGEND_INFERRED_BIOTYPE, LEGEND_MUTS_BAR, NULL,
                     LEGEND_CN_BAR, 
                     nrow = 2, ncol = 2, heights = c(1, 1), widths = c(5, 12))

# w = 6, h = 4
ggarrange(plotlist = list(INFERRED_TILE + theme(legend.position = "none"),
                          KNOWN_TILE + theme(legend.position = "none"),
                          MUTS_BAR + theme(legend.position = "none"),
                          CNA_BAR + theme(legend.position = "none"),
                          LEGENDS),
          heights = c(0.025, 0.025, 0.3, 0.5, 0.1), 
          align = 'v', ncol = 1)