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
# Functions: pre-processing ---------------------------------------------------
#' prepDISCOVERforPlotting
#' @description Modifies results of co-occurence/exclusivity analysis in order 
#' to plot them as tile plot. Only results of analysis of one tumor subtype 
#' should be submitted.
#' @author Maria Litovchenko
#' @param drivesCombatDT result of of co-occurence/exclusivity analysis
#' @param pValCol column containing p-values which will be used for coloring
#' @return data table with added column log_pAdj
prepDISCOVERforPlotting <- function(drivesCombatDT, pValCol = 'p.value') {
  # order by adjusted p-values and select co-oc/excl. which is the most 
  # significant per gene pair
  result <- copy(drivesCombatDT)
  setnames(result, pValCol, 'colWpval')
  result <- result[order(Gene_1, Gene_2, colWpval)]
  result <- result[,.SD[1], by = .(Gene_1, Gene_2)]
  result[, logPvalCol := -log10(colWpval)]
  result[is.infinite(logPvalCol)]$logPvalCol <- max(result$logPvalCol,
                                                    na.rm = T)
  toRevert <- result[mode == 'exclusivity']$logPvalCol
  result[mode == 'exclusivity']$logPvalCol <- -toRevert
  
  geneLvls <- sort(unique(unlist(result[,.(Gene_1, Gene_2)])))
  result[, Gene_1 := factor(Gene_1, geneLvls)]
  result[, Gene_2 := factor(Gene_2, geneLvls)]
  
  setnames(result, c('colWpval', 'logPvalCol'),
           c(pValCol, paste0('log10.', pValCol)))
  colnames(result) <- gsub('_1$', '_Y', colnames(result))
  colnames(result) <- gsub('_2$', '_X', colnames(result))
  
  result
}

#' addGenomicRegionsSpacer
#' @description Adds empty & dummy cells in the DISCOVER results matrix which 
#' will serve as spacers between drivers of different genomic regions in the 
#' future heatmap.
#' @author Maria Litovchenko 
#' @param discoverResultsOneGR data table containing DISCOVER results for one
#' type of genomic region (i.e. CDS only)
#' @param dummyGene string, name of dummy genes
#' @param axisID, string, one of X or Y which designates axis into which 
#' spacers will be inserted
#' @return data table discoverResultsOneGR with added rows which will act as 
#' spacers
addGenomicRegionsSpacer <- function(discoverResultsOneGR, dummyGene, axisID) {
  gr_id <- discoverResultsOneGR[, paste0('gr_id_', axisID), with = F]
  gr_id <- unique(unlist(gr_id))
  
  spacerDT <- data.table('Gene' = paste0(gr_id, ', ', dummyGene),  
                         'gr_id' = gr_id,
                         'gene_name' = dummyGene,
                         'gene_id' = dummyGene)
  colnames(spacerDT) <- paste0(colnames(spacerDT), '_', axisID)
  spacerDT[, p.value := NA]
  spacerDT[, p.adj := NA]
  spacerDT[, log10.p.value := NA]
  if ('supported_by_individual_subtypes' %in% colnames(discoverResultsOneGR)) {
    spacerDT[, supported_by_individual_subtypes := NA]
  } 
  if ('incompatible_due_subtype_specificity' %in% 
      colnames(discoverResultsOneGR)) {
    spacerDT[, incompatible_due_subtype_specificity := NA]
  }
  
  otherAxisCols <- setdiff(colnames(discoverResultsOneGR), colnames(spacerDT))
  otherAxisData <- discoverResultsOneGR[, otherAxisCols, with = F]
  otherAxisData <- unique(otherAxisData)
  spacerDT <- cbind(spacerDT, otherAxisData)
  setcolorder(spacerDT, colnames(discoverResultsOneGR))
  
  return(rbind(discoverResultsOneGR, spacerDT))
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

gradientColors <- c('#125A56', '#00767B', '#238F9D', '#42A7C6', '#60BCE9',
                    '#9DCCEF', '#C6DBED', '#DEE6E7', '#ECEADA', '#F0E6B2',
                    '#F9D576', '#FFB954', '#FD9A44', '#F57634', '#E94C1F',
                    '#D11807', '#A01813')

# Test input arguments --------------------------------------------------------
args <- list(cancer_subtype = 'Panlung', 
             discover_table = "completed_runs/2023_12_25/results/discover/discoverResults-Panlung--hg19.csv",
             comparison = 'all', #all, coding, noncoding
             p_value_min = 0.05)

args <- list(cancer_subtype = 'Adenocarcinoma',
             discover_table = "completed_runs/2023_12_25/results/discover/discoverResults-Adenocarcinoma--hg19.csv",
             comparison = 'all', #all, coding, noncoding
             p_value_min = 0.05)

COLUMN_WITH_RAW_P <- 'p.value'
COLUMN_WITH_ADJ_P <- 'p.adj'

# Read discover table ---------------------------------------------------------
discover <- fread(args$discover_table, header = T, stringsAsFactors = F)
discover <- discover[comparison == args$comparison]
discover[, gene_name_1 := gsub('__and__', '&', gene_name_1)]
discover[, gene_name_2 := gsub('__and__', '&', gene_name_2)]
discover[, gene_id_1 := gsub('__and__', '&', gene_id_1)]
discover[, gene_id_2 := gsub('__and__', '&', gene_id_2)]
message('[', Sys.time(), '] Read --discover_table: ', args$discover_table)

# Prepare discover table for plotting -----------------------------------------
discover[, Gene_1 := paste0(gr_id_1, ', ', gene_name_1)]
discover[, Gene_2 := paste0(gr_id_2, ', ', gene_name_2)]

discover <- prepDISCOVERforPlotting(discover, COLUMN_WITH_RAW_P)
discover <- split(discover, by = 'gr_id_Y')
discover <- lapply(discover, addGenomicRegionsSpacer, 'placeholder', 'Y')
discover <- do.call(rbind, discover)
extraPlaceholderY <- unique(discover$Gene_Y)[length(unique(discover$Gene_Y))]
discover <- discover[Gene_Y != extraPlaceholderY]
discover[, Gene_Y := factor(Gene_Y, unique(Gene_Y))]
if ('tumor_subtype_spec_Y' %in% colnames(discover)) {
  discover[gene_name_Y == 'placeholder']$tumor_subtype_spec_Y <- 'placeholder'
}


discover <- split(discover, by = 'gr_id_X')
discover <- lapply(discover, addGenomicRegionsSpacer, 'placeholder', 'X')
discover <- do.call(rbind, discover)
extraPlaceholderX <- unique(discover$Gene_X)[length(unique(discover$Gene_X))]
discover <- discover[Gene_X != extraPlaceholderX]
discover[, Gene_X := factor(Gene_X, unique(Gene_X))]
if ('tumor_subtype_spec_X' %in% colnames(discover)) {
  discover[gene_name_X == 'placeholder']$tumor_subtype_spec_X <- 'placeholder'
}

# Assign plot labels ----------------------------------------------------------
discover[, cell_label := character()]
p_levels <- c(args$p_value_min, args$p_value_min/10, args$p_value_min/100)
for (p_level in p_levels) {
  labels_to_update_idx <- which(discover[, COLUMN_WITH_ADJ_P, 
                                         with = F] <= p_level)
  discover[labels_to_update_idx]$cell_label <- paste('<', p_level)
}
discover[, cell_label := factor(cell_label, levels = paste('<', p_levels))]

# Format axis tick labels -----------------------------------------------------
# x axis labels & color
if ('tumor_subtype_spec_X' %in% colnames(discover)) {
  xAnnoCol <- 'tumor_subtype_spec_X'
  highlightColors <- c(tumourTypeColorPalette, 'multiple' = 'grey',
                       'placeholder' = 'white')
} else {
  discover[, is_placeholder := gene_name_X != 'placeholder']
  xAnnoCol <- 'is_placeholder'
  highlightColors <- 'white'
}
xAxisLabs <- format_xLabels(discover, 'Gene_X', 'gene_name_X', xAnnoCol, 
                            defaultColor = 'black',
                            highlightColor = highlightColors)
xAxisLabs$dt$gr_id <- gsub(',.*', '', xAxisLabs$dt$Gene_X)

# y axis labels & color
if ('tumor_subtype_spec_Y' %in% colnames(discover)) {
  yAnnoCol <- 'tumor_subtype_spec_Y'
  highlightColors <- c(tumourTypeColorPalette, 'multiple' = 'grey',
                       'placeholder' = 'white')
} else {
  discover[, is_placeholder := gene_name_Y != 'placeholder']
  yAnnoCol <- 'is_placeholder'
  highlightColors <- 'white'
}
yAxisLabs <- format_xLabels(discover, 'Gene_Y', 'gene_name_Y', yAnnoCol, 
                            defaultColor = 'black',
                            highlightColor = highlightColors)
yAxisLabs$dt$gr_id <- gsub(',.*', '', yAxisLabs$dt$Gene_Y)

# Create plot -----------------------------------------------------------------
color_scale_limit <- max(abs(log10(discover$p.value)), na.rm = T)

HEATMAP <- ggplot() + 
  geom_tile(data = discover, 
            mapping = aes_string(x = 'Gene_X', y = 'Gene_Y', 
                                 fill = paste0("log10.", COLUMN_WITH_RAW_P))) +
  geom_point(data = discover[!is.na(cell_label)], size = 1, inherit.aes = F,
             mapping = aes(x = Gene_X, y = Gene_Y, color = cell_label)) + 
  scale_fill_gradientn(colours = gradientColors, na.value = '#FFFFFF', 
                       limits = c(-color_scale_limit, color_scale_limit)) +
  scale_color_manual(values = c('#74ee15', '#f52789', '#e900ff')) +
  scale_x_discrete(position = "top", labels = xAxisLabs$label) + 
  scale_y_discrete(position = "left", labels = yAxisLabs$label) + 
  xlab('') + ylab('') + customGgplot2Theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = -1, size = 6,
                                   color = xAxisLabs$color, face = ),
        axis.text.y = element_text(size = 6, color = yAxisLabs$color),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.line = element_blank(), 
        legend.direction = 'horizontal', legend.position = 'top',
        legend.key.height = unit(0.25, 'cm'), 
        legend.key.width = unit(1, 'cm')) + 
  labs(fill = '-log10(raw p.value)', color = 'Adjusted p-value')

if ('supported_by_individual_subtypes' %in% colnames(discover)) {
  HEATMAP <- HEATMAP + 
    geom_point(data = discover[supported_by_individual_subtypes == T & 
                                 !is.na(p.value)],
               mapping = aes(x = Gene_X, y = Gene_Y), shape = 3, size = 2,
               inherit.aes = F, color = '#555555')
}
if ('incompatible_due_subtype_specificity' %in% colnames(discover)) {
  HEATMAP <- HEATMAP + 
    geom_point(data = discover[incompatible_due_subtype_specificity == T & 
                                 !is.na(p.value)],
               mapping = aes(x = Gene_X, y = Gene_Y), shape = 4, size = 2, 
               inherit.aes = F, color = '#555555') 
}

LEGENDS <- extract_legend(HEATMAP)

TOP_LABELS <- ggplot(data = discover, aes(x = Gene_X, y = '1')) +
  geom_point(color = 'white')
TOP_LABELS <- annotatePlotWithGRlines(basePlot = TOP_LABELS,
                                      plotDT = cbind(discover, dummy = '1',
                                                     gr_id = discover$gr_id_X),
                                      axisLabs = xAxisLabs, 
                                      xOrderCol = 'Gene_X', 
                                      yCol = 'dummy', 
                                      direction = 'vertical',  
                                      vjust = -0.5, size = 2)
TOP_LABELS <- TOP_LABELS + theme_void()

LEFT_LABELS <- ggplot(data = discover, aes(y = Gene_Y, x = '1')) +
  geom_point(color = 'white')
LEFT_LABELS <- annotatePlotWithGRlines(basePlot = LEFT_LABELS,
                                       plotDT = cbind(discover, dummy = '1',
                                                      gr_id = discover$gr_id_Y),
                                       axisLabs = yAxisLabs, 
                                       xOrderCol = 'Gene_Y', 
                                       yCol = 'dummy', 
                                       direction = 'horizontal',  
                                       hjust = 1.1, size = 2)
LEFT_LABELS <- LEFT_LABELS + theme_void()

ANNOTATED_HEATMAP <- ggarrange(LEGENDS,
                               ggarrange(plotlist = list(NULL, TOP_LABELS,
                                                         LEFT_LABELS, 
                                                         HEATMAP + 
                                                           theme(legend.position = 'none')),
                                         ncol = 2, nrow = 2, align = 'hv', 
                                         heights = c(0.15, 1),
                                         widths = c(0.15, 1)),
                               ncol = 1, heights = c(0.05, 1))
