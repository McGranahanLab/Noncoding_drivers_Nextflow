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

# Functions: pre-processing ---------------------------------------------------
#' prepCompatibilityForPlotting
#' @description Modifies results of co-occurence/exclusivity analysis in order 
#' to plot them as tile plot. Only results of analysis of one tumor subtype 
#' should be submitted.
#' @param drivesCombatDT result of of co-occurence/exclusivity analysis
#' @param pValCol column containing p-values which will be used for coloring
#' @return data table with added column log_pAdj
prepCompatibilityForPlotting <- function(drivesCombatDT, pValCol = 'p.value') {
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
             comparison = 'all', #all, coding, noncoding
             drivers = 'completed_runs/2023_12_25/results/tables/drivers/drivers-Panlung--hg19.csv',
             discover_table = "completed_runs/2023_12_25/results/discover/discoverResults-Panlung--hg19.csv",
             p_value_min = 0.05)

COLUMN_WITH_RAW_P <- 'p.value_subtypeBiasRemoved'
COLUMN_WITH_ADJ_P <- 'p.adj_subtypeBiasRemoved'

# Read unfiltered drivers -----------------------------------------------------
drivers <- fread(args$drivers, header = T, stringsAsFactors = F)
drivers <- drivers[!is.na(tier) & FILTER %in% 'PASS']
drivers <- drivers[,.(gr_id, gene_id, gene_name, is_known_cancer)]
drivers[, gene_name := gsub('__and__', '&', gene_name)]
drivers[, gene_id := gsub('__and__', '&', gene_id)]
drivers <- unique(drivers)

message('[', Sys.time(), '] Read --drivers: ', args$drivers)

# Read discover table ---------------------------------------------------------
discover <- fread(args$discover_table, header = T, stringsAsFactors = F)
discover <- discover[comparison == args$comparison]
discover[, gene_name_1 := gsub('__and__', '&', gene_name_1)]
discover[, gene_name_2 := gsub('__and__', '&', gene_name_2)]
discover[, gene_id_1 := gsub('__and__', '&', gene_id_1)]
discover[, gene_id_2 := gsub('__and__', '&', gene_id_2)]
message('[', Sys.time(), '] Read --discover_table: ', args$discover_table)

# add info about if gene is known to be related to cancer
setnames(discover, c('gr_id_1', 'gene_name_1', 'gene_id_1'),
         c('gr_id', 'gene_name', 'gene_id'))
discover <- merge(discover, drivers, by = c('gr_id', 'gene_name', 'gene_id'),
                  all.x = T)
setnames(discover, c('gr_id', 'gene_name', 'gene_id', 'is_known_cancer'),
         c('gr_id_1', 'gene_name_1', 'gene_id_1', 'is_known_cancer_1'))
setnames(discover, c('gr_id_2', 'gene_name_2', 'gene_id_2'),
         c('gr_id', 'gene_name', 'gene_id'))
discover <- merge(discover, drivers, by = c('gr_id', 'gene_name', 'gene_id'),
                  all.x = T)
setnames(discover, c('gr_id', 'gene_name', 'gene_id', 'is_known_cancer'),
         c('gr_id_2', 'gene_name_2', 'gene_id_2', 'is_known_cancer_2'))

discover[, Gene_1 := paste0(gr_id_1, ', ', gene_name_1)]
discover[, Gene_2 := paste0(gr_id_2, ', ', gene_name_2)]
discover <- prepCompatibilityForPlotting(discover, COLUMN_WITH_RAW_P)
discover[, gr_id_1 := factor(gr_id_1, rev(unique(gr_id_1)))]
discover[, gr_id_2 := factor(gr_id_2, unique(gr_id_2))]

# Assign plot labels ----------------------------------------------------------
discover[, cell_label := character()]
p_levels <- c(args$p_value_min, args$p_value_min/10, args$p_value_min/100)
for (p_level in p_levels) {
  labels_to_update_idx <- which(discover[, COLUMN_WITH_ADJ_P, 
                                         with = F] <= p_level)
  discover[labels_to_update_idx]$cell_label <- paste('<', p_level)
}
discover[, cell_label := factor(cell_label, levels = paste('<', p_levels))]

support_by_indivTT

# Create plot -----------------------------------------------------------------
# x axis labels & color
xAxisLabs <- format_xLabels(discover, 'Gene_1', 'gene_name_1', 
                            'is_known_cancer_1') 
# y axis labels & color
yAxisLabs <- format_xLabels(discover, 'Gene_2', 'gene_name_2', 
                            'is_known_cancer_2')

color_scale_limit <- 10
COMPAT_PLOT <- ggplot() + 
  geom_tile(data = discover, 
            mapping = aes_string(x = 'Gene_2', y = 'Gene_1', 
                                 fill = paste0("log10.", COLUMN_WITH_RAW_P))) +
  geom_point(data = discover[!is.na(cell_label)], size = 2,inherit.aes = F,
             mapping = aes(x = Gene_2, y = Gene_1, shape = cell_label)) +
  facet_grid(cols = vars(gr_id_2), rows = vars(gr_id_1), scales = 'free',
             space = 'free', switch = 'y') +
  scale_fill_gradientn(colours = c(tealPalette, neutralColor,
                                   orangePalette),
                       na.value = '#FFFFFF', 
                       limits = c(-color_scale_limit, color_scale_limit)) +
  scale_x_discrete(position = "top", labels = xAxisLabs$label) + 
  scale_y_discrete(position = "left", labels = yAxisLabs$label) + 
  xlab('') + ylab('') + customGgplot2Theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = -1, size = 6),
        axis.text.y = element_text(size = 6),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.direction = 'horizontal', legend.position = 'top',
        axis.line = element_blank(), 
        legend.key.height = unit(0.25, 'cm'), 
        legend.key.width = unit(2, 'cm')) + 
  labs(fill = '-log10(raw p.value)') + guides(color = F) +
  ggtitle(paste0('Co-occurence and exclusivity in ', 
                 unique(discover$tumor_subtype), ', ', args$comparison))
