#!/usr/bin/env Rscript
# FILE:  bin/figures_create_mutual_exclusivity_and_coocurence_plot.R ---------------------------
#
# DESCRIPTION: An R script which creates a plot showing
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla  bin/figures_create_mutual_exclusivity_and_coocurence_plot.R -h
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

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'figures_create_mutual_exclusivity_and_coocurence_plot.R')

subtypeHelp <- paste('A name of the tumor subtype in which drivers were',
                     'detected.')
parser$add_argument("-c", "--cancer_subtype", required = T, 
                    type = 'character', help = subtypeHelp)

discoverHelp <- 'Path to DISCOVER result for the tumor subtype'
parser$add_argument("-d", "--discover_table", required = T, 
                    type = 'character', help = discoverHelp)

comparisonHelp <- paste('A comparison to be plotted. One of "all", "coding",',
                        'or "noncoding".')
parser$add_argument("-comp", "--comparison", required = T, 
                    type = 'character', default = NULL, help = comparisonHelp,
                    choices = c("all", "coding", "noncoding"))

pvalueHelp <- paste("Maximum p-value which is considered statistically",
                    "significant.")
parser$add_argument("-p", "--p_value", required = T, 
                    type = 'numeric', default = 0.05, help = pvalueHelp)

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

args$p_value <- as.numeric(args$p_value)

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

COLUMN_WITH_RAW_P <- 'p.value'
COLUMN_WITH_ADJ_P <- 'p.adj'

# Test input arguments --------------------------------------------------------
# args <- list(cancer_subtype = 'Panlung', 
#              discover_table = "completed_runs/2023_12_25/results/discover/discoverResults-Panlung--hg19.csv",
#              comparison = 'all', #all, coding, noncoding
#              p_value = 0.05, 
#              visuals_json = 'data/visual_parameters.json', 
#              output_type = 'pdf', output = 'test.pdf')

# args <- list(cancer_subtype = 'Adenocarcinoma',
#              discover_table = "completed_runs/2023_12_25/results/discover/discoverResults-Adenocarcinoma--hg19.csv",
#              comparison = 'all', #all, coding, noncoding
#              p_value = 0.05, 
#              visuals_json = 'data/visual_parameters.json', 
#              output_type = 'pdf', output = 'test.pdf')

# Read visuals JSON -----------------------------------------------------------
ESSENTIAL_VISUAL_NAMES <- c('ggplot2_theme', 'colors_tumor_type',
                            'colors_divergent_palette',
                            'excusivity_cooccurrence_plot_width', 
                            'excusivity_cooccurrence_plot_heigth')

visualParams <- readJsonWithVisualParameters(args$visuals_json)
message('[', Sys.time(), '] Read --visuals_json: ', args$visuals_json)

notFoundVisuals <- setdiff(ESSENTIAL_VISUAL_NAMES, names(visualParams))
if (length(notFoundVisuals)) {
  stop('[', Sys.time(), '] Following visuals: ', 
       paste(notFoundVisuals, collapse = ', '), ' not found in JSON.')
}

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
p_levels <- c(args$p_value, args$p_value/10, args$p_value/100)
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
  highlightColors <- c(visualParams$colors_tumor_type, 'multiple' = 'grey',
                       'placeholder' = 'white')
} else {
  discover[, is_placeholder := gene_name_X != 'placeholder']
  xAnnoCol <- 'is_placeholder'
  highlightColors <- 'white'
}
xAxisLabs <- formatAxisLabels(discover, 'Gene_X', 'gene_name_X', xAnnoCol, 
                              defaultColor = 'black',
                              highlightColor = highlightColors)
xAxisLabs$dt$gr_id <- gsub(',.*', '', xAxisLabs$dt$Gene_X)

# y axis labels & color
if ('tumor_subtype_spec_Y' %in% colnames(discover)) {
  yAnnoCol <- 'tumor_subtype_spec_Y'
  highlightColors <- c(visualParams$colors_tumor_type, 'multiple' = 'grey',
                       'placeholder' = 'white')
} else {
  discover[, is_placeholder := gene_name_Y != 'placeholder']
  yAnnoCol <- 'is_placeholder'
  highlightColors <- 'white'
}
yAxisLabs <- formatAxisLabels(discover, 'Gene_Y', 'gene_name_Y', yAnnoCol, 
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
  scale_fill_gradientn(colours = visualParams$colors_divergent_palette, na.value = '#FFFFFF', 
                       limits = c(-color_scale_limit, color_scale_limit)) +
  scale_color_manual(values = c('#74ee15', '#f52789', '#e900ff')) +
  scale_x_discrete(position = "top", labels = xAxisLabs$label) + 
  scale_y_discrete(position = "left", labels = yAxisLabs$label) + 
  xlab('') + ylab('') + visualParams$ggplot2_theme + 
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

# Output plot to file ---------------------------------------------------------
if (args$output_type == 'pdf') {
  pdf(args$output, width = visualParams$excusivity_cooccurrence_plot_width, 
      height = visualParams$excusivity_cooccurrence_plot_heigth)
} else {
  png(args$output, width = visualParams$excusivity_cooccurrence_plot_width, 
      height = visualParams$excusivity_cooccurrence_plot_heigth)
}
ANNOTATED_HEATMAP
dev.off()

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')