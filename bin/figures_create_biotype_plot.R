#!/usr/bin/env Rscript
# FILE:  bin/figures_create_biotype_plot.R ------------------------------------
#
# DESCRIPTION: An R script which creates a plot showing cancer biotype (i.e.
# tumor suppressor (TSG) or oncogene(OG)) of detected drivers.
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla  bin/figures_create_biotype_plot.R -h
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

# Functions: misc -------------------------------------------------------------
#' simplifyBiotype
#' @description Simplifies notation of biotype
#' @author Maria Litovchenko
#' @param biotype_str string, notation of biotype
#' @return string, one of "TSG & OG", "TSG" or "OG"
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

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'figures_create_biotype_plot.R')

subtypeHelp <- paste('A name of the tumor subtype in which drivers were',
                     'detected.')
parser$add_argument("-c", "--cancer_subtype", required = T, 
                    type = 'character', help = subtypeHelp)

driversHelp <- paste('Path to file containing inferred biotypes of drivers')
parser$add_argument("-d", "--drivers_biotyped", required = T,
                    type = 'character', default = NULL, help = driversHelp)

foldHelp <- paste0('Indicates, whether or not splice sites should be',
                   'considered as part of CDS for DISCOVER analysis.',
                   'Default: T.')
parser$add_argument("-f", "--fold_splicesites_in_coding", required = F,
                    default = 'T', choices = c('T', 'F'), type = 'character', 
                    help = foldHelp)

weakTSGhelp <- paste('A munimal cut off on percentage of biallelically',
                     'inactivated patients from the total number of patients',
                     'mutated (SNV and small indels) in a genomic region for',
                     'the region to be assigned as weak tumor suppressor gene',
                     '(TSG). In reverse, for a weak oncogene (OG) the ',
                     'percentage should be below this number. Should range',
                     'between 0 and 1.')
parser$add_argument("-wTSG", "--weak_tsg", required = T, type = 'double', 
                    default = 0.33, help = weakTSGhelp)

tsgHelp <- paste('A munimal cut off on percentage of biallelically',
                 'inactivated patients from the total number of patients',
                 'mutated (SNV and small indels) in a genomic region for',
                 'the region to be assigned as tumor suppressor gene (TSG).',
                 'In reverse, for an oncogene (OG) the percentage should be',
                 'below this number. Should range between 0 and 1.')
parser$add_argument("-TSG", "--tsg", required = T, type = 'double',
                    default = 0.50, help = tsgHelp)

weakOGhelp <- paste('A munimal cut off on percentage of patients with an',
                    'amplification from the total number of patients',
                    'with a copy number variatons in a genomic region for',
                    'the region to be assigned as a weak oncogene (OG).',
                    'In reverse, for a weak tumor supressor gene (TSG) the',
                    'percentage should be below this number. Should range',
                    'between 0 and 1.')
parser$add_argument("-wOG", "--weak_og", required = T, type = 'double',
                    default = 0.33, help = weakOGhelp)

ogHelp <- paste('A munimal cut off on percentage of patients with an',
                'amplification from the total number of patients',
                'with a copy number variatons in a genomic region for',
                'the region to be assigned as an oncogene (OG). In reverse,',
                'for a tumor supressor gene (TSG) the percentage should be',
                'below this number. Should range between 0 and 1.')
parser$add_argument("-OG", "--og", required = T, type = 'double',
                    default = 0.50, help = ogHelp)

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
args$fold_splicesites_in_coding <- as.logical(args$fold_splicesites_in_coding)

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test input arguments --------------------------------------------------------
args <- list(cancer_subtype = 'Panlung',
             drivers_biotyped = "completed_runs/2023_12_25/results/tables/drivers_biotyped/driversBiotyped-Panlung--hg19.csv",
             fold_splicesites_in_coding = T, weak_tsg = 0.33, tsg = 0.5,
             weak_og = 0.33, og = 0.5, visuals_json = '', 
             output_type = 'pdf', output = 'test.pdf')

# Read visuals JSON -----------------------------------------------------------
ESSENTIAL_VISUAL_NAMES <- c('ggplot2_theme', 'color_divider',
                            'colors_biotype', 'colors_divergent_palette',
                            'biotype_plot_width', 'biotype_plot_heigth')

visualParams <- readJsonWithVisualParameters(args$visuals_json)
message('[', Sys.time(), '] Read --visuals_json: ', args$visuals_json)

notFoundVisuals <- setdiff(ESSENTIAL_VISUAL_NAMES, names(visualParams))
if (length(notFoundVisuals)) {
  stop('[', Sys.time(), '] Following visuals: ', 
       paste(notFoundVisuals, collapse = ', '), ' not found in JSON.')
}

# Read driver biotype ---------------------------------------------------------
biotypes <- fread(args$drivers_biotyped, header = T, stringsAsFactors = F)
message('[', Sys.time(), '] Read --biotypes: ', args$biotypes)
biotypes[, gene_plots_id := paste(gene_name, gr_id)]
# add filtering
biotypes[, perc_group := 100 * perc_group]
biotypes[, simplified_known_biotype := sapply(known_cancer_biotype, 
                                              simplifyBiotype)]

# Fold splice sites into coding regions, if requested -------------------------
if (args$fold_splicesites_in_coding) {
  message('[', Sys.time(), '] Will fold splice sites into corresponding ',
          'coding driver genetic elements.')
  
  # assign coding genomic regions - anything which contains CDS as gr_code
  coding_gr_id <- unique(analysisInv[gr_code == 'CDS']$gr_id)
  message('[', Sys.time(), '] Following genomic regions: ', 
          paste0(coding_gr_id, collapse = ', '), ', will be considered as ',
          'coding.')
  # assign splice site genomic regions - anything which contains ss as gr_code
  ss_gr_id <- unique(analysisInv[gr_code == 'ss']$gr_id)
  message('[', Sys.time(), '] Following genomic regions: ', 
          paste0(ss_gr_id, collapse = ', '), ', will be considered as ',
          'containing splice sites.')
  
  if (length(coding_gr_id) != 0 & length(ss_gr_id) != 0) {
    biotypes <- foldSplicSiteDriversIntoCodingDrivers(coding_gr_id, ss_gr_id,
                                                      biotypes)
    message('[', Sys.time(), "] Driver genomic elements considered as ",
            "containing splice sites will be treated as noncoding for genes ",
            "which do not have coding part detected as driver.")
  } else {
    message('[', Sys.time(), '] Did not find either coding or splice site ',
            'regions. Can not fold splice site drivers into coding ones.')
  }
}

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
  geom_tile(data = unique(biotypes[,.(gene_plots_id, 
                                      simplified_known_biotype)]),
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
  geom_hline(yintercept = c(100*args$weak_tsg, 100*args$tsg), linetype = 2,
             linewidth = 0.5, color = '#999999')

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
  geom_hline(yintercept = c(100*args$weak_og, 100*args$og), linetype = 2, 
             linewidth = 0.5, color = '#999999')

# Aggange plots ---------------------------------------------------------------
LEGEND_MUTS_BAR <- extract_legend(MUTS_BAR)
LEGEND_CN_BAR <- extract_legend(CNA_BAR)
LEGEND_INFERRED_BIOTYPE <- extract_legend(INFERRED_TILE) 
LEGENDS <- ggarrange(LEGEND_INFERRED_BIOTYPE, LEGEND_MUTS_BAR, NULL,
                     LEGEND_CN_BAR, 
                     nrow = 2, ncol = 2, heights = c(1, 1), widths = c(5, 12))

# Output plot to file ---------------------------------------------------------
# w = 6, h = 4
if (args$output_type == 'pdf') {
  pdf(args$output, width = visualParams$biotype_plot_width, 
      height = visualParams$biotype_plot_heigth)
} else {
  png(args$output, width = visualParams$biotype_plot_width, 
      height = visualParams$biotype_plot_heigth)
}
ggarrange(plotlist = list(INFERRED_TILE + theme(legend.position = "none"),
                          KNOWN_TILE + theme(legend.position = "none"),
                          MUTS_BAR + theme(legend.position = "none"),
                          CNA_BAR + theme(legend.position = "none"),
                          LEGENDS),
          heights = c(0.025, 0.025, 0.3, 0.5, 0.1), 
          align = 'v', ncol = 1)
dev.off()

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')