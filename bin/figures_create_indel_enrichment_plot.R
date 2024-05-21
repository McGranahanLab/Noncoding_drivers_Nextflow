#!/usr/bin/env Rscript
# FILE:  bin/figures_create_indel_enrichment_plot.R ---------------------------
#
# DESCRIPTION: An R script which creates a plot showing enrichment of drivers
# in certain variant type (i.e. SNPs or INDELs of a specific size) categories.
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla  bin/figures_create_indel_enrichment_plot.R -h
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

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'figures_create_indel_enrichment_plot.R')

subtypeHelp <- paste('A name of the tumor subtype in which drivers were',
                     'detected.')
parser$add_argument("-c", "--cancer_subtype", required = T, 
                    type = 'character', help = subtypeHelp)

driversHelp <- paste('Path to file containing de-novo detected cancer',
                     'drivers in the cancer subtype submitted via',
                     '--cancer_subtype')
parser$add_argument("-d", "--drivers", required = T, type = 'character',
                    default = NULL, help = driversHelp)

filterHelp <- paste("Values of FILTER column which are accepted to be plotted")
parser$add_argument("-afv", "--allowed_filter_values", required = F, 
                    type = 'character', default = NULL, nargs = '+',
                    help = filterHelp)

countsHelp <- paste("Path to file containing information about counts of",
                    "variants of different categories in genomic elements.",
                    "Required columns: gr_id, gene_id, gene_name, var_cat",
                    "Nvars, N_gene, binomPadj.")
parser$add_argument("-vcc", "--variant_categories_counts", required = T, 
                    type = 'character', default = NULL, help = countsHelp)

categoryHelp <- paste("Variant category for which p-values should be marked",
                      "as significant.")
parser$add_argument("-coi", "--category_of_interest", required = T, 
                    type = 'character', default = NULL, help = categoryHelp, 
                    choices = c("SNP", "INDEL, 1bp", "INDEL, 2-5bp", 
                                "INDEL, 6-9bp", "INDEL, 10+bp", "MNP"))

pAdjHelp <- paste("Maximum adjusted for multiple testing p-value which",
                  "is considered statistically significant.")
parser$add_argument("-p", "--p_adj", required = T, type = 'numeric', 
                    default = 0.05, help = pAdjHelp) 

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

if (is.null(args$allowed_filter_values)) {
  args$allowed_filter_values <- list('PASS')
}
args$p_adj <- as.numeric(args$p_adj)

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

VAR_CAT_LEVELS <- rev(c("SNP", "INDEL, 1bp", "INDEL, 2-5bp", "INDEL, 6-9bp",
                        "INDEL, 10+bp", "MNP"))

# Test input arguments --------------------------------------------------------
# args <- list(cancer_subtype = 'Panlung',
#              drivers = "completed_runs/2023_12_25/results/tables/drivers/drivers-Panlung--hg19.csv",
#              allowed_filter_values = list("PASS", "INDEL, 2-5bp"),
#              variant_categories_counts = "completed_runs/2023_12_25/results/mut_rates/varCatEnrich-Panlung--hg19.csv",
#              category_of_interest = "INDEL, 2-5bp", p_adj = 0.05)

# Read visuals JSON -----------------------------------------------------------
ESSENTIAL_VISUAL_NAMES <- c('ggplot2_theme', 'color_divider',
                            'colors_divergent_palette',
                            'indel_enrichment_plot_width', 
                            'indel_enrichment_plot_heigth')

visualParams <- readJsonWithVisualParameters(args$visuals_json)
message('[', Sys.time(), '] Read --visuals_json: ', args$visuals_json)

notFoundVisuals <- setdiff(ESSENTIAL_VISUAL_NAMES, names(visualParams))
if (length(notFoundVisuals)) {
  stop('[', Sys.time(), '] Following visuals: ', 
       paste(notFoundVisuals, collapse = ', '), ' not found in JSON.')
}

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
xAxisLabs <- formatAxisLabels(varCatEnrich, 'gene_plots_id', 'gene_name',
                              'is_known_cancer')
xAxisLabs$labels <- gsub("__and__", "&", xAxisLabs$labels)
# background plot to split different genomic regions from each other
plotBackground <- create_plot_background(varCatEnrich, xAxisLabs, 
                                         'gene_plots_id', 'perc_cat', 
                                         direction = 'vertical')

# Plot bar plot enrichment of variant categories  -----------------------------
neutralColorIdx <- median(1:length(visualParams$colors_divergent_palette))
neutralColor <- visualParams$colors_divergent_palette[neutralColorIdx]
coldPalette <- visualParams$colors_divergent_palette[1:(neutralColorIdx - 1)]
categoriesPallete <- colorRampPalette(coldPalette)
categoriesPallete <- categoriesPallete(length(unique(varCatEnrich$var_cat)))

VAR_CT_ENRICH_BAR <- plotBackground + 
  geom_bar(data = varCatEnrich,
           mapping = aes(x = gene_plots_id, y = perc_cat, fill = var_cat),
           position = "stack", stat = "identity") +
  geom_text(data = varCatEnrich[binomPadj <= args$p_adj & 
                                  var_cat %in% args$category_of_interest],
            mapping = aes(x = gene_plots_id, y = 101, label = "*")) + 
  scale_fill_manual(values = categoriesPallete) +
  scale_x_discrete(drop = F, labels = xAxisLabs$labels) + 
  xlab('driver genomic element') + ylab('% of variants') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  visualParams$ggplot2_theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   family = visualParams$ggplot2_theme[[1]]$axis.text$family,
                                   size = visualParams$ggplot2_theme[[1]]$axis.text$size,
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
# Output plot to file ---------------------------------------------------------
# w = 6, h = 4
if (args$output_type == 'pdf') {
  pdf(args$output, width = visualParams$indel_enrichment_plot_width, 
      height = visualParams$indel_enrichment_plot_heigth)
} else {
  png(args$output, width = visualParams$indel_enrichment_plot_width, 
      height = visualParams$indel_enrichment_plot_heigth)
}
VAR_CT_ENRICH_BAR
dev.off()

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')