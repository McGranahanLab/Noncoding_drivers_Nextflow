#!/usr/bin/env Rscript
# FILE:  figures_create_oncoplot.R --------------------------------------------
#
# DESCRIPTION: An R script which creates a plot showing noncoding driver 
# mutations in patients in which no coding driver mutations were detected.
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla figures_create_oncoplot.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: R v4.1.0, data.table, ggplot2, ggpubr, plyr
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  01.02.2023
# REVISION: 02.06.2024

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

# Test input arguments --------------------------------------------------------
args <- list(cancer_subtype = 'Panlung',
             driver_mutations = 'completed_runs/2023_12_25/results/tables/driver_mutations/driverMutationsNoncodingDriversOnly-Panlung--hg19.csv',
             visuals_json = 'data/visual_parameters.json',
             output_type = 'pdf', output = 'test.pdf')

# Read visuals JSON -----------------------------------------------------------
ESSENTIAL_VISUAL_NAMES <- c('ggplot2_theme', 'colors_confidence_level',
                            'onco_plot_width', 'onco_plot_heigth')

visualParams <- readJsonWithVisualParameters(args$visuals_json)
message('[', Sys.time(), '] Read --visuals_json: ', args$visuals_json)

notFoundVisuals <- setdiff(ESSENTIAL_VISUAL_NAMES, names(visualParams))
if (length(notFoundVisuals)) {
  stop('[', Sys.time(), '] Following visuals: ', 
       paste(notFoundVisuals, collapse = ', '), ' not found in JSON.')
}

# Read driver mutations -------------------------------------------------------
colsToKeep <- c("gr_id", "gene_id", "gene_name", "key", "var_type", 
                "participant_id", "confidenceLvl", "prob_is_driver_mle")
driverMuts <- fread(args$driver_mutations, header = T, select = colsToKeep,
                    stringsAsFactors = F)
driverMuts <- unique(driverMuts)
driverMuts[, is_driver := NULL]
message('[', Sys.time(), '] Read --driver_mutations: ', args$driver_mutations)

if (nrow(driverMuts) == 0) {
  message('[', Sys.time(), '] No mutations found. Exiting.')
  stop_quietly()
}

# Process for plotting --------------------------------------------------------
oncoplotDT <- driverMuts[,.(participant_id, gene_name, gene_id, 
                            gr_id, confidenceLvl)]
oncoplotDT <- oncoplotDT[, confidenceLvl := paste(unique(confidenceLvl), 
                                                  collapse = ';'), 
                         by = .(participant_id, gene_name, gene_id, gr_id)]
oncoplotDT <- unique(oncoplotDT)
oncoplotDT[grepl(';', confidenceLvl)]$confidenceLvl <- 'multihit'

oncoplotDT[, gene_plots_id := paste0(gene_name, '(', gr_id, ')')]
oncoplotDT[, n_tumours := length(unique(participant_id)), by = gene_plots_id]
oncoplotDT <- oncoplotDT[order(-n_tumours)]
oncoplotDT[, gene_plots_id := factor(gene_plots_id, 
                                     rev(unique(gene_plots_id)))]
oncoplotDT[, participant_id := factor(participant_id, unique(participant_id))]

ONCOPLOT <- ggplot(oncoplotDT, 
                   aes(x = participant_id, y = gene_plots_id, 
                       fill = confidenceLvl)) +
  geom_tile() + visualParams$ggplot2_theme + 
  labs(fill = 'confidence level') + ylab('driver gene') + 
  scale_fill_manual(values = visualParams$colors_confidence_level) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position = 'bottom',
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(nrow = 1, byrow = F))

# Output to file --------------------------------------------------------------
if (args$output_type == 'pdf') {
  pdf(args$output, width = visualParams$onco_plot_width, 
      height = visualParams$onco_plot_heigth)
} else {
  png(args$output, width = visualParams$onco_plot_width, 
      height = visualParams$onco_plot_heigth)
}
ONCOPLOT
dev.off()

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')