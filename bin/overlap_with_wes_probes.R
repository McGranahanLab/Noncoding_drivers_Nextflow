#!/usr/bin/env Rscript
# FILE: overlap_with_wes_probes.R ---------------------------------------------
#
# DESCRIPTION: checks the percetage of overlap of genomic elements detected as
# drivers with the WES probes
# USAGE: 
#
# OPTIONS: Run 
#          Rscript --vanilla overlap_with_wes_probes.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: 
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  22.06.2023
# REVISION: 29.05.2024

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
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(plyranges)))

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'overlap_with_wes_probes.R')

args <- parser$parse_args()
timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test inputs -----------------------------------------------------------------
args <- list('probes' = 'data/assets/Agilent_all_human_exome_v4_hg19/S04380110_hs_hg19/S04380110_Covered.bed',
             'genomic_regions' = 'completed_runs/2023_12_25/inputs/inputGR-Panlung-hg19.bed',
             'drivers' = 'completed_runs/2023_12_25/results/tables/drivers/drivers-Panlung--hg19.csv',
             'visuals_json' = 'data/visual_parameters.json',
             'output_type' = 'pdf', 'output' = 'test.pdf')

# Read visuals JSON -----------------------------------------------------------
ESSENTIAL_VISUAL_NAMES <- c('ggplot2_theme', 'color_divider',
                            'indel_enrichment_plot_width', 
                            'indel_enrichment_plot_heigth')

visualParams <- readJsonWithVisualParameters(args$visuals_json)
message('[', Sys.time(), '] Read --visuals_json: ', args$visuals_json)

notFoundVisuals <- setdiff(ESSENTIAL_VISUAL_NAMES, names(visualParams))
if (length(notFoundVisuals)) {
  stop('[', Sys.time(), '] Following visuals: ', 
       paste(notFoundVisuals, collapse = ', '), ' not found in JSON.')
}

# Read probes BED -------------------------------------------------------------
probes <- fread(args$probes, header = F, stringsAsFactors = F)
colnames(probes) <- c('chr', 'start', 'end', 'id')
probes[, chr := gsub('chr', '', chr)]
probes <- makeGRangesFromDataFrame(probes, keep.extra.columns = T)
message('[', Sys.time(), '] Read --probes: ', args$probes)

# Read drivers ----------------------------------------------------------------
drivers <- fread(args$drivers, header = T, stringsAsFactors = F, 
                 select = c('gr_id', 'gene_id', 'gene_name', 'FILTER', 'tier',
                            'is_known_cancer'))
drivers <- unique(drivers)
drivers <- drivers[FILTER == 'PASS' & !is.na(tier)]
drivers <- drivers[, FILTER := NULL]
drivers <- drivers[, tier := NULL]
message('[', Sys.time(), '] Finished reading ', args$drivers)
if (nrow(drivers) == 0) {
  stop('[', Sys.time(), '] no significant (FILTER is PASS and tier is not ',
       'NA) driver genes is found in ', args$drivers, ' table.')
}
message('[', Sys.time(), '] Read --drivers: ', args$drivers)

# Read scanned for drivers genomic regions ------------------------------------
genomicRegions <- readBED12(args$genomic_regions)
message('[', Sys.time(), '] Read --genomic_regions: ', args$genomic_regions)
driverGenomicRegions_idx <- as.data.table(mcols(genomicRegions))
driverGenomicRegions_idx[, idx := 1:nrow(driverGenomicRegions_idx)]
driverGenomicRegions_idx <- merge(driverGenomicRegions_idx, drivers, 
                              by = c('gr_id', 'gene_id', 'gene_name'))
driverGenomicRegions <- genomicRegions[driverGenomicRegions_idx$idx]
driverGenomicRegions <- split(driverGenomicRegions,
                              f = driverGenomicRegions$gr_name)

# Overlap with probes -------------------------------------------------
overlapWithProbes <- lapply(driverGenomicRegions,
                            function(x) join_overlap_intersect(reduce(x),
                                                               probes))
overlapWithProbes_len <- lapply(overlapWithProbes, 
                                function(x) sum(width(x)))
totalLen <- lapply(driverGenomicRegions,
                   function(x) sum(width(reduce(x))))

overlapPercentage <- data.table(gr_name = names(overlapWithProbes), 
                                total = unlist(totalLen[names(overlapWithProbes)]),
                                ovrl = unlist(overlapWithProbes_len[names(overlapWithProbes)]))
overlapPercentage[, perc := 100 * ovrl/total]
overlapPercentage <- merge(overlapPercentage, driverGenomicRegions_idx,
                           by = 'gr_name')
overlapPercentage <- overlapPercentage[order(-perc)]
overlapPercentage[, gr_id := factor(gr_id, unique(gr_id))]
overlapPercentage <- overlapPercentage[order(gr_id, -perc)]
overlapPercentage[, gene_plots_id := paste(gene_name, gr_id)] 
overlapPercentage[, gene_plots_id := factor(gene_plots_id, 
                                            unique(gene_plots_id))]

# Create plot background (rectangles to split genomic regions) ----------------
# color of plot labels
xAxisLabs <- formatAxisLabels(overlapPercentage, 'gene_plots_id', 'gene_name',
                              'is_known_cancer')
xAxisLabs$labels <- gsub("__and__", "&", xAxisLabs$labels)
# background plot to split different genomic regions from each other
plotBackground <- create_plot_background(overlapPercentage, xAxisLabs, 
                                         'gene_plots_id', 'perc', 
                                         direction = 'vertical')

# Plot bar plot of overlap with probes ----------------------------------------
OVERLAP_WITH_PROBES_BAR <- plotBackground + 
  geom_bar(data = overlapPercentage,
           mapping = aes(x = gene_plots_id, y = perc), fill = 'black',
           color = 'black', position = "stack", stat = "identity") +
  scale_x_discrete(drop = F, labels = xAxisLabs$labels) + 
  xlab('driver genomic element') + ylab('% overlap with WES probes') + 
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


# Output plot to file ---------------------------------------------------------
# w = 6, h = 4
if (args$output_type == 'pdf') {
  pdf(args$output, width = visualParams$indel_enrichment_plot_width, 
      height = visualParams$indel_enrichment_plot_heigth)
} else {
  png(args$output, width = visualParams$indel_enrichment_plot_width, 
      height = visualParams$indel_enrichment_plot_heigth)
}
OVERLAP_WITH_PROBES_BAR
dev.off()

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')