#!/usr/bin/env Rscript
# FILE:  figures_create_n_coding_vs_noncoding_driver_mutations_plot.R ---------
#
# DESCRIPTION: An R script which creates a plot showing number patients with a
# specific combination of number of coding and noncoding driver mutations.
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla figures_create_n_coding_vs_noncoding_driver_mutations_plot.R -h
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
parser <- ArgumentParser(prog = 'figures_create_n_coding_vs_noncoding_driver_mutations_plot.R')

subtypeHelp <- paste('A name of the tumor subtype in which drivers were',
                     'detected.')
parser$add_argument("-c", "--cancer_subtype", required = T, 
                    type = 'character', help = subtypeHelp)

mutationsHelp <- paste('Path to file with the counts of the detected driver',
                       'mutations in the cancer subtype. Columns gr_type,',
                       'struct_type, participant_id, n_resolved,',
                       'n_unresolved, N are needed.')
parser$add_argument("-d", "--n_driver_mutations", required = T, 
                    type = 'character', default = NULL, help = mutationsHelp)

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

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Test input arguments --------------------------------------------------------
# args <- list(cancer_subtype = 'Panlung',
#              n_driver_mutations = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutationsCountsPerPatient-Panlung--hg19.csv",
#              visuals_json = 'data/visual_parameters.json',
#              output_type = 'pdf', output = 'test.pdf')

# Read visuals JSON -----------------------------------------------------------
ESSENTIAL_VISUAL_NAMES <- c('ggplot2_theme',
                            'colors_divergent_palette',
                            'n_coding_vs_noncoding_driver_muts_plot_width', 
                            'n_coding_vs_noncoding_driver_muts_plot_heigth')

visualParams <- readJsonWithVisualParameters(args$visuals_json)
message('[', Sys.time(), '] Read --visuals_json: ', args$visuals_json)

notFoundVisuals <- setdiff(ESSENTIAL_VISUAL_NAMES, names(visualParams))
if (length(notFoundVisuals)) {
  stop('[', Sys.time(), '] Following visuals: ', 
       paste(notFoundVisuals, collapse = ', '), ' not found in JSON.')
}

# Read number of driver mutations per patient ---------------------------------
colsToKeep <- c("gr_type", 'struct_type', 'participant_id', 'n_resolved',
                'n_unresolved', "N")
nDriverMuts <- fread(args$n_driver_mutations, header = T, select = colsToKeep,
                     stringsAsFactors = F)
message('[', Sys.time(), '] Read --n_driver_mutations: ',
        args$n_driver_mutations)

# Prepare for plotting --------------------------------------------------------
neutralColorIdx <- median(1:length(visualParams$colors_divergent_palette))
nColors <- visualParams$colors_divergent_palette
hotPalette <- visualParams$colors_divergent_palette[(neutralColorIdx + 1):length(nColors)]

nDriverMutsWide <- nDriverMuts[,.(participant_id, gr_type, N)]
nDriverMutsWide[, N := sum(N), by = .(participant_id, gr_type)]
nDriverMutsWide <- unique(nDriverMutsWide)
nDriverMutsWide <- dcast(nDriverMutsWide, participant_id ~ gr_type, 
                         value.var = "N")
nDriverMutsWide[is.na(nDriverMutsWide)] <- 0

nDriverMutsPlotDT <- nDriverMutsWide[,.(length(unique(participant_id))),
                                     by = .(coding, `non-coding`)]
setnames(nDriverMutsPlotDT, 'V1', 'N')

N_breaks <- c(0, 5, 10)
if (max(nDriverMutsPlotDT$N) > 25) {
  N_breaks <- c(N_breaks, seq(25, round_any(max(nDriverMutsPlotDT$N), 50),
                              by = 25))
} else {
  N_breaks <- c(N_breaks, 20)
}
if (max(N_breaks) < max(nDriverMutsPlotDT$N)) {
  N_breaks <- c(N_breaks, max(nDriverMutsPlotDT$N))
}
nDriverMutsPlotDT[, N_bin := cut(N, breaks = N_breaks)]

nCodVsNC_PLOT <- ggplot(nDriverMutsPlotDT, 
                        aes(x = coding, y = `non-coding`, size = N_bin, 
                            color = N_bin)) +
  geom_text(data = nDriverMutsPlotDT[N >= 5], 
            inherit.aes = F, vjust = -1, size = 3, 
            mapping = aes(x = coding, y = `non-coding`, label = N)) +
  geom_point() + coord_fixed() +
  scale_x_continuous('N. coding driver mutations', 
                     breaks = seq(0, 
                                  max(c(nDriverMutsPlotDT$coding,
                                        nDriverMutsPlotDT$`non-coding`))+1)) +
  scale_y_continuous('N. non-coding driver mutations', 
                     breaks = seq(0, 
                                  max(c(nDriverMutsPlotDT$coding,
                                        nDriverMutsPlotDT$`non-coding`))+1)) +
  scale_color_manual(values = colorRampPalette(hotPalette)(length(N_breaks) - 2)) +
  visualParams$ggplot2_theme + labs(color = 'N. patients', size = 'N. patients') +
  guides(color = guide_legend(ncol = 1), size = guide_legend(nrow = 1)) +
  theme(legend.position = 'right', legend.direction = 'vertical')

# Output plot to file ---------------------------------------------------------
if (args$output_type == 'pdf') {
  pdf(args$output, width = visualParams$n_coding_vs_noncoding_driver_muts_plot_width, 
      height = visualParams$n_coding_vs_noncoding_driver_muts_plot_heigth)
} else {
  png(args$output, width = visualParams$n_coding_vs_noncoding_driver_muts_plot_width, 
      height = visualParams$n_coding_vs_noncoding_driver_muts_plot_heigth)
}
nCodVsNC_PLOT
dev.off()

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')