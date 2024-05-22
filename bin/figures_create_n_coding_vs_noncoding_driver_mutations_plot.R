#!/usr/bin/env Rscript
# FILE:  bin/figures_create_n_coding_vs_noncoding_driver_mutations_plot.R -----
#
# DESCRIPTION: An R script which creates a plot showing number patients with a
# specific combination of number of coding and noncoding driver mutations.
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla  bin/figures_create_n_coding_vs_noncoding_driver_mutations_plot.R -h
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

mutationsHelp <- paste('Path to detected driver mutations file in the ',
                       'cancer subtype. Columns is_driver, gr_id, gene_id,',
                       'gene_name, var_type, key, participant_id, n_total,',
                       'n_driverMut_low, n_driverMut_mle, n_driverMut_high,',
                       'n_tums_w_driver_mut are needed.')
parser$add_argument("-d", "--driver_mutations", required = T, 
                    type = 'character', default = NULL, help = mutationsHelp)

analysisHelp <- paste('Path to inventory table containing details of the',
                      'future analysis to be conducted. Minimal columns:',
                      'tumor_subtype,', 'software,', 'gr_id,', 'gr_code,', 
                      'gr_file,', 'gr_upstr,', 'gr_downstr,', 'gr_genome,', 
                      'gr_excl_id,', 'gr_excl_code,', 'gr_excl_file,',
                      'gr_excl_upstr,', 'gr_excl_downstr,', 'gr_excl_genome,',
                      'blacklisted_codes.')
parser$add_argument("-a", "--inventory_analysis", required = T, 
                    type = 'character', help = analysisHelp)

patientsHelp <- paste('Path to table listing information about patients,',
                      'their cancer types and mutation files. Minimal',
                      'columns: participant_id, tumor_subtype,',
                      'participant_tumor_subtype, somatic_path,',
                      'somatic_genome, cohort_name.')
parser$add_argument("-p", "--inventory_patients", required = T, 
                    type = 'character', help = patientsHelp)

foldHelp <- paste0('Indicates, whether or not splice sites should be',
                   'considered as part of CDS for the counting.',
                   'Default: T.')
parser$add_argument("-f", "--fold_splicesites_in_coding", required = F,
                    default = 'T', choices = c('T', 'F'), type = 'character', 
                    help = foldHelp)

excludeCNVhelp <- paste0('Indicates, whether or not driver copy number',
                         'variantions should be excluded from counting.',
                         'Default: T.')
parser$add_argument("-e", "--exclude_cnv", required = F,
                    default = 'T', choices = c('T', 'F'), type = 'character', 
                    help = excludeCNVhelp)

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
args$exclude_cnv <- as.logical(args$exclude_cnv)

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

CNV_VAR_TYPES <- c('amp', 'hom.del.')
LOW_CONFIDENCE <- c('only n. driver mutations')

# Test input arguments --------------------------------------------------------
# args <- list(inventory_analysis = 'data/inventory/inventory_analysis.csv',
#              inventory_patients = 'data/inventory/inventory_patients.csv',
#              driver_mutations = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Panlung--hg19.csv",
#              cancer_subtype = 'Panlung', exclude_cnv = T, 
#              fold_splicesites_in_coding = T, 
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

# Read in patients inventory --------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients, 1)
patientsInv <- patientsInv[tumor_subtype %in% args$cancer_subtype]
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)

# Read driver mutations -------------------------------------------------------
colsToKeep <- c("is_driver", "gr_id", "gene_id", "gene_name", "key", 
                "var_type", "participant_id", "confidenceLvl",
                "prob_is_driver_mle")
driverMuts <- fread(args$driver_mutations, header = T, select = colsToKeep,
                    stringsAsFactors = F)
driverMuts <- unique(driverMuts)
driverMuts <- driverMuts[is_driver == T]
driverMuts[, is_driver := NULL]
message('[', Sys.time(), '] Read --driver_mutations: ', args$driver_mutations)

# Read in analysis inventory --------------------------------------------------
analysisInv <- readAnalysisInventory(args$inventory_analysis)
message('[', Sys.time(), '] Read --inventory_analysis: ', 
        args$inventory_analysis)
analysisInv <- analysisInv[tumor_subtype == args$cancer_subtype]

# Fold splice sites into coding regions, if requested -------------------------
# assign coding genomic regions - anything which contains CDS as gr_code
coding_gr_id <- unique(analysisInv[gr_code == 'CDS']$gr_id)
message('[', Sys.time(), '] Following genomic regions: ', 
        paste0(coding_gr_id, collapse = ', '), ', will be considered as ',
        'coding.')

if (args$fold_splicesites_in_coding) {
  message('[', Sys.time(), '] Will fold splice sites into corresponding ',
          'coding driver genetic elements.')
  # assign splice site genomic regions - anything which contains ss as gr_code
  ss_gr_id <- unique(analysisInv[gr_code == 'ss']$gr_id)
  if (length(coding_gr_id) != 0 & length(ss_gr_id) != 0) {
    message('[', Sys.time(), '] Following genomic regions: ', 
            paste0(ss_gr_id, collapse = ', '), ', will be considered as ',
            'containing splice sites.')
    driverMuts <- foldSplicSiteDriversIntoCodingDrivers(coding_gr_id, ss_gr_id,
                                                        driverMuts)
  } else {
    message('[', Sys.time(), '] Did not find either coding or splice site ',
            'regions. Can not fold splice site drivers into coding ones.')
  }
}

noncoding_gr_id <- setdiff(unique(analysisInv$gr_id), coding_gr_id)
msg <- paste('[', Sys.time(), '] Following genomic regions:', 
             paste0(noncoding_gr_id, collapse = ', '), ', will be considered',
             'as noncoding.')
if (args$fold_splicesites_in_coding) {
  msg <- paste(msg, "Driver genomic elements considered as containing splice",
               "sites will be treated as noncoding for genes which do not",
               "have coding part detected as driver.")
}
message(msg)

# Set driver genomic element type (coding vs noncoding) -----------------------
driverMuts[, gr_type := ifelse(gr_id %in% coding_gr_id, 'coding', 
                               'non-coding')]

# Count number of SNV/small indel driver mutations per patient ----------------
# count how many RESOLVED driver SNVs/small indels per patient
nDriverSNVresolved <- driverMuts[!var_type %in% CNV_VAR_TYPES]
nDriverSNVresolved <- nDriverSNVresolved[!confidenceLvl %in% LOW_CONFIDENCE & 
                                           confidenceLvl != 'passenger']
nDriverSNVresolved <- nDriverSNVresolved[,.(length(unique(key))), 
                                         by = .(participant_id, gr_type)]
setnames(nDriverSNVresolved, 'V1', 'n_resolved')
nDriverSNVresolved[, struct_type := 'SNVs/small indels']

# next, count how many unresolved driver mutations per patient we have
nDriverSNVunresolved <- driverMuts[!var_type %in% CNV_VAR_TYPES]
nDriverSNVunresolved <- nDriverSNVunresolved[confidenceLvl == LOW_CONFIDENCE]
nDriverSNVunresolved <- nDriverSNVunresolved[,.(participant_id, gr_id, 
                                                gene_name, gene_id, key,
                                                prob_is_driver_mle, gr_type)]
nDriverSNVunresolved <- nDriverSNVunresolved[!is.na(prob_is_driver_mle)]
nDriverSNVunresolved <- nDriverSNVunresolved[,.(sum(prob_is_driver_mle)),
                                             by = .(participant_id, gr_type)]
nDriverSNVunresolved[, V1 := round(V1)]
nDriverSNVunresolved <- nDriverSNVunresolved[V1 > 0]
setnames(nDriverSNVunresolved, 'V1', 'n_unresolved')
nDriverSNVunresolved[, struct_type := 'SNVs/small indels']

nDriverSNV <- merge(nDriverSNVresolved, nDriverSNVunresolved, 
                    by = c('participant_id', 'gr_type', 'struct_type'),
                    all = T)
nDriverSNV[is.na(n_resolved)]$n_resolved <- 0
nDriverSNV[is.na(n_unresolved)]$n_unresolved <- 0

# in case we do have CNA - handle them
# FIX CNA COUNTING! NEED cn_reg_id, put it in key

# Combine SNV/small indel and CNA driver mutations ----------------------------
nDriverMuts <- copy(nDriverSNV)
nDriverMuts[, N := n_resolved + n_unresolved]

# Add patients where no driver mutations were found --------------------------
patientsWithoutDrivers <- setdiff(unique(patientsInv$participant_id),
                                  unique(nDriverMuts$participant_id))

if (length(patientsWithoutDrivers) > 0) {
  patientsWithoutDrivers <- expand.grid(participant_id = patientsWithoutDrivers,
                                        gr_type = c('coding', 'non-coding'),
                                        struct_type = unique(nDriverMuts$struct_type))
  patientsWithoutDrivers <- as.data.table(patientsWithoutDrivers)
  patientsWithoutDrivers <- cbind(patientsWithoutDrivers,
                                  n_resolved = 0, n_unresolved = 0, N = 0)
  nDriverMuts <- rbind(nDriverMuts, patientsWithoutDrivers)
}

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