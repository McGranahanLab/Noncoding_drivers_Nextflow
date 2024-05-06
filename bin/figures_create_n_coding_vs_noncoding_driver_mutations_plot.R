library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)

source('bin/custom_functions.R')

# Functions -------------------------------------------------------------------
summarizeDriverMutationsPerPatient <- function(filePath, exclude_cnv = F) {
  colsToKeep <- c("is_driver", "gr_id", "gene_id", "gene_name", "var_type",
                  "participant_id", "n_driverMut_mle")
  driverMuts <- fread(filePath, header = T, select = colsToKeep,
                      stringsAsFactors = F)
  driverMuts <- unique(driverMuts)
  driverMuts <- driverMuts[is_driver == T]
  
  if (exclude_cnv) {
    driverMuts <- driverMuts[!var_type %in% c('amp', 'hom.del.')]
    message('[', Sys.time(), '] Removed amplifications and homozygous ',
            'deletions from number of driver mutations per patient.')
  }
  
  nDriverMuts <- driverMuts[,.(sum(n_driverMut_mle)), 
                            by = .(participant_id, gr_id)]
  setnames(nDriverMuts, 'V1', 'n_driver_muts')
  nDriverMuts
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
                                          linewidth = 0.25, linetype = 2),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm")))

DIVIDER_COLOR <-  '#eed9c4'
N_PATIENTS_COLOR_SCHEME <- c('#90C987', '#CAE0AB', '#F7F056', '#F6C141',
                             '#F1932D','#E8601C', '#DC050C', '#42150A')

# Test input arguments --------------------------------------------------------
args <- list(driver_mutations = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Panlung--hg19.csv",
             exclude_cnv = T, fold_splice_sites_into_coding = T)

CODING_GR <- c("CDS", "ss")
CNV_VAR_TYPES <- c('amp', 'hom.del.')
LOW_CONFIDENCE <- c('only n. driver mutations')

# Read driver mutations -------------------------------------------------------
colsToKeep <- c("is_driver", "gr_id", "gene_id", "gene_name", "key", 
                "var_type", "participant_id", "confidenceLvl",
                "prob_is_driver_mle")
driverMuts <- fread(args$driver_mutations, header = T, select = colsToKeep,
                    stringsAsFactors = F)
driverMuts <- unique(driverMuts)
driverMuts <- driverMuts[is_driver == T]
driverMuts[, is_driver := NULL]

# Count number of SNV/small indel driver mutations per patient ----------------
# count how many RESOLVED driver SNVs/small indels per patient
nDriverSNVresolved <- driverMuts[!var_type %in% CNV_VAR_TYPES]
nDriverSNVresolved <- nDriverSNVresolved[!confidenceLvl %in% LOW_CONFIDENCE & 
                                           confidenceLvl != 'passenger']
nDriverSNVresolved[, gr_type := ifelse(gr_id == 'CDS', 'coding', 'non-coding')]
# THIS NEEDS TO BE FIXED
if (args$fold_splice_sites_into_coding) {
  nDriverSNVresolved[gr_id == 'ss']$gr_type <- 'coding'
}
nDriverSNVresolved <- nDriverSNVresolved[,.(length(unique(key))), 
                                         by = .(participant_id, gr_type)]
setnames(nDriverSNVresolved, 'V1', 'n_resolved')
nDriverSNVresolved[, struct_type := 'SNVs/small indels']

# next, count how many unresolved driver mutations per patient we have
nDriverSNVunresolved <- driverMuts[!var_type %in% CNV_VAR_TYPES]
nDriverSNVunresolved <- nDriverSNVunresolved[confidenceLvl == LOW_CONFIDENCE]
nDriverSNVunresolved <- nDriverSNVunresolved[,.(participant_id, gr_id, 
                                                gene_name, gene_id, key,
                                                prob_is_driver_mle)]
nDriverSNVunresolved <- nDriverSNVunresolved[!is.na(prob_is_driver_mle)]
nDriverSNVunresolved[, gr_type := ifelse(gr_id == 'CDS', 'coding', 
                                         'non-coding')]
# THIS NEEDS TO BE FIXED
if (args$fold_splice_sites_into_coding) {
  nDriverSNVunresolved[gr_id == 'ss']$gr_type <- 'coding'
}
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
patientsWithoutDrivers <- setdiff(unique(mleMuts$participant_id),
                                  unique(nDriverMuts$participant_id))

if (length(patientsWithoutDrivers) > 0) {
  patientsWithoutDrivers <- expand.grid(participant_id = patientsWithoutDrivers,
                                        gr_type = c('coding', 'non-coding'),
                                        struct_type = unique(n_driver_muts$struct_type))
  patientsWithoutDrivers <- as.data.table(patientsWithoutDrivers)
  patientsWithoutDrivers <- cbind(patientsWithoutDrivers,
                                  n_resolved = 0, n_unresolved = 0)
  nDriverMuts <- rbind(nDriverMuts, patientsWithoutDrivers)
}

# Prepare for plotting --------------------------------------------------------
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
            inherit.aes = F, vjust = -1, size = 5, 
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
  scale_color_manual(values = colorRampPalette(N_PATIENTS_COLOR_SCHEME)(length(N_breaks) - 2)) +
  customGgplot2Theme + labs(color = 'N. patients', size = 'N. patients') +
  guides(color = guide_legend(nrow = 1), size = guide_legend(nrow = 1)) +
  theme(legend.position = 'bottom', legend.direction = 'horizontal')

