library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)

source('bin/custom_functions.R')

# Test input arguments --------------------------------------------------------
args <- list(inventory_analysis = 'data/inventory/inventory_analysis.csv',
             inventory_patients = 'data/inventory/inventory_patients.csv',
             driver_mutations = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Panlung--hg19.csv",
             cancer_subtype = 'Panlung', exclude_cnv = T, 
             fold_splicesites_in_coding = T)

CNV_VAR_TYPES <- c('amp', 'hom.del.')
LOW_CONFIDENCE <- c('only n. driver mutations')

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
driverMuts[, gr_type := ifelse(gr_id %in% coding_gr_id, 'coding', 'non-coding')]

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
  scale_color_manual(values = colorRampPalette(N_PATIENTS_COLOR_SCHEME)(length(N_breaks) - 2)) +
  customGgplot2Theme + labs(color = 'N. patients', size = 'N. patients') +
  guides(color = guide_legend(nrow = 1), size = guide_legend(nrow = 1)) +
  theme(legend.position = 'bottom', legend.direction = 'horizontal')

