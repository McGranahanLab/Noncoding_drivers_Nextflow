library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)

# Visuals setup ---------------------------------------------------------------
font_config <- list('family' = 'Helvetica', 'base' = 6, 'text' = 6, 
                    'axis_title' = 8, 'plot_title' = 10, 'legend_text' = 6, 
                    'legend_title' = 8)

customGgplot2Theme <- list(
  theme_classic(base_size = font_config$base, 
                base_family = font_config$family) +
    theme(axis.line = element_line(colour = 'black', linewidth = 0.5,
                                   linetype = 'solid'),
          axis.ticks = element_line(colour = 'black', linewidth = 0.5,
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

# Test input arguments --------------------------------------------------------
args <- list(cancer_subtype_1 = 'Adenocarcinoma',
             cancer_subtype_2 = 'Squamous_cell',
             driver_mutations_1 = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Adenocarcinoma--hg19.csv",
             driver_mutations_2 = "completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Squamous_cell--hg19.csv",
             subtype_specificity_file = 'completed_runs/2023_12_25/results/tables/subtype_specificity/subtypeSpecificity---hg19.csv')

NON_CNV_VAR_TYPES = c("subs", "mis", "non", "indels")

# Read results of subtype specificity analysis --------------------------------
subtypeSpec <- fread(args$subtype_specificity_file, header = T, 
                     stringsAsFactors = F, 
                     select = c('tumor_subtype_1', 'tumor_subtype_2', 'gr_id',
                                'gene_id', 'gene_name', 'tumor_subtype_spec',
                                'specificity_mode', 'is_driver_ts_1',
                                'is_driver_ts_2'))
subtypeSpec <- subtypeSpec[tumor_subtype_1 %in% args[c('cancer_subtype_1', 
                                                       'cancer_subtype_2')] &
                             tumor_subtype_2 %in% args[c('cancer_subtype_1', 
                                                         'cancer_subtype_2')]]
subtypeSpec <- subtypeSpec[is_driver_ts_1 == T | is_driver_ts_2 == T]

driversInSubtype1 <- subtypeSpec[is_driver_ts_1 == T]
driversInSubtype1 <- driversInSubtype1[,.(gr_id, gene_id, gene_name, 
                                          tumor_subtype_1)]
setnames(driversInSubtype1, 'tumor_subtype_1', 'tumor_subtype')
driversInSubtype2 <- subtypeSpec[is_driver_ts_2 == T]
driversInSubtype2 <- driversInSubtype2[,.(gr_id, gene_id, gene_name, 
                                          tumor_subtype_2)]
setnames(driversInSubtype2, 'tumor_subtype_2', 'tumor_subtype')
driversPerSubtype <- unique(rbind(driversInSubtype1, driversInSubtype2))
driversPerSubtype[, is_driver := T]
rm(driversPerSubtype1, driversPerSubtype2)

subtypeSpec <- subtypeSpec[,.(gr_id, gene_id, gene_name, tumor_subtype_spec,
                              specificity_mode)]
subtypeSpec <- unique(subtypeSpec)

# Read driver mutations -------------------------------------------------------
readAndSumUpDriverMutations <- function(filePath, tumor_subtype) {
  colsToKeep <- c('gr_id', 'gene_id', 'gene_name', "var_type", "n_total", 
                  'n_driverMut_low', 'n_driverMut_mle', 'n_driverMut_high')
  driverMuts <- fread(filePath, header = T, select = colsToKeep,
                      stringsAsFactors = F)
  driverMuts <- driverMuts[var_type %in% NON_CNV_VAR_TYPES]
  driverMuts <- unique(driverMuts)
  driverMuts[, n_total := sum(n_total), by = .(gr_id, gene_id, gene_name)]
  driverMuts[, n_driverMut_low := sum(n_driverMut_low),
               by = .(gr_id, gene_id, gene_name)]
  driverMuts[, n_driverMut_mle := sum(n_driverMut_mle),
               by = .(gr_id, gene_id, gene_name)]
  driverMuts[, n_driverMut_high := sum(n_driverMut_high),
               by = .(gr_id, gene_id, gene_name)]
  driverMuts[, var_type := NULL]
  driverMuts <- unique(driverMuts)
  driverMuts[, tumor_subtype := tumor_subtype]
  driverMuts
}

excessMuts <- rbind(readAndSumUpDriverMutations(args$driver_mutations_1,
                                                args$cancer_subtype_1),
                    readAndSumUpDriverMutations(args$driver_mutations_2,
                                                args$cancer_subtype_2))
excessMuts[, perc_driverMut_mle := 100 * n_driverMut_mle/n_total]
excessMuts[, perc_driverMut_low := 100 * n_driverMut_low/n_total]
excessMuts[, perc_driverMut_high := 100 * n_driverMut_high/n_total]
excessMuts[, n_total := NULL]

excessMuts <- merge(excessMuts, driversPerSubtype, all = T,
                    by = c("tumor_subtype", "gr_id", "gene_id", "gene_name"))

# Prepare data for plotting ---------------------------------------------------
MODE_UPD <- "found as driver in a tumor subtype,\nbut is preferential/specific to the other subtype"

subtypeSpec <- merge(excessMuts, subtypeSpec, all = T,
                     by = c('gr_id', 'gene_id', 'gene_name'))

subtypeSpec[tumor_subtype_spec != tumor_subtype & 
              tumor_subtype_spec %in% 
              c(args$cancer_subtype_1, args$cancer_subtype_2)]$specificity_mode <- MODE_UPD
subtypeSpec[tumor_subtype_spec != tumor_subtype & 
              tumor_subtype_spec %in% 
              c(args$cancer_subtype_1, args$cancer_subtype_2)]$tumor_subtype_spec <- ""
subtypeSpec[, specificity := paste0(tumor_subtype_spec, '-', specificity_mode)]
subtypeSpec[, specificity := gsub('^-', '', specificity)]

# sort so that tumor subtype with more specific/preferential drivers would 
# appear first
subtypeOrder <- subtypeSpec[tumor_subtype_spec %in% c(args$cancer_subtype_1,
                                                      args$cancer_subtype_2)]

if (nrow(subtypeOrder) > 0) {
  subtypeOrder <- subtypeOrder[,.N, by = tumor_subtype_spec]
  subtypeOrder <- subtypeOrder[order(-N)]
  subtypeOrder <- subtypeOrder$tumor_subtype_spec
} else {
  subtypeOrder <- sort(subtypeOrder, descreasing = T)
}
specificity_lvls <- c(paste0(subtypeOrder[2], "-specific"),
                      paste0(subtypeOrder[2], "-preferential"),
                      "common", "not enough evidence",
                      paste0(subtypeOrder[1], "-preferential"),
                      paste0(subtypeOrder[1], "-specific"),
                      MODE_UPD)
subtypeSpec[, specificity := factor(specificity, specificity_lvls, 
                                    ordered = T)]
subtypeSpec[tumor_subtype_spec == ""]$tumor_subtype_spec <- subtypeSpec[tumor_subtype_spec == ""]$specificity_mode
subtypeSpec[, tumor_subtype_spec := factor(tumor_subtype_spec, 
                                           c(subtypeOrder[2], 
                                             "common", 
                                             "not enough evidence",
                                             subtypeOrder[1],
                                             MODE_UPD))]
subtypeSpec[, specificity_mode := NULL]

subtypeSpec <- subtypeSpec[order(-perc_driverMut_mle)]
subtypeSpec[, gr_id := factor(gr_id, rev(unique(gr_id)))]
subtypeSpec <- subtypeSpec[order(gr_id, tumor_subtype_spec, perc_driverMut_mle)]
subtypeSpec[, gene_plots_id := paste(gene_name, gr_id)]
subtypeSpec[, gene_plots_id := factor(gene_plots_id, unique(gene_plots_id))]

# Plot ------------------------------------------------------------------------
leftPlot <- ggplot() +
  geom_bar(data = subtypeSpec[tumor_subtype == subtypeOrder[1]],
           mapping = aes(y = gene_plots_id, x = -perc_driverMut_mle,
                         fill = specificity, color = is_driver,
                         alpha = is_driver), 
           stat = 'identity', linewidth = 0.5) +
  geom_point(data = subtypeSpec[tumor_subtype == subtypeOrder[1]],
             mapping = aes(y = gene_plots_id, x = -perc_driverMut_mle),
             color = 'black', size = 1) +
  geom_errorbar(data = subtypeSpec[tumor_subtype == subtypeOrder[1]],
                mapping = aes(y = gene_plots_id, xmin = -perc_driverMut_low,
                              xmax = -perc_driverMut_high),
                color = 'black', linewidth = 0.5, width = 0.25) +
  scale_color_manual(values = c('FALSE' = '#CCCCCC', 'TRUE' = 'black')) +
  scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 1)) +
  scale_x_continuous(breaks = seq(0, -100, by = -25), 
                     labels = seq(0, 100, by = 25)) +
  scale_y_discrete(drop = F, position = 'right') +
  xlab('% excess of mutations') + ylab('') + 
  labs(color = 'is driver\n in tumor subtype', 
       alpha = 'is driver\n in tumor subtype') +
  customGgplot2Theme + ggtitle(subtypeOrder[1]) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        axis.text.y = element_blank(), axis.line.y = element_blank(),
        plot.title = element_text(size = 8))

rightPlot <- ggplot() +
  geom_bar(data = subtypeSpec[tumor_subtype == subtypeOrder[2]],
           mapping = aes(y = gene_plots_id, x = perc_driverMut_mle,
                         fill = specificity, color = is_driver,
                         alpha = is_driver), 
           stat = 'identity', linewidth = 0.5) +
  geom_point(data = subtypeSpec[tumor_subtype == subtypeOrder[2]],
             mapping = aes(y = gene_plots_id, x = perc_driverMut_mle),
             color = 'black', size = 1) +
  geom_errorbar(data = subtypeSpec[tumor_subtype == subtypeOrder[2]],
                mapping = aes(y = gene_plots_id, xmin = perc_driverMut_low,
                              xmax = perc_driverMut_high),
                color = 'black', linewidth = 0.5, width = 0.25) +
  scale_color_manual(values = c('FALSE' = '#CCCCCC', 'TRUE' = 'black')) +
  scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 1)) +
  scale_x_continuous(breaks = seq(0, 100, by = 25), 
                     labels = seq(0, 100, by = 25)) +
  scale_y_discrete(drop = F) +
  xlab('% excess of mutations') + ylab('') + 
  labs(color = 'is driver\n in tumor subtype', 
       alpha = 'is driver\n in tumor subtype') +
  customGgplot2Theme + ggtitle(subtypeOrder[2]) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 8))

ggarrange(leftPlot + theme(legend.position = 'none'),
          rightPlot + theme(legend.position = 'none'), 
          nrow = 1, align = 'h', widths = c(3, 4))