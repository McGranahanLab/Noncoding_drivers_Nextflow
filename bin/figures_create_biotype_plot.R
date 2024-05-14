library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)

source('bin/custom_functions.R')

# Functions: misc -------------------------------------------------------------
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

# Test input arguments --------------------------------------------------------
args <- list(cancer_subtype = 'Panlung',
             drivers_biotyped = "completed_runs/2023_12_25/results/tables/drivers_biotyped/driversBiotyped-Panlung--hg19.csv",
             min_n_tumors = 10)
# add args for 0.5 and 0.33
# min_n_tumors is not needed? Because it was already done

# Read driver biotype ---------------------------------------------------------
biotypes <- fread(args$drivers_biotyped, header = T, stringsAsFactors = F)
biotypes[, gene_plots_id := paste(gene_name, gr_id)]
# add filtering
biotypes[, perc_group := 100 * perc_group]
biotypes[, simplified_known_biotype := sapply(known_cancer_biotype, 
                                              simplifyBiotype)]

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
  geom_tile(data = unique(biotypes[,.(gene_plots_id, simplified_known_biotype)]),
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
  geom_hline(yintercept = c(33, 50), linetype = 2, linewidth = 0.5, 
             color = '#999999')

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
  geom_hline(yintercept = c(33, 50), linetype = 2, linewidth = 0.5, 
             color = '#999999')

# Aggange plots ---------------------------------------------------------------
LEGEND_MUTS_BAR <- extract_legend(MUTS_BAR)
LEGEND_CN_BAR <- extract_legend(CNA_BAR)
LEGEND_INFERRED_BIOTYPE <- extract_legend(INFERRED_TILE) 
LEGENDS <- ggarrange(LEGEND_INFERRED_BIOTYPE, LEGEND_MUTS_BAR, NULL,
                     LEGEND_CN_BAR, 
                     nrow = 2, ncol = 2, heights = c(1, 1), widths = c(5, 12))

# w = 6, h = 4
ggarrange(plotlist = list(INFERRED_TILE + theme(legend.position = "none"),
                          KNOWN_TILE + theme(legend.position = "none"),
                          MUTS_BAR + theme(legend.position = "none"),
                          CNA_BAR + theme(legend.position = "none"),
                          LEGENDS),
          heights = c(0.025, 0.025, 0.3, 0.5, 0.1), 
          align = 'v', ncol = 1)