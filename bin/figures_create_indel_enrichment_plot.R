library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)

source('bin/custom_functions.R')

# Test input arguments --------------------------------------------------------
args <- list(cancer_subtype = 'Panlung',
             drivers = "completed_runs/2023_12_25/results/tables/drivers/drivers-Panlung--hg19.csv",
             allowed_filter_values = list("PASS", "INDEL, 2-5bp"),
             variant_categories_counts = "completed_runs/2023_12_25/results/mut_rates/varCatEnrich-Panlung--hg19.csv",
             categoty_of_interest = "INDEL, 2-5bp", p_abj = 0.05)

VAR_CAT_LEVELS <- rev(c("SNP", "INDEL, 1bp", "INDEL, 2-5bp", "INDEL, 6-9bp",
                        "INDEL, 10+bp", "MNP"))

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
xAxisLabs <- format_xLabels(varCatEnrich, 'gene_plots_id', 'gene_name',
                            'is_known_cancer')
xAxisLabs$labels <- gsub("__and__", "&", xAxisLabs$labels)
# background plot to split different genomic regions from each other
plotBackground <- create_plot_background(varCatEnrich, xAxisLabs, 
                                         'gene_plots_id', 'perc_cat', 
                                         direction = 'vertical')

# Plot bar plot enrichment of variant categories  -----------------------------
categoriesPallete <- colorRampPalette(tealPalette)
categoriesPallete <- categoriesPallete(length(unique(varCatEnrich$var_cat)))

VAR_CT_ENRICH_BAR <- plotBackground + 
  geom_bar(data = varCatEnrich,
           mapping = aes(x = gene_plots_id, y = perc_cat, fill = var_cat),
           position = "stack", stat = "identity") +
  geom_text(data = varCatEnrich[binomPadj <= args$p_abj & 
                                  var_cat %in% args$categoty_of_interest],
            mapping = aes(x = gene_plots_id, y = 101, label = "*")) + 
  scale_fill_manual(values = categoriesPallete) +
  scale_x_discrete(drop = F, labels = xAxisLabs$labels) + 
  xlab('driver genomic element') + ylab('% of variants') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  customGgplot2Theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   family = customGgplot2Theme[[1]]$axis.text$family,
                                   size = customGgplot2Theme[[1]]$axis.text$size,
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
