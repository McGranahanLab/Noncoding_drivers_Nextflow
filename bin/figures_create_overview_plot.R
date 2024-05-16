#!/usr/bin/env Rscript
# FILE: figures_create_overview_plot.R ----------------------------------------
#
# DESCRIPTION: An R script which creates an overview plot of detected de-novo
# drivers.
# USAGE: 
# OPTIONS: Run 
#          Rscript --vanilla figures_create_overview_plot.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: R v4.1.0, data.table, ggplot2, ggpubr, plyr
# BUGS: --
# NOTES:  
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  01.02.2023
# REVISION: 15.05.2024

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

# Functions : overview of all detected drivers --------------------------------
#' get_bar_transparency
#' @description 
#' @author Maria Litovchenko
#' @param plotDT data table which holds data for plotting. Have to have columns
#'        tumor_subtype, tier, participant_tumor_subtype and a column with the 
#'        name identical to axisOrderCol
#' @param pancanCode string encoding tumor subtype code for pancancer analysis
#' @param axisOrderCol name of the column which contains IDs of genes to plot and
#'        by which they will be ordered on the plot
#' @note Driver status of the gene will control the opacity of the bars. 
#'       1) If gene is a driver in a particular tumor type(s), but not in 
#'       pancancer analysis, the bar for it (them) should be opaque and for the
#'       rest of the tumor types - transparent.
#'       2) If gene is a driver in pancancer analysis, all bars should be 
#'       opaque. However, in the plotDT table now we do not have a breakdown of 
#'       number of patients per tumor type in case gene is a driver in the 
#'       pancancer analysis. Let's add the breakdown by retrieving number of 
#'       patients from plotDT They will be lines which may not have significant
#'       p-values in case a gene was found significant only!) in pancancer 
#'       analysis.
get_bar_transparency <- function(plotDT, pancanCode = NULL, axisOrderCol) {
  result <- copy(plotDT)
  # transparent, if significant in individual tumor subtypes
  result[, transparent := ifelse(tier == 0, 1, 0)]
  # all tumor_subtypes opaque, if significant in pancan analysis
  if (!is.null(pancanCode)) {
    signPancan <- result[tumor_subtype == pancanCode & tier == 1]
    signPancan <- unlist(signPancan[, axisOrderCol, with = F])
    signPancan <- unlist(result[, axisOrderCol, with = F]) %in% signPancan
    result[signPancan]$transparent <- 0
    # no need to keep pancan in the table anymore
    result <- split(result, by = axisOrderCol, drop = T)
    result <- lapply(result, 
                     function(x) x[tumor_subtype != pancanCode | 
                                     (tumor_subtype == pancanCode & 
                                        !participant_tumor_subtype %in% tumor_subtype)])
    result <- do.call(rbind, result)
    result[, tumor_subtype := NULL]
    setnames(result, 'participant_tumor_subtype', 'tumor_subtype')
    if (class(plotDT$tumor_subtype) == 'factor') {
      result[, tumor_subtype := factor(tumor_subtype, 
                                       levels = levels(plotDT$tumor_subtype))]
    }
  }
  result
}

#' create_tumortype_colorpalette
#' @description Creates color palette for tumor subtypes
#' @author Maria Litovchenko
#' @param plotDT data table with data for plotting. Have to have columns 
#'        tumor_subtype
#' @param colorpalette character vector, color palette submitted by user
#' @return named character vector with colors for tumor types. Colors are 
#'.        repeated for the case then driver is significant in that tumor type
#'         (--0) and not (--1)
create_tumortype_colorpalette <- function(plotDT, colorpalette) {
  if (is.null(levels(plotDT$tumor_subtype))) {
    tumor_types_lvl <- as.character(unique(plotDT$tumor_subtype))
  } else {
    tumor_types_lvl <- levels(plotDT$tumor_subtype)
  }
  
  if (is.null(colorpalette)) {
    colorpalette <- palette(rainbow(length(tumor_types_lvl)))
    names(colorpalette) <- tumor_types_lvl
  } else {
    if (is.null(names(colorpalette))) {
      colorpalette <- colorRampPalette(colorpalette)(length(tumor_types_lvl))
      names(colorpalette) <- tumor_types_lvl
    } else {
      if (any(!unique(plotDT$tumor_subtype) %in% names(colorpalette))) {
        stop('[', Sys.time(), '] No colors provided for following tumor ',
             'subtypes: ', paste0(unique(plotDT$tumor_subtype), 
                                  collapse = ', '))
      }
      colorpalette <- colorpalette[tumor_types_lvl]
    }
  }

  result <- rep(colorpalette, 2)
  names(result) <- c(paste0(names(colorpalette), '--0'), 
                     paste0(names(colorpalette), '--1'))
  result
}

#' create_pvalue_colorpalette
#' @description Creates binned color palette for p-values
#' @author Maria Litovchenko
#' @param plotDT data table with columns colorByLog10
#' @param doLog10 whatever -log10 operation on colorBy column should be 
#'                performed
#' @param scaleLims limits of the scale
#' @param nColors number of colors to use in the color palette 
#' @param colorPalette character vector of colors
#' @return list with members colorPalette and colorBreaks
create_pvalue_colorpalette <- function(plotDT, doLog10, scaleLims, nColors,
                                       colorPalette) {
  # color breaks
  colorBreaks <- round(seq(0, max(plotDT$colorByLog10, na.rm = T), 
                           length.out = nColors + 1), 2)
  if (!is.null(scaleLims)) {
    colorBreaks <- round(seq(scaleLims[1], scaleLims[2],
                             length.out = nColors + 1), 2)
  }
  if (doLog10 == T) {
    colorBreaks <- round(colorBreaks)
  }
  # actual colors
  if (is.null(colorPalette)) {
    colorPalette <- colorRampPalette(c("grey0", "grey80"))(nColors)
    if (doLog10 == T) {
      colorPalette <- rev(colorPalette)
    }
  }
  result <- list('palette' = colorPalette, 'colorBreaks' = colorBreaks)
  result
}

#' barplotNpatient
#' @description Creates stacked barplot with number of patients per gene with 
#' histological subtypes on Y axis. In case several regions are given, they are
#' groupped and indicated with horizontal bar.
#' @author Maria Litovchenko
#' @param xLabelCol column name containing label to be on X axis
#' @param colorPalette named vector of color for different tumour type
#' @param xOrderCol column name containing order of X axis
#' @param annoCol column containing annotation, i.e. COSMIC 
#' @param annoColors named vector of colors (in hex) in which X axis labels 
#'                   will be colored
#' @param divider_color hex color which should be used for separating different
#'        genomic regions
#' @param ggplot2Theme ggplot2 custom theme
#' @return ggplot2: barplot with number of patients per gene with histological
#' subtypes on Y axis
barplotNpatient <- function(driversDT, xLabelCol, colorPalette = NULL,
                            pancanCode = NULL, xOrderCol = NULL, 
                            annoCol = NULL, annoColors = 'red',
                            divider_color = '#eed9c4',
                            ggplot2Theme = NULL) {
  if (is.null(xOrderCol)) {
    xOrderCol <- xLabelCol
  }
  colsToGet <- unique(c('tumor_subtype', 'gr_id', 'gene_id', 'gene_name',
                        'tier', 'nParts', 'participant_tumor_subtype', 
                        'FILTER', xOrderCol, xLabelCol, annoCol))
  dt <- driversDT[, colsToGet, with = F]
  # we will not display exact tier information
  dt[, tier := ifelse(is.na(tier), 0, 1)]
  
  # assign transparency values based on significance
  dt <- get_bar_transparency(dt, pancanCode, xOrderCol)
  
  # in order to put first tumor types for which gene was found a driver, and
  # then the other one for which not, we have to modify tumor_subtype by adding 
  # is_driver_in_tt to it.
  dt[, tumor_subtype_upd := paste0(tumor_subtype, '--', transparent)]
  tumor_subtype_upd_lvls <- sort(unique(dt$tumor_subtype))
  if (class(dt$tumor_subtype) == 'factor') {
    tumor_subtype_upd_lvls <- levels(dt$tumor_subtype)
  }
  tumor_subtype_upd_lvls <- c(paste0(rev(tumor_subtype_upd_lvls), '--1'), 
                              paste0(rev(tumor_subtype_upd_lvls), '--0'))
  dt[, tumor_subtype_upd := factor(tumor_subtype_upd, tumor_subtype_upd_lvls)]
  
  # tumor type tier palette
  tumorType_tier_palette <- create_tumortype_colorpalette(dt, colorPalette)
  # x axis labels & color
  xAxisLabs <- formatAxisLabels(dt, xOrderCol, xLabelCol, annoCol, 
                                highlightColor = annoColors)
  # transparency scales
  alphaVals <- c('yes' = 1, 'no' = 0.5)
  if (all(dt$transparent == 1)) {
    alphaVals <- c('no' = 0.5)
  }
  if (all(dt$transparent == 0)) {
    alphaVals <- c('yes' = 1)
  }
  
  # first of all, if there are several genomic regions, create shadowing
  result <- ggplot()
  Y_MAX <- NULL
  if (length(unique(dt$gr_id)) > 1) {
    # compute coordinates for future lines showing genome region 
    grShadeCoords <- get_gr_change_coords(plotDT = dt, labsDT = xAxisLabs$dt,
                                          axisOrderCol = xOrderCol, 
                                          yCol = 'nParts')
    grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                    nrow(grShadeCoords$gr_coords),
                                                    by = 2)]
    result <- result + 
      geom_rect(data = grShadeCoordsCut, fill = divider_color, inherit.aes = F,
                aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, 
                    ymax = Inf), alpha = 0.5) + 
      geom_segment(data = grShadeCoords$gr_coords, inherit.aes = F,
                   aes(x = start, xend = end, 
                       y =  grShadeCoords$lab_coords$y,
                       yend = grShadeCoords$lab_coords$y)) +
      geom_text(mapping = aes(x = x, y = 1.01 * y, label = gr, angle = 90, 
                              hjust = -0.01), inherit.aes = F,
                data = grShadeCoords$lab_coords, size = 2)
    Y_MAX <- 1.03 * unique(grShadeCoords$lab_coords$y)
  }
  # and if there are some filtered out drivers needed to be shown - shadowed 
  # areas for them
  if (any(dt$FILTER != "PASS")) {
    filteredOutShadeCoords <- get_filtered_out_gene_coords(xAxisLabs$dt, xOrderCol)
    filteredOutShadeCoordsCut <- filteredOutShadeCoords[seq(1, nrow(filteredOutShadeCoords), 
                                                            by = 2)]
    result <- result + 
      geom_rect(data = filteredOutShadeCoordsCut, fill = 'red', inherit.aes = F, 
                aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, 
                    ymax = Inf), alpha = 0.15)
  }
  
  result <- result + 
    geom_bar(data = dt, stat = 'identity',
             mapping = aes_string(x = xOrderCol, y = 'nParts', 
                                  alpha = 'transparent',
                                  fill = 'tumor_subtype_upd')) +
    xlab('genomic element') + ylab('N. participants') + 
    scale_fill_manual(values = tumorType_tier_palette, 
                      labels = gsub('--0', '', 
                                    grep('--0', names(tumorType_tier_palette),
                                         value = T)),
                      breaks = grep('--0', names(tumorType_tier_palette),
                                    value = T)) +
    scale_alpha_continuous(breaks = alphaVals, range = alphaVals, 
                           labels = names(alphaVals)) + 
    scale_x_discrete(labels = xAxisLabs$label, expand = c(0, 0)) +
    labs(fill = 'Histological\nsubtype', alpha = 'Found as driver')
  
  if (!is.null(ggplot2Theme)) {
    result <- result + ggplot2Theme
  }
  
  result <- result + 
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 90, size = 6,
                                     color = xAxisLabs$color))
  if (!is.null(Y_MAX)) {
    result <- result + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX))
  } else {
    result <- result + scale_y_continuous(expand = c(0, 0))
  }
  
  result
}

#' tilePlotPvalues
#' @description Plots p-value matrix.
#' @author Maria Litovchenko
#' @param pValsDT data table
#' @param xLabelCol column name containing label to be on X axis
#' @param yLabelCol column name containing label to be on Y axis
#' @param colorBy name of column containing values to be used as color value of
#'                cells
#' @param xOrderCol column name containing order of X axis
#' @param yOrderCol column name containing order of Y axis
#' @param annoCol column containing annotation, i.e. COSMIC
#' @param annoColors named vector of colors (in hex) in which X axis labels 
#'                   will be colored
#' @param divider_color hex color which should be used for separating different
#'        genomic regions
#' @param doLog10 whatever -log10 operation on colorBy column should be 
#'                performed
#' @param scaleLims 
#' @param nColors number of colors to use in the color palette 
#' @param colorPalette character vector of colors
#' @param colorPaletteName name of color legend
#' @param ggplot2Theme ggplot2 theme
tilePlotPvalues <- function(pValsDT, xLabelCol, yLabelCol, colorBy,
                            xOrderCol = NULL, yOrderCol = NULL, 
                            annoCol = NULL, annoColors = 'red', 
                            doLog10 = T, scaleLims = NULL, 
                            nColors = 5, colorPalette = NULL, 
                            divider_color = '#eed9c4',
                            ggplot2Theme = NULL) {
  if (!is.null(colorPalette) & nColors != length(colorPalette)) {
    stop('[', Sys.time(), '] Length of nColors and colorPalette should match')
  }
  dt <- copy(pValsDT)
  if (is.null(xOrderCol)) {
    xOrderCol <- xLabelCol
  }
  if (is.null(yOrderCol)) {
    yOrderCol <- yLabelCol
  }
  dt[, colorByLog10 := unlist(dt[, colorBy, with = F])]
  if (doLog10 == T) {
    dt[, colorByLog10 := -log10(colorByLog10 + 10^(-16))]
  }
  # p value palette
  if (class(dt$colorByLog10) %in% c('numeric', 'integer', 'double')) {
    pval_palette <- create_pvalue_colorpalette(dt, doLog10, scaleLims, nColors,
                                               colorPalette)
  } else {
    pval_palette <- list(palette = colorPalette)
  }
  # x axis labels & color
  xAxisLabs <- formatAxisLabels(dt, xOrderCol, xLabelCol, annoCol, 
                                highlightColor = annoColors)

  # first of all, if there are several genomic regons, create shadowing
  result <- ggplot()
  if (length(unique(dt$gr_id)) > 1) {
    # compute coordinates for future lines showing genome region 
    grShadeCoords <- get_gr_change_coords(plotDT = dt, labsDT = xAxisLabs$dt, 
                                          axisOrderCol = xOrderCol,
                                          yCol = yOrderCol)
    grShadeCoordsCut <- grShadeCoords$gr_coords[seq(2, 
                                                    nrow(grShadeCoords$gr_coords),
                                                    by = 2)]
    result <- result + 
      geom_rect(data = grShadeCoordsCut, fill = divider_color, inherit.aes = F,
                aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, 
                    ymax = Inf), alpha = 0.5)
  }
  # and if there are some filtered out drivers needed to be shown - shadowed 
  # areas for them
  if (any(dt$FILTER != "PASS")) {
    filteredOutShadeCoords <- get_filtered_out_gene_coords(xAxisLabs$dt, xOrderCol)
    filteredOutShadeCoordsCut <- filteredOutShadeCoords[seq(1, 
                                                            nrow(filteredOutShadeCoords),
                                                            by = 2)]
    result <- result + 
      geom_rect(data = filteredOutShadeCoordsCut, fill = 'red', 
                inherit.aes = F, alpha = 0.15,
                aes(xmin = start - 0.5, xmax = end + 0.5, 
                    ymin = -Inf, ymax = Inf))
  } 
  
  result <- result + 
    geom_tile(data = dt,
              mapping = aes_string(x = xOrderCol, y = yOrderCol,
                                   fill = 'colorByLog10')) + 
    scale_x_discrete(labels = xAxisLabs$label)
  if (class(dt$colorByLog10) %in% c('numeric', 'integer', 'double')) {
    result <- result + scale_fill_stepsn(na.value = '#FFFFFF',
                                         breaks = pval_palette$colorBreaks, 
                                         colours = pval_palette$palette)
  } else {
    result <- result + scale_fill_manual(na.value = '#FFFFFF', 
                                         values = pval_palette$palette)
  }
  if (!is.null(ggplot2Theme)) {
    result <- result + ggplot2Theme
  }
  result <- result + xlab('genomic element') + ylab('') +
    scale_y_discrete(expand = c(0, 0)) +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, size = 6, hjust = 1, 
                                     color = xAxisLabs$color))
  result
}

# Test input arguments --------------------------------------------------------
# args <- list(composite_cancer_subtype = 'Panlung',
#              drivers_composite_subtype = "completed_runs/2023_12_25/results/tables/drivers/drivers-Panlung--hg19.csv",
#              inventory_patients = 'data/inventory/inventory_patients.csv',
#              excluded_patients = list("completed_runs/2023_12_25/inputs/hypermutated-Panlung.csv"), 
#              allowed_filter_values = list("PASS", "INDEL, 2-5bp"),
#              drivers_uniform_subtypes = list("completed_runs/2023_12_25/results/tables/drivers/drivers-Adenocarcinoma--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Adenocarcinoma_met--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Adenosquamous--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Carcinoid--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Large_cell--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Mesothelioma--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Neuroendocrine_carcinoma--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Small_cell--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Squamous_cell--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Squamous_cell_met--hg19.csv"),
#              extra_studies = list("data/assets/intogene_detectedCancerGenes.csv",
#                                   "data/assets/mc3_detectedCancerGenes.csv",
#                                   "data/assets/cgc_knownCancerGenes.csv"),
#              extra_studies_names = list("intogen", "mc3", "CGC"),
#              extra_studies_tumorsubtype = list("LNET,LUAD,LUSC,NSCLC,SCLC", 
#                                                "LUAD,LUSC", "nsclc,sclc,lung"),
#              visuals_json = 'visual_parameters.json')

# args <- list(cancer_subtype = NULL, drivers_composite_subtype = NULL,
#              inventory_patients = 'data/inventory/inventory_patients.csv',
#              excluded_patients = list("completed_runs/2023_12_25/inputs/hypermutated-Adenocarcinoma.csv",
#                                       "completed_runs/2023_12_25/inputs/hypermutated-Squamous_cell.csv"), 
#              allowed_filter_values = list("PASS", "INDEL, 2-5bp"),
#              drivers_uniform_subtypes = list("completed_runs/2023_12_25/results/tables/drivers/drivers-Adenocarcinoma--hg19.csv",
#                                              "completed_runs/2023_12_25/results/tables/drivers/drivers-Squamous_cell--hg19.csv"),
#              extra_studies = list("data/assets/intogene_detectedCancerGenes.csv",
#                                   "data/assets/mc3_detectedCancerGenes.csv",
#                                   "data/assets/cgc_knownCancerGenes.csv"),
#              extra_studies_names = list("intogen", "mc3", "CGC"),
#              extra_studies_tumorsubtype = list("LNET,LUAD,LUSC,NSCLC,SCLC", 
#                                                "LUAD,LUSC", "nsclc,sclc,lung"),
#              visuals_json = 'visual_parameters.json')

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = 'figures_create_overview_plot.R')

compositeSubtypeHelp <- paste('A name of composite cancer subtype in which',
                              'de-novo drivers given in',
                              '--drivers_composite_subtype were detected.')
parser$add_argument("-c", "--composite_cancer_subtype", required = F, 
                    type = 'character', help = compositeSubtypeHelp)

compositeDriversHelp <- paste('Path to file containing de-novo detected',
                              'cancer drivers of composite subtype.')
parser$add_argument("-dcs", "--drivers_composite_subtype", required = F,  
                    type = 'character', default = NULL,
                    help = compositeDriversHelp)

patientsHelp <- paste('Path to patientsInv table listing information about',
                      'patients, their cancer types and mutation files.',
                      'Minimal columns: participant_id, tumor_subtype,',
                      'participant_tumor_subtype, somatic_path,',
                      'somatic_genome, cohort_name')
parser$add_argument("-p", "--inventory_patients", required = T, 
                    type = 'character', help = patientsHelp)

excludedHelp <- paste('Path to file(s) containing IDs of excluded patients.')
parser$add_argument("-ep", "--excluded_patients", required = F, 
                    type = 'character', nargs = '*', default = NULL, 
                    help = excludedHelp)

filterHelp <- paste("Values of FILTER column which are accepted to be plotted")
parser$add_argument("-afv", "--allowed_filter_values", required = F, 
                    type = 'character', default = NULL, nargs = '+',
                    help = filterHelp)

uniformDriversHelp <- paste('Path to files containing de-novo detected',
                            'cancer drivers in uniform subtypes. In case',
                            '--drivers_composite_subtype is given the',
                            'uniform subtypes should be part of composite',
                            'subtype.')
parser$add_argument("-dus", "--drivers_uniform_subtypes", required = F, 
                    type = 'character', default = NULL, nargs = '+',
                    help = uniformDriversHelp)

extraNamesHelp <- paste('Name(s) of extra studies to annotate drivers on the',
                        'plot with. One name per file given in',
                        '--extra_studies.')
parser$add_argument("-esn", "--extra_studies_names", required = F, 
                    type = 'character', default = NULL, nargs = '+',
                    help = extraNamesHelp)

extraSubtypesHelp <- paste('Tumor subtype to be used from the extra studies',
                           'to annotate drivers on the plot with. Multiple',
                           'tumor subtypes can be used per study if given',
                           'as string separated by comma, i.e. LUAD,LUSC. In',
                           'case several studies are given, provide several',
                           'such strings.')
parser$add_argument("-est", "--extra_studies_tumorsubtype", required = F, 
                    type = 'character', default = NULL, nargs = '+', 
                    help = extraSubtypesHelp)

extraHelp <- paste('Path to file(s) with extra studies. Must have columns:',
                   'gene_name, gene_id, known_in_tumor_subtype')
parser$add_argument("-es", "--extra_studies", required = F, 
                    type = 'character', default = NULL, nargs = '+',
                    help = extraHelp)

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

if (!is.null(args$drivers_composite_subtype)) {
  if (!file.exists(args$drivers_composite_subtype)) {
    stop('[', Sys.time(), '] ', args$drivers_composite_subtype, ' does not ',
         'exist.')
  }
}

if (!is.null(args$drivers_composite_subtype) & 
    is.null(args$composite_cancer_subtype)) {
  stop('[', Sys.time(), '] --drivers_composite_subtype is given, but ',
       '--composite_cancer_subtype is not.')
}

if (!is.null(args$composite_cancer_subtype) & 
    is.null(args$drivers_composite_subtype)) {
  stop('[', Sys.time(), '] --composite_cancer_subtype is given, but ',
       '--drivers_composite_subtype is not.')
}

if (!file.exists(args$inventory_patients)) {
  stop('[', Sys.time(), '] ', args$inventory_patients, ' does not exist.')
}

if (!is.null(args$excluded_patients)) {
  if (!all(sapply(args$excluded_patients, file.exists))) {
    stop('[', Sys.time(), '] some of the files submitted to ',
         '--excluded_patients do not exist.')
  }
}

if (!is.null(args$allowed_filter_values)) {
  args$allowed_filter_values <- list('PASS')
}

if (!is.null(args$drivers_uniform_subtypes)) {
  if (!all(sapply(args$drivers_uniform_subtypes, file.exists))) {
    stop('[', Sys.time(), '] some of the files submitted to ',
         '--drivers_uniform_subtypes do not exist.')
  }
}

if (is.null(args$drivers_composite_subtype) & 
    is.null(args$drivers_uniform_subtypes)) {
  stop('[', Sys.time(), '] Please submit file to --drivers_composite_subtype,',
       ' --drivers_uniform_subtypes or both')
}

if (!is.null(args$extra_studies)) {
  if (length(args$extra_studies) != length(args$extra_studies_names)) {
    stop('[', Sys.time(), '] Number of items submitted to --extra_studies and',
         ' --extra_studies_names does not match.')
  }
  if (length(args$extra_studies) != length(args$extra_studies_tumorsubtype)) {
    stop('[', Sys.time(), '] Number of items submitted to --extra_studies and',
         ' --extra_studies_tumorsubtype does not match.')
  }
  args$extra_studies_tumorsubtype <- lapply(args$extra_studies_tumorsubtype,
                                            function(x) unlist(strsplit(x, 
                                                                        ',')))
  args$extra_studies_tumorsubtype <- lapply(args$extra_studies_tumorsubtype,
                                            function(x) gsub(' ', '', x))
  args$extra_studies_tumorsubtype <- lapply(args$extra_studies_tumorsubtype,
                                            paste0, collapse = '|')
} else {
  if (!is.null(args$extra_studies_names)) {
    stop('[', Sys.time(), '] no --extra_studies_names is given')
  }
  if (!is.null(args$extra_studies_tumorsubtype)) {
    stop('[', Sys.time(), '] no --extra_studies_tumorsubtype is given.')
  }
}
# create directory for plots

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
printArgs(args)

# Read visuals JSON -----------------------------------------------------------
ESSENTIAL_VISUAL_NAMES <- c('ggplot2_theme', 'color_divider', 
                            'colors_tumor_type', 'overview_plot_width', 
                            'overview_plot_heigth')

visualParams <- readJsonWithVisualParameters(args$visuals_json)
message('[', Sys.time(), '] Read --visuals_json: ', args$visuals_json)

notFoundVisuals <- setdiff(ESSENTIAL_VISUAL_NAMES, names(visualParams))
if (length(notFoundVisuals)) {
  stop('[', Sys.time(), '] Following visuals: ', 
       paste(notFoundVisuals, collapse = ', '), ' not found in JSON.')
}

# Read in patients inventory --------------------------------------------------
patientsInv <- readParticipantInventory(args$inventory_patients, 1)
message('[', Sys.time(), '] Read --inventory_patients: ', 
        args$inventory_patients)

# Read in file with hypermutated patients IDs ---------------------------------
if (!is.null(args$excluded_patients)) {
  hypermutated <- lapply(args$excluded_patients, fread, header = T, 
                         stringsAsFactors = F,
                         select = c('tumor_subtype', 'participant_id'))
  message('[', Sys.time(), '] Read --excluded_patients: ', 
          paste(args$excluded_patients, collapse = ', '))
  hypermutated <- do.call(rbind, hypermutated)
  patientsInv <- patientsInv[!participant_id %in% hypermutated$participant_id]
  message('[', Sys.time(), '] Excluded hypermutated participants ids')
}

# Read composite subtype drivers ----------------------------------------------
if (!is.null(args$drivers_composite_subtype)) {
  driversCompositeSubtypesUnfiltered <- fread(args$drivers_composite_subtype, 
                                              header = T, 
                                              stringsAsFactors = F)
  driversCompositeSubtypesUnfiltered[, gene_name := gsub('__and__', '&', 
                                                         gene_name)]
  driversCompositeSubtypesUnfiltered[, gene_id := gsub('__and__', '&',
                                                       gene_id)]
  message('[', Sys.time(), '] Read --drivers: ', 
          args$drivers_composite_subtype)
}

# Read drivers detected in uniform subtypes -----------------------------------
if (!is.null(args$drivers_uniform_subtypes)) {
  driversUniformSubtypesUnfiltered <- lapply(args$drivers_uniform_subtypes, 
                                             fread, header = T,
                                             stringsAsFactors = F)
  message('[', Sys.time(), '] Read --drivers_uniform_subtypes: ',
          paste0(args$drivers_uniform_subtypes, collapse = ', '))
  driversUniformSubtypesUnfiltered <- do.call(rbind.fill, 
                                              driversUniformSubtypesUnfiltered)
  driversUniformSubtypesUnfiltered <- as.data.table(driversUniformSubtypesUnfiltered)
  driversUniformSubtypesUnfiltered[, gene_name := gsub('__and__', '&', 
                                                       gene_name)]
  driversUniformSubtypesUnfiltered[, gene_id := gsub('__and__', '&', gene_id)]
  
  if (sum(driversUniformSubtypesUnfiltered$tumor_subtype != 
          driversUniformSubtypesUnfiltered$participant_tumor_subtype, 
          na.rm = T) > 1) {
    stop() #exclusion of pan-lung mets here
  }
  
  if (!is.null(args$drivers_composite_subtype)) {
    uniformSubtypes <- unique(driversUniformSubtypesUnfiltered$tumor_subtype)
    partOfComposite <- unique(driversCompositeSubtypesUnfiltered$participant_tumor_subtype)
    notInComposite <- setdiff(uniformSubtypes, partOfComposite)
    if (length(notInComposite) > 0) {
      stop('[', Sys.time(), '] drivers for ', 
           paste(notInComposite, collapse = ', '), ' submitted via ',
           '--drivers_uniform_subtypes is not part of ',
           args$composite_cancer_subtype)
    }
  }
}

# Select genes - drivers in args$tumor_subtype or in uniform subtypes ---------
driverIDs <- data.table(gr_id = character(), gene_id = character(),
                        gene_name = character())

if (!is.null(args$drivers_composite_subtype)) {
  driversCompositeSubtypes <- driversCompositeSubtypesUnfiltered[!is.na(tier) & 
                                                                   FILTER %in% 
                                                                   args$allowed_filter_values]
  if (nrow(driversCompositeSubtypes) == 0) {
    message('[', Sys.time(), '] Did not produce plot because no drivers ',
            'in the composite subtype were detected.')
    stop_quietly()
  }
  driverIDs <- rbind(driverIDs, 
                     driversCompositeSubtypes[,.(gr_id, gene_id, gene_name)])
}

if (!is.null(args$drivers_uniform_subtypes)) {
  driversUniformSubtypes <- driversUniformSubtypesUnfiltered[!is.na(tier) & 
                                                               FILTER %in%
                                                               args$allowed_filter_values]
  if (is.null(args$drivers_composite_subtype) & nrow(driversUniformSubtypes) == 0) {
    message('[', Sys.time(), '] Did not produce plot because no drivers ',
            'in uniform subtypes were detected.')
    stop_quietly()
  }
  driverIDs <- rbind(driverIDs, 
                     driversUniformSubtypes[,.(gr_id, gene_id, gene_name)])
}

driverIDs <- unique(driverIDs)

drivers <- data.table()
if (!is.null(args$drivers_composite_subtype)) {
  drivers <- rbind.fill(drivers,
                        merge(driversCompositeSubtypesUnfiltered, driverIDs,
                              by = c("gr_id", "gene_id", "gene_name")))
}
if (!is.null(args$drivers_uniform_subtypes)) {
  drivers <- rbind.fill(drivers,
                        merge(driversUniformSubtypesUnfiltered, driverIDs,
                              by = c("gr_id", "gene_id", "gene_name")))
}
drivers <- as.data.table(drivers)

if (!is.null(args$drivers_composite_subtype)) {
  # Situations when certain drivers are only found in uniform subtype(s), but
  # not in the composite one, may occur. In such case those drivers will not 
  # have an appropriate FILTER value in composite subtype. Therefore, extra 
  # application of selection on FILTER column is needed.
  drivers <- drivers[FILTER %in% args$allowed_filter_values]
  # We don't want to keep results of uniformal subtypes where a genomic element
  # is not found as drivers
  drivers <- drivers[tumor_subtype == args$composite_cancer_subtype | !is.na(tier)]
} else {
  # in case no composite sybtype was submitted, some modifications to drivers
  # data table are needed to be done in order to be displayed properly. 
  # first, nParts_total needs to reflect total number of patients a driver is
  # mutated in across tumor subtypes
  drivers[, nParts_total := sum(nParts_total, na.rm = T), 
          by = .(gr_id, gene_id, gene_name)]
  # second, if a driver is not shared across all subtypes, it will have values
  # which differ from the accepted values in the FILTER column. Therefore, to
  # preserve such entries (we need them to see in how many patients a driver is
  # mutated) we need to force FILTER value to the same filter value as in the 
  # subtype where genomic element is found significant.
  drivers <- split(drivers, by = c('gr_id', 'gene_id', 'gene_name'), drop = T)
  for (i in 1:length(drivers)) {
    drivers[[i]]$FILTER <- intersect(drivers[[i]]$FILTER, 
                                     args$allowed_filter_values)[1]
  }
  drivers <- do.call(rbind, drivers)
}

# Read extra databases --------------------------------------------------------
if (!is.null(args$extra_studies)) {
  extra_studies <- lapply(args$extra_studies, fread, header = T,
                          stringsAsFactors = F, 
                          select = c('gene_name', 'gene_id', 
                                     'known_in_tumor_subtype'))
  names(extra_studies) <- unlist(args$extra_studies_names)
  message('[', Sys.time(), '] Read --extra_studies: ', 
          paste0(args$extra_studies, collapse = ', '))
  
  for (i in 1:length(extra_studies)) {
    extra_studies[[i]][, known_in_current_tumor_subtype := grepl(tolower(args$extra_studies_tumorsubtype[i]),
                                                                 tolower(known_in_tumor_subtype))]
  }
}

# Calculate cohort sizes ------------------------------------------------------
cohort_sizes <- patientsInv[,.(length(unique(participant_id))), 
                            by = participant_tumor_subtype]
setnames(cohort_sizes, c("participant_tumor_subtype", 'V1'), 
         c("tumor_subtype", 'cohort_size'))
cohort_sizes <- cohort_sizes[order(cohort_size, decreasing = T)]

# Set up levels(order) for genomic regions and tumor types --------------------
# sort tumor types by number of patients
TUMOR_SUBTYPE_ORDER <- cohort_sizes$tumor_subtype

# and genomic regions by maximum number of mutated patients
GR_ORDER <- drivers[FILTER == "PASS"][order(-nParts_total)]
GR_ORDER <- GR_ORDER[,.SD[1], by = gr_id]$gr_id

# in case plotting of 2-5bp filtered out drivers is requested too
GR_FILTERED_OUT_ORDER <- NULL
if (nrow(drivers[FILTER != 'PASS']) > 0) {
  GR_FILTERED_OUT_ORDER <- drivers[FILTER != 'PASS'][order(-nParts_total)]
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
if (nrow(drivers[FILTER != 'PASS' ]) > 0) {
  GENE_FILTERED_OUT_ORDER <- drivers[FILTER != 'PASS']
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
drivers[, participant_tumor_subtype := factor(participant_tumor_subtype, 
                                              TUMOR_SUBTYPE_ORDER)]
drivers[, tumor_subtype := factor(tumor_subtype, 
                                  unique(c(args$composite_cancer_subtype, 
                                           TUMOR_SUBTYPE_ORDER)))]

# Overlap data from extra databases with detected drivers ---------------------
drivers_restricted <- drivers[!is.na(tier) & 
                                FILTER %in% args$allowed_filter_values]
drivers_restricted <- drivers_restricted[,.(gr_id, gene_name, gene_id, 
                                            gene_plots_id, FILTER)]
drivers_restricted <- unique(drivers_restricted)
if (max(drivers_restricted[,.N, by = gene_plots_id]$N) > 1) {
  offenders <- drivers_restricted[,.N, by = gene_plots_id][N > 1]
  offenders <- offenders[,.(gene_name, gene_id)]
  stop('[', Sys.time(), '] Detected some drivers with > 1 FILTER value. ',
       "Can't handle this situation. Offending values: ", 
       paste(apply(offenders, 1, paste, collapse = '-'), collapse = ','))
}

for (i in 1:length(extra_studies)) {
  keys_i <- intersect(c("gene_name", "gene_id"), colnames(extra_studies[[i]]))
  if (!'is_known_cancer' %in% colnames(extra_studies[[i]])) {
    extra_studies[[i]][, is_known_cancer := T]
  }
  
  extra_studies[[i]] <- merge(drivers_restricted, extra_studies[[i]], 
                              by = keys_i, all.x = T)
  extra_studies[[i]][is.na(is_known_cancer)]$is_known_cancer <- F
  extra_studies[[i]][is.na(known_in_current_tumor_subtype)]$known_in_current_tumor_subtype <- F
  extra_studies[[i]] <- extra_studies[[i]][,.(gr_id, gene_name, gene_id, 
                                              gene_plots_id, FILTER, 
                                              is_known_cancer,
                                              known_in_current_tumor_subtype)]
  extra_studies[[i]] <- unique(extra_studies[[i]])
  extra_studies[[i]] <- extra_studies[[i]][order(gr_id, gene_name, gene_id, 
                                                 gene_plots_id, FILTER,
                                                 -is_known_cancer,
                                                 -known_in_current_tumor_subtype)]
  extra_studies[[i]] <- extra_studies[[i]][,.SD[1], 
                                           by = .(gr_id, gene_name, gene_id,
                                                  gene_plots_id, FILTER)]
  extra_studies[[i]][, dummy := names(extra_studies)[i]]
}

# Barplot number of patients -------------------------------------------------
barNpats <- barplotNpatient(driversDT = drivers, 
                            xLabelCol = 'gene_name',
                            xOrderCol = 'gene_plots_id', 
                            annoCol = 'is_known_cancer', 
                            colorPalette = visualParams$colors_tumor_type, 
                            pancanCode = args$composite_cancer_subtype, 
                            divider_color = visualParams$color_divider,
                            ggplot2Theme = visualParams$ggplot2_theme)

barNpats <- barNpats + theme(legend.position = 'bottom')
barNpats <- barNpats + guides(fill = guide_legend(nrow = 2, byrow = T)) +
    theme(legend.box = "horizontal")

# Plot annotation with extra databases ----------------------------------------
if (!is.null(args$extra_studies)) {
  extraDbPlots <- lapply(1:length(extra_studies),
                         function(x) tilePlotPvalues(extra_studies[[x]],
                                                     xLabelCol = 'gene_name',
                                                     yLabelCol = 'dummy', 
                                                     colorBy = 'is_known_cancer', 
                                                     xOrderCol = 'gene_plots_id',
                                                     doLog10 = F,
                                                     annoCol = 'is_known_cancer',
                                                     nColors = 2, 
                                                     divider_color = visualParams$color_divider,
                                                     colorPalette = c('FALSE' = 'white', 
                                                                      'TRUE' = '#8e918f'),
                                                     ggplot2Theme = visualParams$ggplot2_theme) +
                           geom_text(data = extra_studies[[x]][known_in_current_tumor_subtype == T], 
                                     mapping = aes(x = gene_plots_id, y = dummy, 
                                                   label = '*'),
                                     vjust = 0.75) + 
                           theme(legend.position = 'none',
                                 axis.title = element_blank(), 
                                 axis.ticks.y = element_blank()))
}

# Plot tile plot of p-values --------------------------------------------------
pValMatrTiles <- tilePlotPvalues(pValsDT = drivers[!is.na(tier) & 
                                                     FILTER %in% args$allowed_filter_values],
                                 xLabelCol = 'gene_name', 
                                 yLabelCol = 'tumor_subtype', 
                                 colorBy = 'brown.bh_p', 
                                 xOrderCol = 'gene_plots_id',
                                 annoCol = 'is_known_cancer',
                                 doLog10 = T,
                                 ggplot2Theme = visualParams$ggplot2_theme)
pValMatrTiles <- pValMatrTiles + labs(fill = '-log10(Brown p.value)') +
  theme(legend.direction = 'horizontal')

# Assemble --------------------------------------------------------------------
assemblyTheme <- list(theme(legend.position = 'none',
                            axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.line.x = element_blank(),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.grid.minor.y = element_blank())) 
ovrvFigPlotList <- list(barNpats + assemblyTheme)
if (!is.null(args$extra_studies)) {
  extraDbPlots <- lapply(extraDbPlots, function(x) x + assemblyTheme)
  ovrvFigPlotList <- append(ovrvFigPlotList, extraDbPlots)
}
ovrvFigPlotList[[length(ovrvFigPlotList) + 1]] <- pValMatrTiles + 
  theme(legend.position = 'none')

# main plot
plotHeights <- c(1)
if (!is.null(args$extra_studies)) {
  plotHeights <- c(plotHeights, rep(0.025, length(extra_studies)))
}
nTumorSubtypes <- length(unique(drivers[!is.na(tier)]$tumor_subtype))
plotHeights <- c(plotHeights, 0.03 * (nTumorSubtypes + 5))
DRIVERS_OVERVIEW_PLOT <- ggarrange(plotlist = ovrvFigPlotList, align = 'v',
                                   ncol = 1, heights = plotHeights)
DRIVERS_OVERVIEW_LEGENDS <- ggarrange(extract_legend(barNpats), 
                                      extract_legend(pValMatrTiles),
                                      ncol = 1)
DRIVERS_OVERVIEW_PLOT <- ggarrange(DRIVERS_OVERVIEW_PLOT, 
                                   DRIVERS_OVERVIEW_LEGENDS, ncol = 1,
                                   heights = c(1, 0.15))

# Output to file --------------------------------------------------------------
if (args$output_type == 'pdf') {
  pdf(args$output, width = visualParams$overview_plot_width, 
      height = visualParams$overview_plot_heigth)
} else {
  png(args$output, width = visualParams$overview_plot_width, 
      height = visualParams$overview_plot_heigth)
}
DRIVERS_OVERVIEW_PLOT
dev.off()

message("End time of run: ", Sys.time())
message('Total execution time: ', 
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')