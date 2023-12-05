#!/usr/bin/env Rscript

# FILE: 3c_write_regions_for_dndscv.R -----------------------------------------
#
# DESCRIPTION: Formats BED12 file containing genomic regions for one tumor 
#              subtype to dNdScv input format.
#
# USAGE: Rscript --vanilla 3c_write_regions_for_dndscv.R \
#                --gtf [path to GTF file(s)] \
#                --gtf_genomes [genome version(s) of GTF file(s)] \
#                --target_genome_path [path to reference fasta file] \
#                --target_genome_version [final genome version, i.e. hg19] \
#                --chain [path to chain file] \
#                --bed [path to BED12 file with all regions for that subtype] \
#                --cancer_subtype [cancer subtype of interest] \
#                --output [path to file to write files to] \
#
# OPTIONS: Run 
#          Rscript --vanilla 3c_write_regions_for_dndscv.R -h
#          to see the full list of options and their descriptions.
#
# REQUIREMENTS: 
# BUGS: --
# NOTES:
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  03.10.2023
# REVISION: 03.10.2023

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
source(paste0(srcDir, '/custom_functions_preprocessing.R'))
source(paste0(srcDir, '/faster_buildref_dNdScv.R'))

# Libraries -------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dndscv)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicFeatures)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(tools)))
suppressWarnings(suppressPackageStartupMessages(library(utils)))
options(scipen = 999)

# Functions -------------------------------------------------------------------
#' readMultiGtf
#' @description Reads GTF file(s) into GenomicRanges object
#' @author Maria Litovchenko
#' @param gtfPaths vector, path to GTF file(s)
#' @param gtfGenomes vector of strings, genome version of submitted gtf files
#' @param acceptedChrCodes vector with accepted chromosomal codes
#' @param chrStyle character, one of NCBI or UCSC which determine chromosome
#' naming style (1 or chr1 respectively). Final result will have this 
#' chromosome naming style.
#' @param chain chain file for liftover. Leave NULL if liftover is FALSE
#' @return GenomicRanges
readMultiGtf <- function(gtfPaths, gtfGenomes, targetGenome, chainObj = NULL,
                         acceptedChrCodes, chrStyle) {  
  if (length(gtfGenomes) != length(gtfPaths)) {
    stop('[', Sys.time(), '] Length of gtfGenomes should be the same as ',
         'length gtfPaths')
  }
  performLO <- any(!gtfGenomes %in% targetGenome)
  if (performLO & is.null(chainObj)) {
    stop('[', Sys.time(), ']  chainObj is required for liftover')
  }
  
  # import gtf(s) as GRanges and put them into txdb object
  message('[', Sys.time(), '] Started reading ', paste(gtfPaths, sep = ', '))
  result <- lapply(gtfPaths, import)
  
  # perform liftover, if needed
  if (performLO) {
    for (gtfIdx in seq(1, length(result))[gtfGenomes != targetGenome]) {
      loGR <- liftOverGenomicRegion(result[[gtfIdx]], chainObj)
      if ('transcript_id' %in% colnames(mcols(result[[gtfIdx]]))) {
        loGR <- checkLiftOverByTr(result[[gtfIdx]], loGR)
      }
      result[[gtfIdx]] <- loGR
    }
  }
  result <- unlist(GRangesList(result))
  result <- result[!duplicated(as.data.table(result))]
  
  if (!is.null(acceptedChrCodes)) {
    newSeqLevels <- intersect(as.character(unique(seqnames(sort(result)))),
                              acceptedChrCodes)
    result <- keepSeqlevels(result, value = as.factor(newSeqLevels),
                            pruning.mode = "coarse")
    result <- sort(result)
  }
  
  if (any(c(23, '23', 'chr23') %in% seqlevels(result))) {
    message('[', Sys.time(), '] Found that chr X is coded with 23. Changed it',
            ' to X.')
    result <- as.data.table(result)
    result[, seqnames := gsub('23', 'X', seqnames)]
    result[, seqnames := gsub('24', 'Y', seqnames)]
    result[, seqnames := gsub('25', 'M', seqnames)]
    result <- makeGRangesFromDataFrame(result, keep.extra.columns = T)
  }
  message('[', Sys.time(), '] Finished reading ', paste(gtfPaths, sep =', '))
  
  seqlevelsStyle(result) <- intersect(c('UCSC', 'NCBI'), chrStyle)
  
  result
}

#' removeNotTargetGenes
#' @description Removes genes which are not in targetGenes from data table with
#' genomic regions of interest.
#' @author Maria Litovchenko
#' @param gr GRanges object. Must have columns: gene_id, transcript_biotype
#' or transcript_type
#' @param targetGenes vector of gene IDs which length fully passed filtering
#' @return filtered GenomicRegions
removeNotTargetGenes <- function(gr, targetGenes) {
  # summarize transcript biotypes of removed genes
  removedGenes <- gr[!gr$gene_id %in% targetGenes]
  removedGenes <- mcols(removedGenes)
  removedGenes <- removedGenes[, intersect(colnames(removedGenes),
                                           c('gene_id', 'transcript_biotype',
                                             'transcript_type'))]
  removedGenes <- unique(as.data.table(removedGenes))
  setnames(removedGenes, 'transcript_type', 'transcript_biotype',
           skip_absent = T)
  removedGenes <- removedGenes[,.N, by = transcript_biotype][order(-N)]
  
  nbefore <- length(unique(gr$gene_id))
  result <- gr[gr$gene_id %in% targetGenes]
  nafter <- length(unique(result$gene_id))
  ndiff <- length(setdiff(unique(targetGenes), unique(result$gene_id)))
  message('[', Sys.time(), '] Removed ', nbefore - nafter, '(',
          round(100 * (nbefore - nafter) / nbefore, 2), '%) genes due to ',
          'not being present in target ranges. Summary of biotype ',
          'distribution: ')
  message(paste0(capture.output(knitr::kable(removedGenes, 
                                             format = "markdown")),
                 collapse = '\n'))
  message('Number of genes from target ranges not present in GTF files: ', 
          ndiff)
  
  result
}

#' removeGenesNotIndNdScvCovs
#' @description Removes genes for which covariates are not available in dNdScv.
#' @author Maria Litovchenko
#' @param buildRefRes result of function buildref
#' @return filtered result of function buildref
removeGenesNotIndNdScvCovs <- function(buildRefRes) {
  data("covariates_hg19", package = "dndscv")
  
  # process gr_genes
  buildRefRes$gr_genes <- buildRefRes$gr_genes[buildRefRes$gr_genes$names %in% 
                                                 c('CDKN2A', rownames(covs))]
  
  # process RefCDS 
  nbefore <- length(buildRefRes$RefCDS)
  genesInRefcds <- sapply(buildRefRes$RefCDS, function(x) x$gene_name)
  buildRefRes$RefCDS <- buildRefRes$RefCDS[genesInRefcds %in% 
                                             c('CDKN2A', rownames(covs))] 
  nafter <- length(buildRefRes$RefCDS)
  message('[', Sys.time(), '] Removed ', nbefore - nafter, '(',
          round(100 * (nbefore - nafter) / nbefore, 2), '%) genes from ', 
          'GTF due to not being present in covariates object of dNdScv.')
  

  rm(covs)
  buildRefRes
}

#' create_transcript_table
#' Creates custom transcript table from GTF file(s) which can further be used 
#' in dNdScv buildref function.
#' @param gtfPaths vector of paths to GTF files. Although a vector is accepted,
#' use of just one file is recommended. 
#' @param gtfGenomes vector of strings, genome version of submitted GTF files
#' @param targetGenome string, target (final) genome version
#' @param refGenFastaPath character, path to reference genome file. Have to 
#' be the same genome version as target genome.
#' @param acceptChrCodes character vector, accepted chromosomal names
#' @param targetGenes vector of gene IDs which length fully passed filtering
#' @param outputPath string, path to where save the transcript table
#' @return file with transcript table
create_transcript_table <- function(gtfPaths, gtfGenomes, 
                                    targetGenome, refGenFastaPath, 
                                    chainObj = NULL, acceptedChrCodes = NULL,
                                    targetGenes = NULL, outputPath) {
  # import gtf(s) as GRanges and put them into txdb object
  message('[', Sys.time(), '] Started generation of transcript table ',
          'from ', paste(gtfPaths, collapse = ', '))
  gtfGR <- readMultiGtf(gtfPaths, gtfGenomes, targetGenome, chainObj, 
                        acceptedChrCodes, getSeqlevelsStyle(refGenFastaPath))
  # remove gene which are not in targetGenes
  if (!is.null(targetGenes)) {
    gtfGR <- removeNotTargetGenes(gtfGR, targetGenes)
  }
  
  # extract cds for protein coding genes
  gtfCDS <- cdsBy(makeTxDbFromGRanges(gtfGR), 'tx', use.names = T)
  
  protCodTx <- as.data.table(mcols(gtfGR))
  if ('protein_id' %in% colnames(protCodTx)) {
    protCodTx <- protCodTx[!is.na(protein_id)]
  }
  if ('ccdsid' %in% colnames(protCodTx)) {
    protCodTx <- protCodTx[!is.na(ccdsid)]
  }
  if ('transcript_type' %in% colnames(protCodTx)) {
    protCodTx <- protCodTx[transcript_type == 'protein_coding']
  }
  if ('transcript_biotype' %in% colnames(protCodTx)) {
    protCodTx <- protCodTx[transcript_biotype == 'protein_coding']
  }
  protCodTx <- protCodTx[, intersect(c('gene_id', 'gene_name', 'transcript_id',
                                       'protein_id'), colnames(protCodTx)), 
                         with = F]
  if (!'protein_id' %in% colnames(protCodTx) &  
      !'transcript_id' %in% colnames(protCodTx)) {
    stop('[', Sys.time(), '] None of protein_id or transcript_id is present ',
         'in GTF file.')
  }
  
  protCodTx <- unique(protCodTx)
  gtfCDS <- gtfCDS[intersect(protCodTx$transcript_id, names(gtfCDS))]
  gtfCDS <- unlist(gtfCDS)
  mcols(gtfCDS) <- NULL
  
  # annotation columns
  setkey(protCodTx, transcript_id)
  genomeCDSmcols <- protCodTx[names(gtfCDS)]
  setnames(genomeCDSmcols,
           c('gene_id', 'gene_name', 'transcript_id', 'protein_id'), 
           c('gene.id', 'gene.name',  'transcript_id', 'cds.id'),
           skip_absent = T)
  if (!'cds.id' %in% colnames(genomeCDSmcols)) {
    message('[', Sys.time(), '] protein_id is not present in GTF file. Will ',
            'use transcript_id instead.')
    genomeCDSmcols[, cds.id := transcript_id]
  }
  mcols(gtfCDS) <- genomeCDSmcols
  
  # convert to data.table for speed
  gtfCDS <- as.data.table(gtfCDS)
  gtfCDS[, strand := ifelse(strand == '+', 1, -1)]
  gtfCDS[, length := sum(width), by = transcript_id]
  # add cds.start and cds.end
  gtfCDS[, cds.end := cumsum(width), by = transcript_id]
  gtfCDS[, cds.start := c(1, cds.end[-length(cds.end)] + 1),
            by = transcript_id]
  # format for output
  gtfCDS <- gtfCDS[, c('gene.id', 'gene.name', 'cds.id', 'seqnames', 'start',
                       'end', 'cds.start', 'cds.end', 'length', 'strand')]
  setnames(gtfCDS, c('seqnames', 'start', 'end'), 
           c('chr', 'cds.start', 'cds.end'))
 
  write.table(gtfCDS, outputPath, sep = '\t', append = F, quote = F, 
              col.names = T, row.names = F)

  message('[', Sys.time(), '] Finished generation of transcript table ',
          'from ', paste(gtfPaths, collapse = ', '))
}

# Parse input arguments -------------------------------------------------------
# create parser object
parser <- ArgumentParser(prog = '3c_write_regions_for_dndscv.R')

gtfHelp <- 'A path to GTF file(s) for that cancer subtype'
parser$add_argument("-a", "--gtf", required = T, type = 'character', 
                    nargs = '+', default = NULL, help = gtfHelp)

gtfGenomeHelp <- 'Genome version of submitted GTF files'
parser$add_argument("-ag", "--gtf_genomes", required = T, type = 'character', 
                    nargs = '+', default = NULL, help = gtfGenomeHelp)

subtypeHelp <- paste('A cancer subtype to select from patientsInv table. Only',
                     'mutations from patients with that cancer type will be',
                     'selected. In case an analysis of several cancer types',
                     'needed to be performed please run this script ',
                     'separetedly for each cancer type.')
parser$add_argument("-c", "--cancer_subtype", required = T, type = 'character',
                    default = NULL, help = subtypeHelp)

targetGenomePathHelp <- paste('Path to the fasta file, genome version of',
                              'which should be same as ',
                              '--target_genome_version.')
parser$add_argument("-f", "--target_genome_path", required = T, default = '', 
                    type = 'character', help = targetGenomePathHelp)

targetGenomeHelp <- paste("Genome version, i.e. hg19, to which input variants",
                          "files for software should be brought.",
                          "Default: hg19.")
parser$add_argument("-g", "--target_genome_version", required = F,
                    default = 'hg19', type = 'character',
                    help = targetGenomeHelp)

bedHelp <- 'A path to BED12 file with all regions for that cancer subtype'
parser$add_argument("-b", "--bed", required = F, type = 'character', 
                    default = NULL, help = bedHelp)

chainHelp <- paste('Path to chain file in case genome version of mutations is',
                   'not the same as --target_genome_version')
parser$add_argument("-l", "--chain", required = F, default = NULL,
                    type = 'character', help = chainHelp)

coresHelp <- 'How many cores the script should use. Default: 1.'
parser$add_argument("-n", "--cores", required = F, type = 'integer', 
                    default = 1, help = coresHelp)

parser$add_argument("-o", "--output", required = T, type = 'character',
                    help = "Path to the output folder")

args <- parser$parse_args()
check_input_arguments_preproc(args, outputType = 'folder')

timeStart <- Sys.time()
message('[', Sys.time(), '] Start time of run')
args$target_genome_version <- args$target_genome_version[args$target_genome_version != '']
args$target_genome_version <- args$target_genome_version[args$target_genome_version != ' ']
printArgs(args)

# Test inputs -----------------------------------------------------------------
# args <- list(gtf = list('../data/genomic_regions/standard_gtf/Homo_sapiens_assembly19.gencode.gtf'),
#              gtf_genomes = list('hg19'),
#              target_genome_path = '../data/assets/reference_genome/Homo_sapiens_assembly19.fasta',
#              target_genome_version = 'hg19', chain = NULL, 
#              cancer_subtype = 'LUAD',
#              bed = '../TEST/inputs/inputGR-LUAD-hg19.bed', output = '.',
#              cores = 4)

# Create transcript table -----------------------------------------------------
transrTabPath <- paste0(args$cancer_subtype, '_transcript_table.csv')

targetGeneIDs <- NULL
if (!is.null(args$bed)) {
  targetGeneIDs <- unique(readBED12(args$bed)$gene_id)
}

# read in chain file
chain <- NULL
if (any(!args$gtf_genomes %in% args$target_genome_version)) {
  chain <- import.chain(args$chain)
}

create_transcript_table(gtfPaths = args$gtf, gtfGenomes = args$gtf_genomes,
                        targetGenome = args$target_genome_version,
                        refGenFastaPath = args$target_genome_path,
                        chainObj = chain, acceptedChrCodes = acceptedChrNames,
                        targetGenes = targetGeneIDs, 
                        outputPath = transrTabPath)

# Build RefCDS with dNdScv ----------------------------------------------------
refGenStyle <- getSeqlevelsStyle(args$target_genome_path)

message('[', Sys.time(), '] Started building RefRda object for ',
        'dNdScv/DIGdriver')
rda <- buildref_faster(cdsfile = transrTabPath, cores = args$cores,
                       genomefile = args$target_genome_path)

# Save to file RefCDS without restriction to covariates -----------------------
outputFile <- paste0(args$output, '/', args$cancer_subtype, '_', refGenStyle,
                     '.Rda')
save(rda$RefCDS, rda$gr_genes, file = outputFile)

# change refGenStyle to the opposite (NCBI -> UCSC, UCSC -> NCBI) to save 
# Rda in both formats. 
if (refGenStyle == 'NCBI') {
  seqlevelsStyle(rda$gr_genes) <- 'UCSC'
  
  for (RefCDSidx in 1:length(rda$RefCDS)) {
    rda$RefCDS[[RefCDSidx]]$chr <- paste0('chr', rda$RefCDS[[RefCDSidx]]$chr)
  }
  
  outputFile <- paste0(args$output, '/', args$cancer_subtype, '_UCSC.Rda')
}
if (refGenStyle == 'UCSC') {
  seqlevelsStyle(rda$gr_genes) <- 'NCBI'
  
  for (RefCDSidx in 1:length(rda$RefCDS)) {
    rda$RefCDS[[RefCDSidx]]$chr <- gsub('^chr', '', 
                                        rda$RefCDS[[RefCDSidx]]$chr)
  }
   
  outputFile <- paste0(args$output, '/', args$cancer_subtype, '_NCBI.Rda')
}
save(rda$RefCDS, rda$gr_genes, file = outputFile)

# Save to file RefCDS with restriction to covariates --------------------------
RefCDS <- removeGenesNotIndNdScvCovs(RefCDS)

outputFile <- paste0(args$output, '/', args$cancer_subtype, '_', refGenStyle,
                     '_removedGenesNotInCovs.Rda')
save(rda$RefCDS, rda$gr_genes, file = outputFile)

# change refGenStyle to the opposite (NCBI -> UCSC, UCSC -> NCBI) to save 
# Rda in both formats. 
if (refGenStyle == 'NCBI') {
  seqlevelsStyle(rda$gr_genes) <- 'UCSC'
  
  for (RefCDSidx in 1:length(rda$RefCDS)) {
    rda$RefCDS[[RefCDSidx]]$chr <- paste0('chr', rda$RefCDS[[RefCDSidx]]$chr)
  }
  
  outputFile <- paste0(args$output, '/', args$cancer_subtype, 
                       '_UCSC_removedGenesNotInCovs.Rda')
}
if (refGenStyle == 'UCSC') {
  seqlevelsStyle(rda$gr_genes) <- 'NCBI'
  
  for (RefCDSidx in 1:length(rda$RefCDS)) {
    rda$RefCDS[[RefCDSidx]]$chr <- gsub('^chr', '', 
                                        rda$RefCDS[[RefCDSidx]]$chr)
  }
  
  outputFile <- paste0(args$output, '/', args$cancer_subtype, 
                       '_NCBI_removedGenesNotInCovs.Rda')
}
save(rda$RefCDS, rda$gr_genes, file = outputFile)

message('[', Sys.time(), '] Finished building RefRda object for ',
        'dNdScv/DIGdriver')

message("End time of run: ", Sys.time())
message('Total execution time: ',  
        difftime(Sys.time(), timeStart, units = 'mins'), ' mins.')
message('Finished!')