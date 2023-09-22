# FILE: NBR_create_intervals.R ------------------------------------------------
#
# DESCRIPTION: Creates input .region and .txt files containing region of 
#              interest and trinucluotide content of region of interest 
#              respectively for NBR.
#
# USAGE: Rscript --vanilla NBR_create_intervals.R \
#                --region_bed path_to_genomic_region_bed \
#                --genomeFile genome_fasta \
#                --chr chromosomes_to_select
#                
#
# OPTIONS:  
# REQUIREMENTS: argparse, Biostrings, GenomicRanges, Rsamtools
# BUGS: --
# NOTES:  ---
# AUTHOR:  Inigo Martincorena - 2014
# MODERATOR: Maria Litovchenko, m.litovchenko@ucl.ac.uk
# VERSION:  1
# CREATED:  2014
# REVISION: 17.02.2020

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))

# Functions -------------------------------------------------------------------
#' getSeqlevelsStyle
#' @description Detemines a seqlevelstyle (aka chromosome naming style) from
#' the reference genome
#' @author Maria Litovchenko
#' @param fastaPath path to fasta file
#' @return UCSC, in case naming format is chr1, and NCBI otherwise
getSeqlevelsStyle <- function(fastaPath) {
  # read just the first line
  result <- fread(fastaPath, nrows = 1, header = F)$V1
  result <- ifelse(grepl('>chr', result), 'UCSC', 'NCBI')
  result
}

#' buildrefnbr
#' Function to generate reference files for NBR from a reference genome and 
#' interval or bed-like files.
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @details Rheinbay, Muhlig Nielsen, Abascal, et al. (2020) Analyses of 
#'          non-coding somatic drivers in 2,658 cancer whole genomes. Nature. 
#'          578, 102–111.
#' @param bedfile Path to a bed file with the regions to be analysed. Files 
#'                must be 4-column or 12-column standard bed files with names. 
#' @param genomefile Path to the indexed reference genome file.
#' @param onlychrs Vector of valid chromosome names (default: all chromosomes
#'                 will be included)
#' @export
buildrefnbr <- function(bedfile, genomefile, onlychrs = NULL) {
  refGenChrStyle <- getSeqlevelsStyle(genomefile)
  
    # 1. Processing the input bedfile file
    message("[1/3] Processing the input bedfile file... ", Sys.time())
    
    # Loading and processing the input bed file
    bed <- fread(bedfile, header = F, sep = "\t", stringsAsFactors = F)
    message("\t Read in bed file ... ", Sys.time())
    setnames(bed, colnames(bed)[1:4], c("chr","start","end","name"))
    
    if (!is.null(onlychrs)) {
      bed <- bed[chr %in% onlychrs]
    }
    if (nrow(bed) == 0) {
      stop("Empty bed file: if you are using the onlychrs argument, ",
           "ensure that chromosome names match those in the bed file.")
    }
    
    warning('Restricted bed file to standard chromosomes: 1-22, X, Y')
    nBefore <- nrow(bed)
    bed <- bed[chr %in% c(1:22, 'X', 'Y') | 
               chr %in% paste0('chr', c(1:22, 'X', 'Y'))]
    warning('Removed ', nBefore - nrow(bed), 'entries')
    
    if (ncol(bed) == 4) { # bed4
      # Adding a 1bp flank (the start field in a bed file is already 0-based)
      bed[, end := end + 1]
    } else if (ncol(bed) == 12) { # bed12
      bed[, segment_sizes := sapply(unlist(bed[, 11]), 
                                    function(x) strsplit(x, ',')[[1]])]
      bed[, segment_sizes := as.numeric(segment_sizes)]
      bed[, segment_starts := sapply(unlist(bed[, 12]), 
                                     function(x) strsplit(x, ',')[[1]])]
      bed[, segment_starts := as.numeric(segment_starts)]
      
      bed[, start := start + 1 + segment_starts - 1]
      bed[, end := start + 1 + segment_starts + segment_sizes]
      bed <- bed[, c('chr', 'start', 'end', 'name')]
    } else { # Invalid: bed file must have 4 or 12 columns
        stop("Invalid input bed file: please use a bed file with 4 columns, ",
             "with the 4th column being the element name, or a bed12 file ",
             "with 12 columns.")
    }
    
    regionGRlist <- makeGRangesFromDataFrame(bed, keep.extra.columns = T)
    seqlevelsStyle(regionGRlist) <- refGenChrStyle
    regionGRlist <- split(regionGRlist, bed$name)
    message("\t Split bed file ... ", Sys.time())
    regionGRlist <- reduce(regionGRlist)
    message("\t Reduced bed file ... ", Sys.time())
    regionGR <- unlist(regionGRlist)
    message("\t Performed unlist ... ", Sys.time())
    regionGR$bin <- rep(1:length(regionGRlist), sapply(regionGRlist, length))
    message("\t Added bin ID ... ", Sys.time())
    
    ## 2. Calculating the trinucleotide frequencies for each element
    message("[2/3] Calculating the trinucleotide frequencies for each ",
            "element... ", Sys.time())
    # trinucleotide content
    tmpRegionsTriNucl <- trinucleotideFrequency(scanFa(genomefile, regionGR))
    
    regionsTriNucl <- lapply(unique(regionGR$bin), 
                             function(x) which(regionGR$bin == x))
    regionsTriNucl <- lapply(regionsTriNucl,
                             function(x) if(length(x) == 1) {
                               tmpRegionsTriNucl[x, ]
                             } else {
                               colSums(tmpRegionsTriNucl[x, ])
                             })
    regionsTriNucl <- do.call(rbind, regionsTriNucl) 
    # Name and coordinates of each element
    regionslist <- data.frame(region = names(regionGRlist),
                              segments = sapply(regionGRlist, 
                                                function(x) paste(unique(seqnames(x)),
                                                                  paste(start(x) + 1,
                                                                        end(x) + 1,
                                                                        sep = '-'),
                                                                  sep = ':', 
                                                                  collapse = ',')))
    
    out1 <- cbind(regionslist, regionsTriNucl)
    names(regionGR) <- NULL
    out2 <- as.data.frame(regionGR)[, c("bin", "seqnames", "start", "end")]
    out2$start = out2$start + 1 # Removing the trinucleotide flank
    out2$end = out2$end - 1 # Removing the trinucleotide flank
  
    # Writing output files
    message("[3/3] Writing output files ",  Sys.time())
    
    if (!is.null(onlychrs)) {
      selChr <- paste(onlychrs, collapse = '_')
      write.table(out1, col.names = T, row.names = F, sep = "\t", quote = F,
                  file = sprintf(paste0("%s", "_", selChr, ".txt"), bedfile))
      write.table(out2, col.names = T, row.names = F, sep = "\t", quote = F,
                  file = sprintf(paste0("%s", "_", selChr, ".regions"),
                                 bedfile))
      write.table(regionslist, col.names = T, row.names = F, sep = "\t", 
                  quote = F,
                  file = sprintf(paste0("%s", "_", selChr, ".intervals"),
                                 bedfile))
    } else {
      write.table(out1, file = sprintf("%s.txt", bedfile), col.names = T,
                  row.names = F, sep = "\t", quote = F)
      write.table(out2, file = sprintf("%s.regions", bedfile), col.names = T,
                  row.names = F, sep = "\t", quote = F)
      write.table(regionslist, file = sprintf("%s.intervals", bedfile), 
                  col.names = T, row.names = F, sep = "\t", quote = F)
    }
}

# Inputs ----------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument('--region_bed', nargs = 1, required = T,
                    help = 'Path to file with regions of interest, bed format')
parser$add_argument('--genomeFile', nargs = 1, required = T,
                    help = 'Path to fasta file with genome')
parser$add_argument('--chr', nargs = 1, required = F,
                    help = 'Selected chromosome')
args <- parser$parse_args()

# Run buildrefnbr -------------------------------------------------------------
buildrefnbr(args$region_bed, args$genomeFile, args$chr)