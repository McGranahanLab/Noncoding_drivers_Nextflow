#!/usr/bin/env Rscript

# FILE: NBR.R -----------------------------------------------------------------
# DESCRIPTION: Non-coding driver detection with precomputed trinucleotide 
# composition and using the local density of putatively neutral mutations as 
# covariates ("t") in the framework of a negative binomial regression.
#
# USAGE: Rscript --vanilla NBR.R \
#                --mutations_file path_to_mutation_file \
#                --target_regions_path path_to_regions_file \
#                --target_regions_trinucfreqs_path path_to_trinucleotideCont_of_regions_file \
#                --genomeFile genome_fasta \
#                --max_num_muts_perRegion_perSample 2 \
#                --unique_indelsites_FLAG 1 \
#                --gr_drivers driver_regions \
#                --regions_neutralbins_file neutral_regions_file \
#                --trinucfreq_neutralbins_file trinucfreqs_regions_file \
#                --out_prefix results_file_prefix
#
# OPTIONS:  
# REQUIREMENTS: argparse, GenomicRanges, Rsamtools, MASS
# BUGS: --
# NOTES:  ---
# AUTHOR:  Inigo Martincorena - 2014
# MODERATOR: Maria Litovchenko, m.litovchenko@ucl.ac.uk
# VERSION:  1
# CREATED:  2014
# REVISION: 17.02.2020

# Instructions for Mark Cowley et al: --------
# Mutations file must be a tab-separated matrix, with columns: "sampleID","chr","pos","ref","mut"
# Example:
#  sampleID	chr	pos	ref	mut
#  0009b464-b376-4fbc-8a56-da538269a02f	1	1230448	G	A
#  0009b464-b376-4fbc-8a56-da538269a02f	1	1609723	C	T
#  0009b464-b376-4fbc-8a56-da538269a02f	1	1903276	C	T
#  0009b464-b376-4fbc-8a56-da538269a02f	1	2574999	C	T
#  ...
#
# Use chromosome names without "chr", and use build 37 (hg19) coordinates
#
# Modify PATH_TO_DATA_AND_GENOME and genomeFile to match your paths and file names
# PATH_TO_DATA_AND_GENOME: path to GRanges_driver_regions.RData, 
#                                  Trinucfreqs_within_100kb_bins.txt 
#                                  Neutral_regions_within_100kb_bins.txt
#                                  Regions/
#
# Examples:
# miRNA promoters:
#    Rscript find_noncoding_drivers_precomp_correctoverlaps.R "mutations.txt" "Regions/mirna.prom.bed" 
# miRNA precursors:
#    Rscript find_noncoding_drivers_precomp_correctoverlaps.R "mutations.txt" Regions/mirna.pre.bed
# miRNA mature:
#    Rscript find_noncoding_drivers_precomp_correctoverlaps.R "mutations.txt" Regions/mirna.mat.bed
# Protein-coding genes promoters:
#    Rscript find_noncoding_drivers_precomp_correctoverlaps.R "mutations.txt" "Regions/gc19_pc.promCore.bed" 
#
# Input files:
# - mutations.txt: Table of mutations in the 5 column format (sampleID\tchr\tpos\tref\tmut)
# - out_prefix: Region_name\tIntervals_string\tTrinucleotide composition
#
# Test:
#     Rscript find_noncoding_drivers_precomp_correctoverlaps.R  example_Thy-AdenoCa.muts.5cols Regions/gc19_pc.promCore.bed
#   Main output:
#     Regions/gc19_pc.promCore.bed-Selection_output.txt
#       region	chr	start	end	exp_subs	exp_indels	obs_subs	obs_indels	local_t	local_t_indels	pval_subs	pval_indels	pval_both	qval_subs	qval_indels	qval_both	exclude_forfit	cv_predicted_subs	cv_predicted_indels	pval_subs_CV	pval_indels_CV	pval_both_CV	qval_subs_CV	qval_indels_CV	qval_both_CV	qval_both_CV_tiered	obsexp_subs_mle	obsexp_subs_low	obsexp_subs_high	obsexp_indels_mle	obsexp_indels_low	obsexp_indels_high
#       gc19_pc.promCore::gencode::TERT::ENSG00000164362.14	5	1295105	1295362	0.009328714	0.000690961	12	0	1.360598307	0.838612166	7.35E-21	1	0	1.48E-16	1	0	TRUE	0.008465574	0.000666505	6.23E-26	1	0	1.26E-21	1	0	0	1417.505812	721.8452856	2562.228503	0	0	2935.603388
#       gc19_pc.promCore::gencode::PLEKHS1::ENSG00000148735.10	10	115511013	115534874	0.043607871	0.004911716	4	0	0.932850568	1.141215427	2.48E-05	1	0.000287769	0.249859112	1	1	FALSE	0.042695237	0.005474887	2.45E-06	1	3.41E-05	0.024713939	1	0.343976004	0.680680582	93.68726473	28.34270653	229.041039	0	0	357.5185204
#       gc19_pc.promCore::gencode::OR4F5::ENSG00000186092.4	1	68891	69090	0.005602044	0.000535629	0	0	0.721785491	0.552476592	1	1	1	1	1	1	FALSE	0.005694223	0.00045065	1	1	1	1	1	1	1	0	0	343.7822825	0	0	4341.780179
#       gc19_pc.promCore::gencode::AL627309.1::ENSG00000237683.5	1	139310	139579	0.008023746	0.000723099	0	0	0.703778928	0.537202648	1	1	1	1	1	1	FALSE	0.008181887	0.000603953	1	1	1	1	1	1	1	0	0	239.1610734	0	0	3239.682298
#       gc19_pc.promCore::gencode::OR4F29::ENSG00000235249.1	1	367440	367658	0.004864302	0.000586513	0	0	0.727080947	0.509529225	1	1	1	1	1	1	FALSE	0.004939695	0.000483438	1	1	1	1	1	1	1	0	0	396.0928871	0	0	4047.297819
#
#   Columns of interest: mainly "qval_subs_CV","qval_indels_CV","qval_both_CV","qval_both_CV_tiered"
#   The corresponding pvalues can be of interest if you want to do restricted hypotheses testing on a group of
#    genes or regions.
# 
#   If the script is going to be run for different sets of mutations and the same regions, 
#    consider moving the output files to avoid overwritting. Or modify the "write.table" lines
#    to choose your own output naming choices.
#

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(MASS))

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

map_mutations = function(mutations, intervals, chrStyle) {
  rangesREGS = GRanges(intervals[,2], IRanges(as.numeric(intervals[,3]), 
                                              as.numeric(intervals[,4])))
  seqlevelsStyle(rangesREGS) <- chrStyle
  
  mut_start = as.numeric(as.vector(mutations$pos))
  mut_end = mut_start + nchar(as.vector(mutations$ref)) - 1
  rangesMUTS = GRanges(as.vector(mutations$chr), IRanges(mut_start, mut_end))
  seqlevelsStyle(rangesMUTS) <- chrStyle
  
  olmatrix = as.matrix(findOverlaps(rangesMUTS, rangesREGS, type="any",
                                    select="all"))
  return(olmatrix)
}

# Subfunction that subsamples mutations to a maximum of nmax mutations per sample
# per region of interest
cap_mutations = function(m = sample_and_region,
                         nmax = max_num_muts_perRegion_perSample) {
  
  nmut = dim(m)[1]
  mapped_muts = which(m[,2]!="0")
  m = m[mapped_muts,]
  
  rows = paste(m[,1],m[,2],sep=",")
  freqs = table(rows)
  freqs = freqs[freqs>nmax]
  duplrows = strsplit(names(freqs),split=",")
  rmrows = array(0,dim(m)[1])
  
  if (length(freqs)>0) {
    for (j in 1:length(freqs)) {
      vals = duplrows[[j]]
      if (vals[2]!="0") { # Only for mapped mutations we apply the capping
        pos = which(m[,1]==vals[1] & m[,2]==vals[2])
        rmrows[sample(pos, length(pos)-nmax)] = 1
      }
    }
  }
  rmrows_allmuts = array(0,nmut)
  rmrows_allmuts[mapped_muts[rmrows==1]] = 1
  return(rmrows_allmuts)
}

#' fit_global_NBR_model
#' @description Estimates global selection coefficience for NBR
#' @author Maria Litovchenko
#' @param global_l_matrix matrix with column names A, C, G, T and row names
#' equal to all possible trinucleotides (64). A value in the cell = number of 
#' all possible mutation in target regions
#' @param global_l_matrix_neutral matrix with column names A, C, G, T and row 
#' names equal to all possible trinucleotides (64). A value in the cell = 
#' number of possible mutation in neutral regions
#' @param global_n_matrix matrix with column names A, C, G, T and row names
#' equal to all possible trinucleotides (64). A value in the cell = number of 
#' actually detected mutation in target regions
#' @param global_n_matrix_neutral matrix with column names A, C, G, T and row 
#' names equal to all possible trinucleotides (64). A value in the cell =  
#' number of actually detected mutation in neutral regions
#' @return named list with members "par" and "model", where par is data frame,
#' with column names name, mle, cilow, cihigh and model is model object from
#' glm command
fit_global_NBR_model <- function(global_l_matrix, global_l_matrix_neutral,
                                 global_n_matrix, global_n_matrix_neutral) {
  # 1. create substitution model
  trinucl_change <- as.data.table(expand.grid(rownames(global_l_matrix),
                                              colnames(global_l_matrix)))
  colnames(trinucl_change) <- c('refCod', 'altCod')
  trinucl_change[, altCod := sapply(1:nrow(trinucl_change), 
                                    function(x) paste0(substr(refCod[x], 1, 1),
                                                       altCod[x]))]
  trinucl_change[, altCod := sapply(1:nrow(trinucl_change), 
                                    function(x) paste0(altCod[x],
                                                       substr(refCod[x], 3,
                                                              3)))]
  trinucl_change <- apply(trinucl_change, 1, paste0, collapse = '>')
  # it's ok to keep AAA>AAA because they will be NA in global_l_matrix and
  # global_l_matrix_neutral
  substmodel <- cbind(paste0('t*', gsub('TTT>TGT', 't', trinucl_change)), 
                      paste0('t*', gsub('TTT>TGT', 't', trinucl_change)))
  colnames(substmodel) <- c('neutr', 'target')
  rownames(substmodel) <- trinucl_change
  substmodel[, 2] <- gsub('t[*]t', 't', paste0(substmodel[, 2], '*wall'))
  substmodel[, 1] <- gsub('t[*]t', 't', substmodel[, 1])
  # substitution model created
  
  # 2. create L and N matrix, similar to matrices of dNdScv
  L <- cbind(c(global_l_matrix_neutral), c(global_l_matrix))
  rownames(L) <- rownames(substmodel)
  L <- L[sort(rownames(L)), ]
  N <- cbind(c(global_n_matrix_neutral), c(global_n_matrix))
  rownames(N) <- rownames(substmodel)
  N <- N[sort(rownames(N)), ]
  substmodel <- substmodel[sort(rownames(substmodel)), ]
  
  # 3. From here on - fit_substmodel as in dNdScv
  l = c(L)
  n = c(N)
  r = c(substmodel)
  to_keep <- which(!(l == 0 | is.na(l)))
  n = n[to_keep]
  r = r[to_keep] 
  l = l[to_keep] 
  
  params = unique(base::strsplit(x = paste(r, collapse = "*"), 
                                 split = "\\*")[[1]])
  indmat = as.data.frame(array(0, dim = c(length(r), length(params))))
  colnames(indmat) = params
  for (j in 1:length(r)) {
    indmat[j, base::strsplit(r[j], split = "\\*")[[1]]] = 1
  }
  model = glm(formula = n ~ offset(log(l)) + . - 1, data = indmat, 
              family = poisson(link = log))
  mle = exp(coefficients(model))
  ci = exp(confint.default(model))
  par = data.frame(name = gsub("`", "", rownames(ci)), mle = mle[rownames(ci)],
                   cilow = ci[, 1], cihigh = ci[, 2])
  return(list(par = par, model = model))
}


# Subfunction to calculate CI95% for the obs/exp ratios of subs and indels in NBR
ci95nbr = function(n_obs,n_exp,nb_size) {
  wmax = 100000 
  iter = 6
  grid_size = 10
  cutoff = qchisq(p = 0.95, df = 1) # Default params
  wmle = n_obs/n_exp # MLE for w
  ml = dnbinom(x=n_obs, mu=n_exp*wmle, size=nb_size, log=T) # LogLik under MLE
  
  if (!is.na(n_exp)) {
    if (wmle<wmax) {
      # 1. Iterative search of lower bound CI95%
      if (wmle>0) {
        search_range = c(1e-9, wmle)
        for (it in 1:iter) {
          wvec = seq(search_range[1], search_range[2], length.out=grid_size)
          ll = dnbinom(x=n_obs, mu=n_exp*wvec, size=nb_size, log=T)
          lr = 2*(ml-ll) > cutoff
          ind = max(which(wvec<=wmle & lr))
          search_range = c(wvec[ind], wvec[ind+1])
        }
        w_low = wvec[ind]
      } else {
        w_low = 0
      }
      
      # 2. Iterative search of higher bound CI95%
      search_range = c(wmle, wmax)
      llhighbound = dnbinom(x=n_obs, mu=n_exp*wmax, size=nb_size, log=T)
      outofboundaries = !(2*(ml-llhighbound) > cutoff)
      if (!outofboundaries) {
        for (it in 1:iter) {
          wvec = seq(search_range[1], search_range[2],length.out=grid_size)
          ll = dnbinom(x=n_obs, mu=n_exp*wvec, size=nb_size, log=T)
          lr = 2*(ml-ll) > cutoff
          ind = min(which(wvec>=wmle & lr))
          search_range = c(wvec[ind-1], wvec[ind])
        }
        w_high = wvec[ind]
      } else {
        w_high = wmax
      }
    } else {
      wmle = w_low = w_high = wmax # Out of bounds
    }
  } else {
    wmle = w_low = w_high = NA # invalid
  }
  
  return(c(wmle,w_low,w_high))
}

#' returnNULL
#' @description Returns null data table which emulates ActiveDriverWGS results
#' @param c error message
#' @return data table 
returnNULL <- function(c) {
  message('\n[', Sys.time(), '] WARNING: NBR did not converge: ', c)
  nullDT <- suppressWarnings(data.table(region = character(0), 
                                        chr = character(0),
                                        start = integer(0), end = integer(0),
                                        exp_subs = integer(0), 
                                        exp_indels = double(0),
                                        obs_subs = integer(0), 
                                        obs_indels = integer(0),
                                        local_t = double(0), 
                                        local_t_indels = double(0),
                                        pval_subs = double(0),
                                        pval_indels = double(0),
                                        pval_both = double(0), 
                                        qval_subs = double(0),
                                        qval_indels = double(0), 
                                        qval_both = double(0),
                                        exlude_forfit = double(0), 
                                        cv_predicted_subs = double(0),
                                        cv_predicted_indels = double(0), 
                                        pval_subs_CV = double(0),
                                        pval_indels_CV = double(0), 
                                        qval_subs_CV = double(0),
                                        qval_indels_CV = double(0), 
                                        qval_both_CV = double(0),
                                        qval_both_CV_tiered = double(0), 
                                        obsexp_subs_mls = double(0),
                                        obsexp_subs_low = double(0),
                                        obsexp_subs_high = double(0),
                                        obsexp_indels_mle = double(0),
                                        obsexp_indels_low = double(0),
                                        obsexp_indels_high = double(0)))
  return(nullDT)
}

#' stop_quietly
#' @description Stops script without error
#' @return void
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

# Inputs ----------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument('--mutations_file', nargs = 1, required = T,
                    help = 'Path to file with mutations')
parser$add_argument('--target_regions_path', nargs = 1, required = T,
                    help = 'Path to targer regions file (extension .regions)')
parser$add_argument('--target_regions_trinucfreqs_path', nargs = 1, required = T,
                    help = paste('Path to file with trinuclotide content of ',
                                 'target regions (extension .txt)'))
parser$add_argument('--genomeFile', nargs = 1, required = T,
                    help = 'Path to fasta file with genome')
parser$add_argument('--max_num_muts_perRegion_perSample', nargs = 1, 
                    required = F, default = 2,
                    help = 'Maximum number mutations per region per sample')
parser$add_argument('--unique_indelsites_FLAG', nargs = 1, 
                    required = F, default = 1,
                    help = 'Flag for indels')
parser$add_argument('--gr_drivers', nargs = 1, required = T, 
                    help = paste('Path to regions to be excluded for the ', 
                                 'background model fit'))
parser$add_argument('--regions_neutralbins_file', nargs = 1, required = T, 
                    help = paste('Path to the file with neutral regions'))
parser$add_argument('--trinucfreq_neutralbins_file', nargs = 1, required = T, 
                    help = paste('Path to trinucleotide frequency of neutral',
                                'regions in 100kb bins'))
parser$add_argument('--out_prefix', nargs = 1, required = T,
                    help = 'Prefix for output')

args <- parser$parse_args()
mutations_file = args$mutations_file
target_regions_path = args$target_regions_path
target_regions_trinucfreqs_path = args$target_regions_trinucfreqs_path
genomeFile = args$genomeFile
max_num_muts_perRegion_perSample = as.integer(args$max_num_muts_perRegion_perSample)
unique_indelsites_FLAG = as.integer(args$unique_indelsites_FLAG)
regions_neutralbins_file = args$regions_neutralbins_file
trinucfreq_neutralbins_file = args$trinucfreq_neutralbins_file
out_prefix = args$out_prefix

message('Inputs: ')
message('mutations_file: ', mutations_file)
message('target_regions_path: ', target_regions_path)
message('target_regions_trinucfreqs_path: ', target_regions_trinucfreqs_path)
message('genomeFile: ', genomeFile)
message('max_num_muts_perRegion_perSample: ', max_num_muts_perRegion_perSample)
message('unique_indelsites_FLAG: ', unique_indelsites_FLAG)
message('gr_drivers: ', args$gr_drivers)
message('regions_neutralbins_file: ', regions_neutralbins_file)
message('trinucfreq_neutralbins_file: ', trinucfreq_neutralbins_file)
message('out_prefix: ', out_prefix)

message('Started at: ', Sys.time())

# get chromosome naming styles from reference genome
refGenChrStyle <- getSeqlevelsStyle(genomeFile)

if (grepl('rdata', tolower(args$gr_drivers))) {
  load(args$gr_drivers)
} else {
  gr_drivers <- fread(args$gr_drivers, header = T, stringsAsFactors = F)
  colnames(gr_drivers) <- c('chr', 'start', 'end')
  gr_drivers <- makeGRangesFromDataFrame(gr_drivers)
  seqlevelsStyle(gr_drivers) <- refGenChrStyle
}

dir.create(dirname(out_prefix), recursive = T)

# 0. Loading the precomputed region files -------------------------------------
message("Loading target regions... ", Sys.time())
target_regions = read.table(target_regions_path, header = T, sep = "\t")
trinucfreq_table = read.table(target_regions_trinucfreqs_path, header = T, 
                              sep = "\t")
neutral_regions = read.table(regions_neutralbins_file, header = T, 
                             sep = "\t")
trinucfreq_table_neutral = read.table(trinucfreq_neutralbins_file, header = T, 
                                      sep="\t")

# 1. Loading the mutations and extracting the trinucleotide context -----------
message("Loading mutations... ", Sys.time())
chr_list = as.character(c(1:22, "X", "Y", paste0('chr', c(1:22, 'X', 'Y'))))
# Loading the file
mutations = read.table(mutations_file, header = F, sep = "\t", 
                       stringsAsFactors = F) 
mutations = mutations[,1:5]
colnames(mutations) = c("sampleID","chr","pos","ref","mut")
mutations = mutations[as.character(mutations$chr) %in% chr_list, ]
mutations$pos = as.numeric(mutations$pos)

# Extracting the trinucleotide context of each substitution
message("Indels... ", Sys.time())
indels_pos = as.vector(mutations$ref)=="-" | as.vector(mutations$mut)=="-" | 
             nchar(as.vector(mutations$ref))!=1 | 
             nchar(as.vector(mutations$mut))!=1
indels = mutations[indels_pos,]
subs = mutations[!indels_pos,]

message("Trinucleotides... ", Sys.time())
nt = c("A","C","G","T")
base1 = rep(nt, each = 16, times = 1)
base2 = rep(nt, each = 4, times = 4)
base3 = rep(nt, each = 1, times = 16)
trinuc_list = paste(base1,base2,base3, sep="")

subsGR <- GRanges(as.vector(subs$chr), IRanges(subs$pos - 1, subs$pos + 1))
seqlevelsStyle(subsGR) <- refGenChrStyle
seqs = scanFa(genomeFile, subsGR)
rm(subsGR)
muts_trinuc = as.vector(seqs)
subs$trinuc_ref = muts_trinuc

if (nrow(indels) > 0) { indels$trinuc_ref = NA }

# Annotating unique indels: indels with the same coordinates and same ref>mut 
# event are flagged as duplicates
indels$unique_indels = rownames(indels) %in% rownames(unique(indels[,2:5]))
subs$unique_indels = TRUE

mutations = rbind(subs,indels)
mutations = mutations[order(mutations$sampleID, mutations$chr, mutations$pos),]

# Creating the 2 global L_matrix
message("Defining L matrix... ", Sys.time())
L_matrix_ref = array(0, c(64, 4))
rownames(L_matrix_ref) = trinuc_list
colnames(L_matrix_ref) = nt
for (j in 1:64) { L_matrix_ref[j,base2[j]] = NA }

trin_freqs = colSums(trinucfreq_table[,trinuc_list])
L_matrix_global = L_matrix_ref
L_matrix_global[names(trin_freqs),] = L_matrix_global[names(trin_freqs),] + 
                                      array(rep(trin_freqs,4), 
                                            dim = c(length(trin_freqs),4))

trin_freqs = colSums(trinucfreq_table_neutral[,trinuc_list])
L_matrix_global_neutral = L_matrix_ref
L_matrix_global_neutral[names(trin_freqs),] = L_matrix_global_neutral[names(trin_freqs),] + 
                                              array(rep(trin_freqs,4), 
                                                    dim = c(length(trin_freqs),4))

# 2. Mapping the mutations to the regions and to neutral bins------------------
# Subfunction intersecting the mutations with the regions of interest
# a. Mapping the mutations to the regions of interest
message("Mapping mutations... ", Sys.time())
olm = tryCatch(map_mutations(mutations, target_regions, refGenChrStyle),
               error = function(c) {returnNULL(c)}, 
               warning = function(c) {returnNULL(c)},
               message = function(c) {returnNULL(c)})
if (nrow(olm) == 0) {
  write.table(olm, file = paste(out_prefix,"-Selection_output.txt", sep = ""),
              sep = "\t", row.names = F, col.names = T, append = F, 
              quote = F)
  stop_quietly()
}

mutations$target_region = 0
m1 = mutations[olm[,1],] # Duplicating subs if they hit more than one region
m1$target_region = target_regions[olm[,2],1] # Annotating hit region
m2 = mutations[-unique(olm[,1]),] # Mutations not mapping to any element
mutations = rbind(m1,m2)

# b. Mapping the remaining mutations to the neutral regions
message("Mapping mutations (2)... ", Sys.time())
olm = tryCatch(map_mutations(mutations, neutral_regions, refGenChrStyle), 
               error = function(c) {returnNULL(c)}, 
               warning = function(c) {returnNULL(c)},
               message = function(c) {returnNULL(c)})
if (nrow(olm) == 0) {
  write.table(olm, file = paste(out_prefix,"-Selection_output.txt", sep = ""),
              sep = "\t", row.names = F, col.names = T, append = F, 
              quote = F)
  stop_quietly()
}
mutations$neutral_region = 0
m1 = mutations[olm[,1],] # Duplicating subs if they hit more than one region
m1$neutral_region = neutral_regions[olm[,2],1] # Annotating hit region
m2 = mutations[-unique(olm[,1]),] # Mutations not mapping to any element
mutations = rbind(m1,m2)

# c. Masking out duplicate (non-unique) indels (if desired)
if (unique_indelsites_FLAG==1) { # We mask out duplicate indels
    mutations$target_region[mutations$unique_indels==F] = 0
    mutations$neutral_region[mutations$unique_indels==F] = 0
}

# c. Subsampling mutations to a maximum of nmax mutations per sample per region of interest
sample_and_region = cbind(as.vector(mutations$sampleID), 
                          as.character(mutations$target_region))
maskmuts = tryCatch(cap_mutations(sample_and_region, 
                                  max_num_muts_perRegion_perSample), 
                    error = function(c) {returnNULL(c)}, 
                    warning = function(c) {returnNULL(c)},
                    message = function(c) {returnNULL(c)})
if (nrow(maskmuts) == 0) {
  write.table(maskmuts, 
              file = paste(out_prefix,"-Selection_output.txt", sep = ""),
              sep = "\t", row.names = F, col.names = T, append = F, 
              quote = F)
  stop_quietly()
}
mutations$masked_muts_bynmax = maskmuts
mutations$target_region[maskmuts==1] = 0

# 3. Calculate global trinucleotide rates and local density of mutations-------
mutations$trinuc_ref[!(as.vector(mutations$trinuc_ref) %in% trinuc_list)] = NA
n_matrix_global = L_matrix_ref
n_matrix_global_neutral = L_matrix_ref

subsin = (!is.na(mutations$trinuc_ref)) & (mutations$target_region!=0)
aux = cbind(as.vector(mutations$trinuc_ref[subsin]),
            as.vector(mutations$mut[subsin]))
for (j in 1:dim(aux)[1]) {
    n_matrix_global[aux[j,1],aux[j,2]] = n_matrix_global[aux[j,1],aux[j,2]] + 1
}
numindels_target = sum((is.na(mutations$trinuc_ref)) & 
                         (mutations$target_region != 0))

subsin = (!is.na(mutations$trinuc_ref)) & (mutations$neutral_region!=0) & 
  (mutations$target_region==0)
aux = cbind(as.vector(mutations$trinuc_ref[subsin]), 
            as.vector(mutations$mut[subsin]))
for (j in 1:dim(aux)[1]) {
    n_matrix_global_neutral[aux[j,1],aux[j,2]] = n_matrix_global_neutral[aux[j,1],aux[j,2]] + 1
}
numindels_neutral = sum((is.na(mutations$trinuc_ref)) & 
                          (mutations$neutral_region!=0))

# [Maria's edit] estimate global coefficients for substitutions ---------------
global_mle_subs = fit_global_NBR_model(global_l_matrix = L_matrix_global,
                                       global_l_matrix_neutral = L_matrix_global_neutral, 
                                       global_n_matrix = n_matrix_global,
                                       global_n_matrix_neutral = n_matrix_global_neutral)
saveRDS(global_mle_subs, paste0(out_prefix, "-globalRates.Rds"))
write.table(global_mle_subs$par, paste0(out_prefix, "-globalRates.csv"),
            append = F, col.names = T, sep = ',', quote = F, row.names = F)
        
# 3a. Rates -------------------------------------------------------------------
rate_names = array(NA, dim = c(64,4))
for (j in 1:dim(L_matrix_ref)[1]) {
    for (h in 1:dim(L_matrix_ref)[2]) {
        rate_names[j,h] = sprintf("%s>%s%s%s", trinuc_list[j], 
                                  substr(trinuc_list[j], 1, 1), nt[h],
                                  substr(trinuc_list[j], 3, 3))
    }
}

rates_target = c(n_matrix_global / L_matrix_global)
names(rates_target) = c(rate_names)
rates_neutral = c(n_matrix_global_neutral/L_matrix_global_neutral)
names(rates_neutral) = c(rate_names)

indelsrate_target = numindels_target / sum(L_matrix_global,na.rm=T)*3
indelsrate_neutral = numindels_neutral / sum(L_matrix_global_neutral, na.rm=T)*3

# b. Local density of mutations -----------------------------------------------
# Expected number of subs and indels
targetregions_df = data.frame(region=as.vector(trinucfreq_table[,1]))
aux = strsplit(as.vector(trinucfreq_table[,2]), split=":")
targetregions_df$chr = sapply(aux, function(x) x[1])
#aux2 = sapply(aux, function(x) min(suppressWarnings(as.numeric(unlist(strsplit(x[-1], split="-")))), na.rm=T))

targetregions_df$start = sapply(aux, function(x) min(suppressWarnings(as.numeric(unlist(strsplit(x[-1], split="-")))), na.rm=T))
targetregions_df$end = sapply(aux, function(x) max(suppressWarnings(as.numeric(unlist(strsplit(x[-1], split="-")))), na.rm=T))

#targetregions_df$start = sapply(aux2, function(x) min(x,na.rm=T))
#targetregions_df$end = sapply(aux2, function(x) max(x,na.rm=T))
# Shouldn't the previous two lines be replaced by the following two instead? At least on two occasions I have needed to run it like this...
#targetregions_df$start = apply(aux2,2,min)
#targetregions_df$end = apply(aux2,2,max)


tf = as.matrix(trinucfreq_table[,trinuc_list])
targetregions_df$exp_subs = apply(tf, 1, function(x) sum(rep(x,4)*rates_target, na.rm=T) )
targetregions_df$exp_indels = apply(tf, 1, function(x) sum(x)*indelsrate_target )

neutralregions_df = data.frame(region=as.vector(trinucfreq_table_neutral[,1:3]))
tf = as.matrix(trinucfreq_table_neutral[,trinuc_list])
neutralregions_df$exp_subs = apply(tf, 1, function(x) sum(rep(x,4)*rates_neutral, na.rm=T) )
neutralregions_df$exp_indels = apply(tf, 1, function(x) sum(x)*indelsrate_neutral )

# Observed number of subs and indels

targetregions_df$obs_subs = 0
targetregions_df$obs_indels = 0
indel_pos = is.na(mutations$trinuc_ref)
numsubs = table( mutations$target_region[mutations$target_region > 0 & !indel_pos] )
targetregions_df$obs_subs[as.numeric(names(numsubs))] = numsubs
numinds = table( mutations$target_region[mutations$target_region > 0 & indel_pos] )
targetregions_df$obs_indels[as.numeric(names(numinds))] = numinds

neutralregions_df$obs_subs = 0
neutralregions_df$obs_indels = 0
indel_pos = is.na(mutations$trinuc_ref)
numsubs = table( mutations$neutral_region[mutations$neutral_region > 0 & !indel_pos] )
neutralregions_df$obs_subs[as.numeric(names(numsubs))] = numsubs
numinds = table( mutations$neutral_region[mutations$neutral_region > 0 & indel_pos] )
neutralregions_df$obs_indels[as.numeric(names(numinds))] = numinds


# Estimating the neighbourhood "t" for every target region
# We choose as neighbouring neutral regions of a given target region all those that are 
# within "neighbourhood_localrate" distance of the target region. We consider any overlap
# between the segments as valid, which means that the neighbourhood of a target region
# will always be contained within the interval "neighbourhood_localrate + 100kb" around
# the target region (for a 100kb binning of the genome in the neutral reference)
# For example, if neighbourhood_localrate=1e5, only regions less or equal to 200kb away 
# from the ends of the target region are considered in the calculation of the local rate.

neighbourhood_localrate = ceiling(0.001/dim(mutations)[1]*3e9)*100000

rangesTARGET = GRanges(target_regions[,2], IRanges(as.numeric(target_regions[,3])-neighbourhood_localrate, as.numeric(target_regions[,4])+neighbourhood_localrate))
seqlevelsStyle(rangesTARGET) <- refGenChrStyle
rangesNEUTRAL = GRanges(neutralregions_df[,1], IRanges(as.numeric(neutralregions_df[,2]), as.numeric(neutralregions_df[,3])))
seqlevelsStyle(rangesNEUTRAL) <- refGenChrStyle
ol = findOverlaps(rangesTARGET, rangesNEUTRAL, type="any", select="all")
olmatrix = as.matrix(ol)
neighbours = unique(cbind(target_regions[olmatrix[,1],1],olmatrix[,2]))

targetregions_df$local_t = NA
targetregions_df$local_t_indels = NA

for (j in 1:dim(targetregions_df)[1]) {
    neutral_bins = neighbours[neighbours[,1]==j,2]
    targetregions_df$local_t[j] = sum(neutralregions_df$obs_subs[neutral_bins])/
                                  sum(neutralregions_df$exp_subs[neutral_bins])
    targetregions_df$local_t_indels[j] = sum(neutralregions_df$obs_indels[neutral_bins])/
                                         sum(neutralregions_df$exp_indels[neutral_bins])
    if ((j/1000) == round(j/1000)) { 
      print(paste0('Progress: ', 100 * j/dim(targetregions_df)[1], '%'))
    }
}

# 4. Negative binomial regression with local density of mutations and covariates ----
# Masking out regions with ZERO expected values
targetregions_df$exp_subs[targetregions_df$exp_subs == 0] = NA
targetregions_df$exp_indels[targetregions_df$exp_indels == 0] = NA

## 4a. MODEL 1: No use of local mutation rates or covariates ------------------
# Negative binomial regression
model_subs = glm.nb(formula = targetregions_df$obs_subs ~ offset(log(targetregions_df$exp_subs)) -1 )
nb_size_subs = model_subs$theta
if (numindels_target>0) {
    model_indels = glm.nb(formula = targetregions_df$obs_indels ~ offset(log(targetregions_df$exp_indels)) -1 )
    nb_size_indels = model_indels$theta
}

# P-values (Neg Binom)
targetregions_df$pval_subs = NA
targetregions_df$pval_indels = NA
targetregions_df$pval_both = NA

for (j in 1:dim(targetregions_df)[1]) {
    targetregions_df$pval_subs[j] = pnbinom(q=targetregions_df$obs_subs[j]-0.1,
                                            mu=targetregions_df$exp_subs[j], 
                                            size=nb_size_subs, 
                                            lower.tail=F)
    
    if (numindels_target>0) {
        targetregions_df$pval_indels[j] = pnbinom(q=targetregions_df$obs_indels[j]-0.1,
                                                  mu=targetregions_df$exp_indels[j],
                                                  size=nb_size_indels, 
                                                  lower.tail=F)
        # Fisher combined p-value
        p_vec = c(targetregions_df$pval_subs[j], targetregions_df$pval_indels[j])
        targetregions_df$pval_both[j] = 1-pchisq(-2*sum(log(p_vec)),length(p_vec)*2)
    } else {
        # We use only subs
        targetregions_df$pval_both[j] = targetregions_df$pval_subs[j]
    }
    if (round(j/1000)==(j/1000)) { print(j/dim(targetregions_df)[1]) }
}
# Adjusted q-value
targetregions_df$qval_subs = p.adjust(targetregions_df$pval_subs,
                                      method="BH") 
targetregions_df$qval_indels = p.adjust(targetregions_df$pval_indels,
                                        method="BH")
targetregions_df$qval_both = p.adjust(targetregions_df$pval_both,
                                      method="BH")

## 4b. MODEL 2: Using local mutation rates and covariates ---------------------
# Excluding elements overlapping "driver" genomic regions from the negbin
# background fit
gr_elements = GRanges(targetregions_df$chr, 
                      IRanges(targetregions_df$start, targetregions_df$end))
seqlevelsStyle(gr_elements) <- refGenChrStyle
olmatrix = as.matrix(findOverlaps(gr_elements, gr_drivers, type = "any",
                                  select="all"))
exclude_forfit = (1:nrow(targetregions_df)) %in% unique(olmatrix[,1])
targetregions_df$exclude_forfit = exclude_forfit

# Negative binomial regression
model_subs = glm.nb(formula = obs_subs ~ offset(log(exp_subs)) + local_t, 
                    data = targetregions_df[!exclude_forfit,])
nb_size_subs_cv = model_subs$theta
targetregions_df$cv_predicted_subs = exp(predict(model_subs, 
                                                 newdata = targetregions_df))

if (numindels_target>0) {
  ## FOR INDELS THE LOCAL RATE IS A COVARIATE
  model_indels = glm.nb(formula = obs_indels ~ offset(log(exp_indels)) + 
                          local_t_indels,
                        data = targetregions_df[!exclude_forfit,])
  nb_size_indels_cv = model_indels$theta
  targetregions_df$cv_predicted_indels = exp(predict(model_indels, 
                                                     newdata = targetregions_df))
}

# P-values (Neg Binom)
targetregions_df$pval_subs_CV = NA
targetregions_df$pval_indels_CV = NA
targetregions_df$pval_both_CV = NA

for (j in 1:dim(targetregions_df)[1]) {
    targetregions_df$pval_subs_CV[j] = pnbinom(q=targetregions_df$obs_subs[j]-0.1, mu=targetregions_df$cv_predicted_subs[j], size=nb_size_subs_cv, lower.tail=F)
    #targetregions_df$pval_indels_CV[j] = pnbinom(q=targetregions_df$obs_indels[j]-0.1, mu=targetregions_df$cv_predicted_indels[j], size=nb_size_indels_cv, lower.tail=F)
    
    if (numindels_target>0) {
        targetregions_df$pval_indels_CV[j] = pnbinom(q=targetregions_df$obs_indels[j]-0.1, mu=targetregions_df$cv_predicted_indels[j], size=nb_size_indels_cv, lower.tail=F)
        # Fisher combined p-value
        p_vec = c(targetregions_df$pval_subs_CV[j], targetregions_df$pval_indels_CV[j])
        targetregions_df$pval_both_CV[j] = 1-pchisq(-2*sum(log(p_vec)),length(p_vec)*2)
    } else {
        # We use only subs
        targetregions_df$pval_both[j] = targetregions_df$pval_subs[j]
    }
    if (round(j/1000)==(j/1000)) { print(j/dim(targetregions_df)[1]) }
}
targetregions_df$qval_subs_CV = p.adjust(targetregions_df$pval_subs_CV, method="BH") # Adjusted q-value
targetregions_df$qval_indels_CV = p.adjust(targetregions_df$pval_indels_CV, method="BH") # Adjusted q-value
targetregions_df$qval_both_CV = p.adjust(targetregions_df$pval_both_CV, method="BH") # Adjusted q-value


# Tiered FDR correction

targetregions_df$qval_both_CV_tiered = NA
inds = targetregions_df$exclude_forfit
targetregions_df$qval_both_CV_tiered[inds] = p.adjust(targetregions_df$pval_both_CV[inds], method="BH") # Adjusted q-value
targetregions_df$qval_both_CV_tiered[!inds] = p.adjust(targetregions_df$pval_both_CV[!inds], method="BH") # Adjusted q-value
targetregions_df = targetregions_df[order(targetregions_df$qval_both_CV_tiered),]

write.table(targetregions_df, file=paste(out_prefix,"-Selection_output.txt",sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, append=FALSE, quote=FALSE)


# Plot of the underlying gamma distributions of the 4 models (subs and indels with and without covariates)

pdf(paste(out_prefix,"-Underlying_gamma_distributions.pdf",sep=""),height=4,width=4)
xvec = seq(0,4,by=0.0001)
plot(xvec,0.0001*dgamma(x=xvec,shape=nb_size_subs_cv,rate=nb_size_subs_cv),type="l",lty=1,xlab="Relative mutation rate",ylab="",col="cadetblue",las=1,main="PDFs of the underlying Gamma distributions",cex.axis=0.8,cex.main=0.85)
lines(xvec,0.0001*dgamma(x=xvec,shape=nb_size_subs,rate=nb_size_subs),lty=2,col="cadetblue")
if (numindels_target>0) { 
  lines(xvec,0.0001*dgamma(x=xvec,shape=nb_size_indels,rate=nb_size_subs),
        lty=2,col="chocolate")
  lines(xvec,0.0001*dgamma(x=xvec,shape=nb_size_indels_cv,
                           rate=nb_size_subs_cv),lty=1,col="chocolate")
  legend("topright", lty=c(2,1,2,1),
         col=c("cadetblue","cadetblue","chocolate","chocolate"), 
         title="Size parameters", 
         legend=c(sprintf("subs=%0.3g",nb_size_subs),
                  sprintf("subs_cv=%0.3g",nb_size_subs_cv),
                  sprintf("indels=%0.3g",nb_size_indels),
                  sprintf("indels_cv=%0.3g",nb_size_indels_cv)), 
         box.col=NA, cex=0.8)
  
} else {
  legend("topright", lty=c(2,1,2,1),
         col=c("cadetblue","cadetblue"), 
         title="Size parameters", 
         legend=c(sprintf("subs=%0.3g",nb_size_subs),
                  sprintf("subs_cv=%0.3g",nb_size_subs_cv)), 
         box.col=NA, cex=0.8)
}
dev.off()


## Confidence intervals for the obs/exp ratios in NBR (23.12.2016)
calculate_CI95_flag = 1

if (calculate_CI95_flag == 1) {
  ci95_subs = as.matrix(targetregions_df[, c("obs_subs","cv_predicted_subs")])
  ci95_subs = t(apply(ci95_subs, 1, 
                      function(x) tryCatch(ci95nbr(x[1], x[2],
                                                   nb_size_indels_cv),
                                           error=function(err) rep(NA,3)))) 
  targetregions_df$obsexp_subs_mle = ci95_subs[,1]
  targetregions_df$obsexp_subs_low = ci95_subs[,2]
  targetregions_df$obsexp_subs_high = ci95_subs[,3]
    
  if (numindels_target>0) {
    ci95_indels = as.matrix(targetregions_df[,c("obs_indels","cv_predicted_indels")])
    ci95_indels = t(apply(ci95_indels, 1, 
                          function(x) tryCatch(ci95nbr(x[1], x[2],
                                                       nb_size_indels_cv),
                                               error=function(err) rep(NA,3))))
    targetregions_df$obsexp_indels_mle = ci95_indels[,1]
    targetregions_df$obsexp_indels_low = ci95_indels[,2]
    targetregions_df$obsexp_indels_high = ci95_indels[,3]
  } else {
    targetregions_df$obsexp_indels_mle = NA
    targetregions_df$obsexp_indels_low = NA
    targetregions_df$obsexp_indels_high = NA
  }
  write.table(targetregions_df, file=paste(out_prefix,"-Selection_output.txt",sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, append=FALSE, quote=FALSE)
}