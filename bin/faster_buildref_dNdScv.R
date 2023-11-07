#!/usr/bin/env Rscript

# FILE: faster_buildref_dNdScv.R ----------------------------------------------
#
# DESCRIPTION: a modified, faster version of buildref function from dNdScv 
# package. Speed increase is achieved by using mclapply instead of loops and
# by introducing support for multiple cores.
#
# USAGE: source('faster_buildref_dNdScv.R') or box::use(./faster_buildref_dNdScv[...])
#
# OPTIONS: None
#
# REQUIREMENTS: 
# BUGS: --
# NOTES:
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, Cancer Institute, London, the UK
# VERSION:  1
# CREATED:  03.10.2023
# REVISION: 03.10.2023

box::use(Biostrings[...])
box::use(dndscv[...])
box::use(GenomicRanges[...])
box::use(IRanges[...])
box::use(parallel[...])
box::use(Rsamtools[...])
box::use(seqinr[...])
box::use(utils[...])

get_CDSseq = function(gr, strand, genomefile) {
  cdsseq = strsplit(paste(as.vector(Rsamtools::scanFa(genomefile, gr)), 
                          collapse = ""), "")[[1]]
  if (strand == -1) {
    cdsseq = rev(seqinr::comp(cdsseq, forceToLower = F, ambiguous = T))
  }
  return(cdsseq)
}

get_splicesites = function(cds) {
  splpos = numeric(0)
  if (nrow(cds) > 1) {
    if (cds[1, 10] == 1) {
      spl5prime = cds[-nrow(cds), 6]
      spl3prime = cds[-1, 5]
      splpos = unique(sort(c(spl5prime + 1, spl5prime + 2, spl5prime + 5, 
                             spl3prime - 1, spl3prime - 2)))
    }
    else if (cds[1, 10] == -1) {
      spl5prime = cds[-1, 5]
      spl3prime = cds[-nrow(cds), 6]
      splpos = unique(sort(c(spl5prime - 1, spl5prime - 2, spl5prime - 5, 
                             spl3prime + 1, spl3prime + 2)))
    }
  }
  return(splpos)
}

get_spliceseq = function(gr, strand, genomefile) {
  spliceseq = unname(as.vector(Rsamtools::scanFa(genomefile, gr)))
  if (strand == -1) {
    spliceseq = seqinr::comp(spliceseq, forceToLower = F, ambiguous = T)
  }
  return(spliceseq)
}

build_RefCDS_one_gene <- function(gene_cdss, cdsList, genomefile, numcode) {
  result <- list()
  h = keeptrying = 1
  
  while (h <= nrow(gene_cdss) & keeptrying) {
    print(keeptrying)
    pid = gene_cdss[h, 3]
    cds = cdsList[[pid]]
    strand = cds[1, 10]
    chr = cds[1, 4]
    gr = GRanges(chr, IRanges(cds[, 5], cds[, 6]))
    cdsseq = get_CDSseq(gr, strand, genomefile)
    pseq = seqinr::translate(cdsseq, numcode = numcode)
    if (all(pseq[-length(pseq)] != "*") & all(cdsseq != "N")) {
      splpos = get_splicesites(cds)
      if (length(splpos) > 0) {
        gr_spl = GRanges(chr, IRanges(splpos, splpos))
        splseq = get_spliceseq(gr_spl, strand, genomefile)
      }
      if (strand == 1) {
        cdsseq1up = get_CDSseq(GRanges(chr,
                                       IRanges(cds[, 5] - 1, cds[, 6] - 1)),
                               strand, genomefile)
        cdsseq1down = get_CDSseq(GRanges(chr, 
                                         IRanges(cds[, 5] + 1, cds[, 6] + 1)),
                                 strand, genomefile)
        if (length(splpos) > 0) {
          splseq1up = get_spliceseq(GRanges(chr, 
                                            IRanges(splpos - 1, splpos - 1)), 
                                    strand, genomefile)
          splseq1down = get_spliceseq(GRanges(chr, 
                                              IRanges(splpos + 1, splpos + 1)), 
                                      strand, genomefile)
        }
      }
      else if (strand == -1) {
        cdsseq1up = get_CDSseq(GRanges(chr, 
                                       IRanges(cds[, 5] + 1, cds[, 6] + 1)),
                               strand, genomefile)
        cdsseq1down = get_CDSseq(GRanges(chr,
                                         IRanges(cds[, 5] - 1, cds[, 6] - 1)),
                                 strand, genomefile)
        if (length(splpos) > 0) {
          splseq1up = get_spliceseq(GRanges(chr, 
                                            IRanges(splpos + 1, splpos + 1)), 
                                    strand, genomefile)
          splseq1down = get_spliceseq(GRanges(chr, 
                                              IRanges(splpos - 1, splpos - 1)), 
                                      strand, genomefile)
        }
      }
      result$gene_name = gene_cdss[h, 2]
      result$gene_id = gene_cdss[h, 1]
      result$protein_id = gene_cdss[h, 3]
      result$CDS_length = gene_cdss[h, 4]
      result$chr = cds[1, 4]
      result$strand = strand
      result$intervals_cds = unname(as.matrix(cds[, 5:6]))
      result$intervals_splice = splpos
      result$seq_cds = Biostrings::DNAString(paste(cdsseq, collapse = ""))
      result$seq_cds1up = Biostrings::DNAString(paste(cdsseq1up,
                                                      collapse = ""))
      result$seq_cds1down = Biostrings::DNAString(paste(cdsseq1down, 
                                                        collapse = ""))
      if (length(splpos) > 0) {
        result$seq_splice = Biostrings::DNAString(paste(splseq, collapse = ""))
        result$seq_splice1up = Biostrings::DNAString(paste(splseq1up,
                                                           collapse = ""))
        result$seq_splice1down = Biostrings::DNAString(paste(splseq1down, 
                                                             collapse = ""))
      }
      keeptrying = 0
    }
    h = h + 1
  }
  result
}

calculate_impact <- function(oneGeneRefCDS, nt, impMatrix, trinucSubsInd, 
                             trinucInd) {
  L = array(0, dim = c(192, 4))
  
  cdsseq = as.character(as.vector(oneGeneRefCDS$seq_cds))
  cdsseq1up = as.character(as.vector(oneGeneRefCDS$seq_cds1up))
  cdsseq1down = as.character(as.vector(oneGeneRefCDS$seq_cds1down))
  
  ind = rep(1:length(cdsseq), each = 3)
  old_trinuc = paste(cdsseq1up[ind], cdsseq[ind], cdsseq1down[ind], sep = "")
  new_base = c(sapply(cdsseq, function(x) nt[nt != x]))
  new_trinuc = paste(cdsseq1up[ind], new_base, cdsseq1down[ind], sep = "")
  codon_start = rep(seq(1, length(cdsseq), by = 3), each = 9)
  old_codon = paste(cdsseq[codon_start], cdsseq[codon_start + 1], 
                    cdsseq[codon_start + 2], sep = "")
  pos_in_codon = rep(rep(1:3, each = 3), length.out = length(old_codon))
  aux = strsplit(old_codon, "")
  
  new_codon = sapply(1:length(old_codon), function(x) {
    new_codonx = aux[[x]]
    new_codonx[pos_in_codon[x]] = new_base[x]
    return(new_codonx)
  })
  new_codon = paste(new_codon[1, ], new_codon[2, ], new_codon[3, ], sep = "")
  imp = impMatrix[(trinucInd[new_codon] - 1) * 64 + trinucInd[old_codon]]
  matrind = trinucSubsInd[paste(old_trinuc, new_trinuc, sep = ">")]
  matrix_ind = table(matrind[which(imp == 1)])
  L[as.numeric(names(matrix_ind)), 1] = matrix_ind
  matrix_ind = table(matrind[which(imp == 2)])
  L[as.numeric(names(matrix_ind)), 2] = matrix_ind
  matrix_ind = table(matrind[which(imp == 3)])
  L[as.numeric(names(matrix_ind)), 3] = matrix_ind
  if (length(oneGeneRefCDS$intervals_splice) > 0) {
    splseq = as.character(as.vector(oneGeneRefCDS$seq_splice))
    splseq1up = as.character(as.vector(oneGeneRefCDS$seq_splice1up))
    splseq1down = as.character(as.vector(oneGeneRefCDS$seq_splice1down))
    old_trinuc = rep(paste(splseq1up, splseq, splseq1down, sep = ""), each = 3)
    new_trinuc = paste(rep(splseq1up, each = 3), 
                       c(sapply(splseq, function(x) nt[nt != x])), 
                       rep(splseq1down, each = 3), sep = "")
    matrind = trinucSubsInd[paste(old_trinuc, new_trinuc, sep = ">")]
    matrix_ind = table(matrind)
    L[as.numeric(names(matrix_ind)), 4] = matrix_ind
  }
  oneGeneRefCDS$L <- L
  oneGeneRefCDS
}

buildref_faster <- function(cdsfile, genomefile, numcode = 1, 
                            excludechrs = NULL, onlychrs = NULL,
                            useids = F, cores = 1) {
  message(Sys.time(), " [1/3] Preparing the environment...")
  reftable = read.table(cdsfile, header = 1, sep = "\t", stringsAsFactors = F, 
                        quote = "\"", na.strings = "-", fill = T)
  colnames(reftable) = c("gene.id", "gene.name", "cds.id", "chr",
                         "chr.coding.start", "chr.coding.end", "cds.start", 
                         "cds.end", "length", "strand")
  reftable[, 5:10] = suppressWarnings(lapply(reftable[, 5:10], as.numeric))
  longname = paste(reftable$gene.id, reftable$gene.name, sep = ":")
  if (useids == T) {
    reftable$gene.name = longname
  }
  if (length(unique(reftable$gene.name)) < length(unique(longname))) {
    warning(sprintf("%0.0f unique gene IDs (column 1) found. %0.0f unique gene names (column 2) found. Consider combining gene names and gene IDs or replacing gene names by gene IDs to avoid losing genes (see useids argument in ? buildref)", 
                    length(unique(reftable$gene.id)), 
                    length(unique(reftable$gene.name))))
  }
  validchrs = as.character(GenomicRanges::seqnames(Rsamtools::scanFaIndex(genomefile)))
  validchrs = setdiff(validchrs, excludechrs)
  if (length(onlychrs) > 0) {
    validchrs = validchrs[validchrs %in% onlychrs]
  }
  if (any(validchrs %in% unique(reftable$chr))) {
    validchrs = validchrs[validchrs %in% unique(reftable$chr)]
  } else {
    reftable$chr = paste("chr", reftable$chr, sep = "")
    validchrs = validchrs[validchrs %in% unique(reftable$chr)]
    if (length(validchrs) == 0) {
      stop("No chromosome names in common between the genome file and the CDS table")
    }
  }
  reftable = reftable[reftable[, 1] != "" & reftable[, 2] != "" & 
                        reftable[, 3] != "" & !is.na(reftable[, 5]) & 
                        !is.na(reftable[, 6]), ]
  reftable = reftable[which(reftable$chr %in% validchrs), ]
  transc_gr = GRanges(reftable$chr, IRanges(reftable$chr.coding.start, 
                                             reftable$chr.coding.end))
  chrs_gr = Rsamtools::scanFaIndex(genomefile)
  ol = as.data.frame(GenomicRanges::findOverlaps(transc_gr, chrs_gr, 
                                                 type = "within", 
                                                 select = "all"))
  if (length(unique(ol[, 1])) < nrow(reftable)) {
    stop(sprintf("Aborting buildref. %0.0f rows in cdsfile have coordinates that fall outside of the corresponding chromosome length. Please ensure that you are using the same assembly for the cdsfile and genomefile", 
                 nrow(reftable) - length(unique(ol[, 1]))))
  }
  reftable = reftable[unique(ol[, 1]), ]
  fullcds = intersect(reftable$cds.id[reftable$cds.start == 1],
                      reftable$cds.id[reftable$cds.end == reftable$length])
  ol_start = as.data.frame(GenomicRanges::findOverlaps(transc_gr, chrs_gr, 
                                                       type = "start", 
                                                       select = "all"))[, 1]
  if (any(ol_start)) {
    reftable[ol_start, "chr.coding.start"] = reftable[ol_start, 
                                                      "chr.coding.start"] + 3
    reftable[ol_start, "cds.start"] = reftable[ol_start, "cds.start"] + 3
  }
  ol_end = as.data.frame(GenomicRanges::findOverlaps(transc_gr, chrs_gr, 
                                                     type = "end", 
                                                     select = "all"))[, 1]
  if (any(ol_end)) {
    reftable[ol_end, "chr.coding.end"] = reftable[ol_end, "chr.coding.end"] - 3
    reftable[ol_end, "cds.end"] = reftable[ol_end, "cds.end"] - 3
  }
  if (any(c(ol_start, ol_end))) {
    warning(sprintf("The following genes were found to start or end at the first or last base of their contig. Since dndscv needs trinucleotide contexts for all coding bases, codons overlapping ends of contigs have been trimmed. Affected genes: %s.", 
                    paste(reftable[unique(c(ol_start, ol_end)), "gene.name"], 
                          collapse = ", ")))
  }
  cds_table = unique(reftable[, c(1:3, 9)])
  cds_table = cds_table[order(cds_table$gene.name, -cds_table$length), ]
  cds_table = cds_table[(cds_table$length%%3) == 0, ]
  cds_table = cds_table[cds_table$cds.id %in% fullcds, ]
  reftable = reftable[reftable$cds.id %in% fullcds, ]
  gene_list = unique(cds_table$gene.name)
  reftable = reftable[order(reftable$chr, reftable$chr.coding.start), ]
  cds_split = split(reftable, f = reftable$cds.id)
  gene_split = split(cds_table, f = cds_table$gene.name)
  
  message(Sys.time(), " [2/3] Building the RefCDS object...")
  RefCDS <- mclapply(1:length(gene_split),
                     function(idx) build_RefCDS_one_gene(gene_split[[idx]],
                                                         cds_split, 
                                                         genomefile, 
                                                         numcode),
                       mc.cores = cores)
  RefCDS <- RefCDS[sapply(RefCDS, length) != 0]
  names(RefCDS) <- sapply(RefCDS, function(x) x$gene_name)
  
  message(Sys.time(), " [3/3] Calculating the impact of all possible coding changes...")
  nt = c("A", "C", "G", "T")
  trinuc_list = paste(rep(nt, each = 16, times = 1), 
                      rep(nt, each = 4, times = 4), 
                      rep(nt, each = 1, times = 16), sep = "")
  trinuc_ind = structure(1:64, names = trinuc_list)
  trinuc_subs = NULL
  for (j in 1:length(trinuc_list)) {
    trinuc_subs = c(trinuc_subs, 
                    paste(trinuc_list[j], 
                          paste(substr(trinuc_list[j], 1, 1), 
                                setdiff(nt, substr(trinuc_list[j], 2, 2)), 
                                substr(trinuc_list[j], 3, 3), sep = ""), 
                          sep = ">"))
  }
  trinuc_subsind = structure(1:192, names = trinuc_subs)
  impact_matrix = array(NA, dim = c(64, 64))
  colnames(impact_matrix) = rownames(impact_matrix) = trinuc_list
  for (j in 1:64) {
    for (h in 1:64) {
      from_aa = seqinr::translate(strsplit(trinuc_list[j], "")[[1]],
                                  numcode = numcode)
      to_aa = seqinr::translate(strsplit(trinuc_list[h], "")[[1]],
                                numcode = numcode)
      if (to_aa == from_aa) {
        impact_matrix[j, h] = 1
      } else if (to_aa == "*") {
        impact_matrix[j, h] = 3
      } else if ((to_aa != "*") & (from_aa != "*") & (to_aa != from_aa)) {
        impact_matrix[j, h] = 2
      } else if (from_aa == "*") {
        impact_matrix[j, h] = NA
      }
    }
  }
  saveRDS(RefCDS, 'RefCDS_before.rds')
  RefCDS <- mclapply(1:length(RefCDS),
                     function(idx) calculate_impact(RefCDS[[idx]], nt, 
                                                    impact_matrix,
                                                    trinuc_subsind,
                                                    trinuc_ind),
                     mc.cores = cores)
  saveRDS(RefCDS, 'RefCDS_after.rds')
  names(RefCDS) <- sapply(RefCDS, function(x) x$gene_name)
  aux = unlist(sapply(1:length(RefCDS), 
                      function(x) t(cbind(x, 
                                          rbind(RefCDS[[x]]$intervals_cds, 
                                                cbind(RefCDS[[x]]$intervals_splice,
                                                      RefCDS[[x]]$intervals_splice))))))
  df_genes = as.data.frame(t(array(aux, dim = c(3, length(aux)/3))))
  colnames(df_genes) = c("ind", "start", "end")
  df_genes$chr = unlist(sapply(1:length(RefCDS), 
                               function(x) rep(RefCDS[[x]]$chr, 
                                               nrow(RefCDS[[x]]$intervals_cds) +
                                                 length(RefCDS[[x]]$intervals_splice))))
  df_genes$gene = sapply(RefCDS, function(x) x$gene_name)[df_genes$ind]
  gr_genes = GRanges(df_genes$chr, IRanges(df_genes$start, df_genes$end))
  GenomicRanges::mcols(gr_genes)$names = df_genes$gene
  
  list('RefCDS' = RefCDS, 'gr_genes' = gr_genes)
}