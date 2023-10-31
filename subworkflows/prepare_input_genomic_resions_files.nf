/* 
    A workflow to create input genomic regions files for DIGdriver, dNdScv 
    (Rda), NBR and OncodriveFML.
*/

include { FILTER_GENOMIC_REGIONS } from '../modules/prepare_genomic_regions.nf'
include { CREATE_RDA_FOR_DNDSCV_DIGDRIVER } from '../modules/prepare_genomic_regions.nf'
include { WRITE_REGIONS_FOR_DIGDRIVER } from '../modules/prepare_genomic_regions.nf'
include { WRITE_REGIONS_FOR_NBR } from '../modules/prepare_genomic_regions.nf'
include { WRITE_REGIONS_FOR_ONCODRIVEFML } from '../modules/prepare_genomic_regions.nf'

workflow PREPARE_INPUT_GENOMIC_REGIONS_FILES {
    take:
    analysis_inv
    blacklist_inv
    target_genome_fasta
    target_genome_chr_len
    chain
    inventories_pass

    main:
    // pre-process regions
    filtered_regions = FILTER_GENOMIC_REGIONS(inventories_pass.combine(analysis_inv)
                                                              .combine(blacklist_inv)
                                                              .combine(target_genome_fasta)
                                                              .combine(target_genome_chr_len)
                                                              .combine(chain)).bed.flatten()
    filtered_regions.map{ it ->
                           def tumor_subtype = it.name.toString().tokenize('-').get(1)
                           return tuple(tumor_subtype, it)
                        }.set{filtered_regions}
 
    tumor_subtypes_and_gr = analysis_inv.splitCsv(header: true)
                                        .map{row -> tuple(row.tumor_subtype, row.software, row.gr_id)}
                                        .unique()
                                        .groupTuple(by: [0, 1])
                                        .combine(filtered_regions, by: 0)

    // write regions in software format
    digdriver_regions = WRITE_REGIONS_FOR_DIGDRIVER(tumor_subtypes_and_gr)
    nbr_regions = WRITE_REGIONS_FOR_NBR(tumor_subtypes_and_gr)
    oncodrivefml_regions = WRITE_REGIONS_FOR_ONCODRIVEFML(tumor_subtypes_and_gr)

    // create Rda files with reference genome information for dNdScv and DIGdriver
    gtf_for_rda = analysis_inv.splitCsv(header: true)
                              .map{row -> tuple(row.tumor_subtype, row.software, 
                                                row.gr_code, row.gr_file, row.gr_genome)}
                              .unique()
                              .map { it ->
                                  if ((it[1] == 'digdriver' || it[1] == 'dndscv') && it[2] == 'CDS') {
                                      return it
                                  }
                              }
                              .map { it ->
                                  return tuple(it[0], it[3], it[4])
                              }
                              .unique()
                              .groupTuple(by: [0])
                              .combine(filtered_regions, by: [0])   
    ch_dndscv_digdriver_rda = CREATE_RDA_FOR_DNDSCV_DIGDRIVER(gtf_for_rda.combine(target_genome_fasta)
                                                                         .combine(chain)).rda

    emit:
    digdriver = digdriver_regions
    dndscv_digdriver_rda = ch_dndscv_digdriver_rda
    nbr = nbr_regions
    oncodrivefml = oncodrivefml_regions
}