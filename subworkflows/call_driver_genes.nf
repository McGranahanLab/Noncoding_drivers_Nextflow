def infer_tumor_subtype (filePath) {
    result = filePath.name.toString().tokenize('-').get(1)
    return result
}

def infer_genomic_region(filePath) {
    result = filePath.name.toString().tokenize('-').get(3)
    return result
}

/* 
    A workflow to create input mutation files for CHASM+, DIGdriver, dNdScv, 
    MutPanning, NBR and OncodriveFML.
*/

include { CHASMplus } from '../modules/run_software.nf'
include { DIGDRIVER } from '../modules/run_software.nf'
include { DNDSCV } from '../modules/run_software.nf'
include { MUTPANNING } from '../modules/run_software.nf'
include { NBR } from '../modules/run_software.nf'
include { ONCODRIVEFML } from '../modules/run_software.nf'

workflow RUN_CHASMplus {
    take:
    analysis_inv_path
    chasmplus_mutations
    chusmplusAnno_inv_path
    target_genome_version

    main:
    chasm_analysis_inv = analysis_inv_path.splitCsv(header: true)
                                          .map {row -> 
                                              return tuple(row.tumor_subtype,
                                                           row.gr_id,
                                                           row.software)
                                          }
                                          .unique()
    chusmplusAnno_inv = chusmplusAnno_inv_path.splitCsv(header: true)
                                              .map {row -> 
                                                  return tuple(row.tumor_subtype,
                                                               row.chasm_annotator)
                                              }
    chasmplus_results = CHASMplus(chasm_analysis_inv.combine(chasmplus_mutations, by: [0])
                                                    .filter{ it[2] == 'chasmplus' }
                                                    .combine(chusmplusAnno_inv, by: [0])).csv
    emit:
    results = chasmplus_results
}

workflow RUN_DIGDRIVER {
    take:
    digdriver_inv
    digdriver_mutations
    digdriver_regions
    dndscv_digdriver_rda
    digdriver_elements
    target_genome_fasta

    main:
	digdriver_inv = digdriver_inv.splitCsv(header: true)
                                 .map { it ->
                                     return(tuple(it.tumor_subtype, it.model_file))
                                 }
    digdriver_results = DIGDRIVER(digdriver_regions.flatten()
                                                   .map { it ->
                                                            return tuple(infer_tumor_subtype(it),
                                                                         infer_genomic_region(it),
                                                                         'digdriver', it)
                                                    }
                                                    .combine(dndscv_digdriver_rda, by: [0])
                                                    .combine(digdriver_mutations, by: [0])
                                                    .combine(digdriver_inv, by: [0])
                                                    .combine(digdriver_elements)
                                                    .combine(target_genome_fasta)).csv

    emit:
    results = digdriver_results
}

workflow RUN_DNDSCV {
    take:
    analysis_inv
    dndscv_mutations
    dndscv_digdriver_rda

    main:
    dndscv_results = DNDSCV(analysis_inv.splitCsv(header: true)
                                        .map{row -> tuple(row.tumor_subtype, row.software, row.gr_id)}
                                        .unique()
                                        .map { it ->
                                            if (it[1] == 'dndscv') {
                                                return tuple(it[0], it[2], 'dndscv')
                                            }
                                        }
                                        .combine(dndscv_digdriver_rda, by: [0])
                                        .combine(dndscv_mutations, by: [0])).csv
    emit:
    results = dndscv_results
}

workflow RUN_MUTPANNING {
    take:
    analysis_inv
    mutpanning_mutations

    main:
    mutpanning_results = MUTPANNING(analysis_inv.splitCsv(header: true)
                                                .map{row -> tuple(row.tumor_subtype, row.software, row.gr_id)}
                                                .unique()
                                                .map { it ->
                                                    if (it[1] == 'mutpanning') {
                                                         return tuple(it[0], it[2], 'mutpanning')
                                                    }
                                                }
                                                .combine(mutpanning_mutations, by: [0])).csv
    emit:
    results = mutpanning_results
}

workflow RUN_NBR {
    take:
    nbr_mutations
    nbr_regions
    target_genome_fasta
    nbr_neutral_bins
    nbr_neutral_trinucfreq
    nbr_driver_regs

    main:
    nbr_results = NBR(nbr_regions.flatten()
                                 .map { it ->
                                     return tuple(infer_tumor_subtype(it),
                                                  infer_genomic_region(it), 'nbr',
                                                  it)
                                 }
                                 .combine(nbr_mutations, by: [0])
                                 .combine(target_genome_fasta)
                                 .combine(nbr_neutral_bins)
                                 .combine(nbr_neutral_trinucfreq)
                                 .combine(nbr_driver_regs)).csv
    emit:
    results = nbr_results
}

workflow RUN_ONCODRIVEFML {
    take:
    oncodrivefml_mutations
    oncodrivefml_regions
    oncodrivefml_config

    main:
    oncodrivefml_results = ONCODRIVEFML(oncodrivefml_regions.flatten()
                                                            .map { it ->
                                                               return tuple(infer_tumor_subtype(it),
                                                                            infer_genomic_region(it),
                                                                            'oncodrivefml', it)}
                                                            .combine(oncodrivefml_mutations, by: [0])
                                                            .combine(oncodrivefml_config)).csv
    emit:
    results = oncodrivefml_results
}