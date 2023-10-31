/* 
    A workflow to create input mutation files for CHASM+, DIGdriver, dNdScv, 
    MutPanning, NBR and OncodriveFML.
*/

include { FILTER_INPUT_MUTATIONS } from '../modules/prepare_mutations.nf'
include { WRITE_MUTATIONS_FOR_CHASMPLUS } from '../modules/prepare_mutations.nf'
include { WRITE_MUTATIONS_FOR_DIGDRIVER } from '../modules/prepare_mutations.nf'
include { WRITE_MUTATIONS_FOR_DNDSCV } from '../modules/prepare_mutations.nf'
include { WRITE_MUTATIONS_FOR_MUTPANNING } from '../modules/prepare_mutations.nf'
include { WRITE_MUTATIONS_FOR_NBR } from '../modules/prepare_mutations.nf'
include { WRITE_MUTATIONS_FOR_ONCODRIVEFML } from '../modules/prepare_mutations.nf'

workflow PREPARE_INPUT_MUTATION_FILES {
    take:
    analysis_inv
    patients_inv
    blacklist_inv
    target_genome_fasta
    target_genome_chr_len
    chain
    inventories_pass

    main:
	// filter mutations
    tumor_subtypes = analysis_inv.splitCsv(header: true)
                                 .map{row -> row.tumor_subtype}
                                 .unique()
    filtered_mutations = FILTER_INPUT_MUTATIONS(inventories_pass.combine(tumor_subtypes)
                                                                .combine(patients_inv)
                                                                .combine(analysis_inv)
                                                                .combine(blacklist_inv)
                                                                .combine(target_genome_fasta)
                                                                .combine(target_genome_chr_len)
                                                                .combine(chain))

    // output them in a format suitable to run driver calling tools on
    tumor_subtypes = analysis_inv.splitCsv(header: true)
                                 .map{row -> tuple(row.tumor_subtype, row.software)}
                                 .unique()
                                 .combine(filtered_mutations.mutations, by: [0])
    chasmplus_mutations = WRITE_MUTATIONS_FOR_CHASMPLUS(tumor_subtypes)
    digdriver_mutations = WRITE_MUTATIONS_FOR_DIGDRIVER(tumor_subtypes)
    dndscv_mutations = WRITE_MUTATIONS_FOR_DNDSCV(tumor_subtypes)
    mutpanning_mutations = WRITE_MUTATIONS_FOR_MUTPANNING(tumor_subtypes)
    nbr_mutations = WRITE_MUTATIONS_FOR_NBR(tumor_subtypes)
    oncodrivefml_mutations = WRITE_MUTATIONS_FOR_ONCODRIVEFML(tumor_subtypes)

    emit:
    chasmplus = chasmplus_mutations
    digdriver = digdriver_mutations
    dndscv = dndscv_mutations
    mutpanning = mutpanning_mutations
    nbr = nbr_mutations
    oncodrivefml = oncodrivefml_mutations
}