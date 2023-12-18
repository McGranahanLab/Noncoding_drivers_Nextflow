#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
    Organization: UCL Cancer Institute
    Laboratory: Cancer Genome Evolution
    Authors: Maria Litovchenko
    Purpose: 
    Results:
    Notes:
*/

/* ----------------------------------------------------------------------------
* Include modules & workflows
*----------------------------------------------------------------------------*/
include { check_inventories } from './modules/check_input_files.nf'
include { check_fasta_is_uscs } from './modules/check_input_files.nf'
include { check_digdriver_files } from './modules/check_input_files.nf'
include { check_nbr_files } from './modules/check_input_files.nf'
include { check_oncodrivefml_files } from './modules/check_input_files.nf'
include { PREPARE_INPUT_GENOMIC_REGIONS_FILES } from './subworkflows/prepare_input_genomic_resions_files.nf'
include { PREPARE_INPUT_MUTATION_FILES } from './subworkflows/prepare_input_mutations_files.nf'
include { CALCULATE_MUTATION_RATES } from './modules/mutation_rates.nf'
include { RUN_DIGDRIVER } from './subworkflows/call_driver_genes.nf'
include { RUN_DNDSCV } from './subworkflows/call_driver_genes.nf'
include { RUN_MUTPANNING } from './subworkflows/call_driver_genes.nf'
include { RUN_NBR } from './subworkflows/call_driver_genes.nf'
include { RUN_ONCODRIVEFML } from './subworkflows/call_driver_genes.nf'
include { COMBINE_P_VALS_AND_ANNOTATE } from './modules/postprocessing.nf'
include { ASSIGN_TIER } from './modules/postprocessing.nf'
include { FILTER_TIERED_DRIVERS } from './modules/postprocessing.nf'
include { ANNOTATE_GENOMICRANGES_WITH_CN } from './modules/drivers_characterization.nf'
include { ANNOTATE_MUTATIONS_WITH_MULTIPLICITY } from './modules/drivers_characterization.nf'

/* ----------------------------------------------------------------------------
* Custom functions
*----------------------------------------------------------------------------*/
def channel_from_params_path (staticPath) {
    // create an empty dummy file
    def no_file = new File(".NO_FILE")
    no_file.createNewFile()

    if (staticPath == '') {
        result = Channel.fromPath(".NO_FILE")
    } else {
        result = Channel.fromPath(staticPath, checkIfExists: true)
                        .ifEmpty { exit 1, "[ERROR]: " + staticPath + "file not found"}
    } 
    return result
}

def create_result_file_tuple (inventoryRow, outputDir, genomeVersion) {
    mutMapPath = file(outputDir + '/results/mut_rates/mutMapToGR-' + 
                      inventoryRow.tumor_subtype + '--' + genomeVersion + 
                      '.csv', checkIfExists: false)
    mutRatePath = file(outputDir + '/results/mut_rates/meanMutRatePerGR-' + 
                       inventoryRow.tumor_subtype + '--' + genomeVersion + 
                       '.csv', checkIfExists: false)
    scannedGRpath = file(outputDir + '/inputs/inputGR-' + 
                         inventoryRow.tumor_subtype + '-' + 
                         genomeVersion + '.bed', checkIfExists: false)
    resultsPath = file(outputDir + '/results/' + inventoryRow.software + '/' +
                       inventoryRow.software + 'Results-' + 
                       inventoryRow.tumor_subtype + '-' + inventoryRow.gr_id +
                       '-' + genomeVersion + '.csv', checkIfExists: false)
    return tuple(inventoryRow.tumor_subtype, inventoryRow.gr_id, 
                 mutMapPath, mutRatePath, scannedGRpath, inventoryRow.software,
                 resultsPath)
}

/* ----------------------------------------------------------------------------
* Main workflows
*----------------------------------------------------------------------------*/
workflow CALL_DE_NOVO_CANCER_DRIVERS {
    /*
        Step 1a: check that inventories have all the needed columns and values
                 are acceptable
    */
    patients_inv     = Channel.fromPath(params.patients_inventory, 
                                        checkIfExists: true)
                              .ifEmpty { exit 1, "[ERROR]: patients inventory file not found" }
    analysis_inv     = Channel.fromPath(params.analysis_inventory,
                                        checkIfExists: true)
                              .ifEmpty { exit 1, "[ERROR]: analysis inventory file not found" }
    blacklist_inv    = channel_from_params_path(params.blacklist_inventory)
    digdriver_inv    = channel_from_params_path(params.digdriver_models_inventory)
    chasmplus_inv    = channel_from_params_path(params.chasmplus_annotators_inventory)
    // perform the checks
    inventories_pass = check_inventories(patients_inv, analysis_inv, 
                                         blacklist_inv, digdriver_inv, 
                                         chasmplus_inv)
    inventories_pass = inventories_pass.collect()

    /*
        Step 1b: check that given reference genome is in UCSC format
    */
    target_genome_fasta   = Channel.fromPath(params.target_genome_path,
                                             checkIfExists: true)
                                   .ifEmpty { exit 1, 
                                              "[ERROR]: target genome fasta file not found" }
    target_genome_chr_len = Channel.fromPath(params.target_genome_chr_len,
                                             checkIfExists: true)
                                   .ifEmpty { exit 1, 
                                            "[ERROR]: chromosomal lengths of target genome file not found" }
    chain                 = channel_from_params_path(params.chain)
    // perform the check for reference genome format (only UCSC)
    target_genome_pass    = check_fasta_is_uscs(target_genome_fasta)

    /*
        Step 1c: check that essential files for various softwares are present
    */
    software               = analysis_inv.splitCsv(header: true).map{row -> row.software}
                                         .unique()
    digdriver_elements     = channel_from_params_path(params.digdriver_elements)
    nbr_neutral_bins       = channel_from_params_path(params.nbr_regions_neutralbins_file)
    nbr_neutral_trinucfreq = channel_from_params_path(params.nbr_trinucfreq_neutralbins_file)
    nbr_driver_regs        = channel_from_params_path(params.nbr_driver_regs_file)
    oncodrivefml_config    = channel_from_params_path(params.oncodrivefml_config)
    essential_files        = software.combine(digdriver_elements)
                                     .combine(nbr_neutral_bins)
                                     .combine(nbr_neutral_trinucfreq)
                                     .combine(nbr_driver_regs)
                                     .combine(oncodrivefml_config)
    // perform the checks
    digdriver_files_pass = check_digdriver_files(essential_files)
    nbr_files_pass = check_nbr_files(essential_files)
    oncodrivefml_files_pass = check_oncodrivefml_files(essential_files)

    // combine all the checks together. Will work if and only if all checks are
    // passed.
    inventories_pass = inventories_pass.combine(target_genome_pass)
                                       .combine(digdriver_files_pass)
                                       .combine(nbr_files_pass)
                                       .combine(oncodrivefml_files_pass)
     
    /* 
        Step 2: create input mutations files
    */
    PREPARE_INPUT_MUTATION_FILES (analysis_inv, patients_inv, blacklist_inv,
                                  target_genome_fasta, target_genome_chr_len,
                                  chain, inventories_pass)
    /* 
        Step 3: create input genomic regions files
    */
    PREPARE_INPUT_GENOMIC_REGIONS_FILES (analysis_inv, blacklist_inv,
                                         target_genome_fasta,
                                         target_genome_chr_len, chain,
                                         inventories_pass)
    /* 
        Step 4: calculate mutation rates
    */
    gene_name_synonyms = channel_from_params_path(params.gene_name_synonyms)
    varanno_conversion_table = channel_from_params_path(params.varanno_conversion_table)
    coding_gr_id = analysis_inv.splitCsv(header: true)
                               .map { row -> 
                                        tuple(row.tumor_subtype, row.gr_code,
                                              row.gr_id)
                               }
                               .filter{ it[1] == 'CDS' }
                               .map { it ->
                                         return(tuple(it[0], it[2]))
                               }
                               .unique()
                               .groupTuple(by: [0, 1], remainder: true)
    CALCULATE_MUTATION_RATES (PREPARE_INPUT_GENOMIC_REGIONS_FILES.out.bed
                                  .combine(PREPARE_INPUT_MUTATION_FILES.out.maf, by: [0])
                                  .combine(target_genome_chr_len)
                                  .combine(gene_name_synonyms)
                                  .combine(varanno_conversion_table)
                                  .combine(coding_gr_id, by: [0]))

    /* 
        Step 5a: run CHASMplus
    
    RUN_CHASMplus (analysis_inv, PREPARE_INPUT_MUTATION_FILES.out.chasmplus, 
                   chasmplus_inv)
    */

    /* 
        Step 5b: run DIGdriver
    */
    RUN_DIGDRIVER (digdriver_inv, 
                   PREPARE_INPUT_MUTATION_FILES.out.digdriver,
                   PREPARE_INPUT_GENOMIC_REGIONS_FILES.out.digdriver,
                   PREPARE_INPUT_GENOMIC_REGIONS_FILES.out.dndscv_digdriver_rda,
                   digdriver_elements, target_genome_fasta)
    /* 
        Step 5c: run dNdScv
    */
    RUN_DNDSCV (analysis_inv,
                PREPARE_INPUT_MUTATION_FILES.out.dndscv,
                PREPARE_INPUT_GENOMIC_REGIONS_FILES.out.dndscv_digdriver_rda)
    /* 
        Step 5d: run MutPanning
    */
    RUN_MUTPANNING (analysis_inv, 
                    PREPARE_INPUT_MUTATION_FILES.out.mutpanning)
    /* 
        Step 5e: run NBR
    */
    RUN_NBR (PREPARE_INPUT_MUTATION_FILES.out.nbr,
             PREPARE_INPUT_GENOMIC_REGIONS_FILES.out.nbr,
             target_genome_fasta, nbr_neutral_bins,
             nbr_neutral_trinucfreq, nbr_driver_regs)
    /* 
        Step 5f: run OncodriveFML
    */
    RUN_ONCODRIVEFML (PREPARE_INPUT_MUTATION_FILES.out.oncodrivefml,
                      PREPARE_INPUT_GENOMIC_REGIONS_FILES.out.oncodrivefml,
                      oncodrivefml_config)
 }

workflow POSTPROCESSING {
    /*
        Step 1: check inventory which defines tiers, inventories with 
                expression data and check rawP cap
    */

    /*
        Step 2: load analyis inventory & check that all software run & produced
                 result files
    */
    analysis_inv = Channel.fromPath(params.analysis_inventory,
                                    checkIfExists: true)
                          .ifEmpty { exit 1, "[ERROR]: analysis inventory file not found" }

    // DO NOT FORGET TO TURN checkIfExists TO true IN FUNCTION create_path_to_result_file

    analysis_inv  = analysis_inv.splitCsv(header: true)
                                .map { row -> 
                                           return create_result_file_tuple(row,
                                                                           params.outdir,
                                                                           params.target_genome_version)
                                     }
                                 .unique()
                                 .groupTuple(by: [0, 1, 2, 3, 4], remainder: true)
    // load tier definition table
    tier_inventory = Channel.fromPath(params.tier_inventory, 
                                      checkIfExists: true)
                            .ifEmpty { exit 1, "[ERROR]: tier inventory file not found" }
    // load patients inventory
    patients_inv   = Channel.fromPath(params.patients_inventory,
                                      checkIfExists: true)
                            .ifEmpty { exit 1, "[ERROR]: patients inventory file not found" }
    chain          = channel_from_params_path(params.chain)

    /*
        Step 3: load data used for driver annotation (i.e. known cancer gene
                status, olfactory genes, expression, etc)
    */
    gene_name_synonyms = channel_from_params_path(params.gene_name_synonyms)
    known_cancer_genes = channel_from_params_path(params.known_cancer_genes)
    olfactory_genes    = channel_from_params_path(params.olfactory_genes)
    gtex_inventory     = channel_from_params_path(params.gtex_inventory) 
    gtex_expression    = channel_from_params_path(params.gtex_expression)
    tcga_inventory     = channel_from_params_path(params.tcga_inventory) 
    tcga_expression    = channel_from_params_path(params.tcga_expression)

    combined_pvals = COMBINE_P_VALS_AND_ANNOTATE (analysis_inv.combine(Channel.from(params.rawP_cap))
                                                              .combine(gene_name_synonyms)
                                                              .combine(known_cancer_genes)
                                                              .combine(olfactory_genes)
                                                              .combine(gtex_inventory)
                                                              .combine(gtex_expression)
                                                              .combine(tcga_inventory)
                                                              .combine(tcga_expression)).csv
    tiered_pvals   = ASSIGN_TIER (combined_pvals.groupTuple(by: [0], remainder: true)
                                                .combine(tier_inventory)
                                                .combine(Channel.of(params.combine_p_method))).csv
 
    drivers        = FILTER_TIERED_DRIVERS (tiered_pvals.combine(Channel.fromPath(params.analysis_inventory,
                                                                                  checkIfExists: true))).csv
    // only tumor subtypes which do have drivers detected will be processed 
    // further
    drivers        = drivers.filter { it[2]  == 'yes' }
                            .map { it -> return(tuple(it[0], it[1])) }

    /* 
        Step 4a: annotate drivers with copy number (if provided)
    */
    // INSERT HERE A BREAK TO STOP PROCESS FROM HAPPENING IN CASE NO CN WAS PROVIDED
    // DO IT VIA MOVING CN AND MUTATION MULTIPLICITY FILES INTO SEPARATE INVENTORY
    // THIS WILL ALSO MEAN THAT THE 1st PIPELINE WON'T BE TRIGGERED UPON CHANGE
    // IN INVENTORY
    driver_gr_annotated_with_cn = ANNOTATE_GENOMICRANGES_WITH_CN(analysis_inv.map { it -> return(tuple(it[0], it[4])) }
                                                                             .unique()
                                                                             .combine(drivers, by: [0])
                                                                             .combine(patients_inv)
                                                                             .combine(chain))
    mutsinDrives_annotated_with_mults = ANNOTATE_MUTATIONS_WITH_MULTIPLICITY(analysis_inv.map { it -> return(tuple(it[0], it[2]))}
                                                                                         .unique()
                                                                                         .combine(drivers, by: [0])
                                                                                         .combine(patients_inv))
}


// inform about completition
workflow.onComplete {
    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}