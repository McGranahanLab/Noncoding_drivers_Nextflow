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
* Include workflows
*----------------------------------------------------------------------------*/
include { PREPARE_INPUT_GENOMIC_REGIONS_FILES } from './subworkflows/prepare_input_genomic_resions_files.nf'
include { PREPARE_INPUT_MUTATION_FILES } from './subworkflows/prepare_input_mutations_files.nf'

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

def infer_tumor_subtype (filePath) {
    result = filePath.name.toString().tokenize('-').get(1)
    return result
}

def infer_genomic_region(filePath) {
    result = filePath.name.toString().tokenize('-').get(3)
    return result
}

/* ----------------------------------------------------------------------------
* Processes
*----------------------------------------------------------------------------*/
process check_inventories {
    input:
    path patients_inventory_path
    path analysis_inventory_path
    path blacklist_inventory_path

    output:
    stdout emit: inventories_pass

    script:
    """
    1_check_inventories.R --inventory_patients ${patients_inventory_path} \
                          --inventory_analysis ${analysis_inventory_path} \
                          --target_genome_version ${params.target_genome_version} \
                          --inventory_blacklisted ${blacklist_inventory_path}
    """
}

/* ----------------------------------------------------------------------------
* Workflows
*----------------------------------------------------------------------------*/
workflow {
    /*
        Step 1: convert inventories & reference genome files to channels
    */
    // create channels to all inventories
    patients_inv = Channel.fromPath(params.patients_inventory, 
                                    checkIfExists: true)
                          .ifEmpty { exit 1, "[ERROR]: patients inventory file not found" }
    analysis_inv = Channel.fromPath(params.analysis_inventory,
                                    checkIfExists: true)
                          .ifEmpty { exit 1, "[ERROR]: analysis inventory file not found" }
    blacklist_inv = channel_from_params_path(params.blacklist_inventory)

    // create channels to target genome verion and to chain file for liftover
    target_genome_fasta = Channel.fromPath(params.target_genome_path,
                                           checkIfExists: true)
                                 .ifEmpty { exit 1, 
                                            "[ERROR]: target genome fasta file not found" }
    target_genome_chr_len = Channel.fromPath(params.target_genome_chr_len,
                                           checkIfExists: true)
                                   .ifEmpty { exit 1, 
                                            "[ERROR]: chromosomal lengths of target genome file not found" }
    chain = channel_from_params_path(params.chain)

    // create channels to gene name synonyms (needed for mutation rate calc)
    gene_name_synonyms = channel_from_params_path(params.gene_name_synonyms)
    varanno_conversion_table = channel_from_params_path(params.varanno_conversion_table)

    /*
        Step 2: convert additional files for the software (loaded only if 
        corresponding software is requested) to channels
    */
    software = analysis_inv.splitCsv(header: true).map{row -> row.software}
                           .unique()
    /*
    This code makes nextflow never finish
    check that all needed files are present
    software.map { it ->
        if (it == 'digdriver') {
            Channel.fromPath(params.digdriver_models_inventory, checkIfExists: true)
                           .ifEmpty { exit 1, "[ERROR]: digdriver models inventory file not found" }
            Channel.fromPath(params.digdriver_elements, checkIfExists: true)
                           .ifEmpty { exit 1, "[ERROR]: digdriver elements file not found" }
        }
        if (it == 'nbr') {
            nbr_neutral_bins = Channel.fromPath(params.nbr_regions_neutralbins_file,
                                                checkIfExists: true)
                                      .ifEmpty { exit 1, "[ERROR]: nbr_regions_neutralbins_file not found" }
            nbr_neutral_trinucfreq = Channel.fromPath(params.nbr_trinucfreq_neutralbins_file, 
                                                      checkIfExists: true)
                                            .ifEmpty { exit 1, "[ERROR]: nbr_trinucfreq_neutralbins_file not found" }
            nbr_driver_regs = Channel.fromPath(params.nbr_driver_regs_file, 
                                               checkIfExists: true)
                                     .ifEmpty { exit 1, "[ERROR]: nbr_driver_regs_file not found" }
        }
        if (it == 'oncodrivefml') {
            oncodrivefml_config = Channel.fromPath(params.oncodrivefml_config,
                                                   checkIfExists: true)
                                         .ifEmpty { exit 1, "[ERROR]: oncodrivefml_config not found" }
        }
    }*/
     
    // now all the files needed for individual software could be loaded via
    // channel_from_params_path function. If file does not exist, it will just 
    // create an empty file & continue 
    digdriver_inv = channel_from_params_path(params.digdriver_models_inventory)
    digdriver_elements = channel_from_params_path(params.digdriver_elements)
    nbr_neutral_bins = channel_from_params_path(params.nbr_regions_neutralbins_file)
    nbr_neutral_trinucfreq = channel_from_params_path(params.nbr_trinucfreq_neutralbins_file)
    nbr_driver_regs = channel_from_params_path(params.nbr_driver_regs_file)
    oncodrivefml_config = channel_from_params_path(params.oncodrivefml_config)

    /*
        Step 1: check that inventories have all the needed columns and values
                are acceptable 
    */
    inventories_pass = check_inventories(patients_inv, analysis_inv, blacklist_inv)
    inventories_pass = inventories_pass.collect()

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