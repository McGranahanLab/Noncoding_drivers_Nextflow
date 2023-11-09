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
include { CALCULATE_MUTATION_RATES } from './modules/mutation_rates.nf'
include { RUN_DIGDRIVER } from './subworkflows/call_driver_genes.nf'
include { RUN_DNDSCV } from './subworkflows/call_driver_genes.nf'
include { RUN_MUTPANNING } from './subworkflows/call_driver_genes.nf'
include { RUN_NBR } from './subworkflows/call_driver_genes.nf'
include { RUN_ONCODRIVEFML } from './subworkflows/call_driver_genes.nf'

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

// same function is also defined in call_driver_genes.nf
def infer_tumor_subtype (filePath) {
    result = filePath.name.toString().tokenize('-').get(1)
    return result
}

// same function is also defined in call_driver_genes.nf
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
    path digdriverModels_inventory_path

    output:
    stdout emit: inventories_pass

    script:
    """
    1_check_inventories.R --inventory_patients ${patients_inventory_path} \
                          --inventory_analysis ${analysis_inventory_path} \
                          --target_genome_version ${params.target_genome_version} \
                          --inventory_blacklisted ${blacklist_inventory_path} \
                          --inventory_digdriver ${digdriverModels_inventory_path}
    """
}

process check_fasta_is_uscs {
    input:
    path fasta_file

    output:
    stdout emit: fasta_pass

    script:
    """
    grep '>' ${fasta_file} > chrs.txt
    nChr=`wc -l chrs.txt | sed 's/ chrs.txt//g'`
    nChrUCSC=`grep '>chr' chrs.txt | wc -l`
    if [[ \$nChr -eq \$nChrUCSC ]]
    then
        echo 'FASTA is UCSC'
    else
        echo ${fasta_file} is not in UCSC format
        exit 1
    fi
    """
}

/* ----------------------------------------------------------------------------
* Workflows
*----------------------------------------------------------------------------*/
workflow CALL_DE_NOVO_CANCER_DRIVERS {
    /*
        Step 1a: check that inventories have all the needed columns and values
                 are acceptable
    */
    // create channels to all inventories
    patients_inv = Channel.fromPath(params.patients_inventory, 
                                    checkIfExists: true)
                          .ifEmpty { exit 1, "[ERROR]: patients inventory file not found" }
    analysis_inv = Channel.fromPath(params.analysis_inventory,
                                    checkIfExists: true)
                          .ifEmpty { exit 1, "[ERROR]: analysis inventory file not found" }
    blacklist_inv = channel_from_params_path(params.blacklist_inventory)
    digdriver_inv = channel_from_params_path(params.digdriver_models_inventory)
    // perform the checks
    inventories_pass = check_inventories(patients_inv, analysis_inv, 
                                         blacklist_inv, digdriver_inv)
    inventories_pass = inventories_pass.collect()

    /*
        Step 1b: check that given reference genome is in UCSC format
    */
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
    // perform the checks
    target_genome_pass = check_fasta_is_uscs(target_genome_fasta)



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
    mut_rates = CALCULATE_MUTATION_RATES (PREPARE_INPUT_GENOMIC_REGIONS_FILES.out.bed
                                          .combine(PREPARE_INPUT_MUTATION_FILES.out.maf, by: [0])
                                          .combine(target_genome_chr_len)
                                          .combine(gene_name_synonyms)
                                          .combine(varanno_conversion_table))
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
    nbr_results = RUN_NBR.out.results
                         .map { it ->
                                    return tuple(infer_tumor_subtype(it),
                                                 infer_genomic_region(it), 
                                                 'nbr', it)
                              }
    /* 
        Step 5f: run OncodriveFML
    */
    RUN_ONCODRIVEFML (PREPARE_INPUT_MUTATION_FILES.out.oncodrivefml,
                      PREPARE_INPUT_GENOMIC_REGIONS_FILES.out.oncodrivefml,
                      oncodrivefml_config)
    oncodrivefml_results = RUN_ONCODRIVEFML.out.results
                                           .map { it ->
                                                    return tuple(infer_tumor_subtype(it),
                                                                 infer_genomic_region(it), 
                                                                 'oncodrivefml', it)
                                                }

    /* PERFORM CHECK FOR THE */
    nbr_results.join(oncodrivefml_results, by: [0, 1], remainder: true).view()
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