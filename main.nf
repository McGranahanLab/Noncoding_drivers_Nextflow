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
    // filter mutations
    tumor_subtypes = analysis_inv.splitCsv(header: true)
                                 .map{row -> row.tumor_subtype}
                                 .unique()
    filtered_mutations = filter_input_mutations(inventories_pass.combine(tumor_subtypes)
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
    chasmplus_mutations = write_mutations_for_chasmplus(tumor_subtypes)
    digdriver_mutations = write_mutations_for_digdriver(tumor_subtypes)
    dndscv_mutations = write_mutations_for_dndscv(tumor_subtypes)
    mutpanning_mutations = write_mutations_for_mutpanning(tumor_subtypes)
    nbr_mutations = write_mutations_for_nbr(tumor_subtypes)
    oncodrivefml_mutations = write_mutations_for_oncodrivefml(tumor_subtypes)

    /* 
        Step 3: create input genomic regions files
    */
    // pre-process regions
    filtered_regions = filter_genomic_regions(inventories_pass.combine(analysis_inv)
                                                              .combine(blacklist_inv)
                                                              .combine(target_genome_fasta)
                                                              .combine(target_genome_chr_len)
                                                              .combine(chain)).bed.flatten()
    filtered_regions.map { it ->
                           def tumor_subtype = it.name.toString().tokenize('-').get(1)
                           return tuple(tumor_subtype, it)
                         }.set{filtered_regions}
    /* in case set of the same regions is analysed for a lot of tumor subtypes
       it is faster to create region file for one tumor subtype and then copy 
       it for the other tumor subtypes (for one software of course). This is 
       what is done here */
    tumor_subtypes_and_gr = analysis_inv.splitCsv(header: true)
                                        .map{row -> tuple(row.tumor_subtype, row.software, row.gr_id)}
                                        .unique()
                                        .groupTuple(by: [0, 1])
                                        .combine(filtered_regions, by: 0)
    digdriver_regions = write_regions_for_digdriver(tumor_subtypes_and_gr)
    nbr_regions = write_regions_for_nbr(tumor_subtypes_and_gr)
    oncodrivefml_regions = write_regions_for_oncodrivefml(tumor_subtypes_and_gr)
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
    dndscv_digdriver_rda = create_rda_for_dndscv_digdriver(gtf_for_rda.combine(target_genome_fasta)
                                                                      .combine(chain))
    /* 
        Step 4b: run DIGdriver
    */
    digdriver_inv = digdriver_inv.splitCsv(header: true)
                                 .map { it ->
                                     return(tuple(it.tumor_subtype, it.model_file))
                                 }
    digdriver_results = digdriver(digdriver_regions.flatten()
                                                   .map { it ->
                                                            return tuple(infer_tumor_subtype(it),
                                                                         infer_genomic_region(it),
                                                                         'digdriver', it)
                                                    }
                                                    .combine(dndscv_digdriver_rda.rda, by: [0])
                                                    .combine(digdriver_mutations, by: [0])
                                                    .combine(digdriver_inv, by: [0])
                                                    .combine(digdriver_elements)
                                                    .combine(target_genome_fasta))
    /* 
        Step 4c: run dNdScv
    */
    dndscv_results = dndscv(analysis_inv.splitCsv(header: true)
                                        .map{row -> tuple(row.tumor_subtype, row.software, row.gr_id)}
                                        .unique()
                                        .map { it ->
                                            if (it[1] == 'dndscv') {
                                                return tuple(it[0], it[2], 'dndscv')
                                            }
                                        }
                                        .combine(dndscv_digdriver_rda.rda, by: [0])
                                        .combine(dndscv_mutations, by: [0]))
    /* 
        Step 4d: run MutPanning
    */
    mutpanning_results = mutpanning(analysis_inv.splitCsv(header: true)
                                                .map{row -> tuple(row.tumor_subtype, row.software, row.gr_id)}
                                                .unique()
                                                .map { it ->
                                                    if (it[1] == 'mutpanning') {
                                                         return tuple(it[0], it[2], 'mutpanning')
                                                    }
                                                }
                                                .combine(mutpanning_mutations, by: [0]))
    /* 
        Step 4e: run NBR
    */
    nbr_results = nbr(nbr_regions.flatten()
                                 .map { it ->
                                     return tuple(infer_tumor_subtype(it),
                                                  infer_genomic_region(it), 'nbr',
                                                  it)
                                 }
                                 .combine(nbr_mutations, by: [0])
                                 .combine(target_genome_fasta)
                                 .combine(nbr_neutral_bins)
                                 .combine(nbr_neutral_trinucfreq)
                                 .combine(nbr_driver_regs))
    /* 
        Step 4f: run OncodriveFML
    */
    oncodrivefml_results = oncodrivefml(oncodrivefml_regions.flatten()
                                                            .map { it ->
                                                               return tuple(infer_tumor_subtype(it),
                                                                            infer_genomic_region(it),
                                                                            'oncodrivefml', it)}
                                                            .combine(oncodrivefml_mutations, by: [0])
                                                            .combine(oncodrivefml_config))
}

// inform about completition
workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Script finished SUCCESSFULLY after $workflow.duration"
    } else {
      log.info "[$workflow.complete] >> Script finished with ERRORS after $workflow.duration"
    }
}