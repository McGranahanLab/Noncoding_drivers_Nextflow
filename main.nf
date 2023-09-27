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

process create_input_mutation_files {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), val(inventory_check_res), 
          path(patients_inventory_path), path(analysis_inventory_path),
          path(blacklist_inventory_path), path(target_genome_fasta),
          path(chain)

    output:
    path('inputs/*inputMutations*')

    script:
    """
    software_flatten=`echo ${software} | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/,//g'`
    2_create_input_mutation_files.R --inventory_patients ${patients_inventory_path} \
                                    --inventory_analysis ${analysis_inventory_path} \
                                    --inventory_blacklisted ${blacklist_inventory_path} \
                                    --cancer_subtype ${tumor_subtype} \
                                    --software \$software_flatten \
                                    --min_depth ${params.min_depth} \
                                    --min_tumor_vac ${params.min_tumor_vac} \
                                    --min_tumor_vaf ${params.min_tumor_vaf} \
                                    --max_germline_vaf ${params.max_germline_vaf} \
                                    --max_germline_vac ${params.max_germline_vac} \
                                    --max_n_vars ${params.max_n_vars} \
                                    --target_genome_path ${target_genome_fasta} \
                                    --target_genome_version ${params.target_genome_version} \
                                    --chain ${chain} \
                                    --output 'inputs/' \
                                    --cores ${task.cpus}
    
    """
}

process create_input_genomic_regions_files {
    input:
    tuple val(inventory_check_res), path(analysis_inventory_path),
          path(blacklist_inventory_path), path(target_genome_fasta),
          path(chain)

    output:
    path('inputs/*inputGR*')

    script:
    """
    3_create_input_genomic_regions_files.R --inventory_analysis ${analysis_inventory_path} \
                                           --inventory_blacklisted ${blacklist_inventory_path} \
                                           --target_genome_path ${target_genome_fasta} \
                                           --target_genome_version ${params.target_genome_version} \
                                           --chain ${chain} \
                                           --ignore_strand ${params.ignore_strand} \
                                           --min_reg_len ${params.min_reg_len} \
                                           --output 'inputs/' \
                                           --cores ${task.cpus}
    """
}

process calculate_mutation_rates {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(bed_path), path(maf_path), 
          path(target_genome_chr_len), path(gene_name_synonyms),
          path(varanno_conversion_table)

    output:
    tuple val(tumor_subtype), path('*meanMutRatePerGR.csv'), 
          path('*mutMapToGR.csv'), path('*varCatEnrich.csv')

    script:
    """
    4_calculate_mutation_rates.R --tumor_subtype ${tumor_subtype} \
                                 --variants ${maf_path} \
                                 --genomic_regions ${bed_path} \
                                 --gene_name_synonyms ${gene_name_synonyms} \
                                 --bin_len ${params.bin_len} \
                                 --target_genome_version ${params.target_genome_version} \
                                 --target_genome_chr_len ${target_genome_chr_len} \
                                 --calc_synonymous ${params.calc_synonymous} \
                                 --cdsAcceptedClass ${params.cdsAcceptedClass} \
                                 --ncAcceptedClass ${params.ncAcceptedClass} \
                                 --varanno_conversion_table ${varanno_conversion_table} \
                                 --annotation_failed_code ${params.annotation_failed_code} \
                                 --output 'mutRate-'${tumor_subtype}'-'
    """
}

process mutpanning {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), val(gr_id), path(mutations),
          path(mutpan_inv)

    output:
    path "${software}-results-${tumor_subtype}-${gr_id}-${params.target_genome_version}*"

    script:
    """
    OUT_FILE=${software}"-results-"${tumor_subtype}"-"${gr_id}'-'${params.target_genome_version}'.csv'

    mpPath=/bin/MutPanningV2/commons-math3-3.6.1.jar:/bin/MutPanningV2/jdistlib-0.4.5-bin.jar:/bin/MutPanningV2
    java -Xmx${params.mutpanning_java_memory} -classpath \$mpPath MutPanning "." \
        $mutations $mutpan_inv "/bin/MutPanningV2/Hg19/"

    cpFrom=`find SignificanceRaw | grep -v Uniform | grep '.txt\$'`
    if [[ -f "\$cpFrom" ]]; then
        cp \$cpFrom \$OUT_FILE
    else
        echo 'MutPanning run on '${tumor_subtype}" "${gr_id}' did not yield results.'
        echo 'It can happen if your cohort is too small. Will create an empty file.'
        echo -e "Name\tTargetSize\tTargetSizeSyn\tCount\tCountSyn\tSignificanceSyn\tFDRSyn\tSignificance\tFDR" > \$OUT_FILE
    fi
    """
}

process nbr {
    tag "$tumor_subtype-$gr_id"

    input:
    tuple val(tumor_subtype), val(software), val(gr_id), path(mutations),
          path(regions), path(target_genome_fasta), path(nbr_neutral_bins),
          path(nbr_neutral_trinucfreq), path(nbr_driver_regs)
    
    output:
    path "${software}-results-${tumor_subtype}-${gr_id}-${params.target_genome_version}*"

    script:
    """
    OUT_FILE=${software}"-results-"${tumor_subtype}"-"${gr_id}'-'${params.target_genome_version}'.csv'

    NBR_create_intervals.R --region_bed ${regions} --genomeFile ${target_genome_fasta}
    regionsFile=`echo ${regions}| sed 's/.csv\$/.csv.regions/g'`
    triNuclregionsFile=`echo ${regions} | sed 's/.csv\$/.csv.txt/g'`

    NBR.R --mutations_file ${mutations} --target_regions_path \$regionsFile \
          --target_regions_trinucfreqs_path \$triNuclregionsFile \
          --genomeFile ${target_genome_fasta} \
          --max_num_muts_perRegion_perSample 2 --unique_indelsites_FLAG 1 \
          --gr_drivers ${nbr_driver_regs} \
          --regions_neutralbins_file ${nbr_neutral_bins} \
          --trinucfreq_neutralbins_file ${nbr_neutral_trinucfreq} \
          --out_prefix \$OUT_FILE
    mv \$OUT_FILE"-Selection_output.txt" \$OUT_FILE 
    OUT_FILE_BASE=`echo \$OUT_FILE | sed 's/.csv//g'`
    # mv \$OUT_FILE"-global_mle_subs.Rds" \$OUT_FILE_BASE"-global_mle_subs.Rds" 
    # mv \$OUT_FILE"-globalRates.csv" \$OUT_FILE_BASE"-globalRates.csv"
    """
}

process oncodrivefml {
    tag "$tumor_subtype"-"$gr_id"

    input:
    tuple val(tumor_subtype), val(software), val(gr_id), path(mutations), 
          path(regions), path(oncodrivefml_config)

    script:
    """
    OUT_FILE=${software}"-results-"${tumor_subtype}"-"${gr_id}'-'${params.target_genome_version}'.csv'

    oncodrivefml -i ${mutations} -e ${regions} --sequencing wgs \
                 --configuration ${oncodrivefml_config} \
                 --output '.'  --cores ${task.cpus}
    mv "oncodrivefml-inputMutations-"$tumor_subtype"-"${params.target_genome_version}"-oncodrivefml.tsv" \$OUT_FILE
    """
}

/* ----------------------------------------------------------------------------
* Workflows
*----------------------------------------------------------------------------*/
workflow {
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
    target_genome_chr_len = channel_from_params_path(params.target_genome_chr_len)
    chain = channel_from_params_path(params.chain)

    // create channels to gene name synonyms (needed for mutation rate calc)
    gene_name_synonyms = channel_from_params_path(params.gene_name_synonyms)
    varanno_conversion_table = channel_from_params_path(params.varanno_conversion_table)

    // NBR specific files
    nbr_neutral_bins = Channel.fromPath(params.nbr_regions_neutralbins_file, 
                                        checkIfExists: true)
                              .ifEmpty { exit 1, "[ERROR]: nbr_regions_neutralbins_file not found" }
    nbr_neutral_trinucfreq = Channel.fromPath(params.nbr_trinucfreq_neutralbins_file, 
                                              checkIfExists: true)
                                    .ifEmpty { exit 1, "[ERROR]: nbr_trinucfreq_neutralbins_file not found" }
    nbr_driver_regs = Channel.fromPath(params.nbr_driver_regs_file, 
                                               checkIfExists: true)
                                     .ifEmpty { exit 1, "[ERROR]: nbr_driver_regs_file not found" }
    //Oncodrivefml specific files
    oncodrivefml_config = Channel.fromPath(params.oncodrivefml_config,
                                           checkIfExists: true)
                                 .ifEmpty { exit 1, "[ERROR]: oncodrivefml_config not found" }

    /*
        Step 1: check that inventories have all the needed columns and values
                are acceptable 
    */
    inventories_pass = check_inventories(patients_inv, analysis_inv, blacklist_inv)
    inventories_pass = inventories_pass.collect()

    /* 
        Step 2: create input mutations files for all requested software, 
                parallelise by cancer subtype. This will not be run if 
                inventories do not pass the check.
    */
    tumor_subtypes = analysis_inv.splitCsv(header: true)
                                 .map{row -> tuple(row.tumor_subtype, row.software)}
                                 .unique().groupTuple()
    input_mutations =  create_input_mutation_files(tumor_subtypes
                                                   .combine(inventories_pass)
                                                   .combine(patients_inv)
                                                   .combine(analysis_inv)
                                                   .combine(blacklist_inv)
                                                   .combine(target_genome_fasta)
                                                   .combine(chain))
                                                   .flatten()

     // split on csv which will go to be processed by driver calling software and bed files
    // which will go to estimation of mutation rates
    input_mutations.branch {
        csv: it.name.toString().endsWith('csv')
            return it
        maf: true
            return it
    }.set{ input_mutations_split }


    /*
        Step 3: create input mutations genomic region files for all requested
                software, parallelize by cancer subtype. This will not be run if
                inventories do not pass the check.
    */
    input_genomic_regions = create_input_genomic_regions_files(inventories_pass
                                                               .combine(analysis_inv)
                                                               .combine(blacklist_inv)
                                                               .combine(target_genome_fasta)
                                                               .combine(chain))
                                                               .flatten()
    // split on csv which will go to be processed by driver calling software and bed files
    // which will go to estimation of mutation rates
    input_genomic_regions.branch {
        csv: it.toString().endsWith('csv')
             return it
        bed: true
             return it
    }.set{ input_genomic_regions_split }

    /*
        Step 4: calculate local mutation rates (used for filtering of drivers)
    input_genomic_regions_split.bed.map { file ->
        def splitted_name = file.name.toString().tokenize('-')
        def tumor_subtype = splitted_name.get(1)
        return tuple(tumor_subtype, file)
    }.set{regions_bed_files}
    input_mutations_split.maf.map { file ->
        def splitted_name = file.name.toString().tokenize('-')
        def tumor_subtype = splitted_name.get(1)
        return tuple(tumor_subtype, file)
    }.set{muts_maf_files}

    mut_rates = calculate_mutation_rates(regions_bed_files.join(muts_maf_files, by: [0])
                                                          .combine(target_genome_chr_len)
                                                          .combine(gene_name_synonyms)
                                                          .combine(varanno_conversion_table)) */

    /*
        Step 5: prepare channel input_to_soft which will contain pairs of 
                genomic regions and mutations for all the software
    */
    input_genomic_regions_split.csv.map { file ->
        def splitted_name = file.name.toString().tokenize('-')
        def tumor_subtype = splitted_name.get(2)
        def gr_id = splitted_name.get(3)
        def software = splitted_name.get(0)
        return tuple(tumor_subtype, software, gr_id, file)
     }.set{ input_genomic_regions_to_soft }

    input_mutations_split.csv.map { file ->
        def splitted_name = file.name.toString().tokenize('-')
        def tumor_subtype = splitted_name.get(2)
        def software = splitted_name.get(0)
        return tuple(tumor_subtype, software, file)
    }.set{ input_mutations_split_to_soft }

    input_mutations_split_to_soft = analysis_inv.splitCsv(header: true)
                                                .map{row -> tuple(row.tumor_subtype, row.software, 
                                                                  row.gr_id)}
                                                .unique()
                                                .combine(input_mutations_split_to_soft, 
                                                         by: [0, 1])
    input_to_soft = input_mutations_split_to_soft.combine(input_genomic_regions_to_soft,
                                                          by: [0, 1, 2])

    input_to_soft.branch {
        digdriver: it[1].toString() == 'digdriver'
            return it
        dndscv: it[1].toString() == 'dndscv'
            return it
        mutpanning: it[1].toString() == 'mutpanning'
            return it
        nbr: it[1].toString() == 'nbr'
            return it
        oncodrivefml: it[1].toString() == 'oncodrivefml'
            return it
    }.set{ input_to_soft_split }

    /*
        Step 6: run mutpanning
    
    mutpanning_results = mutpanning(input_to_soft_split.mutpanning
                                                       .map { it ->
                                                                def mutpan_inv = it[3].toString().replace("inputMutations",
                                                                                                          "patientsInv")
                                                                return tuple(it[0], it[1], it[2], it[3], mutpan_inv)
                                                            })*/

    /*
        Step 7: run NBR
    
    nbr_results = nbr(input_to_soft_split.nbr.combine(target_genome_fasta)
                                             .combine(nbr_neutral_bins)
                                             .combine(nbr_neutral_trinucfreq)
                                             .combine(nbr_driver_regs))*/


    /*
        Step 8: run oncodrivefml
    */
    input_to_soft_split.oncodrivefml.combine(oncodrivefml_config).view()

}

// inform about completition
workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Script finished SUCCESSFULLY after $workflow.duration"
    } else {
      log.info "[$workflow.complete] >> Script finished with ERRORS after $workflow.duration"
    }
}