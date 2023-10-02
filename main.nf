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

process filter_input_mutations {
    tag "$tumor_subtype"

    input:
    tuple val(inventory_check_res), val(tumor_subtype),
          path(patients_inventory_path), path(analysis_inventory_path),
          path(blacklist_inventory_path), path(target_genome_fasta),
          path(target_genome_chr_len), path(chain)

    output:
    tuple val(tumor_subtype), path('*.maf'), emit: mutations
    tuple val(tumor_subtype), path('hypermutated*'), optional: true, emit: hypermutant
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    """ 
    OUT_FILE='inputMutations-'${tumor_subtype}'-'${params.target_genome_version}'.maf'
    2_filter_input_mutations.R --inventory_patients ${patients_inventory_path} \
                               --inventory_analysis ${analysis_inventory_path} \
                               --inventory_blacklisted ${blacklist_inventory_path} \
                               --cancer_subtype ${tumor_subtype} \
                               --min_depth ${params.min_depth} \
                               --min_tumor_vac ${params.min_tumor_vac} \
                               --min_tumor_vaf ${params.min_tumor_vaf} \
                               --max_germline_vaf ${params.max_germline_vaf} \
                               --max_germline_vac ${params.max_germline_vac} \
                               --max_n_vars ${params.max_n_vars} \
                               --target_genome_version ${params.target_genome_version} \
                               --target_genome_path ${target_genome_fasta} \
                               --target_genome_chr_len ${target_genome_chr_len} \
                               --chain ${chain} --output \$OUT_FILE  \
                               --cores ${task.cpus} \
                               1>filter_mutations_${tumor_subtype}.out \
                               2>filter_mutations_${tumor_subtype}.err
    
    """
}

process write_mutations_for_chasmplus {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'chasmplus'

    script:
    """
    OUT_FILE='inputMutations-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    2a_write_mutations_for_chasmplus.R --maf ${maf} \
                                       --cancer_subtype ${tumor_subtype} \
                                       --output \$OUT_FILE
    """
}

process write_mutations_for_digdriver {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'digdriver'

    script:
    """
    OUT_FILE='inputMutations-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    2b_write_mutations_for_digdriver.R --maf ${maf} \
                                       --cancer_subtype ${tumor_subtype} \
                                       --output \$OUT_FILE
    """
}

process write_mutations_for_dndscv {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'dndscv'

    script:
    """
    OUT_FILE='inputMutations-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    2c_write_mutations_for_dndscv.R --maf ${maf} \
                                    --cancer_subtype ${tumor_subtype} \
                                    --output \$OUT_FILE
    """
}


process write_mutations_for_mutpanning {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*'), 
          path('mutpanning-patientsInv*')

    when:
    software == 'mutpanning'

    script:
    """
    OUT_FILE='inputMutations-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    2d_write_mutations_for_mutpanning.R --maf ${maf} \
                                        --cancer_subtype ${tumor_subtype} \
                                        --output \$OUT_FILE
    """
}

process write_mutations_for_nbr {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'nbr'

    script:
    """
    OUT_FILE='inputMutations-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    2e_write_mutations_for_nbr.R --maf ${maf} \
                                 --cancer_subtype ${tumor_subtype} \
                                 --output \$OUT_FILE
    """
}

process write_mutations_for_oncodrivefml {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'oncodrivefml'

    script:
    """
    OUT_FILE='inputMutations-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    2f_write_mutations_for_oncodrivefml.R --maf ${maf} \
                                          --cancer_subtype ${tumor_subtype} \
                                          --output \$OUT_FILE
    """
}

process filter_genomic_regions {
    input:
    tuple val(inventory_check_res), path(analysis_inventory_path),
          path(blacklist_inventory_path), path(target_genome_fasta),
          path(target_genome_chr_len), path(chain)

    output:
    path('*.bed'), emit: bed
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    """
    3_filter_input_genomic_regions.R --inventory_analysis ${analysis_inventory_path} \
                                     --inventory_blacklisted ${blacklist_inventory_path} \
                                     --target_genome_path ${target_genome_fasta} \
                                     --target_genome_chr_len ${target_genome_chr_len} \
                                     --target_genome_version ${params.target_genome_version} \
                                     --chain ${chain} \
                                     --ignore_strand ${params.ignore_strand} \
                                     --min_reg_len ${params.min_reg_len} \
                                     --output '.' --cores ${task.cpus} \
                                     1>filter_genomic_regions.out \
                                     2>filter_genomic_regions.err
    """
}

process write_regions_for_digdriver {
    tag "$tumor_subtype-$software"

    input:
    tuple val(software), val(gr), val(tumor_subtype), path(bed)
    val target_genome_version

    output:
    path('*.csv')

    when:
    software == 'digdriver'

    script:
    """
    IFS='--' read -r -a gr_arr <<< "$gr"
    IFS='--' read -r -a subtypes_arr <<< "$tumor_subtype"

    tumor_subtype=\${subtypes_arr[0]}
    for i in "\${!gr_arr[@]}"; do
        if [[ -n "\${gr_arr[\$i]}" ]]; then
            oneGR="\${gr_arr[\$i]}"
            outFileOneGR='inputGR-'\$tumor_subtype'-digdriver-'\$oneGR'-'${target_genome_version}'.csv'
            grep -w \$oneGR ${bed} > \$outFileOneGR
        fi
    done

    # in case set of the same regions is analysed for a lot of tumor subtypes
    # it is faster to create region file for one tumor subtype and then copy 
    # it for the other tumor subtypes (for one software of course). This is 
    # what is done here                 
    if [ "\${#subtypes_arr[@]}" -gt 1 ]; then
        for i in "\${!subtypes_arr[@]}"; do
            if [[ "\$i" -gt 0 ]] && [[ -n "\${subtypes_arr[\$i]}" ]]; then
                cpFrom=(\$(find . | grep \${subtypes_arr[0]}))
                for j in "\${!cpFrom[@]}"; do
                    cpTo=`echo "\${cpFrom[\$j]}" | sed "s/\${subtypes_arr[0]}/\${subtypes_arr[i]}/g"`
                    cp "\${cpFrom[\$j]}" \$cpTo
                    
                done
            fi
        done
    fi
    """
}

process write_regions_for_nbr {
    tag "$tumor_subtype-$software"

    input:
    tuple val(software), val(gr), val(tumor_subtype), path(bed)

    output:
    path('*.csv')

    when:
    software == 'nbr'

    script:
    """
    gr_arr=`echo ${gr} | sed 's/--/ /g'`
    IFS='--' read -r -a subtypes_arr <<< "$tumor_subtype"

    3e_write_regions_for_nbr.R --bed ${bed} \
                               --cancer_subtype \${subtypes_arr[0]} \
                               --gr_id \$gr_arr --output '.'

    # in case set of the same regions is analysed for a lot of tumor subtypes
    # it is faster to create region file for one tumor subtype and then copy 
    # it for the other tumor subtypes (for one software of course). This is 
    # what is done here                 
    if [ "\${#subtypes_arr[@]}" -gt 1 ]; then
        for i in "\${!subtypes_arr[@]}"; do
            if [[ "\$i" -gt 0 ]] && [[ -n "\${subtypes_arr[\$i]}" ]]; then
                cpFrom=(\$(find . | grep \${subtypes_arr[0]}))
                for j in "\${!cpFrom[@]}"; do
                    cpTo=`echo "\${cpFrom[\$j]}" | sed "s/\${subtypes_arr[0]}/\${subtypes_arr[i]}/g"`
                    cp "\${cpFrom[\$j]}" \$cpTo
                    
                done
            fi
        done
    fi
    """
}

process write_regions_for_oncodrivefml {
    tag "$tumor_subtype-$software"

    input:
    tuple val(software), val(gr), val(tumor_subtype), path(bed)

    output:
    path('*.csv')

    when:
    software == 'oncodrivefml'

    script:
    """
    gr_arr=`echo ${gr} | sed 's/--/ /g'`
    IFS='--' read -r -a subtypes_arr <<< "$tumor_subtype"

    3f_write_regions_for_oncodrivefml.R --bed ${bed} \
                                        --cancer_subtype \${subtypes_arr[0]} \
                                        --gr_id \$gr_arr --output '.'

    # in case set of the same regions is analysed for a lot of tumor subtypes
    # it is faster to create region file for one tumor subtype and then copy 
    # it for the other tumor subtypes (for one software of course). This is 
    # what is done here                 
    if [ "\${#subtypes_arr[@]}" -gt 1 ]; then
        for i in "\${!subtypes_arr[@]}"; do
            if [[ "\$i" -gt 0 ]] && [[ -n "\${subtypes_arr[\$i]}" ]]; then
                cpFrom=(\$(find . | grep \${subtypes_arr[0]}))
                for j in "\${!cpFrom[@]}"; do
                    cpTo=`echo "\${cpFrom[\$j]}" | sed "s/\${subtypes_arr[0]}/\${subtypes_arr[i]}/g"`
                    cp "\${cpFrom[\$j]}" \$cpTo
                    
                done
            fi
        done
    fi
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

process dndscv {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), val(gr_id), path(mutations),
          path(regions)

    output:
    path "${software}-results-${tumor_subtype}-${gr_id}-${params.target_genome_version}*"

    script:
    """
    OUT_FILE=${software}"-results-"${tumor_subtype}"-"${gr_id}'-'${params.target_genome_version}'.csv'

    # covariates
    rdaWithCovs=`echo ${regions} | tr ' ' '\n' | grep withCovs.rda\$`

    run_dndscv.R --with_covariates T --refRda \$rdaWithCovs \
                 --variants ${mutations} --computeCI T --output \$OUT_FILE
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
    tag "$tumor_subtype-$gr_id"

    input:
    tuple val(tumor_subtype), val(software), val(gr_id), path(mutations), 
          path(regions), path(oncodrivefml_config)

    output:
    path "${software}-results-${tumor_subtype}-${gr_id}-${params.target_genome_version}*"

    script:
    """
    OUT_FILE=${software}"-results-"${tumor_subtype}"-"${gr_id}'-'${params.target_genome_version}'.csv'

    oncodrivefml -i ${mutations} -e ${regions} --sequencing wgs \
                 --configuration ${oncodrivefml_config} \
                 --output '.'  --cores ${task.cpus}
    mv "oncodrivefml-inputMutations-"${tumor_subtype}"-"${params.target_genome_version}"-oncodrivefml.tsv" \$OUT_FILE
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
    target_genome_chr_len = Channel.fromPath(params.target_genome_chr_len,
                                           checkIfExists: true)
                                 .ifEmpty { exit 1, 
                                            "[ERROR]: chromosomal lengths of target genome file not found" }
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
    tumor_subtypes_and_gr =  analysis_inv.splitCsv(header: true)
                                         .map{row -> tuple(row.tumor_subtype, row.software, row.gr_id)}
                                         .unique().combine(filtered_regions, by: 0)
                                         .groupTuple(by: [0, 1])
                                         .map { it ->
                                             def gr_id_combo = it[2].sort(mutate = false).join('--')
                                             def bed_file = it[3].unique().flatten()
                                             return tuple(it[1], gr_id_combo, it[0], bed_file)
                                         }
                                         .groupTuple(by: [0, 1])
                                         .map { it ->
                                             def tum_st_combo = it[2].sort(mutate = false).join('--')
                                             return tuple(it[0], it[1], tum_st_combo, it[3].flatten()[0])
                                         }
    nbr_regions = write_regions_for_nbr(tumor_subtypes_and_gr)
    oncodrivefml_regions = write_regions_for_oncodrivefml(tumor_subtypes_and_gr)
    digdriver_regions = write_regions_for_digdriver(tumor_subtypes_and_gr, params.target_genome_version)
}

// inform about completition
workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Script finished SUCCESSFULLY after $workflow.duration"
    } else {
      log.info "[$workflow.complete] >> Script finished with ERRORS after $workflow.duration"
    }
}