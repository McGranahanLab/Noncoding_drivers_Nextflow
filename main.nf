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
	                                --cores ${params.cores}
	
	"""
}

process create_input_genomic_regions_files {
    input:
    val inventory_check_res
    path analysis_inventory_path
    path blacklist_inventory_path
    path target_genome_fasta
    path chain

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
                                           --cores ${params.cores}
    """
}

process mutpanning {
    tag "$tumor_subtype"-"gr_id"

	input:
	tuple val(tumor_subtype), val(software), val(gr_id), path(mutations), path(mutpan_inv)

	script:
	"""
	mpPath=/bin/MutPanningV2/commons-math3-3.6.1.jar:/bin/MutPanningV2/jdistlib-0.4.5-bin.jar:/bin/MutPanningV2
        java -Xmx2G -classpath \$mpPath MutPanning
	"""
}

process nbr {
    tag "$tumor_subtype"-"gr_id"

    input:
    tuple val(tumor_subtype), val(software), val(gr_id), path(mutations), path(regions)
    path target_genome_fasta
    
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
          --gr_drivers ${params.assets_dir}"/NBR/GRanges_driver_regions_hg19.txt" \
          --regions_neutralbins_file ${params.assets_dir}"/NBR/Neutral_regions_within_100kb_bins_hg19.txt" \
          --trinucfreq_neutralbins_file ${params.assets_dir}"/NBR/Trinucfreqs_within_100kb_bins_hg19.txt" \
          --out_prefix \$OUT_FILE
    mv $OUT_FILE"-Selection_output.txt" \$OUT_FILE 
    OUT_FILE_BASE=`echo \$OUT_FILE | sed 's/.csv//g'`
    mv \$OUT_FILE"-global_mle_subs.Rds" \$OUT_FILE_BASE"-global_mle_subs.Rds" 
    mv \$OUT_FILE"-globalRates.csv" \$OUT_FILE_BASE"-globalRates.csv"
    """
}

process oncodrivefml {
    tag "$tumor_subtype"-"gr_id"

    input:
    tuple val(tumor_subtype), val(software), val(gr_id), path(mutations), path(regions)

    script:
    """
    echo ${mutations}
    echo ${regions} 

    oncodrivefml -i ${mutations} -e ${regions} --sequencing wgs \
                 --configuration 'conf/oncodrivefml_'${params.target_genome_version}'.config' \
                 --output '.'  --cores $CORES_TO_USE
    oncodrBaseName="oncodrivefml-inputMutations-"$TUMOR"-"$TARGET_GENOME
    mv $oncodrBaseName"-oncodrivefml.tsv" $OUT_FILE
    OUT_FILE_BASE=`echo $OUT_FILE | sed 's/.csv//g'`
    mv $oncodrBaseName"-oncodrivefml.png" $OUT_FILE_BASE".png"
    mv $oncodrBaseName"-oncodrivefml.html" $OUT_FILE_BASE".html"
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
    if (params.blacklist_inventory == '') {
        def no_file = new File(".NO_FILE")
        no_file.createNewFile()
        blacklist_inv = Channel.fromPath(".NO_FILE")
    } else {
        blacklist_inv = Channel.fromPath(params.blacklist_inventory, checkIfExists: true)
                               .ifEmpty { exit 1, "[ERROR]: black&white lists inventory file not found" }
    } 

    // create channels to target genome verion and to chain file for liftover
    target_genome_fasta = Channel.fromPath(params.target_genome_path,
                                           checkIfExists: true)
                                 .ifEmpty { exit 1, 
                                            "[ERROR]: target genome fasta file not found" }
    if (params.chain == '') {
        def no_file = new File(".NO_FILE")
        no_file.createNewFile()
        chain = Channel.fromPath(".NO_FILE")
    } else {
        chain = Channel.fromPath(params.chain, checkIfExists: true)
                       .ifEmpty { exit 1, "[ERROR]: chain file not found" }
    }

    /*
        Step 1: check that inventories have all the needed columns and values
                are acceptable 
    */
    inventories_pass =  check_inventories(patients_inv, analysis_inv, blacklist_inv)
    inventories_pass = inventories_pass.collect()

    /* 
        Step 2: create input mutations files for all requested software, 
                parallelise by cancer subtype. This will not be run if 
                inventories do not pass the check.
    */
    tumor_subtypes = analysis_inv.splitCsv(header: true)
                                 .map{row -> tuple(row.tumor_subtype, row.software)}
                                 .unique().groupTuple()
    input_mutations =  create_input_mutation_files(tumor_subtypes.combine(inventories_pass)
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
    input_genomic_regions = create_input_genomic_regions_files(inventories_pass, analysis_inv,
                                                               blacklist_inv, target_genome_fasta,
                                                               chain)
                                                               .flatten()
    // split on csv which will go to be processed by driver calling software and bed files
    // which will go to estimation of mutation rates
    input_genomic_regions.branch {
        csv: it.toString().endsWith('csv')
             return it
        bed: true
             return it
    }.set{ input_genomic_regions_split }

    input_genomic_regions_split.bed.view()
    input_mutations_split.maf.view()

    /*
        Step 4: prepare channel input_to_soft which will contain pairs of 
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
    input_to_soft = input_mutations_split_to_soft.join(input_genomic_regions_to_soft,
                                                       by: [0, 1, 2], remainder: true)
    input_to_soft.branch {
        mutpanning: it[1].toString() == 'mutpanning'
            return it
        oncodrivefml: it[1].toString() == 'oncodrivefml'
            return it
    }.set{ input_to_soft_split }

    /*input_to_soft_split.mutpanning.map { it ->
        def mutpan_inv = it[3].toString().replace("inputMutations", "patientsInv")
        return tuple(it[0], it[1], it[2], it[3], mutpan_inv)
    } | mutpanning*/
}
