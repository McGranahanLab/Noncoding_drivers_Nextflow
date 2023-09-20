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
	debug true

	input:
	tuple val(patients_inventory_path), val(analysis_inventory_path), 
	      val(blacklist_inventory_path)

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
	debug true

	input:
	tuple val(patients_inventory_path), val(analysis_inventory_path), 
	      val(blacklist_inventory_path), val(target_genome_fasta),
	      val(chain), val(tumor_subtype), val(software), 
	      val(inventory_check_res)

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
	                                --output ${params.output}'inputs/' \
	                                --cores ${params.cores}
	
	"""
}

process create_input_genomic_regions_files {
	debug true

	input:
	tuple val(analysis_inventory_path), val(blacklist_inventory_path), 
	      val(target_genome_fasta), val(chain), val(inventory_check_res)

	script:
	"""
	3_create_input_genomic_regions_files.R --inventory_analysis ${analysis_inventory_path} \
	                                       --inventory_blacklisted ${blacklist_inventory_path} \
	                                       --target_genome_path ${target_genome_fasta} \
	                                	   --target_genome_version ${params.target_genome_version} \
	                                       --chain ${chain} \
	                                       --ignore_strand ${params.ignore_strand} \
	                                       --min_reg_len ${params.min_reg_len} \
	                                       --output ${params.output}'inputs/' \
	                                       --cores ${params.cores}
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
                          .ifEmpty { exit 1, "[ERROR]: target genome fasta file not found" }
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
    inventories = patients_inv.combine(analysis_inv)
                              .combine(blacklist_inv)
    inventories_pass = inventories | check_inventories
    inventories_pass = inventories_pass.collect()


    /* 
       Step 2: create input mutations files for all requested software, 
       		   parallelize by cancer subtype. This will not be run if 
       		   inventories do not pass the check.
    */
    tumor_subtypes = analysis_inv.splitCsv(header: true)
    	                         .map(row -> tuple(row.tumor_subtype, row.software))
    	                         .unique().groupTuple()
    inventories.combine(target_genome_fasta)
    	       .combine(chain)
    	       .combine(tumor_subtypes)
    	       .combine(inventories_pass) | create_input_mutation_files

	/* 
       Step 3: create input mutations genomic region files for all requested
       		   software, parallelize by cancer subtype. This will not be run if
       		   inventories do not pass the check.
    */
    analysis_inv.combine(blacklist_inv)
    			.combine(target_genome_fasta)
    	        .combine(chain).combine(inventories_pass) | create_input_genomic_regions_files


    
}
