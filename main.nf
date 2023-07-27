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

/*
patients_inventory = Channel.fromPath(params.patients_inventory, checkIfExists: true)
                            .ifEmpty { exit 1, "[ERROR]: patients inventory file not found" }
                            .splitCsv(header: true)
                            .unique()
                            .view()

analysis_inventory = Channel.fromPath(params.analysis_inventory, checkIfExists: true)
                            .ifEmpty { exit 1, "[ERROR]: analysis inventory file not found" }
                            .splitCsv(header: true)
                            .view()
*/



/* ----------------------------------------------------------------------------
* Processes
*----------------------------------------------------------------------------*/
process check_inventories {
	debug true
	
	input:
	tuple val(patients_inventory_path), val(analysis_inventory_path), 
	      val(blacklist_inventory_path)

	script:
	"""
	if echo "${blacklist_inventory_path}" | grep -q "EMPTY_FILE.txt"
	then
	    1_check_inventories.R --inventory_patients ${patients_inventory_path} \
	                          --inventory_analysis ${analysis_inventory_path} \
	                          --target_genome_version ${params.target_genome_version}
    else
        1_check_inventories.R --inventory_patients ${patients_inventory_path} \
	                          --inventory_analysis ${analysis_inventory_path} \
	                          --target_genome_version ${params.target_genome_version} \
	                          --inventory_blacklisted ${blacklist_inventory_path} 
    fi
	
	"""
}

/* ----------------------------------------------------------------------------
* Workflows
*----------------------------------------------------------------------------*/
workflow {
    patients_inv = Channel.fromPath(params.patients_inventory)
    analysis_inv = Channel.fromPath(params.analysis_inventory) 
    blacklist_inv = Channel.fromPath(params.blacklist_inventory) 

    patients_inv.combine(analysis_inv)
                .combine(blacklist_inv) | 
                check_inventories
}