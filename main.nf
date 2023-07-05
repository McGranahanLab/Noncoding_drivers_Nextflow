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
	input:
	tuple val(patients_inventory_path), val(analysis_inventory_path), 
	      val(blacklist_inventory_path)

	script:
	"""
	if [ -z "${blacklist_inventory_path}" ]
	then
	    1_check_patients_inventory.R --inventory_patients ${patients_inventory_path} \
	                                 --inventory_analysis ${analysis_inventory_path}
    else
        1_check_patients_inventory.R --inventory_patients ${patients_inventory_path} \
	                             --inventory_analysis ${analysis_inventory_path} \
	                             --inventory_blacklisted ${blacklist_inventory_path} 
    fi

	
	"""
}

workflow {
    patients_inv = Channel.of(params.patients_inventory)
    analysis_inv = Channel.of(params.analysis_inventory) 
    blacklist_inv = Channel.of(params.blacklist_inventory) 

    patients_inv.combine(analysis_inv)
                .combine(blacklist_inv) | 
                check_inventories
}