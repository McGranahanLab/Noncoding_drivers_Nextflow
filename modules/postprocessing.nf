process COMBINE_P_VALS_AND_ANNOTATE {
    tag "$tumor_subtype-$gr_id"

	input:
    tuple val(tumor_subtype), val(gr_id), val(softwares), 
          path(result_files), path(mut_rate), path(gene_name_synonyms), 
          path(known_cancer_genes), path(olfactory_genes), val(rawP_cap),
          path(expression_inventories), path(expression)

    script:
    """
    # a unique ID to use in all further commands
    RUN_CODE=${tumor_subtype}'-'${gr_id}'-'${params.target_genome_version}

    MSG_FILE="combineAndAnno-"\$RUN_CODE'.out'
    ERR_FILE="combineAndAnno-"\$RUN_CODE'.err'
    OUT_FILE="combinedP-"\$RUN_CODE'.csv'

    software_parsed=`echo \${softwares} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    result_files_parsed=`echo \${result_files} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    expression_inventories_parsed=`echo \${expression_inventories} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    expression_parsed=`echo \${expression} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`

    6_combine_results_and_annotate_with_metadata.R --cancer_subtype ${tumor_subtype} \
        --software \$software_parsed --run_results \$result_files_parsed \
        --mut_rate ${mut_rate} --gene_name_synonyms ${gene_name_synonyms} \
        --known_cancer_genes ${known_cancer_genes} \
        --olfactory_genes ${olfactory_genes} --rawP_cap ${rawP_cap} \
        --expression_inventories \$expression_inventories_parsed \
        --expression \$expression_parsed --output $OUT_FILE \
        1>$MSG_FILE 2>ERR_FILE
    """
}