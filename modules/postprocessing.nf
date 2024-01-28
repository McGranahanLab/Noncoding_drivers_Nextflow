process COMBINE_P_VALS_AND_ANNOTATE {
    tag "$tumor_subtype-$gr_id"

	input:
    tuple val(tumor_subtype), val(gr_id), path(mut_map), path(mut_rate), 
          path(scanned_gr), path(var_cat_enr), val(softwares), 
          path(result_files),
          path(gene_name_synonyms), path(known_cancer_genes),
          path(olfactory_genes), 
          path(inventory_gtex), path(expression_gtex),
          path(inventory_tcga), path(expression_tcga)

    output:
    tuple val(tumor_subtype), path("combinedP-${tumor_subtype}-${gr_id}-${params.target_genome_version}.csv"), emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    def synonyms = gene_name_synonyms.name != '.NO_FILE' ? "--gene_name_synonyms $gene_name_synonyms" : ''
    def known_genes = known_cancer_genes.name != '.NO_FILE' ? "--known_cancer_genes $known_cancer_genes" : ''
    def olfactory = olfactory_genes.name != '.NO_FILE' ? "--olfactory_genes $olfactory_genes" : ''
    """
    # a unique ID to use in all further commands
    RUN_CODE=$tumor_subtype'-'$gr_id'-'$params.target_genome_version
    MSG_FILE="combineAndAnno-"\$RUN_CODE'.out'
    ERR_FILE="combineAndAnno-"\$RUN_CODE'.err'
    OUT_FILE="combinedP-"\$RUN_CODE'.csv'

    software_parsed=`echo $softwares | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    result_files_parsed=`echo $result_files | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`

    6_combine_results_and_annotate_with_metadata.R --cancer_subtype $tumor_subtype \
        --gr_id $gr_id --software \$software_parsed \
        --run_results \$result_files_parsed --mut_rate $mut_rate \
        --scanned_gr $scanned_gr \
        --rawP_cap $params.rawP_cap \
        --var_cat_enrich $var_cat_enr \
        --var_cat_enrich_padj $params.padj_2_5bp_enrich \
        --expression_inventories $inventory_gtex $inventory_tcga \
        --expression $expression_gtex $expression_tcga \
        --output \$OUT_FILE \
        $synonyms $known_genes $olfactory \
        1>\$MSG_FILE 2>\$ERR_FILE
    """
}

process ASSIGN_TIER {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(combined_p_values_tabs), 
          path(inventory_tier), val(combine_p_values_with)

    output:
    tuple val(tumor_subtype), path("assignedTier-${tumor_subtype}--${params.target_genome_version}.csv"), emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    """
    # a unique ID to use in all further commands
    RUN_CODE=$tumor_subtype'--'$params.target_genome_version
    MSG_FILE="assign_tier-"\$RUN_CODE'.out'
    ERR_FILE="assign_tier-"\$RUN_CODE'.err'
    OUT_FILE="assignedTier-"\$RUN_CODE'.csv'

    7_assign_tier.R --combined_p_values_tables $combined_p_values_tabs \
                    --inventory_tier $inventory_tier \
                    --combine_p_value_method $combine_p_values_with \
                    --output \$OUT_FILE 1>\$MSG_FILE 2>\$ERR_FILE
    """
}

process FILTER_TIERED_DRIVERS {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(tiered_values),
          path(analysis_inventory_path)

    output:
    tuple val(tumor_subtype), path("drivers-${tumor_subtype}--${params.target_genome_version}.csv"),
          env(n_drivers), emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    """
    # a unique ID to use in all further commands
    RUN_CODE=$tumor_subtype'--'$params.target_genome_version
    MSG_FILE="filtering-"\$RUN_CODE'.out'
    ERR_FILE="filtering-"\$RUN_CODE'.err'
    OUT_FILE="drivers-"\$RUN_CODE'.csv'

    8_filter_tiered_drivers.R --inventory_analysis $analysis_inventory_path \
                              --cancer_subtype $tumor_subtype \
                              --tiered_pvals $tiered_values \
                              --min_n_soft_cod $params.min_n_soft_cod \
                              --min_n_soft_noncod $params.min_n_soft_noncod \
                              --min_n_muts $params.min_n_muts \
                              --min_n_patients $params.min_n_patients \
                              --max_local_mut_rate_q $params.max_local_mut_rate_q \
                              --max_gr_syn_mut_rate_q $params.max_gr_syn_mut_rate_q \
                              --max_gr_mut_rate_q $params.max_gr_mut_rate_q \
                              --max_gr_len_q $params.max_gr_len_q \
                              --remove_olfactory $params.remove_olfactory \
                              --remove_2_5bp_enrich $params.remove_2_5bp_enrich \
                              --output \$OUT_FILE 1>\$MSG_FILE 2>\$ERR_FILE

    # check, that there are indeed drivers detected by looking at column FILTER
    # and column tier. FILTER should be PASS and tier should not be NA.
    tier_col_idx=`head -1 \$OUT_FILE | tr '\t' '\n' | cat -n | grep -w "tier" | grep -v '[.]' | cut -f1 | sed 's/ //g'`
    filter_col_idx=`head -1 \$OUT_FILE | tr '\t' '\n' | cat -n | grep -w "FILTER" | cut -f1 | sed 's/ //g'`
    n_drivers=`cut -f\$tier_col_idx,\$filter_col_idx \$OUT_FILE | grep PASS | grep -v NA | wc -l`
    if [ "\$n_drivers" -eq "0" ]; then
        n_drivers='no'
    else
        n_drivers='yes'
    fi
    """
}