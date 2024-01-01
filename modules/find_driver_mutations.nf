process FIND_DRIVER_MUTATIONS {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(gr_id), path(muts_to_gr), path(run_result),
          path(subtype_drivers), path(chasmplus), path(cn_of_drivers), 
          path(inferred_biotype), path(analysis_inventory_path), 
          path(patients_inventory_path), path(known_driver_mutations)

    output:
    tuple val(tumor_subtype), path("driverMutations-${tumor_subtype}--${params.target_genome_version}.csv"), emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    def knownMuts = known_driver_mutations.name != '.NO_FILE' ? "--known_driver_mutations $known_driver_mutations" : ''
    """
    RUN_CODE=$tumor_subtype'--'$params.target_genome_version
    MSG_FILE="findDriverMuts-"\$RUN_CODE'.out'
    ERR_FILE="findDriverMuts-"\$RUN_CODE'.err'
    OUT_FILE="driverMutations-"\$RUN_CODE'.csv'

    #def chasm = 
    #def inferred_biotype =
    #def cn_of_drivers =

    excludeSilent=''
    if [ "$params.exclude_silent_from_biotyping" = "T" ]; then
        excludeSilent="--synAcceptedClass $params.synAcceptedClass"
    fi

    12_find_driver_mutations.R \
        --inventory_analysis $analysis_inventory_path \
        --inventory_patients $patients_inventory_path \
        --cancer_subtype $tumor_subtype \
        --drivers $subtype_drivers  \
        --gr_id $gr_id \
        --run_result $run_result \
        --muts_to_gr $muts_to_gr \
        --chasmplus $chasmplus \
        --inferred_biotype $inferred_biotype \
        --cn_of_drivers $cn_of_drivers \
        --chasm_score_min $params.chasm_score_min \
        --chasm_padj $params.chasm_padj \
        --max_fp $params.max_fp $knownMuts \$excludeSilent \
        --output \$OUT_FILE 1>\$MSG_FILE 2>\$ERR_FILE
    """
}
