process ANNOTATE_GENOMICRANGES_WITH_CN {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(bed_path), path(subtype_drivers),
          path(patients_inventory_path), path(chain)

    output:
    tuple val(tumor_subtype), path("driversAnnotatedWitCN-${tumor_subtype}--${params.target_genome_version}.csv"), emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    def chain = chain.name != '.NO_FILE' ? "--chain $chain" : ''
    """
    RUN_CODE=$tumor_subtype'--'$params.target_genome_version
    MSG_FILE="annotateWithCN-"\$RUN_CODE'.out'
    ERR_FILE="annotateWithCN-"\$RUN_CODE'.err'
    OUT_FILE="driversAnnotatedWitCN-"\$RUN_CODE'.csv'

    9_match_CN_and_scanned_genomic_regions.R \
            --inventory_patients $patients_inventory_path \
            --genomic_regions $bed_path \
            --cancer_subtype $tumor_subtype \
            --target_genome_version $params.target_genome_version \
            --drivers $subtype_drivers \
            --amp $params.amp --gain $params.gain --loss $params.loss \
            --min_width $params.min_width --min_gapwidth $params.min_gapwidt \
            $chain \
            --output \$OUT_FILE 1>\$MSG_FILE 2>\$ERR_FILE
    """
}

process ANNOTATE_MUTATIONS_WITH_MULTIPLICITY {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(muts_to_gr), path(subtype_drivers),
          path(patients_inventory_path)
    
    output:
    tuple val(tumor_subtype), path("driversAnnotatedWitMultipl-${tumor_subtype}--${params.target_genome_version}.csv"), emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    """
    RUN_CODE=$tumor_subtype'--'$params.target_genome_version
    MSG_FILE="annotateWithMult-"\$RUN_CODE'.out'
    ERR_FILE="annotateWithMult-"\$RUN_CODE'.err'
    OUT_FILE="driversAnnotatedWitMultipl-"\$RUN_CODE'.csv'

    excludeSilent=''
    if [ "$params.exclude_silent_from_biotyping" = "T" ]; then
        excludeSilent="--synAcceptedClass $params.synAcceptedClass"
    fi

    10_match_mutmulti_and_mutations.R \
        --inventory_patients $patients_inventory_path \
        --cancer_subtype $tumor_subtype --muts_to_gr $muts_to_gr \
        --drivers $subtype_drivers \$excludeSilent \
        --output \$OUT_FILE 1>\$MSG_FILE 2>\$ERR_FILE
    """
}