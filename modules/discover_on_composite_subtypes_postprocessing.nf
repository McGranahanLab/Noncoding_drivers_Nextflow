process CHECK_DISCOVER_ON_COMPOSITE_SUBTYPES {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(discover_composite_subtype),
          val(uniformal_tumor_subtypes), path(discover_uniformal_subtypes)
          
    output:
    tuple val(tumor_subtype), path("discoverResults-${tumor_subtype}--${params.target_genome_version}.csv"), emit: csv
    tuple path('*.out'), path('*.err'), emit: logs
    
    script:
    """
    RUN_CODE=$tumor_subtype'--'$params.target_genome_version
    MSG_FILE="discoverCheckCompositeSubtypes-"\$RUN_CODE'.out'
    ERR_FILE="discoverCheckCompositeSubtypes-"\$RUN_CODE'.err'
    OUT_FILE="discoverResults-"\$RUN_CODE'.csv'
    
    uniformal_tumor_subtypes_parsed=`echo ${uniformal_tumor_subtypes} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    discover_uniformal_subtypes_parsed=`echo ${discover_uniformal_subtypes} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    
    16_cooccurrence_and_exclusivity_supported_by_individual_subtypes.R \
      --cancer_subtype $tumor_subtype \
      --discover_composite_subtype $discover_composite_subtype \
      --uniformal_tumor_subtypes \uniformal_tumor_subtypes_parsed \
      --discover_uniformal_subtypes \$discover_uniformal_subtypes_parsed \
      --p_max_unif_subtype $params.p_max_uniformal_subtypes_discover \
      --output \$OUT_FILE 1>\$MSG_FILE 2>\$ERR_FILE
    """
}