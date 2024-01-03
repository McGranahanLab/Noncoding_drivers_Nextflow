process CALCULATE_SELECTION_RATES {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(gr_id), path(result_files), path(drivers)

    output:
    tuple val(tumor_subtype), path("selectionRates-${tumor_subtype}--${params.target_genome_version}.csv"), emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    """
    RUN_CODE=$tumor_subtype'--'$params.target_genome_version
    MSG_FILE="calcSelectionRates-"\$RUN_CODE'.out'
    ERR_FILE="calcSelectionRates-"\$RUN_CODE'.err'
    OUT_FILE="selectionRates-"\$RUN_CODE'.csv'

    gr_id_parsed=`echo $gr_id | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`

    13_calculate_selection_rates.R \
        --cancer_subtype $tumor_subtype \
        --gr_id \$gr_id_parsed \
        --run_result $result_files \
        --drivers $drivers \
        --output \$OUT_FILE 1>\$MSG_FILE 2>\$ERR_FILE
    """
}