process PLOT_OVERVIEW_OF_DRIVERS {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(patients_inventory_path),
          path(drivers_composite_subtype), path(drivers_uniform_subtypes),
          path(excluded_patients), path(visual_json), path(extra_studies)

    output:
    tuple val(tumor_subtype), path("driversOverview-${tumor_subtype}--${params.target_genome_version}.${params.plot_output_type}"), emit: plot, optional: true
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    def composite_cancer_subtype   = drivers_composite_subtype.name != '.NO_FILE' ? "--composite_cancer_subtype $tumor_subtype" : ''
    def drivers_composite_subtype  = drivers_composite_subtype.name != '.NO_FILE' ? "--drivers_composite_subtype $drivers_composite_subtype" : ''
    def excluded_patients          = !excluded_patients.name.startsWith('.NO_FILE') ? "--excluded_patients $excluded_patients" : ''
    def allowed_filter_values      = params.allowed_filter_values.size != 0 ? "--allowed_filter_values "+params.allowed_filter_values.collect{ '"' + it + '"'}.join(" ") : ''
    def extra_studies_names        = params.extra_studies.size != 0 ? "--extra_studies_names "+params.extra_studies_names.collect{ '"' + it + '"'}.join(" ") : ''
    def extra_studies_tumorsubtype = params.extra_studies.size != 0 ? "--extra_studies_tumorsubtype "+params.extra_studies_tumorsubtype.collect{ '"' + it + '"'}.join(" ") : ''
    """
    # a unique ID to use in all further commands
    RUN_CODE=$tumor_subtype'--'$params.target_genome_version
    MSG_FILE="plotOverview-"\$RUN_CODE'.out'
    ERR_FILE="plotOverview-"\$RUN_CODE'.err'
    OUT_FILE="driversOverview-"\$RUN_CODE'.'$params.plot_output_type

    figures_create_overview_plot.R   $composite_cancer_subtype \
                                     $drivers_composite_subtype \
        --inventory_patients         $patients_inventory_path \
                                     $excluded_patients \
                                     $allowed_filter_values \
        --drivers_uniform_subtypes   $drivers_uniform_subtypes \
        --extra_studies              $extra_studies \
                                     $extra_studies_names \
                                     $extra_studies_tumorsubtype \
        --visuals_json               $visual_json \
        --output_type                $params.plot_output_type \
        --output                     \$OUT_FILE \
        1>\$MSG_FILE 2>\$ERR_FILE
    
    """
}