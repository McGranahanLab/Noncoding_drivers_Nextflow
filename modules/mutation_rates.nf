process CALCULATE_MUTATION_RATES {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(bed_path), path(maf_path), 
          path(target_genome_chr_len), path(gene_name_synonyms),
          path(varanno_conversion_table)

    output:
    tuple val(tumor_subtype), path('*meanMutRatePerGR.csv'), path('*mutMapToGR.csv'), path('*varCatEnrich.csv'), emit: mutrates
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    """
    4_calculate_mutation_rates.R --tumor_subtype ${tumor_subtype} \
                                 --variants ${maf_path} \
                                 --genomic_regions ${bed_path} \
                                 --gene_name_synonyms ${gene_name_synonyms} \
                                 --bin_len ${params.bin_len} \
                                 --target_genome_version ${params.target_genome_version} \
                                 --target_genome_chr_len ${target_genome_chr_len} \
                                 --calc_synonymous ${params.calc_synonymous} \
                                 --cdsAcceptedClass ${params.cdsAcceptedClass} \
                                 --ncAcceptedClass ${params.ncAcceptedClass} \
                                 --varanno_conversion_table ${varanno_conversion_table} \
                                 --annotation_failed_code ${params.annotation_failed_code} \
                                 --output 'mutRate-'${tumor_subtype}'-' \
                                 1>mut_rates_${tumor_subtype}.out \
                                 2>mut_rates_${tumor_subtype}.err
    """
}