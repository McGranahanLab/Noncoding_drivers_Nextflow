process CALCULATE_MUTATION_RATES {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(bed_path), path(maf_path), 
          path(target_genome_chr_len), path(gene_name_synonyms),
          path(varanno_conversion_table)

    output:
    tuple val(tumor_subtype), path('meanMutRatePerGR*'), path('mutMapToGR*'), path('varCatEnrich*'), emit: mutrates
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    def gene_name_synonyms = gene_name_synonyms.name != '.NO_FILE' ? "--gene_name_synonyms $gene_name_synonyms" : ''
    def target_genome_chr_len = target_genome_chr_len.name != '.NO_FILE' ? "--target_genome_chr_len $target_genome_chr_len" : ''
    def varanno_conversion_table = varanno_conversion_table.name != '.NO_FILE' ? "--varanno_conversion_table $varanno_conversion_table" : ''
    """
    4_calculate_mutation_rates.R --cancer_subtype $tumor_subtype \
                                 --variants $maf_path \
                                 --genomic_regions $bed_path \
                                 --bin_len $params.bin_len \
                                 --target_genome_version $params.target_genome_version \
                                 --calc_synonymous $params.calc_synonymous \
                                 --cdsAcceptedClass $params.cdsAcceptedClass \
                                 --synAcceptedClass $params.synAcceptedClass \
                                 --ncAcceptedClass $params.ncAcceptedClass \
                                 --annotation_failed_code $params.annotation_failed_code \
                                 --output '.' $gene_name_synonyms $target_genome_chr_len \
                                 $varanno_conversion_table \
                                 1>mut_rates_${tumor_subtype}.out \
                                 2>mut_rates_${tumor_subtype}.err
    """
}
