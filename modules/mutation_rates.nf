process CALCULATE_MUTATION_RATES {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(bed_path), path(maf_path), 
          path(target_genome_chr_len), path(gene_name_synonyms),
          path(varanno_conversion_table), val(coding_gr_id)

    output:
    tuple val(tumor_subtype), path('meanMutRatePerGR*'), path('mutMapToGR*'), path('varCatEnrich*'), emit: mutrates
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    def gene_name_synonyms       = !gene_name_synonyms.name.startsWith(params.empty_file_prefix) ? "--gene_name_synonyms $gene_name_synonyms" : ''
    def target_genome_chr_len    = !target_genome_chr_len.name.startsWith(params.empty_file_prefix) ? "--target_genome_chr_len $target_genome_chr_len" : ''
    def varanno_conversion_table = !varanno_conversion_table.name.startsWith(params.empty_file_prefix) ? "--varanno_conversion_table $varanno_conversion_table" : ''
    """
    coding_gr_id_parsed=`echo $coding_gr_id | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'` 
    4_calculate_mutation_rates.R --cancer_subtype $tumor_subtype \
                                 --variants $maf_path \
                                 --genomic_regions $bed_path \
                                 --bin_len $params.bin_len \
                                 --target_genome_version $params.target_genome_version \
                                 --calc_synonymous $params.calc_synonymous \
                                 --remove_synonymous_from_coding $params.remove_synon_from_coding \
                                 --cdsAcceptedClass $params.cdsAcceptedClass \
                                 --ncAcceptedClass $params.ncAcceptedClass \
                                 --coding_gr_id \$coding_gr_id_parsed \
                                 --synAcceptedClass $params.synAcceptedClass \
                                 --annotation_failed_code $params.annotation_failed_code \
                                 --output '.' $gene_name_synonyms $target_genome_chr_len \
                                 $varanno_conversion_table \
                                 1>mut_rates_${tumor_subtype}.out \
                                 2>mut_rates_${tumor_subtype}.err
    """
}
