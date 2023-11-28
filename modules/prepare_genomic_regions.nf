process FILTER_GENOMIC_REGIONS {
    input:
    tuple path(analysis_inventory_path), path(blacklist_inventory_path),
          path(target_genome_fasta), path(target_genome_chr_len), path(chain)

    output:
    path('*.bed', emit: bed)
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    def inventory_blacklisted = blacklist_inventory_path.name != '.NO_FILE' ? "--inventory_blacklisted $blacklist_inventory_path" : ''
    def target_genome_chr_len = target_genome_chr_len.name != '.NO_FILE' ? "--target_genome_chr_len $target_genome_chr_len" : ''
    def chain = chain.name != '.NO_FILE' ? "--chain $chain" : ''
    """
    3_filter_input_genomic_regions.R --inventory_analysis $analysis_inventory_path \
                                     --target_genome_path $target_genome_fasta \
                                     --target_genome_version $params.target_genome_version \
                                     --chain $chain \
                                     --ignore_strand $params.ignore_strand \
                                     --min_reg_len $params.min_reg_len \
                                     --output '.' --cores ${task.cpus} \
                                     $inventory_blacklisted \
                                     $target_genome_chr_len $chain\
                                     1>filter_genomic_regions.out \
                                     2>filter_genomic_regions.err
    """
}

process CREATE_RDA_FOR_DNDSCV_DIGDRIVER {
    tag "$tumor_subtype-dndscv-digdriver-rda"

    input:
    tuple val(tumor_subtype), path(gtf), val(gtf_genome_version), path(bed), 
          path(target_genome_fasta), path(chain)

    output:
    tuple val(tumor_subtype), path("*NCBI.Rda"), path("*UCSC.Rda"), emit: rda
    tuple path('*.out'), path('*.err'), emit: logs
    
    script:
    def chain = chain.name != '.NO_FILE' ? "--chain $chain" : ''
    """
    gtf_gv=`echo ${gtf_genome_version} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    3c_write_regions_for_dndscv.R --gtf $gtf --gtf_genomes \$gtf_gv \
                                  --cancer_subtype $tumor_subtype \
                                  --target_genome_path $target_genome_fasta \
                                  --target_genome_version $params.target_genome_version \
                                  --cores ${task.cpus} --output '.' $chain \
                                  1>rda_for_dndscv_digdriver_${tumor_subtype}.out \
                                  2>rda_for_dndscv_digdriver_${tumor_subtype}.err
    """
}

process WRITE_REGIONS_FOR_DIGDRIVER {
    tag "$tumor_subtype-$software"

    input:
    tuple val(tumor_subtype), val(software), val(gr), path(bed)

    output:
    path('*.csv')

    when:
    software == 'digdriver'

    script:
    """
    gr_parsed=`echo ${gr} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    IFS=' ' read -r -a gr_parsed <<< "\$gr_parsed"

    for i in "\${!gr_parsed[@]}"; do
        oneGR="\${gr_parsed[\$i]}"
        outFileOneGR='inputGR-'$tumor_subtype'-digdriver-'\$oneGR'-'$params.target_genome_version'.csv'
        grep -w \$oneGR $bed > \$outFileOneGR
        # digdriver requires NCBI region format
        sed -i 's/^chr//g' \$outFileOneGR
    done
    """
}

process WRITE_REGIONS_FOR_NBR {
    tag "$tumor_subtype-$software"

    input:
    tuple val(tumor_subtype), val(software), val(gr), path(bed)

    output:
    path('*.csv')

    when:
    software == 'nbr'

    script:
    """
    gr_parsed=`echo ${gr} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    3e_write_regions_for_nbr.R --bed $bed --gr_id \$gr_parsed \
                               --cancer_subtype $tumor_subtype \
                               --output '.'
    """
}

process WRITE_REGIONS_FOR_ONCODRIVEFML {
    tag "$tumor_subtype-$software"

    input:
    tuple val(tumor_subtype), val(software), val(gr), path(bed)

    output:
    path('*.csv')

    when:
    software == 'oncodrivefml'

    script:
    """
    gr_parsed=`echo ${gr} | sed 's/\\[//g' | sed 's/,//g' | sed 's/\\]//g'`
    3f_write_regions_for_oncodrivefml.R --bed $bed --gr_id \$gr_parsed \
                                        --cancer_subtype $tumor_subtype \
                                        --output '.'
    """
}