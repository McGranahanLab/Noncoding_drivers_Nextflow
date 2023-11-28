process FILTER_INPUT_MUTATIONS {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(patients_inventory_path),
          path(analysis_inventory_path), path(blacklist_inventory_path),
          path(target_genome_fasta), path(target_genome_chr_len), path(chain)

    output:
    tuple val(tumor_subtype), path('*.maf'), emit: mutations
    tuple val(tumor_subtype), path('hypermutated*'), optional: true, emit: hypermutant
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    def inventory_blacklisted = blacklist_inventory_path.name != '.NO_FILE' ? "--inventory_blacklisted $blacklist_inventory_path" : ''
    def target_genome_chr_len = target_genome_chr_len.name != '.NO_FILE' ? "--target_genome_chr_len $target_genome_chr_len" : ''
    def chain = chain.name != '.NO_FILE' ? "--chain $chain" : ''
    """ 
    OUT_FILE='inputMutations-'$tumor_subtype'-'$params.target_genome_version'.maf'
    2_filter_input_mutations.R --inventory_patients $patients_inventory_path \
                               --inventory_analysis $analysis_inventory_path \
                               --cancer_subtype $tumor_subtype \
                               --min_depth $params.min_depth \
                               --min_tumor_vac $params.min_tumor_vac \
                               --min_tumor_vaf $params.min_tumor_vaf \
                               --max_germline_vaf $params.max_germline_vaf \
                               --max_germline_vac $params.max_germline_vac \
                               --max_n_vars $params.max_n_vars \
                               --target_genome_version $params.target_genome_version \
                               --target_genome_path $target_genome_fasta \
                               --output \$OUT_FILE  --cores ${task.cpus} \
                               $inventory_blacklisted $target_genome_chr_len \
                               $chain \
                               1>filter_mutations_${tumor_subtype}.out \
                               2>filter_mutations_${tumor_subtype}.err
    
    """
}

process WRITE_MUTATIONS_FOR_CHASMPLUS {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'chasmplus'

    script:
    """
    OUT_FILE='inputMutations-'$tumor_subtype'-'$software'-'$params.target_genome_version'.csv'
    2a_write_mutations_for_chasmplus.R --maf $maf \
                                       --cancer_subtype $tumor_subtype \
                                       --output \$OUT_FILE
    """
}

process WRITE_MUTATIONS_FOR_DIGDRIVER {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'digdriver'

    script:
    """
    OUT_FILE='inputMutations-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    2b_write_mutations_for_digdriver.R --maf ${maf} \
                                       --cancer_subtype ${tumor_subtype} \
                                       --output \$OUT_FILE
    """
}

process WRITE_MUTATIONS_FOR_DNDSCV {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'dndscv'

    script:
    """
    OUT_FILE='inputMutations-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    2c_write_mutations_for_dndscv.R --maf ${maf} \
                                    --cancer_subtype ${tumor_subtype} \
                                    --output \$OUT_FILE
    """
}


process WRITE_MUTATIONS_FOR_MUTPANNING {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*'), 
          path('patientsInv*')

    when:
    software == 'mutpanning'

    script:
    """
    OUT_FILE='inputMutations-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    OUT_INV='patientsInv-'${tumor_subtype}'-'$software'-'${params.target_genome_version}'.csv'
    2d_write_mutations_for_mutpanning.R --maf ${maf} \
                                        --cancer_subtype ${tumor_subtype} \
                                        --output \$OUT_FILE \
                                        --output_inventory \$OUT_INV
    """
}

process WRITE_MUTATIONS_FOR_NBR {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'nbr'

    script:
    """
    OUT_FILE='inputMutations-'$tumor_subtype'-'$software'-'$params.target_genome_version'.csv'
    2e_write_mutations_for_nbr.R --maf $maf \
                                 --cancer_subtype $tumor_subtype \
                                 --output \$OUT_FILE
    """
}

process WRITE_MUTATIONS_FOR_ONCODRIVEFML {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), val(software), path(maf)

    output:
    tuple val(tumor_subtype), path('inputMutations*')

    when:
    software == 'oncodrivefml'

    script:
    """
    OUT_FILE='inputMutations-'$tumor_subtype'-'$software'-'$params.target_genome_version'.csv'
    2f_write_mutations_for_oncodrivefml.R --maf $maf \
                                          --cancer_subtype $tumor_subtype \
                                          --output \$OUT_FILE
    """
}