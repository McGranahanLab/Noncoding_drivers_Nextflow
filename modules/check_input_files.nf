process check_inventories {
    input:
    path patients_inventory_path
    path analysis_inventory_path
    path blacklist_inventory_path
    path digdriverModels_inventory_path
    path chasmplusAnno_inventory_path

    output:
    stdout emit: inventories_pass

    script:
    """
    1_check_inventories.R --inventory_patients ${patients_inventory_path} \
                          --inventory_analysis ${analysis_inventory_path} \
                          --target_genome_version ${params.target_genome_version} \
                          --inventory_blacklisted ${blacklist_inventory_path} \
                          --inventory_digdriver ${digdriverModels_inventory_path} \
                          --inventory_chasmplus ${chasmplusAnno_inventory_path}
    """
}

process check_fasta_is_uscs {
    input:
    path fasta_file

    output:
    stdout emit: fasta_pass

    script:
    """
    grep '>' ${fasta_file} > chrs.txt
    nChr=`wc -l chrs.txt | sed 's/ chrs.txt//g'`
    nChrUCSC=`grep '>chr' chrs.txt | wc -l`
    if [[ \$nChr -eq \$nChrUCSC ]]
    then
        echo 'FASTA is UCSC'
    else
        echo ${fasta_file} is not in UCSC format
        exit 1
    fi
    """
}

process check_digdriver_files {
    input:
    tuple val(software), path(digdriverElements), path(nbrNeutralBins),
          path(nbrNeutralTrinucfreq), path(nbrDriverRegs),
          path(oncodrivefmlConfig)

    when:
    software == 'digdriver'

    output:
    stdout emit: digdriver_pass

    script:
    """
    if [ ! -f ${digdriverElements} ]
    then
        echo DIGDRIVER was requested, but ${digdriverElements} does not exist
        exit 1
    fi
    echo "DIGDRIVER files are OK"
    """
}

process check_nbr_files {
    input:
    tuple val(software), path(digdriverElements), path(nbrNeutralBins),
          path(nbrNeutralTrinucfreq), path(nbrDriverRegs),
          path(oncodrivefmlConfig)

    when:
    software == 'nbr'

    output:
    stdout emit: nbr_pass

    script:
    """
    if [ ! -f ${nbrNeutralBins} ]
    then
        echo NBR was requested, but ${nbrNeutralBins} does not exist
        exit 1
    fi
    if [ ! -f ${nbrNeutralTrinucfreq} ]
    then
        echo NBR was requested, but ${nbrNeutralTrinucfreq} does not exist
        exit 1
    fi
    if [ ! -f ${nbrDriverRegs} ]
    then
        echo NBR was requested, but ${nbrDriverRegs} does not exist
        exit 1
    fi

    echo "NBR files are OK"
    """
}

process check_oncodrivefml_files {
    input:
    tuple val(software), path(digdriverElements), path(nbrNeutralBins),
          path(nbrNeutralTrinucfreq), path(nbrDriverRegs),
          path(oncodrivefmlConfig)

    when:
    software == 'oncodrivefml'

    output:
    stdout emit: oncodrivefml_pass

    script:
    """
    if [ ! -f ${oncodrivefmlConfig} ]
    then
        echo OncodriveFML was requested, but ${oncodrivefmlConfig} does not exist
        exit 1
    fi
    echo "OncodriveFML files are OK"
    """
}
