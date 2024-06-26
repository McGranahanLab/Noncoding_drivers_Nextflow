process CHASMplus {
    tag "$tumor_subtype-$gr_id"

    input:
    tuple val(tumor_subtype), val(gr_id), val(software), path(mutations), 
          val(annotator)

    output:
    path "${software}Results-${tumor_subtype}-${gr_id}-${params.target_genome_version}.csv", emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    when:
    software == 'chasmplus'

    script:
    """
    MSG_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.out'
    ERR_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.err'

    # a unique ID to use in all further commands
    RUN_CODE=$software'Results-'$tumor_subtype'-'$gr_id'-'$params.target_genome_version
    OUT_FILE=\$RUN_CODE'.csv'

    oc run $mutations -a $annotator -t text -l $params.target_genome_version \
       --mp ${task.cpus}

    # move error and log files for their final location
    mv $mutations'.err' \$ERR_FILE
    mv $mutations'.log' \$MSG_FILE
    mv $mutations'.sqlite' \$RUN_CODE'.sqlite'

    # unfortunately, chasm plus puts results in the same folder as input files
    # It also collates 4 tables together, so we need to disentangle them
    chasmOutTsv=$mutations'.tsv'
    reportLevels=(\$(grep '^#Report level' \$chasmOutTsv | sed 's/.* //g'))
    tabsStartIdx=(\$(grep -n '^#Report level' \$chasmOutTsv | sed 's/:.*//g'))
    tabsStartIdx+=(`wc -l \$chasmOutTsv | sed 's/ .*//g'`)
    iterateTo="\${#tabsStartIdx[@]}"
    iterateTo=\$((iterateTo - 2)) # -2 due to counting from 0
    for tsi in `seq 0 \$iterateTo`; do
        startRead="\${tabsStartIdx[\$tsi]}"
        startRead=\$((startRead + 3))
        endRead=\$((tsi + 1))
        endRead="\${tabsStartIdx[\$endRead]}"

        readRegion=`echo \$startRead','\$endRead'p'`
        reportLvl="\${reportLevels[\$tsi]}"
        sed -n "\$readRegion" \$chasmOutTsv | grep -v '^#' > \$RUN_CODE'_'\$reportLvl
    done

    awk -F'\t' '\$17!=""' \$RUN_CODE'_variant' > \$OUT_FILE
    
    """
}

process DIGDRIVER {
    tag "$tumor_subtype-$gr_id"

    input:
    tuple val(tumor_subtype), val(gr_id), val(software), path(regions), 
          path(rda_ncbi), path(rda_ucsc), 
          path(rda_ncbi_restrictToCovs), path(rda_ucsc_restrictToCovs), 
          path(mutations), path(digdriver_model), 
          path(digdriver_elements), path(target_genome_fasta)

    output:
    path "${software}Results-${tumor_subtype}-${gr_id}-${params.target_genome_version}.csv", emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    when:
    software == 'digdriver'

    script:
    """
    MSG_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.out'
    ERR_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.err'

    # a unique ID to use in all further commands
    RUN_CODE=$tumor_subtype'-'$gr_id'-'$params.target_genome_version
    OUT_FILE=$software"Results-"\$RUN_CODE'.csv'

    # copy model as it will be modified during the run
    DIG_MODEL='digdriver-model-'\$RUN_CODE'.h5'
    cp $digdriver_model \$DIG_MODEL
    # make it writable
    chmod 700 \$DIG_MODEL

    # path to element_data.h5 referred in DIGdriver tutorials 
    DIG_ELEMENTS='digdriver-elements-'\$RUN_CODE'.h5'
    cp $digdriver_elements \$DIG_ELEMENTS
    chmod 700 \$DIG_ELEMENTS

    # STEP 1: annotate the mutation file 
    ANNOT_VARIANTS=$mutations'-DIGannotated.csv'
    touch \$ANNOT_VARIANTS
    DigPreprocess.py annotMutationFile --n-procs ${task.cpus} \
        $mutations $rda_ncbi $target_genome_fasta \$ANNOT_VARIANTS \
        1>\$MSG_FILE 2>\$ERR_FILE

    # STEP 2: preprocess the nucleotide contexts of the annotations
    DigPreprocess.py preprocess_element_model \$DIG_ELEMENTS \$DIG_MODEL \
        $target_genome_fasta \$RUN_CODE --f-bed $regions \
        1>>\$MSG_FILE 2>>\$ERR_FILE

    # STEP 3: pretrain a neutral mutation model 
    DigPretrain.py elementModel \$DIG_MODEL \$DIG_ELEMENTS \$RUN_CODE \
        1>>\$MSG_FILE 2>>\$ERR_FILE

    # geneDriver is not used even for CDS, because it can handle only SNP
    DigDriver.py elementDriver \$ANNOT_VARIANTS \$DIG_MODEL \$RUN_CODE \
        --f-bed $regions --outdir '.' --outpfx \$RUN_CODE \
        1>>\$MSG_FILE 2>>\$ERR_FILE
    mv \$RUN_CODE'.results.txt' \$OUT_FILE
    """
}

process DNDSCV {
    tag "$tumor_subtype-$gr_id"

    input:
    tuple val(tumor_subtype), val(gr_id), val(software), 
          path(rda_ncbi), path(rda_ucsc), 
          path(rda_ncbi_restrictToCovs), path(rda_ucsc_restrictToCovs),
          path(mutations)

    output:
    path "${software}Results-${tumor_subtype}-${gr_id}-${params.target_genome_version}.csv", emit: csv
    path "*global.csv", emit: global
    path "*Rds", emit: rds
    tuple path('*.out'), path('*.err'), emit: logs

    when:
    software == 'dndscv'

    script:
    """
    OUT_FILE=$software"Results-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.csv'
    OUT_GLOBAL_FILE=$software"Results-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'_global.csv'
    OUT_RDS=$software"Results-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.Rds'
    MSG_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.out'
    ERR_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.err'

    rda=$rda_ncbi_restrictToCovs
    mut_format=`head -2 $mutations | tail -n 1 | cut -f2`
    if [[ \$mut_format == chr* ]];
    then
        rda=$rda_ucsc_restrictToCovs
    fi

    run_dndscv.R --with_covariates T --genomic_regions \$rda \
                 --variants $mutations --computeCI T --output \$OUT_FILE \
                 --outputGlobal \$OUT_GLOBAL_FILE --outputRds \$OUT_RDS \
                 1>\$MSG_FILE 2>\$ERR_FILE
    """
}

process MUTPANNING {
    tag "$tumor_subtype-$gr_id"

    input:
    tuple val(tumor_subtype), val(gr_id), val(software), path(mutations),
          path(mutpan_inv)

    output:
    path "${software}Results-${tumor_subtype}-${gr_id}-${params.target_genome_version}.csv", emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    when:
    software == 'mutpanning'

    script:
    """
    OUT_FILE=$software"Results-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.csv'
    MSG_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.out'
    ERR_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.err'

    mpPath=/bin/MutPanningV2/commons-math3-3.6.1.jar:/bin/MutPanningV2/jdistlib-0.4.5-bin.jar:/bin/MutPanningV2
    java -Xmx${params.mutpanning_java_memory} -classpath \$mpPath MutPanning "." \
        $mutations $mutpan_inv "/bin/MutPanningV2/Hg19/" \
        1>\$MSG_FILE 2>\$ERR_FILE

    nResultFiles=\$(ls -l SignificanceRaw | wc -l)
    if [ \$nResultFiles -ne 1 ]; then
        cp SignificanceRaw/Significancecustom.txt \$OUT_FILE
    else 
        echo 'MutPanning run on '$tumor_subtype" "$gr_id' did not yield results.' >> \$MSG_FILE
        echo 'It can happen if your cohort is too small. Will create an empty file.' >> \$MSG_FILE
        echo -e "Name\tTargetSize\tTargetSizeSyn\tCount\tCountSyn\tSignificanceSyn\tFDRSyn\tSignificance\tFDR" > \$OUT_FILE
    fi
    """
}

process NBR {
    tag "$tumor_subtype-$gr_id"

    input:
    tuple val(tumor_subtype), val(gr_id), val(software), path(regions), 
          path(mutations), path(target_genome_fasta), path(nbr_neutral_bins),
          path(nbr_neutral_trinucfreq), path(nbr_driver_regs)
    
    output:
    path "${software}Results-${tumor_subtype}-${gr_id}-${params.target_genome_version}.csv", emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    when:
    software == 'nbr'

    script:
    """
    OUT_FILE=$software"Results-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.csv'
    MSG_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.out'
    ERR_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.err'

    NBR_create_intervals.R --region_bed $regions \
                           --genomeFile $target_genome_fasta \
                           1>\$MSG_FILE 2>\$ERR_FILE

    regionsFile=`echo $regions| sed 's/.csv\$/.csv.regions/g'`
    triNuclregionsFile=`echo $regions | sed 's/.csv\$/.csv.txt/g'`

    NBR.R --mutations_file $mutations --target_regions_path \$regionsFile \
          --target_regions_trinucfreqs_path \$triNuclregionsFile \
          --genomeFile $target_genome_fasta \
          --max_num_muts_perRegion_perSample 2 --unique_indelsites_FLAG 1 \
          --gr_drivers $nbr_driver_regs \
          --regions_neutralbins_file $nbr_neutral_bins \
          --trinucfreq_neutralbins_file $nbr_neutral_trinucfreq \
          --out_prefix \$OUT_FILE \
          1>>\$MSG_FILE 2>>\$ERR_FILE

    mv \$OUT_FILE"-Selection_output.txt" \$OUT_FILE 
    OUT_FILE_BASE=`echo \$OUT_FILE | sed 's/.csv//g'`
    # mv \$OUT_FILE"-global_mle_subs.Rds" \$OUT_FILE_BASE"-global_mle_subs.Rds" 
    # mv \$OUT_FILE"-globalRates.csv" \$OUT_FILE_BASE"-globalRates.csv"
    """
}

process ONCODRIVEFML {
    tag "$tumor_subtype-$gr_id"

    input:
    tuple val(tumor_subtype), val(gr_id), val(software), path(regions),
          path(mutations), path(oncodrivefml_config)

    output:
    path "${software}Results-${tumor_subtype}-${gr_id}-${params.target_genome_version}.csv", emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    when:
    software == 'oncodrivefml'

    script:
    """
    OUT_FILE=$software"Results-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.csv'
    MSG_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.out'
    ERR_FILE=$software"-"$tumor_subtype"-"$gr_id'-'$params.target_genome_version'.err'

    oncodrivefml -i $mutations -e $regions --sequencing wgs \
                 --configuration $oncodrivefml_config \
                 --output '.'  --cores ${task.cpus} \
                 1>\$MSG_FILE 2>\$ERR_FILE
    mv "inputMutations-"$tumor_subtype"-oncodrivefml-"$params.target_genome_version"-oncodrivefml.tsv" \$OUT_FILE
    """
}


process RUN_DISCOVER {
    tag "$tumor_subtype"

    input:
    tuple val(tumor_subtype), path(muts_to_gr), path(subtype_drivers),
          path(analysis_inventory_path), path(patients_inventory_path),
          val(software), path(subtype_specificity)

    output:
    tuple val(tumor_subtype), path("discoverResults-${tumor_subtype}--${params.target_genome_version}.csv"), emit: csv
    tuple path('*.out'), path('*.err'), emit: logs

    script:
    """
    RUN_CODE=$tumor_subtype'--'$params.target_genome_version
    MSG_FILE="discover-"\$RUN_CODE'.out'
    ERR_FILE="discover-"\$RUN_CODE'.err'
    OUT_FILE="discoverResults-"\$RUN_CODE'.csv'

    15_cooccurrence_and_exclusivity.R \
            --cancer_subtype $tumor_subtype \
            --inventory_patients $patients_inventory_path \
            --inventory_analysis $analysis_inventory_path \
            --drivers $subtype_drivers --muts_to_gr $muts_to_gr \
            --synAcceptedClass $params.synAcceptedClass \
            --fold_splicesites_in_coding $params.fold_splicesites_in_coding \
            --min_patients_discover $params.min_patients_discover \
            --subtype_specificity $subtype_specificity \
            --specificity_mode $params.specificity_mode \
            --output \$OUT_FILE 1>\$MSG_FILE 2>\$ERR_FILE
    """
}
