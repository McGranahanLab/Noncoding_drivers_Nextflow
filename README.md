# Welcome

> [!IMPORTANT]
> THE DOCUMENTATION IS UNDER DEVELOPMENT. PLEASE BEAR WITH US

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

This [Nextflow](https://www.nextflow.io/) pipeline is designed for the *de novo* detection of coding and noncoding somatic driver genomic elements based on single nucleotide variations (SNVs) and small insertions and deletions (indels) in cancer patient cohorts. It currently integrates five advanced calling algorithms: [DIGDriver](https://github.com/maxwellsh/DIGDriver), [dNdScv](https://github.com/im3sanger/dndscv/tree/master), NBR, [MutPanning](https://www.genepattern.org/modules/docs/MutPanning#gsc.tab=0), and [OncodriveFML](https://bbglab.irbbarcelona.org/oncodrivefml/home). [DIGDriver](https://github.com/maxwellsh/DIGDriver), NBR, and OncodriveFML are capable of detecting both coding and noncoding driver genetic elements, whereas [dNdScv](https://github.com/im3sanger/dndscv/tree/master) and advanced calling algorithms: [DIGDriver](https://github.com/maxwellsh/DIGDriver), [dNdScv](https://github.com/im3sanger/dndscv/tree/master), NBR, [MutPanning](https://www.genepattern.org/modules/docs/MutPanning#gsc.tab=0) focus solely on detecting coding drivers. The source code for NBR was provided by Dr. [Inigo Martincorena](https://github.com/im3sanger).

> CHASMplus

Overall, the pipeline can be divided into 3 steps: 1) application of the *de novo* cancer driver detection software to a patient cohort(s) and region(s) of interest to obtain raw p-values 2) postprocessing of the cancer driver detection software output, i.e. combination of raw p-values using Brown or methods  3) plotting of the results.

It is highly recommended to include [dNdScv](https://github.com/im3sanger/dndscv/tree/master) for coding and NBR for noncoding regions as those software are essential in post-processing to estimate the percentage of driver mutations in the discovered driver genomic regions and to pinpoint individual driver mutations shall it be possible. As [dNdScv](https://github.com/im3sanger/dndscv/tree/master) and NBR share concepts behind, it is not recommended to run NBR on CDS at the same time as [dNdScv](https://github.com/im3sanger/dndscv/tree/master).

This documentation provides comprehensive instructions on setting up, configuring, and running the pipeline along with detailed descriptions of the outputs.

## Table of content

* [Software requirements](#software-requirements)
* [Supported genome versions](#supported-genome-versions)
* [Supported NGS types](#supported-ngs-types)
* [Inputs](#inputs)
  * [Genomic variants files (mutations)](#genomic-variants-files-mutations)
    * [Annovar](#annovar)
    * [MAF](#maf)
  * [Genomic regions of interest](#genomic-regions-of-interest)
    * [GTF](#gtf)
    * [BED](#bed)
  * [Mutations multiplicity](#mutations-multiplicity)
  * [Inventory tables](#inventory-tables)
    * [Patients inventory table](#patients-inventory-table)
    * [Analysis inventory table](#analysis-inventory-table)
      * [Example of analysis table entries for CDS](#example-of-analysis-table-entries-for-cds)
      * [Example of analysis table entries for splice sites](#example-of-analysis-table-entries-for-splice-sites)
      * [Example of analysis table entries for 5'UTRs](#example-of-analysis-table-entries-for-5utrs)
      * [Example of analysis table entries for 3'UTRs](#example-of-analysis-table-entries-for-3utrs)
      * [Example of analysis table entries for shortRNA](#example-of-analysis-table-entries-for-shortrna)
    * [Black or whitelisted regions inventory table](#black-or-whitelisted-regions-inventory-table)
    * [DIGDriver models inventory table](#digdriver-models-inventory-table)
    * [CHASMplus annotators inventory table](#chasmplus-annotators-inventory-table)
    * [Expression inventory table](#expression-inventory-table)
    * [Tier definition table (inventory)](#tier-definition-table-inventory)
* [Parameters](#parameters)
  * [General](#general)
    * [Target genome version](#target-genome-version)
    * [Inventories](#inventories)
    * [Output directory](#output-directory)
  * [CHASMplus - specific files](#chasmplus---specific-files)
  * [DIGDriver - specific files](#digdriver---specific-files)
  * [NBR - specific files](#nbr---specific-files)
  * [OncodriveFML - specific files](#oncodrivefml---specific-files)
  * [Containers](#containers)
  * [Mutations filtering parameters](#mutations-filtering-parameters)
  * [Genomic regions filtering parameters](#genomic-regions-filtering-parameters)
  * [Assignment of mutations to regions and mutation rate calculations parameters](#assignment-of-mutations-to-regions-and-mutation-rate-calculations-parameters)
  * [Postprocessing](#postprocessing)
    * [Filtering out olfactory genes](#filtering-out-olfactory-genes)
    * [Filtering out not expressed genes](#filtering-out-not-expressed-genes)
    * [Filtering out hypermutated genomic regions](#filtering-out-hypermutated-genomic-regions)
    * [Biotyping](#biotyping)
    * [CHASMplus](#chasmplus)
    * [Tumor subtype specificity](#tumor-subtype-specificity)
    * [Mutual exclusivity and co-occurence](#mutual-exclusivity-and-co-occurence)
  * [Plotting](#plotting)
* [Profiles](#profiles)
* [Pipeline's execution](#pipelines-execution)
* [Outputs](#outputs)

## Software requirements

* **Nextflow**: The pipeline is written in DSL2 and requires [Nextflow](https://www.nextflow.io/docs/latest/install.html) version 23.04.2 or higher.
* **Singularity**: All software used in the pipeline is containerized. Interactions with containers are executed via [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html). The pipeline has been tested with Singularity version 3.8.3.

## Supported genome versions

Ideally, all input files should be in `hg19` coordinates. However, if this is not the case, avoid performing the liftover as it is already implemented in the pipeline. This approach minimizes the potential inconsistencies introduced by the liftover procedure.

## Supported NGS types

NBR can not run on WES.

## Inputs

Two inputs essential for the *de-novo* detection of cancer driver genomic elements are genetic alterations (mutations) and genomic regions of interest (i.e., set of coordinates which define coding regions, promoter, 5'UTRs, *etc*). To ensure that a signal of positive selection to be detected from the data is not distorted by a lower ability to perform sequencing in some genomic regions, it is also recommended to provide coordinates of [black-/white- listed regions](#black-or-whitelisted-regions-inventory-table).

[Inventory tables](#inventory-tables) are used to tie different input files, input types and software which will be applied to them together.

### Genomic variants files (mutations)

The pipeline can handle genomic variants files of two formats: The annovar-like table and the MAF-like table. Each file should contain genomic alterations (SNVs and small indels) for one patient only. The sections below provide an example of the input tables of two types.

#### Annovar

The table below demonstrates an example of a genomic variant file in Annovar-like format.

| **chr** | **start** | **stop** | **ref** | **var** | **Gene.refGene** | **Func.refGene** | **ExonicFunc.refGene** | **GeneDetail.refGene** | **AAChange.refGene** | **t_depth** | **t_ref_count** | **t_alt_count** | **n_depth** | **n_ref_count** | **n_alt_count** |
|:-------:|:---------:|:--------:|:-------:|:-------:|:----------------:|:----------------:|:----------------------:|:----------------------:|:--------------------:|:------------------:|:---------------:|:---------------:|:-----------:|:---------------:|:---------------:|
| 1 | 67705958 | 67705958 | G | A | IL23R | exonic | nonsynonymous SNV | . | IL23R:NM_144701:exon9:c.G1142A:p.R381Q | 25 | 15 | 9 | 42 | 42 | 0 |
| 2 | 234183368 | 234183368 | A | G | ATG16L1 | exonic | nonsynonymous SNV | . | ATG16L1:NM_198890:exon5:c.A409G:p.T137A,ATG16L1:NM_017974:exon8:c.A841G:p.T281A,ATG16L1:NM_001190266:exon9:c.A646G:p.T216A,ATG16L1:NM_001190267:exon9:c.A550G:p.T184A,ATG16L1:NM_030803:exon9:c.A898G:p.T300A |  57 | 29 | 27 | 114 | 114 | 0 |
| 16 | 50745926 | 50745926 | C | T | NOD2 | exonic | nonsynonymous SNV | . | NOD2:NM_001293557:exon3:c.C2023T:p.R675W,NOD2:NM_022162:exon4:c.C2104T:p.R702W | 38 | 30 | 8 | 38 | 38 | 0 |
| 13 | 20797176 | 21105944 | 0 | - | CRYL1;GJB6 | exonic | frameshift deletion | . | GJB6:NM_001110220:wholegene,GJB6:NM_001110221:wholegene,GJB6:NM_006783:wholegene,GJB6:NM_001110219:wholegene,CRYL1:NM_015974:wholegene | 49 | 24 | 25 | 29 | 29 | 0 |
| 8 | 8887543 | 8887543 | A | T | ERI1 | exonic | stoploss | . | ERI1:NM_153332:exon7:c.A1049T:p.X350L | 32 | 25 | 7 | 42 | 42 | 0 |

where

* **chr** *[essential]*: a chromosome where genomic variant was detected
* **start** *[essential]*: a start position of a genomic variant
* **stop** *[essential]*: an end position of a genomic variant
* **ref** *[essential]*: reference allele
* **var** *[essential]*: alternative allele
* **Gene.refGene** *[essential]*: a name of a gene to which a mutation was mapped to
* **Func.refGene** *[essential]*:
* **ExonicFunc.refGene** *[essential]*:
* **GeneDetail.refGene** *[essential]*:
* **AAChange.refGene** *[essential]*: aminoacid change which a mutation had induced
* **t_depth** *[optional]*: depth of a tumour sample at this position
* **t_ref_count** *[optional]*: number of reads with reference allele at this position in tumour sample
* **t_alt_count** *[optional]*: number of reads with an alternative allele at this position in the tumour sample
* **n_depth** *[optional]*: depth of a normal sample at this position
* **n_ref_count** *[optional]*: number of reads with reference allele at this position in the normal sample
* **n_alt_count** *[optional]*:  number of reads with an alternative allele at this position in the normal sample

#### MAF

The table below demonstrates an example of a genomic variant file in MAF-like format.

| **tumour_Sample_Barcode** | **Chromosome** | **Start_Position** | **End_Position** | **Reference_Allele** | **tumour_Seq_Allele2** | **Gene** | **Variant_Classification** | **Amino_acids** | **t_depth** | **t_ref_count** | **t_alt_count** | **n_depth** | **n_ref_count** | **n_alt_count** |
|:------------------------:|:--------------:|:------------------:|:----------------:|:--------------------:|:---------------------:|:--------:|:--------------------------------:|:---------------:|:-----------:|:---------------:|:---------------:|:-----------:|:---------------:|:---------------:|
| participant_1 | 10 | 96828976 | 96828977 | C | A | CYP2C8 | Intron | . | 25 | 15 | 9 | 42 | 42 | 0 |
| participant_1 | 11 | 118343898 | 118343899 | C | T | KMT2A | Missense_Mutation | S/L | 57 | 29 | 27 | 114 | 114 | 0 |
| participant_1 | 12 | 8074198 | 8074199 | G | A | SLC2A3 | Silent | I | 38 | 30 | 8 | 38 | 38 | 0 |
| participant_1 | 13 | 23909421 | 23909422 | T | C | SACS | Missense_Mutation | H/R | 49 | 24 | 25 | 29 | 29 | 0 |
| participant_1 | 1 | 17720542 | 17720543 | C | A | PADI6 | RNA | . | 32 | 25 | 7 | 42 | 42 | 0 |

where

* **tumour_Sample_Barcode** *[essential]*: The unique ID of a patient, e.g., `participant_1`.
* **Chromosome** *[essential]*: a chromosome where genomic variant was detected
* **Start_Position** *[essential]*: a start position of a genomic variant
* **End_Position** *[essential]*: an end position of a genomic variant
* **Reference_Allele** *[essential]*: reference allele
* **tumour_Seq_Allele2** *[essential]*: alternative allele
* **Gene** *[essential]*: a name of a gene to which a mutation was mapped to
* **Variant_Classification**: mutation's impact on a genomic element to which mutation was mapped. One of following values `Frame_Shift_Del`, `Frame_Shift_Ins`, `In_Frame_Del`, `In_Frame_Ins`, `Missense_Mutation`, `Nonsense_Mutation`, `Silent`, `Translation_Start_Site`, `Nonstop_Mutation`, `De_novo_Start_InFrame`, `De_novo_Start_OutOfFrame`, `Unknown`, `3'UTR`, `5'UTR`, `3'Flank`, `5'Flank`, `IGR`, `Intron`, `RNA`, `Splice_Site`.
* **Amino_acids** *[essential]*: aminoacid change which a mutation had induced
* **t_depth** *[optional]*: depth of a tumour sample at this position
* **t_ref_count** *[optional]*: number of reads with reference allele at this position in tumour sample
* **t_alt_count** *[optional]*: number of reads with an alternative allele at this position in the tumour sample
* **n_depth** *[optional]*: depth of a normal sample at this position
* **n_ref_count** *[optional]*: number of reads with reference allele at this position in the normal sample
* **n_alt_count** *[optional]*:  number of reads with an alternative allele at this position in the normal sample

### Genomic regions of interest

Genomic regions of interest can be provided via files in
[`gtf`](https://www.ensembl.org/info/website/upload/gff.html) or
[`bed`](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. It is
customary that more well established regions of a genome with known
biological function, i.e. CDS, promoters, lncRNA, miRNA, are derived from `gtf`
files and custom genome regions, i.e. regions with yet unknown functionality or
not yet fully experimentally validated ones, such as enhancers, are supplied
via `bed` files.

#### GTF

> [!NOTE]
> Due to inability of some *de-novo* driver calling software to support
> different genome versions (see
> [supported genome versions](#supported-genome-versions) section) it is
> highly recommended to use `gtf` genome annotations files for `hg19`
> genome. Specifically, necessary files for
> [MutPanning](https://www.genepattern.org/modules/docs/MutPanning#gsc.tab=0)
> and [DIGDriver](https://github.com/maxwellsh/DIGDriver) execution are solely
> available on `hg19` coordinates. While it is possible to liftover genomic
> coordinates to `hg19` from the other genome versions, such procedure may
> lead to shattering of regions into smaller ones or disruption of regions'
> internal structure, i.e. codons of CDS.

Standard `gtf` files for `hg19` genome can be downloaded from the
[UCSC web browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/).
It is recommended to use `refGene`, `ncbiRefSeq` or `ensGene` versions of the
`gtf` as these files have both `gene_name` and `gene_id` annotation. However,
all of these file formats are lacking `transcript_biotype` field which is
**essential** for distinguishing between protein coding genes and other types
of genes, i.e. `lncRNA` with the internal structure annotated as `exons`. Here
is an example of `gtf` file without `transcript_biotype` field:

```text
chr1 refGene transcript 11874 14409 . + . gene_id "DDX11L1"; transcript_id "NR_046018";  gene_name "DDX11L1";
chr1 refGene exon 11874 12227 . + . gene_id "DDX11L1"; transcript_id "NR_046018"; exon_number "1"; exon_id "NR_046018.1"; gene_name "DDX11L1";
chr1 refGene exon 12613 12721 . + . gene_id "DDX11L1"; transcript_id "NR_046018"; exon_number "2"; exon_id "NR_046018.2"; gene_name "DDX11L1";
chr1 refGene exon 13221 14409 . + . gene_id "DDX11L1"; transcript_id "NR_046018"; exon_number "3"; exon_id "NR_046018.3"; gene_name "DDX11L1";
```

and with `transcript_biotype` field:

```text
chr1 refGene exon 11874 12227 . + . gene_id "DDX11L1"; transcript_id "NR_046018"; gene_name "DDX11L1"; exon_number "1"; transcript_biotype "misc_RNA"; transcript_id_m "NR_046018";
chr1 refGene transcript 11874 14409 . + . gene_id "DDX11L1"; transcript_id "NR_046018"; gene_name "DDX11L1"; transcript_biotype "misc_RNA"; transcript_id_m "NR_046018";
chr1 refGene exon 12613 12721 . + . gene_id "DDX11L1"; transcript_id "NR_046018"; gene_name "DDX11L1"; exon_number "2"; transcript_biotype "misc_RNA"; transcript_id_m "NR_046018";
chr1 refGene exon 13221 14409 . + . gene_id "DDX11L1"; transcript_id "NR_046018"; gene_name "DDX11L1"; exon_number "3"; transcript_biotype "misc_RNA"; transcript_id_m "NR_046018";
```

`transcript_biotype` field can be added by intersecting `refGene gtf` with
`ensembl gtf` file. The most recent `ensembl` annotation `gtf` for
`hg19`/`GRCh37.75.gtf.gz` can be downloaded
[here](https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/). A script
to perform the intersection is available in
[preprocessing_scripts](preprocessing_scripts/) folder.

#### BED

One of the most convenient ways to define custom genomic regions, i.e.
enhancers, is though the `bed`-like formatted file. Such file must be tab
separated and contain 7 columns, for example:

| **chr** | **start** | **end** | **strand** | **gene_id** | **gene_name** | **rCode** |
|:-------:|:---------:|:-------:|:----------:|:-----------:|:-------------:|:---------:|
|chr1|858255|858648|*|SAMD11|SAMD11|enhancer|
|chr1|893086|897162|*|KLHL17|KLHL17|enhancer|
|chr1|893086|897162|*|NOC2L|NOC2L|enhancer|
|chr1|901067|902970|*|PLEKHN1|PLEKHN1|enhancer|
|chr1|910650|915364|*|PERM1|PERM1|enhancer|

where

* **chr** *[essential]*: chromosome
* **start** *[essential]*: start of the region
* **end** *[essential]*: end of the region
* **strand** *[essential]*: strand on which regions is located, i.e. `+`, `-` or
`*` (no strand)
* **gene_id** *[essential]*: ID of the gene to which region belongs. It can be
any string.
* **gene_name** *[essential]*: gene name to which region belongs. It can be
any string.
* **rCode** *[essential]*: string, region biotype. The `rCode` can have any
string. **The values from the `rCode` column should be used as an entry for**
**`gr_code` column of the [analysis table](#analysis-inventory-table).**

### Mutations multiplicity

### Inventory tables

#### Patients inventory table

The patient inventory table is a comma-separated file that contains detailed information about all participants (patients) in the cohorts, including ID, tumour subtype, path to the mutation table, and other relevant data. This table is also used to define specific cohorts of participants for further analysis, i.e. adenocarcinomas, squamous cell carcinomas, pan-cancer, etc. As parallelisation is ensured by the pipeline architecture as well as by Nextflow itself, there is no need to have separate patient inventory tables for each tumour subtype.

The table below provides an example of a patient inventory table.

| **tumour_subtype** | **participant_id** | **participant_tumour_subtype** | **somatic_genome** | **somatic_path** | **mutmultiplicity_path** | **cn_segments_genome** | **cn_segments_path** | **cohort_name** |
|:-----------------:|:------------------:|:-----------------------------:|:------------------:|:----------------:|:------------------------:|:----------------------:|:--------------------:|:---------------:|
| Adenocarcinoma    | participant_1      | LUAD                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Adenocarcinoma    | participant_2      | LUAD                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Adenocarcinoma    | participant_3      | LUAD                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Squamous_cell    | participant_41      | LUSC                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Squamous_cell    | participant_42      | LUSC                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Squamous_cell    | participant_43      | LUSC                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Panlung    | participant_1      | LUAD                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Panlung    | participant_2      | LUAD                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Panlung    | participant_3      | LUAD                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Panlung    | participant_41      | LUSC                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Panlung    | participant_42      | LUSC                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |
| Panlung    | participant_43      | LUSC                | hg38               | full_path_to_file| full_path_to_file        | hg38                   | full_path_to_file    | GEL            |

where

* **tumour_subtype** *[essential]*: The name of the tumour cohort to be analyzed. For example, all patients with lung adenocarcinomas may be grouped in a cohort named `Adenocarcinoma`. This column must not contain values which are numbers, i.e. "adenocarcinama_1" is an allowed value, but "1234" is not. The values in this column must not contain a "-" character.
* **participant_id** *[essential]*: The unique ID of a patient, e.g., `participant_1`. Each value of the `participant_id` column must be linked to one and only one value of `participant_tumour_subtype` column.
* **participant_tumour_subtype** *[essential]*: The histological subtype of a tumour found in the corresponding participant, e.g., `LUAD` (**Lu**ng **Ad**enocarcinoma). This column must not contain values which are numbers,  i.e. "LUAD_1" is an allowed value, but "78" is not.
* **somatic_genome** *[essential]*: The version of the genome in which the coordinates of mutations are specified, e.g., `hg38`. This column must not contain values which are numbers, i.e. "hg38" is an allowed value, but "38" is not. The genome version must be the same for all files listed in the **somatic_path** column. Ideally, all input files should be in `hg19` coordinates. However, if this is not the case, avoid performing the liftover as it is already implemented in the pipeline. This approach minimizes the potential inconsistencies introduced by the liftover procedure.
* **somatic_path** *[essential]*: The full path on your system (computer/HPC/*etc*) to the file containing **somatic** genetic mutations detected in the tumour of the corresponding individual. For the format of this file, see the section above. The existence of the files in this column will be checked before pipeline execution.
* **mutmultiplicity_path** *[optional]*: The full path on your system (computer/HPC/*etc*) to the file containing mutation multiplicities for the genomic variants defined in the `somatic_path` column. For the format of this file, see the section above. The existence of the files in this column will be checked before pipeline execution (if the column is present in the table).
* **cn_segments_genome** *[optional]*: The version of the genome in which the coordinates of copy number segments are specified, e.g., `hg38`. This column must not contain values which are numbers,  i.e. "hg38" is an allowed value, but "38" is not. The genome version must be the same for all files listed in the **cn_segments_path** column.
* **cn_segments_path** *[optional]*: The full path on your system (computer/HPC/*etc*) to the file containing copy number segments of the tumour genome found in the corresponding participant. For the format of this file, see the section above. The existence of the files in this column will be checked before pipeline execution (if the column is present in the table).
* **cohort_name**  *[essential]*: . This column must not contain values which are numbers, i.e. "GEL_1" is an allowed value, but "009" is not.

Cancer cohorts that include multiple histological subtypes (for example, a `pan-lung` cancer cohort may include tumour samples from adenocarcinomas, squamous cell carcinomas, mesotheliomas, neuroendocrine carcinomas, *etc.*) can be defined as shown in lines 7-12 of the table. It is preferable that the **participant_tumour_subtype** column contains the actual histological subtype of the tumour, rather than a "pan-lung" substitute.

#### Analysis inventory table

The analysis inventory table is a comma-separated file that defines the genomic regions to be scanned for potential cancer driver elements. It also links tumour subtypes defined in the patient inventory table to these genomic regions of interest. Additionally, the table specifies the software to be used for scanning each genomic region. As parallelism is ensured by the pipeline architecture as well as by Nextflow itself, there is no need to have separate analysis inventory tables for each tumour subtype.

The table below provides an example of an analysis inventory table.

| **tumour_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|
| Adenocarcinoma    | dndscv       |  coding   |   CDS       | full_path_to_gtf |   0          |      0         | hg38          | NA             | NA | NA | NA | NA | NA | CRG;DAC;DUKE|

where

* **tumour_subtype** *[essential]*: The name of the tumour cohort to be analyzed. It should match one of the cohort names listed in the `tumour_subtype` column of the [patient inventory table](#patients-inventory-table). The values in this column must not contain a "-" character.
* **software** *[essential]*: The name of the software to be applied to the cohort listed in the `tumour_subtype` column. Permitted values are: `DIGDriver`, `dndscv`, `mutpanning`, `chasmplus`, `nbr`, and `oncodrivefml`.
* **gr_id** *[essential]*: The name of the genomic region(s) set(s) to be analyzed. This can be any user-defined string, such as "coding" for CDS. The values in this column must not contain a "-" character.
* **gr_code** *[essential]*: A code (string) defining the biotypes of component parts of the genomic region(s) set(s) to be analyzed. The biotypes of the component parts are defined via files listed in `gr_file` column. Permitted values are: `3primeUTR`, `5primeUTR`, `CDS`, `lincRNA`, `lincRNA_promoter`, `lincRNA_ss`, `miRNA`, `misc_RNA`, `promoter`, `rRNA`, `snoRNA`, `snRNA`, and `ss` (splice sites). Genomic regions annotated with the corresponding biological function (e.g., all 3'UTRs in the case of `3primeUTR` or long non-coding RNAs in the case of `lincRNA`) will be extracted from the files specified in `gr_file`, combined gene-wise (collapsed across transcripts in the case of UTRs), and included in the set defined in the `gr_id` column. Each set of genomic regions of interest (defined uniquely by its ID in the `gr_id` column) can comprise one or more `gr_code`s. For example, a `gr_id` named "coding_and_UTRs" can have `CDS`, `3primeUTR`, and `5primeUTR` in the `gr_code` column. Please see below for illustrative examples. For [dNdScv](https://github.com/im3sanger/dndscv/tree/master) and `MutPanning` only `CDS` is accepted in this column.
* **gr_file** *[essential]*: The full path on your system (computer/HPC/*etc*) to a GTF or BED file from which genomic elements of interest should be extracted. The BED file must have the following columns: `chr`, `start`, `end`, `strand`, `gene_id`, `gene_name`, and `rCode`, where the `rCode` column contains the biotype of the region as a string, matching one of the values in the `gr_code` column. The GTF file must have the following fields: `gene_name`, `gene_id`, `gene_type`, `gene_biotype`, `transcript_id`, `transcript_type`, `transcript_biotype`. Please refer to the section [genomic regions of interest](#genomic-regions-of-interest) for more details. The existence of the files listed in this column will be checked prior to the pipeline's execution.
* **gr_upstr** *[essential]*: The number of bases upstream of the genomic region of interest to be included. For example, 5'UTRs can be extended upstream by an additional 50bp.
* **gr_downstr** *[essential]*: The number of bases downstream of the genomic region of interest to be included.
* **gr_genome** *[essential]*: The version of the genome in which the genomic coordinates of the regions are presented in the file listed in the `gr_file` column, e.g., `hg38`. This column must not contain values which are numbers, i.e. "hg38" is an allowed value, but "38" is not.
* **blacklisted_codes** *[essential]*: ; separated
* **union_percentage**
* **intersect_percentage**

In many cases, it is necessary to exclude certain genomic regions from the regions of interest to ensure clarity in the analysis. For example, regions overlapping with the cumulative set of CDS coordinates are typically excluded from splice site regions to prevent contamination of the splice site signal with the signal from CDS. The columns `gr_excl_id`, `gr_excl_code`, `gr_excl_file`, `gr_excl_upstr`, `gr_excl_downstr`, and `gr_excl_genome` define these exclusion regions in the same manner as the columns described above define regions of interest.

##### Example of analysis table entries for CDS

| **tumour_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** | **union_percentage** | **intersect_percentage** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|:---------------------:|:---------------------:|
| Adenocarcinoma    | dndscv       |  coding   |   CDS       | full_path_to_gtf |   0          |      0         | hg38          | NA             | NA | NA | NA | NA | NA | CRG;DAC;DUKE| NA | NA|

Typically, no upstream or downstream extensions of the CDS regions are considered during the coding region analysis. Therefore, the values in the `gr_upstr` and `gr_downstr` columns are set to `0`. Additionally, no genomic regions are excluded from the set of coding regions, so all columns related to exclusion regions (`gr_excl_id`, `gr_excl_code`, `gr_excl_file`, `gr_excl_upstr`, `gr_excl_downstr`, and `gr_excl_genome`) are set to `NA`. Furthermore, the `union_percentage` and `intersect_percentage` columns are set to `NA` as neither overlapping nor intersecting the coding genomic regions is performed to avoid disrupting the codon gene architecture.

##### Example of analysis table entries for splice sites

| **tumour_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** | **union_percentage** | **intersect_percentage** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|:---------------------:|:---------------------:|
| Adenocarcinoma    | oncodrivefml       |  splice_sites   |   ss       | full_path_to_gtf |   20          |      6         | hg38          | coding             | CDS | full_path_to_gtf | 0 | 0 | hg38 | CRG;DAC;DUKE| 50 | 80|

Splice sites are defined as intronic regions extending `20`bp from the donor site and `6`bp from the acceptor site. Therefore, values of `gr_upstr` and `gr_downstr` columns were set to `20` and `6` respectively. In order to prevent any spillover of the signal from coding regions of the genome to splice sites and keep them strictly intronic, regions of the genome annotated both as splice sites and coding are excluded from the consideration under splice sites set. To request this, columns `gr_excl_id`, `gr_excl_code`, `gr_excl_file`, `gr_excl_upstr`, `gr_excl_downstr`, and `gr_excl_genome` are set `coding`, `CDS`, `full_path_to_gtf`, `0`, `0`, and `hg38` respectively which matches the definition of `coding` regions from the section above. Shall you wish to define different coding regions to be excluded from the splice sites it is also possible via changes in `gr_excl_id`, `gr_excl_code`, `gr_excl_file`, `gr_excl_upstr`, `gr_excl_downstr`, and `gr_excl_genome` columns. That change will not affect the definition of coding regions to consider outlined in section [Example of analysis table entries for CDS](#example-of-analysis-table-entries-for-cds). To prevent the scoring of the same genomic regions multiple times, a union of splice site regions will be taken if they overlap > 50% (column `union_percentage`) and an intersection if they overlap > 80% (column `intersect_percentage`).

> full_path_to_gtf or bed is ok too?

##### Example of analysis table entries for 5'UTRs

| **tumour_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** | **union_percentage** | **intersect_percentage** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|:---------------------:|:---------------------:|
| Adenocarcinoma    | oncodrivefml       |  5_prime_UTR   |   5primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | coding             | CDS | full_path_to_gtf | 0 | 0 | hg38 | CRG;DAC;DUKE| 50 | 80|
| Adenocarcinoma    | oncodrivefml       |  5_prime_UTR   |   5primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | splice_sites             | ss | full_path_to_gtf | 20 | 6 | hg38 | CRG;DAC;DUKE| 50 | 80|

The table above shows an example definition table for 5'UTRs. Bases overlapping coding and splice site regions are excluded from the set of 5'UTRs.

##### Example of analysis table entries for 3'UTRs

| **tumour_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** | **union_percentage** | **intersect_percentage** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|:---------------------:|:---------------------:|
| Adenocarcinoma    | oncodrivefml       |  3_prime_UTR   |   3primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | coding             | CDS | full_path_to_gtf | 0 | 0 | hg38 | CRG;DAC;DUKE| 50 | 80|
| Adenocarcinoma    | oncodrivefml       |  3_prime_UTR   |   3primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | splice_sites             | ss | full_path_to_gtf | 20 | 6 | hg38 | CRG;DAC;DUKE| 50 | 80|
| Adenocarcinoma    | oncodrivefml       |  3_prime_UTR   |   3primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | 5_prime_UTR             | 5primeUTR | full_path_to_gtf | 0 | 0 | hg38 | CRG;DAC;DUKE| 50 | 80|

The table above shows an example definition table for 3'UTRs. Bases overlapping coding, splice site and 5'UTR regions are excluded from the set of 3'UTRs.

##### Example of analysis table entries for shortRNA

#### Black or whitelisted regions inventory table

The black-/while- listed regions inventory table is a comma-separated file that
defines genomic regions of inclusion (white) and exclusion (black). Base pairs
constituting genomic regions of interest defined in the
[analysis inventory table](#analysis-inventory-table) above will be excluded
from the analysis if they overlap black-listed regions and will be included if
and only if they overlap white-listed regions. This inventory table is
**optional**. If it is not provided, all values in the column
`blacklisted_codes` of
[analysis inventory table](#analysis-inventory-table) should be set to `NA`.

The table below provides an example of black-/white- listed regions inventory
table. The examples show commonly used black-/white- lists in the genomics
studies: CRG alignability for 100mers, DAC blacklisted regions, and Duke
uniqueness. More information about these tracks can be found
[here](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeMapability&db=hg19).

|list_name | file_path                                          | file_genome | file_type |
|:--------:|:--------------------------------------------------:|:-----------:|:---------:|
|CRG       | wgEncodeCrgMapabilityAlign100mer.bigWig            | hg19        | white     |
|DAC       | wgEncodeDacMapabilityConsensusExcludable.bed       | hg19        | black     |
|DUKE      | wgEncodeDukeMapabilityUniqueness35bp_processed.bed | hg19        | white     |

where

* **list_name** *[essential]*: A name of the black-/white- list. This name will
be used in `blacklisted_codes` columns in the
[analysis inventory table](#analysis-inventory-table). Each list name must
appear once and only once in the table and must not be a number.
* **file_path** *[essential]*: The full path on your system
(computer/HPC/*etc*) to a
[`bigWig`](https://genome.ucsc.edu/goldenPath/help/bigWig.html) or
[`bed`](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file containing
black- or white-listed regions. Each file path should appear once and only
once in the table.
* **file_type** *[essential]*: A type of regions: `black` (bases overlapping
them will be excluded) or `white`(only bases overlapping them will be included)
* **file_genome** *[essential]*: The version of the genome in which the genomic
coordinates of the regions are presented in the file listed in the `file_path`
column, e.g., `hg19`. This column must not contain values which are numbers,
i.e. "hg19" is allowed value, but "19" is not.

> [!WARNING]
> Due to possible split of the genomic region into numerous smaller regions and
> the subsequent increase in the computation requirements (i.e. RAM) required
> for the pipeline, black- and white-listed regions must be provided based on
> the same genome version as the target genome set by `target_genome_version`
> parameter.

#### DIGDriver models inventory table

The DIGDriver inventory table is a comma-separated file that defines
relationship between the analysed tumour subtypes and models which will be used
for them during the DIGDriver run. The complete list of the available models can be
found on [DIGDriver data portal](https://cb.csail.mit.edu/cb/DIG/downloads/).

The table below provides an example of a DIGDriver models inventory table.

| tumour_subtype  | model_file                                                |
|:--------------:|:---------------------------------------------------------:|
| Adenocarcinoma | DIGDriver_models/Lung-AdenoCA_SNV_MNV_INDEL.Pretrained.h5 |
| Squamous_cell  | DIGDriver_models/Lung-SCC_SNV_MNV_INDEL.Pretrained.h5     |
| Panlung        | DIGDriver_models/Lung_tumours_SNV_MNV_INDEL.Pretrained.h5  |

where

* **tumour_subtype** *[essential]*: The name of the tumour cohort to be
analyzed. It should match one of the cohort names listed in the `tumour_subtype`
column of the [patient inventory table](#patients-inventory-table). All tumour
subtypes for which analysis with DIGDriver was requested in the
[analysis inventory table](#analysis-inventory-table) must have a DIGDriver
model assigned.
* **model_file** *[essential]*: A full path to one of
[DIGDriver models](https://cb.csail.mit.edu/cb/DIG/downloads/). The same models
could be used for the different tumour subtypes.

#### CHASMplus annotators inventory table

The CHASMplus annotators inventory table is a comma-separated file that defines
the relationship between the analysed tumour subtypes and annotators which will
be used for them during the CHASMplus run. The complete list of the available
annotators can be found
[here](https://chasmplus.readthedocs.io/en/latest/models.html).

The table below provides an example of a CHASMplus annotators inventory table.

| tumour_subtype | chasm_annotator |
|:-------------:|:---------------:|
| Adenocarcinoma| chasmplus_LUAD  |
| Squamous_cell | chasmplus_LUSC  |
| Panlung       | chasmplus       |

where

* **tumour_subtype** *[essential]*: The name of the tumour cohort to be
analyzed. It should match one of the cohort names listed in the `tumour_subtype`
column of the [patient inventory table](#patients-inventory-table). All tumour
subtypes for which analysis with CHASMplus was requested in the
[analysis inventory table](#analysis-inventory-table) must have a CHASMplus
model assigned.
* **chasm_annotator** *[essential]*: A name of a
[CHASMplus annotator](https://chasmplus.readthedocs.io/en/latest/models.html).
The same annotator could be used for the different tumour subtypes.

#### Expression inventory table

#### Tier definition table (inventory)

The tier definition (inventory) table is a comma-separated file that defines
tiers of the future detected driver genomic elements. Tiers permit to use
various combinations of the cut offs on individual (reported by a driver
discovery software) raw, individual multiple test corrected (FDR) and merged
(see available methods for p-value merging [below](#postprocessing)) p-values.
There is no preference of one tier other the other.

The table below provides an example of a tier definition (inventory) table:

| tier | indivRaw_cutoff | nIndivRawSoft_sign | indivFDR_cutoff | nIndivFDRsoft_sign | restrictToKnownCancer | mergedFDR_cutoff |
| ---- | --------------- | ------------------ | --------------- | ------------------ | --------------------- | ---------------- |
| 1    | 0.1             | 2                  | 0.01            | 2                  | F                     | 0.01             |
| 2    | 0.1             | 2                  | 1               | 0                  | T                     | 0.05             |
| 3    | 0.1             | 2                  | 1               | 0                  | F                     | 0.01             |

where

* **tier** *[essential]*: an ID/name for the tier
* **indivRaw_cutoff** *[essential]*: a cut off on raw p-value of an individual
driver detecting software. `indivRaw_cutoff` is a number > 0 and < 1.
* **nIndivRawSoft_sign** *[essential]*: number of individual driver detecting
software which should have a raw p-value < `indivRaw_cutoff`.
`nIndivRawSoft_sign` is a number > 0 and < 1.
* **indivFDR_cutoff** *[essential]*: a cut off on FRD-corrected p-value of an
individual driver detecting software. `indivFDR_cutoff` is a number > 0 and < 1.
* **nIndivFDRsoft_sign** *[essential]*: number of individual driver detecting
software which should have a FRD-corrected p-value < `indivFDR_cutoff`.
`nIndivFDRsoft_sign` is a number > 0 and < 1.
* **restrictToKnownCancer** *[essential]*: boolean, indicating wherether or not
genomic regions should be restricted to the regions associated to the known
cancer driver genes. The list of known cancer driver genomic elements is
supplied via [`known_cancer_genes`](#postprocessing) parameter.
* **mergedFDR_cutoff** *[essential]*: a cut off on FRD-corrected merged p-value.
`mergedFDR_cutoff` is a number > 0 and < 1. The method of p-value merging is
controlled via [`combine_p_method`](#postprocessing) parameter.

For example, a genomic element will be assigned a tier 1 if all three following
conditions are true:

* there are at least 2 driver discovery software (i.e. dNdScv and OncodriveFML)
which raw p-values are < 0.1
* there are 2 driver discovery software for which FDR-corrected p-values are
< 0.01
* merged (Brown method by default) FRD corrected p-value is < 0.01

The first condition helps to ameliorate the effect of individual softwares'
technical biases on the results.

## Parameters

All the parameters described below are set in
[nextflow.config file](nextflow.config). Parameters related to the reference
genome, such as path to the reference fasta (`target_genome_path`), chromosomal
lengths (`target_genome_chr_len`) and chain file for liftover (`chain`) are set
up in configuration files in [conf](conf/) folder.

### General

#### Target genome version

* `target_genome_version`: a genome version to which the data should be brought
to before the driver detection software is applied. Currently, only `hg19` is
accepted. Please see [supported genomes version](#supported-genome-versions)
section for the reasoning.

#### Inventories

* `patients_inventory`: a path to the inventory file, i.e.
`'data/inventory/inventory_patients.csv'` provides detailed information about
all participants (patients) in the cohort(s). See
[patients inventory table](#patients-inventory-table) section for more details.
* `analysis_inventory`: a path to the inventory file, i.e.
`'data/inventory/inventory_analysis.csv'`. The table links together cohorts of
tumour subtypes, genomic regions of interest and software to be applied to the
regions. See
[analysis inventory table](#analysis-inventory-table) section for more details.
* `blacklist_inventory`: a path to the inventory file, i.e.
`'data/inventory/inventory_blacklist.csv'` providing detailed information about
black-/white- listed genomic regions. See
[black or whitelisted regions inventory table](#black-or-whitelisted-regions-inventory-table)
section for more details. This inventory is optional. If no genomic regions are
black-/white- listed, set the value of this parameter to '', i.e.
`blacklist_inventory = ''`.

#### Output directory

* `outdir` path to the output directory where results will be stored, i.e.
`completed_runs/`.

### CHASMplus - specific files

If analysis of data with [CHASMplus](https://chasmplus.readthedocs.io/en/latest/)
is requested, an inventory linking together cohorts of tumour subtypes and
[CHASMplus](https://chasmplus.readthedocs.io/en/latest/) annotators should be
provided via the `chasmplus_annotators_inventory` parameter, i.e.
`chasmplus_annotators_inventory = 'data/inventory/inventory_chasmplus_annotator.csv'`.
If analysis with [CHASMplus](https://chasmplus.readthedocs.io/en/latest/) is
not requested, set this parameter to `''`, i.e.
`chasmplus_annotators_inventory = ''`.

### DIGDriver - specific files

If analysis of genomic regions with
[DIGDriver](https://github.com/maxwellsh/DIGDriver)
is requested, then a file matching
[DIGDriver](https://github.com/maxwellsh/DIGDriver)
models to the tumour subtypes under considereations and a `element_data.h5` file
needed for [DIGDriver](https://github.com/maxwellsh/DIGDriver) training.

* `DIGDriver_models_inventory`: a path to the inventory file, i.e.
`'data/inventory/inventory_DIGDriver_models.csv'`. The table links together
cohorts of tumour subtypes and
[DIGDriver](https://github.com/maxwellsh/DIGDriver)
models to be used for their analysis.
* `DIGDriver_elements`: a path to `element_data.h5` used internally by
[DIGDriver](https://github.com/maxwellsh/DIGDriver). The file can be downloaded
[here](https://cb.csail.mit.edu/cb/DIG/downloads//dig_data_files/).

If analysis with [DIGDriver](https://github.com/maxwellsh/DIGDriver) is not
requested, set these parameters to `''`, i.e. `DIGDriver_models_inventory = ''`
and `DIGDriver_elements = ''`.

### NBR - specific files

If analysis of genomic regions with NBR is requested, then three additional
files are needed to be provided. All of these files can be found in the
[NBR](data/assets/NBR.zip) folder of this GitHub repository.

* `nbr_regions_neutralbins_file`: a path to a file defining neutral regions,
i.e. `'data/assets/NBR/Neutral_regions_within_100kb_bins_hg19.txt'`.
* `nbr_trinucfreq_neutralbins_file`: a path to a file listing trinucleotide
content within 100kb bins, i.e.
`data/assets/NBR/Trinucfreqs_within_100kb_bins_hg19.txt`
* `nbr_driver_regs_file`: a path to a file listing genomic regions containing
known driver genomic elements, i.e.
`data/assets/NBR/GRanges_driver_regions_hg19.txt`

### OncodriveFML - specific files

If analysis of genomic regions with OncodriveFML is requested, then several
additional files are needed to be provided. First of all, a configuration file
containing all OncodriveFML - specific perameters. Example of such file can be
found in [conf/oncodrivefml_hg19.config](conf/oncodrivefml_hg19.config).
Descriptions of the parameters is available
[here](https://oncodrivefml.readthedocs.io/en/latest/configuration.html). To
let the pipeline know which file should be used as OncodriveFML configuration
file, set up `oncodrivefml_config` parameter as path to the provided file,
i.e. `conf/oncodrivefml_hg19.config`. If analysis of genomic regions with
OncodriveFML is not requested, `oncodrivefml_config` parameter can be set to
`''`.

> [!WARNING]
> OncodriveFML requires a large files (17Gb) with genomic scores to be
> executed. It is not possible unfortunately to put such large files in the
> container. Therefore, the files need to be pre-downloaded in your computing
> environment. The download can be completed by executing the
> [example OncodriveFML run](https://oncodrivefml.readthedocs.io/en/latest/includes.html#run-the-example).
> Upon the run completion all the needed files will be stored in your system
> in `~/bgdata` folder.

### Containers

The following set of parameters defines containers to be used during all steps
of the pipeline execution. All containers can be viewed at
[Docker hub](https://hub.docker.com/r/marialitovchenko/noncoding_driver_pipeline/tags).
Recipes for container re-creation can be found in
[container_recipes folder](container_recipes).

* `chasmplus_container`: a container to be used for CHASMplus execution.
Default: `marialitovchenko/noncoding_driver_pipeline:chasmplus`
* `DIGDriver_container`: a container to be used for
[DIGDriver](https://github.com/maxwellsh/DIGDriver) execution.
Default: `marialitovchenko/noncoding_driver_pipeline:DIGDriver`
* `mutpanning_container`: a container to be used for MutPanning execution.
Default: `marialitovchenko/noncoding_driver_pipeline:mutpanning`
* `oncodrivefml_container`: a container to be used for OncodriveFML execution.
Default: `marialitovchenko/noncoding_driver_pipeline:oncodrivefml`
* `r_container`: a container to be used for the creation of input files for all
software, as well as [dNdScv](https://github.com/im3sanger/dndscv/tree/master),
NBR and postprocessing execution.
Default: `marialitovchenko/noncoding_driver_pipeline:r_packages`

### Mutations filtering parameters

* `min_depth`: a minimal depth of coverage of a mutation in a tumour sample for
the mutation to be considered for the *de novo* cancer driver discovery.
Recommended value: `30`.
* `min_tumour_vac`: a minimal number of reads with the alternative (mutated)
allele (also known as variant allele count (VAC)) in the tumour sample for the
mutation to be considered for the *de novo* cancer driver discovery. Recommended
value: `10`.
* `min_tumour_vaf`: a minimal percentage from the total reads with the
alternative (mutated) allele (also known as variant allele fraction (VAF)) in
the tumour sample for the mutation to be considered for the *de novo* cancer
driver discovery. Recommended value: `5` (percent).
* `max_germline_vac`: a maximum number of reads with the alternative (mutated)
allele in the germline sample for the mutation to be considered for the
*de novo* cancer driver discovery. Recommended value: `5`.
* `max_germline_vaf`: a maximum percentage from the total reads with the
alternative (mutated) allele in the germline sample for the mutation to be
considered for the *de novo* cancer driver discovery. Recommended value: `1`
(percent).
* `max_n_vars`: a maximum number of SNVs and small indels discovered in a
participant's tumour so that the sample is *not* recognised as hypermutated. If
mutations' number exceeds `max_n_vars`, then a sample will be removed from the
analysis. Recommended value: `90000`.

### Genomic regions filtering parameters

* `ignore_strand`: boolean, indicating, whatever or not strand information
should be taken into account while removing regions to exclude from the target
regions. For example, if 3'UTR of a gene_A is located on '+' strand and
overlaps CDS of a gene_B located on the '-' strand, then bases of gene_A's
3'UTR will be removed if `ignore_strand` set to `F` and kept otherwide.
Default: `'T'`. Accepted values: `'T'` and `'F'`.
* `min_reg_len`: during removal of unwanted regions from target regions (i.e.
CDS from the set of 3'UTRs) some regions may become very small. The
`min_reg_len` parameter sets the smallest length of a region which still will
be considered. Default: `5`.

### Assignment of mutations to regions and mutation rate calculations parameters

In the course of pipeline execution three types of mutation rate based on
mutation mapping to genomic region can be computed:

1. **genomic region specific mutation rate** then number of mutations within
each specific genomic region, i.e. TP53 CDS, is computed and divided by the
length of the corresponding region,
2. **local mutation rate** computed via bining the genome on consequtive bins
of certain size (for example 50kb), counting mutations inside them and dividing
by the bin size
3. **synonumous mutation rate** computed as number of synonymous mutations in a
certain genomic region; the synonumous mutation rate is only computed for CDS
regions.

These three mutation rates can be later used in postprocessing for filtering
out hypermutated regions of a genome. The genomic region specific mutation rate
will always be computed. Calculations of the local mutation rate can be
switched off by setting `bin_len` parameter to `-1`. Similarly, calculations of
synonumous mutation rate are controlled via `calc_synonymous` parameter.

* `remove_synonymous_from_coding` boolean, indicating whether synonymous
mutations should be removed from calculations of genomic region specific
mutation rate in coding genomic regions. Default: `T`. Accepted values: `T` and
`F`.
* `bin_len` parameter sets the length of the bins (in bp) for the local
mutation rate computations. Set it to `-1` to switch off local mutation rate
calculations. Reccomended value: `50000` (50kb).
* `calc_synonymous` boolean, indicating whether mutation rate based on
synonymous mutations should be computed (CDS regions only). Reccomended value:
`T`. Accepted values: `T` and `F`.

During calculations of the mutation rates genomic variants will be mapped to
the genomic regions of interest. The mapping will be checked for consistency by
ensuring that coding mutations are indeed mapped to the coding regions of a
genome and noncoding to the noncoding ones respectively.

* `gene_name_synonyms` a path to file which contains gene name synonyms, i.e.
`'data/assets/hgnc_complete_set_processed.csv'`. Set this parameter to `''` if
no such file can be provided. This file is used to match mutation to the
corresponding genes even if different gene names synonyms were used. For
example, `RASA1` gene is also
[known](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000145715;r=5:87267883-87391931)
under `CM-AVM`, `GAP`, `P120`, `P120GAP`, `P120RASGAP` and `RASA` names and
mutations assigned to any of these synonums should be mapped to genomic regions
associated with `RASA1`. The file should have two columns: `idx` and
`gene_name`. Gene names synonyms should have the same `idx`, i.e.:
  | idx | gene_name  |
  |:---:|:----------:|
  | 1234| CM-AVM     |  
  | 1234| GAP        |
  | 1234| P120       |
  | 1234| P120GAP    |
  | 1234| P120RASGAP |
  | 1234| RASA       |
  | 1234| RASA1      |
  | 12  | TP53       |
  | 12  | TP53       |
  
  `idx` can be any number.
* `cdsAcceptedClass` variant classifications from MAF format which are
acceptable as annotation of coding variants. Should a mutation not annotated
with the variant classifications set in this parameter be found in the coding
genomic regions, it will be flagged and re-annotated. Reccomended value:
`Frame_Shift_Del Frame_Shift_Ins In_Frame_Del In_Frame_Ins Missense_Mutation Nonsense_Mutation Silent Translation_Start_Site Nonstop_Mutation De_novo_Start_InFrame De_novo_Start_OutOfFrame Unknown`.
It is higly reccomended to include `Unknown` in the `cdsAcceptedClass` to
handle MNPs.
* `synAcceptedClass` variant classifications from MAF format which are
acceptable as annotation of synonymous variants. Reccomended value: `Silent`.
* `ncAcceptedClass` variant classifications from MAF format which are
acceptable as annotation of noncoding variants. Should a mutation not annotated
with the variant classifications set in this parameter be found in the coding
genomic regions, it will be flagged and re-annotated. Reccomended value:
`3primeUTR 3primeFlank 5primeUTR 5primeFlank IGR Intron RNA Targeted_Region Splice_Site Unknown`

As mentioned above in the description of `cdsAcceptedClass` and
`ncAcceptedClass` parameters, should a mismatch between genomic region class
and mutation class be detected, a variant re-annotation will be conducted by
means of [VariantAnnotation](https://www.bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
package. It is very likely that the initial annotation of the genetic variants
was not performed by the VariantAnnotation package, and therefore the new and
the original annotations should be harmonized. Parameter
`varanno_conversion_table` serves exactly this purpose.

* `varanno_conversion_table` a path to file which contains table for
harmonizing annotations created by VariantAnnotation package with the original
ones. The file should contain three columns: `variantAnnotation_anno`,
`var_type` and `var_class`, where
  * `variantAnnotation_anno` denotes a term used by VariantAnnotation package,
   i.e. `synonymous`
  * `var_type` denotes variant structural type, and should be one of `SNP`,
  `DEL`, `INS` or `MNP`.
  * `var_class` denotes a term used by the tool used to produce original
  annotation, i.e. `Silent`.

  Please note that a combination of values from `variantAnnotation_anno` and
  `var_type` columns should be uniquely matched to a value of `var_class`
  column. For example, `synonymous` and `SNP` uniquely maps to `Silent` value
  of `var_class` column.

  | variantAnnotation_anno | var_type  | var_class|
  |:---:|:----------:|:---:|
  | synonymous | SNP | Silent |
  | synonymous | DEL | Unknown|
  | synonymous | INS | Unknown|
  | synonymous | MNP | Silent|
  | nonsynonymous | SNP | Missense_Mutation|
  | nonsynonymous | DEL | In_Frame_Del|
  | nonsynonymous | INS | Frame_Shift_Ins|
  | nonsynonymous | MNP | Missense_Mutation|
  | nonsense | SNP | Nonsense_Mutation|
  | nonsense | DEL | Frame_Shift_Del|
  | nonsense | INS | Frame_Shift_Ins|
  | nonsense | MNP | Nonsense_Mutation|
  | frameshift | SNP | Unknown|
  | frameshift | DEL | Frame_Shift_Del|
  | frameshift | INS | Frame_Shift_Ins|
  | frameshift | MNP | Unknown|
* `annotation_failed_code`: string which should be used if `VariantAnnotation`
package fails to re-annotate a variant. Reccomended value: `Unknown`.

> [!NOTE]
> While a re-annotation will be attempted, it is expected that the amount of
> cases needing such a procedure will be low.

### Postprocessing

* `known_cancer_genes` `'data/assets/cgc_knownCancerGenes.csv'`

* `rawP_cap`: cap on raw p-values. If a raw p-value will be less than `rawP_cap`
it will be replaces with `rawP_cap`. For example, if `rawP_cap` is `1e-8`, all
p values below `1e-8` (i.e. `1e-16`, `1e-15`, `0`) will be set to `1e-8`.
Reccomended value: `'1e-8'`.
* `combine_p_method`: method to use for combining p-values. Note: all available
methods for combining p-values will be computed. The `combine_p_method`
parameter will just choose primary one to be used for tier assignment.
Reccomended value: `'brown'`. Accepted values: `'brown'`, `'fisher'`,
`'stouffer'`, `'harmonic'`.
* `tier_inventory`: a path to
[tier definition table](#tier-definition-table-inventory).

#### Filtering out olfactory genes

* `olfactory_genes`: path to table containing information about known olfactory
genes. The file must have columns `gene_id` and `gene_name`. For example:
  | gene_id | gene_name       |
  |:-------:|:---------------:|
  | OR4G4P  | ENSG00000268020 |
  | OR4G11P | ENSG00000240361 |
  | OR4F5   | ENSG00000186092 |
  | OR4F29  | ENSG00000284733 |
  | OR4F16  | ENSG00000284662 |

  The `data/assets/olfactory_barnes_2020.csv` table shows a list of olfactory
  genes derived from
  [Barnes et al 2020](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6583-3).
  A genomic regions will be considered as olfactory if either `gene_id` or
  `gene_name` is found in the `olfactory_genes` table.
* `remove_olfactory`: boolean, indicating whether or not olfactory genes should
be removed from the further processing. Reccomended value: `T`. Put to `F` to
disable. Accepted values: `T` and `F`.

#### Filtering out not expressed genes

* `gtex_inventory` `'data/inventory/inventory_expression_gtex.csv'`
* `tcga_inventory` `'data/inventory/inventory_expression_tcga.csv'`
* `gtex_expression` `'data/assets/GTEx_expression.csv'`
* `tcga_expression` `'data/assets/TCGA_expression.csv'`

#### Filtering out hypermutated genomic regions

* `min_n_soft_cod`: minimal number of `de-novo` driver detecting software
successfully executed on a coding genomic region in order for that region to be
considered for consequitive p-value merging. Genomic regions with the number of
successfully executed software below `min_n_soft_cod` will be excluded from the
further processing. Reccomended value: `3`. Put to `0` to disable.
* `min_n_soft_noncod`: minimal number of `de-novo` driver detecting software
successfully executed on a noncoding genomic region in order for that region to
be considered for consequitive p-value merging. Genomic regions with the number
of successfully executed software below `min_n_soft_noncod` will be excluded
from the further processing. Reccomended value: `2`. Put to `0` to disable.
* `min_n_muts`: minimal number of somatic mutations detected in a genomic
region for it to be considered for further processing. Genomic regions with the
number of mutations below `min_n_muts` will be excluded from the further
processing. Reccomended value: `3`. Put to `0` to disable.
* `min_n_patients`: minimal number of patients with a somatic mutation in a
genomic region for it to be considered for further processing. Genomic regions
with the number of mutated patients below `min_n_patients` will be excluded
from the further processing. Reccomended value: `3`. Put to `0` to disable.
* `max_local_mut_rate_q`: maximum quantile of a local mutation rate. Genomic
regions for which local mutation rate falls into a quantile exceeding
`max_local_mut_rate_q` will be excluded from the further processing.
Reccomended value: `0.95`. Put to `1` to disable. Accepted values: from `0` to
`1`.
* `max_gr_mut_rate_q`: maximum quantile of a genomic region - specific
mutation rate. Genomic regions for which genomic region - specific mutation
rate falls into a quantile exceeding `max_gr_mut_rate_q` will be excluded
from the further processing. Reccomended value: `1`. Put to `1` to disable.
Accepted values: from `0` to `1`.
* `max_gr_syn_mut_rate_q`: maximum quantile of a synonymous mutation rate.
Genomic regions for which synonymous mutation rate falls into a quantile
exceeding `max_gr_syn_mut_rate_q` will be excluded from the further
processing. Reccomended value: `1`. Put to `1` to disable. Accepted values:
from `0` to `1`.
* `max_gr_len_q`: maximum quantile of a genomic region length. Genomic regions
for which total length falls into a quantile exceeding `max_gr_len_q` will be
excluded from the further processing. This parameter allows to exclude extra
long genes. Reccomended value: `0.99`. Put to `1` to disable. Accepted values:
from `0` to `1`.
* `remove_2_5bp_enrich`: boolean, indicating whether or not genomic regions
enriched in 2-5bp deletions should be removed. Recently
[Reijns et al](https://www.nature.com/articles/s41586-022-04403-y) demonstrated
that such enrichment is a sign of TOP1 transcription-induced mutagenesis. We
hypothesise that genes with high expression level would be especially affected
by TOP1 mutagenesis due to increased transcription. Reccomended value: `T`. Put
to `F` to disable. Accepted values: `T` and `F`.
* `padj_2_5bp_enrich`: a maximum andjusted for multiple testing p-value of 2-5bp
enrichment in genomic regions test. If an andjusted p-values for the test is
smaller than `padj_2_5bp_enrich`, a genomic region will be considered as
enriched in 2-5bp deletions and therefore removed from the further
consideration. Reccomended value: `'5e-2'` = 0.05. Put to `1` to disable.
Accepted values: `T` and `F`.

#### Biotyping

* `min_gapwidt`                   `1000`
Minimum length of a gap between lifted over regions that will prevent them from being merged together. Default = Inf (no reduce will be performed)

* `min_width`                     `1000`
Minimum width of lifted over genomic regions. Default = 0 (no regions will be filtered out).

* `amp`
Cut off (in log 2 scale) on copy number to determine amplification. Default: log2(4/2) = `1`

* `gain`
Cut off (in log 2 scale) on copy number to determine gain. The gain cut off should be < amplification cut off Default: log2(2.5/2) = `0.3219281`.
* `loss` Cut off (in log 2 scale) on copy number to determine loss. Default: log2(1.5/2) = `-0.4150375`

* `exclude_silent_from_biotyping` `'T'`

* `min_biotype_muts_patients` `5`
Minimal number of patients in which a mutation (SNV or small indel) in a genomic element should be found so that mutations would be considered for biotyping.

* `min_biotype_cna_patients`  `10`
Minimal number of patients in which a amplification/gain/loss in a genomic element should be found so that copy number would be considered for biotyping.

* `weak_tsg`                  `0.33`
* `tsg`                       `0.50`
* `weak_og`                   `0.33`
* `og`                        `0.50`

#### CHASMplus

* `chasm_score_min`           `0.5`
* `chasm_padj`                `0.05`
* `max_fp`                    `0.05`
* // extra exclude_silent_from_biotyping?
* `known_driver_mutations`    `'data/assets/CancerGenomeInterpreter_lung_hg19.tsv'`

#### Tumor subtype specificity

* `subtype_spec_pval`         `0.05`

#### Mutual exclusivity and co-occurence

* `fold_splicesites_in_coding`        `'T'`
* `min_patients_discover`             `10`
* `specificity_mode`                  `"specific"`
* `p_max_uniformal_subtypes_discover` `0.1`
* `p_adj_discover`                    `0.05`

### Plotting

* `allowed_filter_values`      ["PASS", "INDEL, 2-5bp"]
* `extra_studies`              ["data/assets/intogene_detectedCancerGenes.csv", "data/assets/mc3_detectedCancerGenes.csv", "data/assets/cgc_knownCancerGenes.csv"]
* `extra_studies_names`        ["intogen", "mc3", "CGC"]
* `extra_studies_tumoursubtype` ["LNET,LUAD,LUSC,NSCLC,SCLC", "LUAD,LUSC", "nsclc,sclc,lung"]
* `plot_output_type`          "pdf" // or "png"
* `visual_json`                "data/visual_parameters.json"

## Profiles

## Pipeline's execution

```bash
```

```bash
```

## Outputs
