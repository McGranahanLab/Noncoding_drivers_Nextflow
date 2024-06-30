# Welcome
**THE DOCUMENTATION IS UNDER DEVELOPMENT. PLEASE BEAR WITH US**

This [Nextflow](https://www.nextflow.io/) pipeline is designed for the *de novo* detection of coding and noncoding somatic driver genomic elements based on single nucleotide variations (SNVs) and small insertions and deletions (indels) in cancer patient cohorts. It currently integrates five advanced calling algorithms: [DIGdriver](https://github.com/maxwellsh/DIGDriver), [dNdScv](https://github.com/im3sanger/dndscv/tree/master), NBR, [MutPanning](https://www.genepattern.org/modules/docs/MutPanning#gsc.tab=0), and [OncodriveFML](https://bbglab.irbbarcelona.org/oncodrivefml/home). [DIGdriver](https://github.com/maxwellsh/DIGDriver), NBR, and OncodriveFML are capable of detecting both coding and noncoding driver genetic elements, whereas [dNdScv](https://github.com/im3sanger/dndscv/tree/master) and advanced calling algorithms: [DIGdriver](https://github.com/maxwellsh/DIGDriver), [dNdScv](https://github.com/im3sanger/dndscv/tree/master), NBR, [MutPanning](https://www.genepattern.org/modules/docs/MutPanning#gsc.tab=0) focus solely on detecting coding drivers. Source code for NBR was provided by Dr. [Inigo Martincorena](https://github.com/im3sanger). 

> CHASMplus

Overall, the pipeline can be divided into 3 steps: 1) application of the *de novo* cancer driver detection software to a patient cohort(s) and region(s) of interest to obtain raw p-values 2) postprocessing of the cancer driver detection software output, i.e. combination of raw p-values using Brown or methods  3) plotting of the results.

It is higly recommended to include [dNdScv](https://github.com/im3sanger/dndscv/tree/master) for coding and NBR for noncoding regions as those software are essential in post processing to estimate percentage of driver mutations in the discovered driver genomic regions and to pinpoint individual driver mutations shall it be possible. As [dNdScv](https://github.com/im3sanger/dndscv/tree/master) and NBR share concepts behind, it is not recommenced to run NBR on CDS at the same time as [dNdScv](https://github.com/im3sanger/dndscv/tree/master).

This documentation provides comprehensive instructions on setting up, configuring, and running the pipeline, along with detailed descriptions of the outputs.

# Table of content
- [Software requirements](#requirements)
- [Supported genome versions](#supported-genome-versions)
- [Inputs](#inputs)
  - [Genomic variants files (mutations)](#genomic-variants-files)
  - [Genomic regions of interest](#genomic-regions-of-interest)
  - [Mutations multiplicity](#mutations-multiplicity)
  - [Inventory tables](#inventory-tables)
    - [Patients inventory table](#patients-inventory-table)
    - [Analysis inventory table](#analysis-inventory-table)
      - [Example of analysis table entries for CDS](#)
      - [Example of analysis table entries for splice sites](#)
      - [Example of analysis table entries for 5'UTRs](#)
      - [Example of analysis table entries for 3'UTRs](#)
      - [Example of analysis table entries for shortRNA](#)
    - [Black or whitelisted regions inventory table](#black-or-whitelisted-regions-inventory-table)
    - [DIGdriver models inventory table](#digdriver-models-inventory-table)
    - [CHASMplus annotators inventory table](#CHASMplus annotators inventory table)
    - [Expression inventory table](#Expression inventory table)
- [Parameters](#Parameters)
  - [General: path to inventories]()
  - [Containers](#Containers)
  - [Mutations filtering parameters](#Mutations-filtering-parameters)
- [Profiles](#Profiles)
- [Pipeline's execution](#Pipeline-s-execution)
- [Outputs](#Outputs)

# Software requirements
* **Nextflow**: The pipeline is written in DSL2 and requires [Nextflow](https://www.nextflow.io/docs/latest/install.html) version 23.04.2 or higher.
* **Singularity**: All software used in the pipeline is containerized. Interactions with containers are executed via [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html). The pipeline has been tested with Singularity version 3.8.3.

# Supported genome versions
Ideally, all input files should be in `hg19` coordinates. However, if this is not the case, avoid performing the liftover as it is already implemented in the pipeline. This approach minimizes the potential inconsistencies introduced by the liftover procedure.

# Supported NGS types
NBR can not run on WES.

# Inputs
Two types of inputs which are absolutely essential for the *de-novo* detection of cancer driver genomic elements are genetic alterations (mutations) and genomic regions of interest (i.e. set of coordinates which defines coding regions, promoter, 5'UTRs, *etc*). To ensure that a signal of positive selection to be detected from the data is not distorted by lower ability to perform sequencing in some genomic regions, it is recommended to also provide coordinates of [black-/white- listed regions](#black-or-whitelisted-regions-inventory-table).

[Inventory tables](#inventory-tables) are used to tie different input files, input types and software which will be applied to them together.

## Genomic variants files (mutations)
The pipeline can handle genomic variants files of two formats: Annovar-like table and MAF-like table. Each file should contain genomic alterations (SNVs and small indels) for one patient only. The sections below provide an example of the input tables of two types. 

### Annovar
The table below demonstrates an example of genomic variant file in Annovar-like format.

| **chr** | **start** | **stop** | **ref** | **var** | **Gene.refGene** | **Func.refGene** | **ExonicFunc.refGene** | **GeneDetail.refGene** | **AAChange.refGene** | **t_depth** | **t_ref_count** | **t_alt_count** | **n_depth** | **n_ref_count** | **n_alt_count** |
|:-------:|:---------:|:--------:|:-------:|:-------:|:----------------:|:----------------:|:----------------------:|:-------------------------------:|:--------------------:|:--------------------:|:----------:|:---------------:|:---------------:|:------------------:|:---------------:|:---------------:|
| 1 | 67705958 | 67705958 | G | A | exonic | IL23R | IL23R:NM_144701:exon9:c.G1142A:p.R381Q | 0 | nonsynonymous SNV | 25 | 15 | 9 | 42 | 42 | 0 |
| 2 | 234183368 | 234183368 | A | G | exonic | ATG16L1 | ATG16L1:NM_198890:exon5:c.A409G:p.T137A,ATG16L1:NM_017974:exon8:c.A841G:p.T281A,ATG16L1:NM_001190266:exon9:c.A646G:p.T216A,ATG16L1:NM_001190267:exon9:c.A550G:p.T184A,ATG16L1:NM_030803:exon9:c.A898G:p.T300A | 0 | nonsynonymous SNV | 57 | 29 | 27 | 114 | 114 | 0 |
| 16 | 50745926 | 50745926 | C | T | exonic | NOD2 | NOD2:NM_001293557:exon3:c.C2023T:p.R675W,NOD2:NM_022162:exon4:c.C2104T:p.R702W | 0 | nonsynonymous SNV | 38 | 30 | 8 | 38 | 38 | 0 |
| 13 | 20797176 | 21105944 | 0 | -0 | exonic | CRYL1;GJB6 | GJB6:NM_001110220:wholegene,GJB6:NM_001110221:wholegene,GJB6:NM_006783:wholegene,GJB6:NM_001110219:wholegene,CRYL1:NM_015974:wholegene | 0 | frameshift deletion | 49 | 24 | 25 | 29 | 29 | 0 |
| 8 | 8887543 | 8887543 | A | T | exonic | ERI1 | ERI1:NM_153332:exon7:c.A1049T:p.X350L | 0 | stoploss | 32 | 25 | 7 | 42 | 42 | 0 |

where

- **chr** *[essential]*: a chromosome where genomic variant was detected
- **start** *[essential]*: a start position of a genomic variant
- **stop** *[essential]*: an end position of a genomic variant
- **ref** *[essential]*: reference allele
- **var** *[essential]*: alternative allele
- **Gene.refGene** *[essential]*: a name of a gene to which a mutation was mapped to
- **Func.refGene** *[essential]*: 
- **ExonicFunc.refGene** *[essential]*:
- **GeneDetail.refGene** *[essential]*:
- **AAChange.refGene** *[essential]*: aminoacid change which a mutation had induced
- **t_depth** *[optional]*: depth of a tumor sample at this position
- **t_ref_count** *[optional]*: number of reads with reference allele at this position in tumor sample
- **t_alt_count** *[optional]*: number of reads with alternative allele at this position in tumor sample
- **n_depth** *[optional]*: depth of a normal sample at this position
- **n_ref_count** *[optional]*: number of reads with reference allele at this position in normal sample
- **n_alt_count** *[optional]*:  number of reads with alternative allele at this position in normal sample
                    
### MAF
The table below demonstrates an example of genomic variant file in MAF-like format.

| **Tumor_Sample_Barcode** | **Chromosome** | **Start_Position** | **End_Position** | **Reference_Allele** | **Tumor_Seq_Allele2** | **Gene** | **Variant_Classification** | **Amino_acids** | **t_depth** | **t_ref_count** | **t_alt_count** | **n_depth** | **n_ref_count** | **n_alt_count** |
|:------------------------:|:--------------:|:------------------:|:----------------:|:--------------------:|:---------------------:|:--------:|:--------------------------------:|:---------------:|:-----------:|:---------------:|:---------------:|:-----------:|:---------------:|:---------------:|
| participant_1 | 10 | 96828976 | 96828977 | C | A | CYP2C8 | Intron | . | 25 | 15 | 9 | 42 | 42 | 0 |
| participant_1 | 11 | 118343898 | 118343899 | C | T | KMT2A | Missense_Mutation | S/L | 57 | 29 | 27 | 114 | 114 | 0 |
| participant_1 | 12 | 8074198 | 8074199 | G | A | SLC2A3 | Silent | I | 38 | 30 | 8 | 38 | 38 | 0 |
| participant_1	| 13 | 23909421 | 23909422 | T | C | SACS | Missense_Mutation | H/R | 49 | 24 | 25 | 29 | 29 | 0 |
| participant_1 | 1 | 17720542 | 17720543 | C | A | PADI6 | RNA | . | 32 | 25 | 7 | 42 | 42 | 0 | 

where

- **Tumor_Sample_Barcode** *[essential]*: The unique ID of a patient, e.g., `participant_1`.
- **Chromosome** *[essential]*: a chromosome where genomic variant was detected
- **Start_Position** *[essential]*: a start position of a genomic variant
- **End_Position** *[essential]*: an end position of a genomic variant
- **Reference_Allele** *[essential]*: reference allele
- **Tumor_Seq_Allele2** *[essential]*: alternative allele
- **Gene** *[essential]*: a name of a gene to which a mutation was mapped to
- **Variant_Classification**: mutation's impact on a genomic element to which mutation was mapped. One of following values `Frame_Shift_Del`, `Frame_Shift_Ins`, `In_Frame_Del`, `In_Frame_Ins`, `Missense_Mutation`, `Nonsense_Mutation`, `Silent`, `Translation_Start_Site`, `Nonstop_Mutation`, `De_novo_Start_InFrame`, `De_novo_Start_OutOfFrame`, `Unknown`, `3'UTR`, `5'UTR`, `3'Flank`, `5'Flank`, `IGR`, `Intron`, `RNA`, `Splice_Site`.
- **Amino_acids** *[essential]*: aminoacid change which a mutation had induced
- **t_depth** *[optional]*: depth of a tumor sample at this position
- **t_ref_count** *[optional]*: number of reads with reference allele at this position in tumor sample
- **t_alt_count** *[optional]*: number of reads with alternative allele at this position in tumor sample
- **n_depth** *[optional]*: depth of a normal sample at this position
- **n_ref_count** *[optional]*: number of reads with reference allele at this position in normal sample
- **n_alt_count** *[optional]*:  number of reads with alternative allele at this position in normal sample

## Genomic regions of interest
## Mutations multiplicity

## Inventory tables
### Patients inventory table
The patient inventory table is a comma-separated file that contains detailed information about all participants (patients) in the cohort, including ID, tumor subtype, path to the mutation table, and other relevant data. This table is also used to define specific cohorts of participants for further analysis, i.e. adenocarcinomas, squamous cell carcinomas, pan-cancer, etc. As parallelisation is ensured by the pipeline architecture as well as by Nextflow itself, there is no need to have separate patient inventory tables for each tumor subtype. 

The table below provides an example of a patient inventory table.

| **tumor_subtype** | **participant_id** | **participant_tumor_subtype** | **somatic_genome** | **somatic_path** | **mutmultiplicity_path** | **cn_segments_genome** | **cn_segments_path** | **cohort_name** |
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

- **tumor_subtype** *[essential]*: The name of the tumor cohort to be analyzed. For example, all patients with lung adenocarcinomas may be grouped in a cohort named `Adenocarcinoma`. This column must not contain values which are numbers, i.e. "adenocarcinama_1" is allowed value, but "1234" is not. The values in this column must not contain "-" character.
- **participant_id** *[essential]*: The unique ID of a patient, e.g., `participant_1`. Each value of `participant_id` column must be linked to one and only one value of `participant_tumor_subtype` column.
- **participant_tumor_subtype** *[essential]*: The histological subtype of a tumor found in the corresponding participant, e.g., `LUAD` (**Lu**ng **Ad**enocarcinoma). This column must not contain values which are numbers,  i.e. "LUAD_1" is allowed value, but "78" is not.
- **somatic_genome** *[essential]*: The version of the genome in which the coordinates of mutations are specified, e.g., `hg38`. This column must not contain values which are numbers, i.e. "hg38" is allowed value, but "38" is not. The genome version must be the same for all files listed in **somatic_path** column. Ideally, all input files should be in `hg19` coordinates. However, if this is not the case, avoid performing the liftover as it is already implemented in the pipeline. This approach minimizes the potential inconsistencies introduced by the liftover procedure.
- **somatic_path** *[essential]*: The full path on your system (computer/HPC/*etc*) to the file containing **somatic** genetic mutations detected in the tumor of the corresponding individual. For the format of this file, see the section above. Existence of the files in this column will be checked on prior to pipeline execution.
- **mutmultiplicity_path** *[optional]*: The full path on your system (computer/HPC/*etc*) to the file containing mutation multiplicities for the genomic variants defined in the `somatic_path` column. For the format of this file, see the section above. Existence of the files in this column will be checked on prior to pipeline execution (if column is present in the table). 
- **cn_segments_genome** *[optional]*: The version of the genome in which the coordinates of copy number segments are specified, e.g., `hg38`. This column must not contain values which are numbers,  i.e. "hg38" is allowed value, but "38" is not. The genome version must be the same for all files listed in **cn_segments_path** column.
- **cn_segments_path** *[optional]*: The full path on your system (computer/HPC/*etc*) to the file containing copy number segments of the tumor genome found in the corresponding participant. For the format of this file, see the section above. Existence of the files in this column will be checked on prior to pipeline execution (if column is present in the table).
- **cohort_name**  *[essential]*: . This column must not contain values which are numbers, i.e. "GEL_1" is allowed value, but "009" is not.

Cancer cohorts that include multiple histological subtypes (for example, a `pan-lung` cancer cohort may include tumor samples from adenocarcinomas, squamous cell carcinomas, mesotheliomas, neuroendocrine carcinomas, *etc.*) can be defined as shown in lines 7-12 of the table. It is preferable that the **participant_tumor_subtype** column contains the actual histological subtype of the tumor, rather than a "pan-lung" substitute.
 
### Analysis inventory table
The analysis inventory table is a comma-separated file that defines the genomic regions to be scanned for potential cancer driver elements. It also links tumor subtypes defined in the patient inventory table to these genomic regions of interest. Additionally, the table specifies the software to be used for scanning each genomic region. As parallelism is ensured by the pipeline architecture as well as by Nextflow itself, there is no need to have separate analysis inventory tables for each tumor subtype. 

The table below provides an example of an analysis inventory table.

| **tumor_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|
| Adenocarcinoma    | dndscv       |  coding   |   CDS       | full_path_to_gtf |   0          |      0         | hg38          | NA             | NA | NA | NA | NA | NA | CRG;DAC;DUKE|

where

- **tumor_subtype** *[essential]*: The name of the tumor cohort to be analyzed. It should match one of the cohort names listed in the `tumor_subtype` column of the [patient inventory table](#Patients-inventory-table). The values in this column must not contain "-" character.
- **software** *[essential]*: The name of the software to be applied to the cohort listed in the `tumor_subtype` column. Permitted values are: `digdriver`, `dndscv`, `mutpanning`, `chasmplus`, `nbr`, and `oncodrivefml`.
- **gr_id** *[essential]*: The name of the genomic region(s) set(s) to be analyzed. This can be any user-defined string, such as "coding" for CDS. The values in this column must not contain "-" character.
- **gr_code** *[essential]*: A code (string) defining the biotypes of component parts of the genomic region(s) set(s) to be analyzed. The biotypes of the component parts are defined via files listed in `gr_file` column. Permitted values are: `3primeUTR`, `5primeUTR`, `CDS`, `lincRNA`, `lincRNA_promoter`, `lincRNA_ss`, `miRNA`, `misc_RNA`, `promoter`, `rRNA`, `snoRNA`, `snRNA`, and `ss` (splice sites). Genomic regions annotated with the corresponding biological function (e.g., all 3'UTRs in the case of `3primeUTR` or long non-coding RNAs in the case of `lincRNA`) will be extracted from the files specified in `gr_file`, combined gene-wise (collapsed across transcripts in the case of UTRs), and included in the set defined in the `gr_id` column. Each set of genomic regions of interest (defined uniquely by its ID in the `gr_id` column) can comprise one or more `gr_code`s. For example, a `gr_id` named "coding_and_UTRs" can have `CDS`, `3primeUTR`, and `5primeUTR` in the `gr_code` column. Please see below for illustrative examples. For [dNdScv](https://github.com/im3sanger/dndscv/tree/master) and `MutPanning` only `CDS` is accepted in this column.
- **gr_file** *[essential]*: The full path on your system (computer/HPC/*etc*) to a GTF or BED file from which genomic elements of interest should be extracted. The BED file must have the following columns: `chr`, `start`, `end`, `strand`, `gene_id`, `gene_name`, `rCode`, where the `rCode` column contains the biotype of the region as a string, matching one of the values in the `gr_code` column. The GTF file must have the following fields: `gene_name`, `gene_id`, `gene_type`, `gene_biotype`, `transcript_id`, `transcript_type`, `transcript_biotype`. Please refer to the section [genomic regions of interest](#Genomic-regions-of-interest) for more details. Existence of the files listed in this column will be checked prior to pipeline's execution.
- **gr_upstr** *[essential]*: The number of bases upstream of the genomic region of interest to be included. For example, 5'UTRs can be extended upstream by an additional 50bp.
- **gr_downstr** *[essential]*: The number of bases downstream of the genomic region of interest to be included.
- **gr_genome** *[essential]*: The version of the genome in which the genomic coordinates of the regions are presented in the file listed in the `gr_file` column, e.g., `hg38`. This column must not contain values which are numbers, i.e. "hg38" is allowed value, but "38" is not.
- **blacklisted_codes** *[essential]*: ; separated
- **union_percentage**
- **intersect_percentage**

In many cases, it is necessary to exclude certain genomic regions from the regions of interest to ensure clarity in the analysis. For example, regions overlapping with the cumulative set of CDS coordinates are typically excluded from splice site regions to prevent contamination of the splice site signal with the signal from CDS. The columns `gr_excl_id`, `gr_excl_code`, `gr_excl_file`, `gr_excl_upstr`, `gr_excl_downstr`, and `gr_excl_genome` define these exclusion regions in the same manner as the columns described above define regions of interest.

#### Example of analysis table entries for CDS

| **tumor_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** | **union_percentage** | **intersect_percentage** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|:---------------------:|:---------------------:|
| Adenocarcinoma    | dndscv       |  coding   |   CDS       | full_path_to_gtf |   0          |      0         | hg38          | NA             | NA | NA | NA | NA | NA | CRG;DAC;DUKE| NA | NA|

Typically, no upstream or downstream extensions of the CDS regions are considered during the coding region analysis. Therefore, the values in the `gr_upstr` and `gr_downstr` columns are set to `0`. Additionally, no genomic regions are excluded from the set of coding regions, so all columns related to exclusion regions (`gr_excl_id`, `gr_excl_code`, `gr_excl_file`, `gr_excl_upstr`, `gr_excl_downstr`, and `gr_excl_genome`) are set to `NA`. Furthermore, the `union_percentage` and `intersect_percentage` columns are set to `NA` as neither overlapping nor intersecting the coding genomic regions is performed to avoid disrupting the codon gene architecture.

#### Example of analysis table entries for splice sites

| **tumor_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** | **union_percentage** | **intersect_percentage** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|:---------------------:|:---------------------:|
| Adenocarcinoma    | oncodrivefml       |  splice_sites   |   ss       | full_path_to_gtf |   20          |      6         | hg38          | coding             | CDS | full_path_to_gtf | 0 | 0 | hg38 | CRG;DAC;DUKE| 50 | 80|

Splice sites are defined as intronic regions extending `20`bp from the donor site and `6`bp from the acceptor site. Therefore, values of `gr_upstr` and `gr_downstr` columns were set to `20` and `6` respectively. In order to prevent any spill over of the signal from coding regions of the genome to splice sites and keep them strictly intronic, regions of the genome annotated both as splice sites and coding are excluded from the consideration under splice sites set. To request this, columns `gr_excl_id`, `gr_excl_code`, `gr_excl_file`, `gr_excl_upstr`, `gr_excl_downstr`, and `gr_excl_genome` are set `coding`, `CDS`, `full_path_to_gtf`, `0`, `0`, and `hg38` respectively which matches the definition of `coding` regions from the section above. Shall you wish to define differently coding regions to be excluded from the splice sites it is also possible via change in `gr_excl_id`, `gr_excl_code`, `gr_excl_file`, `gr_excl_upstr`, `gr_excl_downstr`, and `gr_excl_genome` columns. That change will not affect the definition of coding regions to consider outlined in section [Example of analysis table entries for CDS](#Example-of-analysis-table-entries-for-CDS). To prevent scoring of the same genomic regions multiple times, a union of splice site regions will be taken if they overlap > 50% (column `union_percentage`) and an intersection if they overlap > 80% (column `intersect_percentage`).

> full_path_to_gtf or bed is ok too?

#### Example of analysis table entries for 5'UTRs

| **tumor_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** | **union_percentage** | **intersect_percentage** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|:---------------------:|:---------------------:|
| Adenocarcinoma    | oncodrivefml       |  5_prime_UTR   |   5primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | coding             | CDS | full_path_to_gtf | 0 | 0 | hg38 | CRG;DAC;DUKE| 50 | 80|
| Adenocarcinoma    | oncodrivefml       |  5_prime_UTR   |   5primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | splice_sites             | ss | full_path_to_gtf | 20 | 6 | hg38 | CRG;DAC;DUKE| 50 | 80|

The table above shows an example definition table for 5'UTRs. Bases overlapping coding and splice site regions are excluded from the set of 5'UTRs. 

#### Example of analysis table entries for 3'UTRs

| **tumor_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** | **union_percentage** | **intersect_percentage** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|:---------------------:|:---------------------:|
| Adenocarcinoma    | oncodrivefml       |  3_prime_UTR   |   3primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | coding             | CDS | full_path_to_gtf | 0 | 0 | hg38 | CRG;DAC;DUKE| 50 | 80|
| Adenocarcinoma    | oncodrivefml       |  3_prime_UTR   |   3primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | splice_sites             | ss | full_path_to_gtf | 20 | 6 | hg38 | CRG;DAC;DUKE| 50 | 80|
| Adenocarcinoma    | oncodrivefml       |  3_prime_UTR   |   3primeUTR       | full_path_to_gtf |   0          |      0         | hg38          | 5_prime_UTR             | 5primeUTR | full_path_to_gtf | 0 | 0 | hg38 | CRG;DAC;DUKE| 50 | 80|

The table above shows an example definition table for 3'UTRs. Bases overlapping coding, splice site and 5'UTR regions are excluded from the set of 3'UTRs.

#### Example of analysis table entries for shortRNA

### Black or whitelisted regions inventory table
The black-/while- listed regions inventory table is a comma-separated file that
defines genomic regions of inclusion (white) and exclusion (black). Base pairs 
constituting genomic regions of interest defined in the 
[analysis inventory table](#Analysis-inventory-table) above will be excluded 
from the analysis if they overlap black-listed regions and will be included if
and only if the overlap white-listed regions. This inventory table is 
**optional**. If it is not provided, all values in the column 
`blacklisted_codes` of 
[analysis inventory table](#Analysis-inventory-table) should be set to `NA`.

The table below provides an example of black-/white- listed regions inventory
table. The examples shows commonly used black-/white- lists in the genomics 
studies: CRG alignability for 100mers, DAC blacklisted regions, and Duke 
uniqueness. More information about these tracks can be found 
[here](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeMapability&db=hg19).

|list_name | file_path                                          | file_genome | file_type |
|:--------:|:--------------------------------------------------:|:-----------:|:---------:|
|CRG       | wgEncodeCrgMapabilityAlign100mer.bigWig            | hg19        | white     |
|DAC       | wgEncodeDacMapabilityConsensusExcludable.bed       | hg19        | black     |
|DUKE      | wgEncodeDukeMapabilityUniqueness35bp_processed.bed | hg19        | white     |

where

- **list_name** *[essential]*: A name of the black-/white- list. This name will
be used in `blacklisted_codes` columns in the 
[analysis inventory table](#Analysis-inventory-table). Each list name must 
appear once and only once in the table and must not be a number.
- **file_path** *[essential]*: The full path on your system 
(computer/HPC/*etc*) to a 
[`bigWig`](https://genome.ucsc.edu/goldenPath/help/bigWig.html) or 
[`bed`](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file containing 
black- or white- listed regions. Each file path should appear once and only
once in the table.
- **file_type** *[essential]*: A type of regions: `black` (bases overlapping 
them will be excluded) or `white`(only bases overlapping them will be included)
- **file_genome** *[essential]*: The version of the genome in which the genomic
coordinates of the regions are presented in the file listed in the `file_path`
column, e.g., `hg19`. This column must not contain values which are numbers, 
i.e. "hg19" is allowed value, but "19" is not.

> Due possible split of genomic region into the numerous smaller regions and 
subsequent increase in the computation requeirements (i.e. RAM) requered for
the pipeline, black- and white- listed regions must be provided based on the
same genome version as the target genome set by `target_genome_version` 
parameter.

### DIGdriver models inventory table
The DIGdriver inventory table is a comma-separated file that defines 
relationship between the analysed tumor subtypes and models which will be used
for them during DIGdriver run. The complete list of the available models can be
found on [DIGdriver data portal](https://cb.csail.mit.edu/cb/DIG/downloads/).

The table below provides an example of a DIGdriver models inventory table.

| tumor_subtype  | model_file                                                |
|:--------------:|:---------------------------------------------------------:|
| Adenocarcinoma | DIGdriver_models/Lung-AdenoCA_SNV_MNV_INDEL.Pretrained.h5 |
| Squamous_cell  | DIGdriver_models/Lung-SCC_SNV_MNV_INDEL.Pretrained.h5     |
| Panlung        | DIGdriver_models/Lung_tumors_SNV_MNV_INDEL.Pretrained.h5  |

where

- **tumor_subtype** *[essential]*: The name of the tumor cohort to be analyzed.
It should match one of the cohort names listed in the `tumor_subtype` column of
the [patient inventory table](#Patients-inventory-table). All tumor subtypes 
for which analysis with DIGdriver was requested in the 
[analysis inventory table](#Analysis-inventory-table) must have a DIGdriver 
model assigned.
- **model_file** *[essential]*: A full path to one of [DIGdriver models](https://cb.csail.mit.edu/cb/DIG/downloads/)
Same models could be used for the different tumor subtypes.

### CHASMplus annotators inventory table
The CHASMplus annotators inventory table is a comma-separated file that defines
relationship between the analysed tumor subtypes and annotators which will be 
used for them during CHASMplus run. The complete list of the available annotators
can be found on [here](https://chasmplus.readthedocs.io/en/latest/models.html).

The table below provides an example of a CHASMplus annotators inventory table.

| tumor_subtype | chasm_annotator |
|:-------------:|:---------------:|
| Adenocarcinoma| chasmplus_LUAD  | 
| Squamous_cell | chasmplus_LUSC  | 
| Panlung       | chasmplus       |


where

- **tumor_subtype** *[essential]*: The name of the tumor cohort to be analyzed.
It should match one of the cohort names listed in the `tumor_subtype` column of
the [patient inventory table](#Patients-inventory-table). All tumor subtypes 
for which analysis with CHASMplus was requested in the 
[analysis inventory table](#Analysis-inventory-table) must have a CHASMplus 
model assigned.
- **chasm_annotator** *[essential]*: A name of a 
[CHASMplus annotator](https://chasmplus.readthedocs.io/en/latest/models.html).
Same annotator could be used for the different tumor subtypes.

### Expression inventory table

# Parameters
All the parameters described below are set in [nextflow.config file](nextflow.config).

### General

**Target genome version**
- `target_genome_version` 'hg19'
path to refernce genome fasta are in profiles

**Inventories**
- `patients_inventory`             'data/inventory/inventory_patients.csv' 
- `analysis_inventory`             'data/inventory/inventory_analysis.csv' 
- `blacklist_inventory`            'data/inventory/inventory_blacklist.csv' optional blacklist_inventory '' 

required if [DIGdriver](https://github.com/maxwellsh/DIGDriver) run is requested, otherwise, set to ''

- `digdriver_models_inventory`     'data/inventory/inventory_digdriver_models.csv' 
- `digdriver_elements` 'data/assets/DIGdriver/element_data.h5'

required if CHAMSMplus run is requested, otherwise, set to ''

- `chasmplus_annotators_inventory` 'data/inventory/inventory_chasmplus_annotator.csv' 

NBR 
- `nbr_regions_neutralbins_file`    'data/assets/NBR/Neutral_regions_within_100kb_bins_hg19.txt'
- `nbr_trinucfreq_neutralbins_file` 'data/assets/NBR/Trinucfreqs_within_100kb_bins_hg19.txt'
- `nbr_driver_regs_file`            'data/assets/NBR/GRanges_driver_regions_hg19.txt'

required if CHAMSMplus run is requested, otherwise, set to ''

- `oncodrivefml_config` 'conf/oncodrivefml_hg19.config'

- `outdir` 'completed_runs/2023_12_25/'

### Containers
The following set of parameters define containers to be used during all steps of the pipeline execution. All containers can be viewed at [Docker hub](https://hub.docker.com/r/marialitovchenko/noncoding_driver_pipeline/tags). Recipes for container re-creation can be found in [container_recipes folder](container_recipes).

- `chasmplus_container`: a container to be used for CHASMplus execution. Default: `marialitovchenko/noncoding_driver_pipeline:chasmplus`
- `digdriver_container`: a container to be used for [DIGdriver](https://github.com/maxwellsh/DIGDriver) execution. Default: `marialitovchenko/noncoding_driver_pipeline:digdriver`
- `mutpanning_container`: a container to be used for MutPanning execution. Default: `marialitovchenko/noncoding_driver_pipeline:mutpanning`
- `oncodrivefml_container`: a container to be used for OncodriveFML execution. Default: `marialitovchenko/noncoding_driver_pipeline:oncodrivefml`
- `r_container`: a container to be used for creation of input files for all software, as well as [dNdScv](https://github.com/im3sanger/dndscv/tree/master), NBR and postprocessing execution. Default: `marialitovchenko/noncoding_driver_pipeline:r_packages`

### Mutations filtering parameters

- `min_depth`: a minimal depth of coverage of a mutation in a tumor sample for the mutation to be considered for the *de novo* cancer driver discovery. Recommended value: `30`.
- `min_tumor_vac`: a minimal number of reads with the alternative (mutated) allele (also known as variant allele count (VAC)) in the tumor sample for the mutation to be considered for the *de novo* cancer driver discovery. Recommended value: `10`.
- `min_tumor_vaf`: a minimal percentage from the total reads with the alternative (mutated) allele (also known as variant allele fraction (VAF)) in the tumor sample for the mutation to be considered for the *de novo* cancer driver discovery. Recommended value: `5` (percent).
- `max_germline_vac`: a maximum number of reads with the alternative (mutated) allele in the germline sample for the mutation to be considered for the *de novo* cancer driver discovery. Recommended value: `5`.
- `max_germline_vaf`: a maximum percentage from the total reads with the alternative (mutated) allele in the germline sample for the mutation to be considered for the *de novo* cancer driver discovery. Recommended value: `1` (percent).
- `max_n_vars` `90000`
- `ignore_strand` `'T'`
- `min_reg_len` `5`
- `gene_name_synonyms` `'data/assets/hgnc_complete_set_processed.csv'`   optional
- `bin_len` `50000`  optional, set it to `-1` to switch off filtering based on local mutation rate
- `calc_synonymous` `'T'` optional, set it to `'F'` to switch off filtering based on synonymous mutation rate (CDS only)
- `remove_synonymous_from_coding` `'T'` put to `F` to disable
- `cdsAcceptedClass` `Frame_Shift_Del Frame_Shift_Ins In_Frame_Del In_Frame_Ins Missense_Mutation Nonsense_Mutation Silent Translation_Start_Site Nonstop_Mutation De_novo_Start_InFrame De_novo_Start_OutOfFrame Unknown`
- `synAcceptedClass` `Silent`
- `ncAcceptedClass`  `3primeUTR 3primeFlank 5primeUTR 5primeFlank IGR Intron RNA Targeted_Region Splice_Site Unknown`
- `varanno_conversion_table` `data/assets/variantAnnotation_to_annovar_conversion.txt` optional, set to '' to not use
- `annotation_failed_code` `Unknown`

### Postprocessing
- `known_cancer_genes` `'data/assets/cgc_knownCancerGenes.csv'`
- `olfactory_genes` `'data/assets/olfactory_barnes_2020.csv'`
- `rawP_cap` `'1e-8'`
- `gtex_inventory` `'data/inventory/inventory_expression_gtex.csv'`
- `tcga_inventory` `'data/inventory/inventory_expression_tcga.csv'`
- `gtex_expression` `'data/assets/GTEx_expression.csv'`
- `tcga_expression` `'data/assets/TCGA_expression.csv'`
- `combine_p_method` `'brown'`
- `tier_inventory` `'data/inventory/inventory_tier_definition.csv'`

- `min_n_soft_cod`          3    // put to 0 to disable
- `min_n_soft_noncod`       2    // put to 0 to disable
- `min_n_muts`              3    // put to 0 to disable
- `min_n_patients`          3    // put to 0 to disable 
- `max_local_mut_rate_q`    0.95 // put to 1 to disable
    //max_gr_mut_rate_q       0.99    // put to 1 to disable
    `max_gr_mut_rate_q`       1    // put to 1 to disable
    //max_gr_syn_mut_rate_q   0.99    // put to 1 to disable
    `max_gr_syn_mut_rate_q`   1    // put to 1 to disable
    `max_gr_len_q`            0.99 // put to 1 to disable
    `remove_olfactory`        'T'  // put to 'F' to disable
    `remove_2_5bp_enrich`     'T' // put to 'F' to disable
- `padj_2_5bp_enrich`       `'5e-2'`

- `min_gapwidt`                   1000
- `min_width`                     1000
- `amp`                           1
- `gain`                          0.3219281
- `loss`                          -0.4150375
- `exclude_silent_from_biotyping` 'T'

- `min_biotype_muts_patients` 5
- `min_biotype_cna_patients`  10
- `weak_tsg`                  0.33
- `tsg`                       0.50
- `weak_og`                   0.33
- `og`                        0.50

- `chasm_score_min`           0.5
- `chasm_padj`                0.05
- `max_fp`                    0.05
- // extra exclude_silent_from_biotyping?
- `known_driver_mutations`    `'data/assets/CancerGenomeInterpreter_lung_hg19.tsv'`

- `subtype_spec_pval`         `0.05` 

- `fold_splicesites_in_coding`        `'T'`
- `min_patients_discover`             `10`
- `specificity_mode`                  `"specific"`
- `p_max_uniformal_subtypes_discover` `0.1`
- `p_adj_discover`                    `0.05`
    
### Plotting
- allowed_filter_values      ["PASS", "INDEL, 2-5bp"]
- extra_studies              ["data/assets/intogene_detectedCancerGenes.csv", "data/assets/mc3_detectedCancerGenes.csv", "data/assets/cgc_knownCancerGenes.csv"]
- extra_studies_names        ["intogen", "mc3", "CGC"]
- extra_studies_tumorsubtype ["LNET,LUAD,LUSC,NSCLC,SCLC", "LUAD,LUSC", "nsclc,sclc,lung"]
- plot_output_type           "pdf" // or "png"
- visual_json                "data/visual_parameters.json"

# Profiles

# Pipeline's execution

# Outputs

