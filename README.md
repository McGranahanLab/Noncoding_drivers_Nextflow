## Welcome
**THE DOCUMENTATION IS UNDER DEVELOPMENT. PLEASE BEAR WITH US**

This [Nextflow](https://www.nextflow.io/) pipeline is designed for the de-novo detection of coding and noncoding somatic driver genomic elements in cancer patient cohorts. It currently integrates five advanced calling algorithms: [DIGdriver](https://github.com/maxwellsh/DIGDriver), [dNdScv](https://github.com/im3sanger/dndscv/tree/master), NBR, [MutPanning](https://www.genepattern.org/modules/docs/MutPanning#gsc.tab=0), and [OncodriveFML](https://bbglab.irbbarcelona.org/oncodrivefml/home). DIGdriver, NBR, and OncodriveFML are capable of detecting both coding and noncoding driver genetic elements, whereas dNdScv and MutPanning focus solely on detecting coding drivers. Source code for NBR was provided by Dr. 
[Inigo Martincorena](https://github.com/im3sanger). 

This documentation provides comprehensive instructions on setting up, configuring, and running the pipeline, along with detailed descriptions of the outputs.

## Requirements
* **Nextflow**: The pipeline is written in DSL2 and requires [Nextflow](https://www.nextflow.io/docs/latest/install.html) version 23.04.2 or higher.
* **Singularity**: All software used in the pipeline is containerized. Interactions with containers are executed via [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html). The pipeline has been tested with Singularity version 3.8.3.

## Input formats
### Inventory tables
#### Patients inventory table
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
 
#### Analysis inventory table
The analysis inventory table is a comma-separated file that defines the genomic regions to be scanned for potential cancer driver elements. It also links tumor subtypes defined in the patient inventory table to these genomic regions of interest. Additionally, the table specifies the software to be used for scanning each genomic region. As parallelisation is ensured by the pipeline architecture as well as by Nextflow itself, there is no need to have separate analysis inventory tables for each tumor subtype. 

The table below provides an example of an analysis inventory table.

| **tumor_subtype** | **software** | **gr_id** | **gr_code** | **gr_file**      | **gr_upstr** | **gr_downstr** | **gr_genome** | **gr_excl_id** | **gr_excl_code** | **gr_excl_file** | **gr_excl_upstr** | **gr_excl_downstr** | **gr_excl_genome** | **blacklisted_codes** |
|:-----------------:|:------------:|:---------:|:-----------:|:----------------:|:------------:|:--------------:|:-------------:|:---------------:|:----------------------------|:----------------:|:-----------------:|:-------------------:|:------------------:|:---------------------:|
| Adenocarcinoma    | dndscv       |  coding   |   CDS       | full_path_to_gtf |   0          |      0         | hg38          | NA             | NA | NA | NA | NA | NA | DUKE, DAC|

where

- **tumor_subtype** *[essential]*: The name of the tumor cohort to be analyzed. It should match one of the cohort names listed in the `tumor_subtype` column of the [patient inventory table](#Patients-inventory-table). The values in this column must not contain "-" character.
- **software** *[essential]*: The name of the software to be applied to the cohort listed in the `tumor_subtype` column. Permitted values are: `digdriver`, `dndscv`, `mutpanning`, `chasmplus`, `nbr`, and `oncodrivefml`.
- **gr_id** *[essential]*: The name of the genomic region(s) set(s) to be analyzed. This can be any user-defined string, such as "coding" for CDS. The values in this column must not contain "-" character.
- **gr_code** *[essential]*: A code (string) defining the biotypes of component parts of the genomic region(s) set(s) to be analyzed. The biotypes of the component parts are defined via files listed in `gr_file` column. Permitted values are: `3primeUTR`, `5primeUTR`, `CDS`, `lincRNA`, `lincRNA_promoter`, `lincRNA_ss`, `miRNA`, `misc_RNA`, `promoter`, `rRNA`, `snoRNA`, `snRNA`, and `ss` (splice sites). Genomic regions annotated with the corresponding biological function (e.g., all 3'UTRs in the case of `3primeUTR` or long non-coding RNAs in the case of `lincRNA`) will be extracted from the files specified in `gr_file`, combined gene-wise (collapsed across transcripts in the case of UTRs), and included in the set defined in the `gr_id` column. Each set of genomic regions of interest (defined uniquely by its ID in the `gr_id` column) can comprise one or more `gr_code`s. For example, a `gr_id` named "coding_and_UTRs" can have `CDS`, `3primeUTR`, and `5primeUTR` in the `gr_code` column. Please see below for illustrative examples. For `dNdScv` and `MutPanning` only `CDS` is accepted in this column.
- **gr_file** *[essential]*: The path to a GTF or BED file from which genomic elements of interest should be extracted. The BED file must have the following columns: `chr`, `start`, `end`, `strand`, `gene_id`, `gene_name`, `rCode`, where the `rCode` column contains the biotype of the region as a string, matching one of the values in the `gr_code` column. The GTF file must have the following fields: `gene_name`, `gene_id`, `gene_type`, `gene_biotype`, `transcript_id`, `transcript_type`, `transcript_biotype`. Please refer to the section [genomic regions of interest](#Genomic-regions-of-interest) for more details. Existence of the files listed in this column will be checked prior to pipeline's execution.
- **gr_upstr** *[essential]*: The number of bases upstream of the genomic region of interest to be included. For example, 5'UTRs can be extended upstream by an additional 50bp.
- **gr_downstr** *[essential]*: The number of bases downstream of the genomic region of interest to be included.
- **gr_genome** *[essential]*: The version of the genome in which the genomic coordinates of the regions are presented in the file listed in the `gr_file` column, e.g., `hg38`.
- **blacklisted_codes** *[essential]*: 

In many cases, it is necessary to exclude certain genomic regions from the regions of interest to ensure clarity in the analysis. For example, regions overlapping with the cumulative set of CDS coordinates are typically excluded from splice site regions to prevent contamination of the splice site signal with the signal from CDS. The columns `gr_excl_id`, `gr_excl_code`, `gr_excl_file`, `gr_excl_upstr`, `gr_excl_downstr`, and `gr_excl_genome` define these exclusion regions in the same manner as the columns described above define regions of interest.

