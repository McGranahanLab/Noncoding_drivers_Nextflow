## Welcome
This [Nextflow](https://www.nextflow.io/) pipeline is designed for the de-novo detection of coding and noncoding somatic driver genomic elements in cancer patient cohorts. It currently integrates five advanced calling algorithms: [DIGdriver](https://github.com/maxwellsh/DIGDriver), [dNdScv](https://github.com/im3sanger/dndscv/tree/master), NBR, [MutPanning](https://www.genepattern.org/modules/docs/MutPanning#gsc.tab=0), and [OncodriveFML](https://bbglab.irbbarcelona.org/oncodrivefml/home). DIGdriver, NBR, and OncodriveFML are capable of detecting both coding and noncoding driver genetic elements, whereas dNdScv and MutPanning focus solely on detecting coding drivers. Source code for NBR was provided by Dr. 
[Inigo Martincorena](https://github.com/im3sanger). 

This documentation provides comprehensive instructions on setting up, configuring, and running the pipeline, along with detailed descriptions of the outputs.

## Input formats
### Inventory tables
#### Patients inventory table
The patient inventory table is a comma-separated file that contains detailed information about all participants (patients) in the cohort, including ID, tumor subtype, path to the mutation table, and other relevant data. This table is also used to define specific cohorts of participants for further analysis, i.e. adenocarcinomas, squamous cell carcinomas, pan-cancer, etc. 

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
- **somatic_genome** *[essential]*: The version of the genome in which the coordinates of mutations are specified, e.g., `hg38`. This column must not contain values which are numbers,  i.e. "hg38" is allowed value, but "38" is not. The genome version must be the same for all files listed in **somatic_path** column.
- **somatic_path** *[essential]*: The full path on your system (computer/HPC/*etc*) to the file containing **somatic** genetic mutations detected in the tumor of the corresponding individual. For the format of this file, see the section above. Existence of the files in this column will be checked on prior to pipeline execution.
- **mutmultiplicity_path** *[optional]*: The full path on your system (computer/HPC/*etc*) to the file containing mutation multiplicities for the genomic variants defined in the `somatic_path` column. For the format of this file, see the section above. Existence of the files in this column will be checked on prior to pipeline execution (if column is present in the table). 
- **cn_segments_genome** *[optional]*: The version of the genome in which the coordinates of copy number segments are specified, e.g., `hg38`. This column must not contain values which are numbers,  i.e. "hg38" is allowed value, but "38" is not. The genome version must be the same for all files listed in **cn_segments_path** column.
- **cn_segments_path** *[optional]*: The full path on your system (computer/HPC/*etc*) to the file containing copy number segments of the tumor genome found in the corresponding participant. For the format of this file, see the section above. Existence of the files in this column will be checked on prior to pipeline execution (if column is present in the table).
- **cohort_name**  *[essential]*: . This column must not contain values which are numbers, i.e. "GEL_1" is allowed value, but "009" is not.

Cancer cohorts that include multiple histological subtypes (for example, a `pan-lung` cancer cohort may include tumor samples from adenocarcinomas, squamous cell carcinomas, mesotheliomas, neuroendocrine carcinomas, *etc.*) can be defined as shown in lines 7-12 of the table. It is preferable that the **participant_tumor_subtype** column contains the actual histological subtype of the tumor, rather than a "pan-lung" substitute.
