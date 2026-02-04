# fastq2bam_museum

This is a workflow for museum samples that automates processing of raw Illumina (short-reads) sequencing data from FASTQ to analysis-ready BAM files. This workflow was created following the original [fastq2bam](https://github.com/JGoudetGroup/fastq2bam?tab=readme-ov-file#) workflow developed by [Marianne Bachmann](https://github.com/m-bachmann), implementing a few modifications relevant to museum samples.


This workflow evaluates and trims raw illumina reads, identifying paired- an unpaired-R1 and -R2 reads, as well as merging overlapping-reads (a common outcome of short fragment sizes occurring in degraded DNA). Trimmed reads are mapped to a reference genome with proper read group tags added. Bam files obtained from the same library are are merged generating a single `_mapped.bam` file. Finally, duplicates are marked resulting in a sorted analysis-ready `_marked.bam` file per library.

The pipeline also allows for quality control at each processing step from raw read assessment to the final duplicated-marked bam files, generating comprehensive QC reports and aggregating them via MultiQC.

## Running the pipeline

### 1. Environment setup

Before executing the pipeline, ensure your environment is properly configured:

- Nextflow depends on java and it needs to previously load java.

```bash
module load openjdk/17.0.8.1_1
```
 Note the software version as the default java in cluster is not compatible.

- You can access the centralized Python environment at:

```bash
/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/venv/biopython
```

The `biopython` environment includes `pandas`, `Biopython` and `multiqc`. If you don't have access to the shared environment or are running elsewhere:

- Use the `bin/biopython_requirements.txt` file included in this repository.
- Create a local environment with it:

```bash
python3 -m venv $ENVNAME
source $ENVNAME/bin/activate
# Add required packages
pip install -r biopython_requirements.txt
```
Ensure you specify the path to this environment's `bin/activate` in the `params.json`.


### 2. Input files

The pipeline must be run with individuals sequenced on the **same instrument** (e.g. Novaseq) and with the same **library protocol** (e.g. Nextera or IDTxGen), this will determine the correct adapters to be trimmed as well as the assignment of read groups. This information should be specified in the `params.json` file and applies to all libraries equally. Memory, CPU and time allocation should be adjusted accordingly in the `nextflow.config` file. This workflow has been run in three different sets of museum samples and specifications on the individual changes for each set of samples are specified in their corresponding folder_scripts.

#### 2.1 FASTQ files

The FASTQ files are provided through **a comma separated file** containing:

|source_organism_id |specimen_id |lane | R1           | R2            |
|-------------------|------------|-----|--------------|---------------|
|Sample_ID             |Library name|lane |Path to read 1|Path to read 2 |

> While the source_organism_id is not needed for the naming of the file, it is required for the `SAMPLE`(`SM`) read group attribution. The library name does not always follow the same format, so it is not very flexible to add an automated extraction of the SampleID from the library name inside the pipeline.

#### 2.2 Reference genome

Path to the reference genome's `FASTA` file, where all the corresponding indexes are available too (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`).

#### 2.3 Adapter file

File containing a list of library adapters to remove using `trimmomatic`. See the [trimming](#52-trimming) section for more details.

#### 2.4 Scaffolds file

The [scaffolds bed file](lists/lg_scaffolds.bed) provided to `qualimap` consists of scaffolds from the 40 linkage groups identified in [Topaloudis *et al*., 2025](https://doi.org/10.1093/genetics/iyae190). This BED file enables assessment of coverage distribution across these genomic segments.

### 3. Parameters

#### 3.1 `museum_fastq2bam_params.json`

Below you can find the description of the parameters required in the [`museum_fastq2bam_params.json`](museum_fastq2bam_Germany_4th.samples_params.json) file:

|Parameter          | Description                                                                |Example                                   |
|-------------------|----------------------------------------------------------------------------|------------------------------------------|
|`input`            |Path to the CSV file listing all `FASTQ` files.                       | [example_fastq.csv](lists/museum_short.reads_2025_4th.samples.Germany.csv)|
|`bin_dir`          |Path to auxiliary scripts.                                                  |`"${projectDir}/bin"`                       |
|`result_dir`       |Path to the results directory.                                              |`"/work/user/project/1_Reads2Bam"`         |
|`scratch_dir`      |Path for intermediate or temporary files.                                   |`"/scratch/user/project/fastq2bam/scratch"` |
|`work_dir`         |Path to NF's `work` directory. Set to `/scratch/` partition for faster I/O. |`"/scratch/user/project/fastq2bam/work"`|
|`report_dir`       |Path to publish pipeline reports.                                           |`"/work/user/project/1_Reads2Bam/fastq2bam_report"`|
|`enable_reports`   |If `true`, generate NF report, trace and dag image.                         |`true` |
|`publish`          |The file publishing mode (NF parameter). Usually set to `copy`.             |`"copy"`            |
|`run_id`       |Sequencing run identifier.                                                  |`"Germany_2025_4th"` |
|`python`           |Select version in DCSR's software stack.                                    |`"python/3.12.1"`|
|`trimmomatic`      |Select version in DCSR's software stack.                                    |`"trimmomatic/0.39"`|
|`bwa`              |Select version in DCSR's software stack.                                    |`bwa/0.7.17`|
|`samtools`         |Select version in DCSR's software stack.                                    |`"samtools/1.19.2"`|
|`openjdk`          |Select version in DCSR's software stack.                                    |`"openjdk/17.0.8.1"`|
|`fastqc`           |Select version in DCSR's software stack.                                    |`"fastqc/0.12.1"`|
|`gcc`              |Select version in DCSR's software stack.                                    |`"gcc/12.3.0"`|
| `picard`          |Select version in DCSR's software stack.                                    |`"picard/3.1.1"`|
|`qualimap`         |Select version in DCSR's software stack.                                    |`qualimap/2.2.1"`|
|`bcftools`         |Select version in DCSR's software stack.                                    |`"bcftools/1.21"`|
|`htslib`         |Select version in DCSR's software stack.                                     |`"htslib/1.21"`|
|`python_venv`      |Path to the pipeline's Python virtual environment.                         |`"/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/venv/biopython/bin/activate"`|
|`phred`            |`trimmomatic`: base quality encoding.                                      |`"phred33"`|
|`illumina_adapters`|`trimmomatic`: path to Illumina adapters.                                  |`"/work/FAC/FBM/DEE/jgoudet/barn_owl/apulido/TOOLS/adapters/TruSeq3_NOVOGENE.fa"`|
|`illumina_clipping`|`trimmomatic`: [clipping parameters](#trimming).                           |`"2:30:10:3:70:true"`|
|`min_read_length`  |`trimmomatic`: minimum read length.                                          |`"70"`|
|`ref_genome`       |`bwa`: path to the reference genome FASTA.                                 |`"/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta"`|
|`ref_genome_version` |`sex_assignment`: version of the genome assembly selected                                  | `"2020"`|
|`RG_PL`            |`samtools addreplacerg`: sequencing platform.                                |`"ILLUMINA"`|
|`RG_CN`            |`samtools addreplacerg`: sequencing center.                                |`"NOVOGENE"`|
|`RG_PM`            |`samtools addreplacerg`: sequencing platform model.                          |`"NovaSeqXPlus"`|
|`qualimap_windows` |`qualimap bamqc -nw`: number of windows.                                   |`"1000"`|
|`qualimap_bed`     |`qualimap bamqc -gff`: path to feature file.                               |[lg_scaffolds.bed](lists/lg_scaffolds.bed)|
|`account`          |Cluster account name (`jgoudet_barn_owl`).                                       |`"jgoudet_barn_owl"`|
|`email`            |Mail for SLURM notifications.                                              |`"angelica.pulido@unil.ch"`|

#### 3.2 `nextflow.config`

Resource allocation is primarily configured in `fastq2bam_params.json`. The pipeline implements dynamic resource scaling - processes failing due to insufficient resources (exit codes 137, 140, or 143) will automatically retry with increased CPU, memory, and time allocations.

- Modify initial resource allocations for specific processes in the process-specific section:

  ```config
  // ==========================
  // PROCESS-SPECIFIC PARAMETERS
  // ==========================


  process {
      withName: getRG {
          cpus = 1
          memory = '2 GB'
          time = '10m'
          array = 100
      }
      // all other processes
  }
  ```

- Configure the execution environmet and error strategies in the *profile definition* section:

  ```config
  // ====================
  // PROFILE DEFINITION
  // ====================

  profiles {
      // for local execution - only during development
      standard {
          process.executor = 'local'
          errorStrategy = 'finish'
      }
      // for SLURM execution
      cluster {
          process.executor = 'slurm'
          // how to handle different error types
          errorStrategy = { task.exitStatus in [137, 143, 140] ? 'retry' : 'terminate' }
          // number of automatic retry attempts
          maxRetries = 3
          // SLURM-specific parameters
          process.clusterOptions = "--partition=cpu --account=${params.account}"
      }
  }
  ```


### 4. Nextflow execution

To run the script on the cluster using SLURM use the `-profile` flag with `cluster` as argument:

```bash
./nextflow run museum_fastq2bam.nf \
  -params-file museum_fastq2bam_params.json \
  -c museum_fastq2bam.config \
  -profile cluster
```

If the pipeline gets interrupted, resume from the last completed step, using the flag `-resume`:

```bash
./nextflow run museum_fastq2bam.nf \
  -params-file museum_fastq2bam_params.json \
  -c museum_fastq2bam.config \
  -profile cluster \
  -resume
```

### 5. Expecifications on the pipeline steps

#### 5.1. Read group extraction

This initial step extracts sequencing metadata for read group assignment. Each Illumina `FASTQ` header contains:

```bash
@LH00289:<runID>:<flowcellID>:<flowcellLane>:2220:14684:26984 1:N:0:<index>`
```

The `rg_parser.py` extracts the following from the header:

- sequencer run ID
- flow cell ID
- flow cell lane
- index

And combines it with information provided by the user in the `params.json` to write a `<specimen>_<run>_<lane>_RG_info.txt` file containing:

```text
<specimen>|<seqrun>.<flowcell_id>.<lane>|<index>|<sample>
```

This file is provided to the `samtools addreplacerg` process. The `sample` identifier represents the individual (`source_organism_id`), enabling integration of sequencing data from the same individual across different libraries, sequencing runs, or lanes.

#### 5.2. Merge overlapping reads
We use [`PEAR`](https://www.h-its.org/software/pear-paired-end-read-merger/), a paired-end read merger, to identify and assemble illumina paired-end reads that overlap due to the DNA fragment size being smaller than twice the length of the reads. This is a feature common on museum DNA, given their large degradation profile.

#### 5.3. Trimming

We use [`trimmomatic`](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for adapter removal. 

This step is divided in two processes `trim_PE` and `trim_SE`. In the former `paired` (R1 and R2) as well as `unpaired` (R1 and R2) reads are identified and trimmed. A second process (`trim_SE`) has as input the already assembled overlaping reads identified on the previous process (Merge overlapping reads).

The Illumina clipping parameters have been separated into logical components in the JSON parameter file for clarity.

- Illumina clip definition in `trimmomatic`:

```bash
`ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>`
```

- Illumina clip definition in the pipelines `params.json`:

```json
{
"illumina_adapters":"<fastaWithAdaptersEtc>",
"illumina_clipping":"<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>"
}
```

This structure provides flexibility to easily switch between different adapter sets (TruSeq, Nextera) without modifying the core pipeline code. All trimming logs are automatically captured and can be summarized using `multiqc` for quality assessment. Select the following parameters according to the library preparation kit used:

- **TruSeq**
  - "illumina_adapters": `/software/UHTS/Analysis/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa`
  - "illumina_clipping": `2:30:10:3:"true"`
- **Nextera**
  - "illumina_adapters": `/software/UHTS/Analysis/trimmomatic/0.36/adapters/NexteraPE-PE.fa`
  - "illumina_clipping": `2:30:10:3:"true"`

After trimming, 3 different set of reads are output:
  - `_{R1,R2}_paired.fastq.gz`: Refering to the forward and reverse PAIRED reads.
  - `_{R1,R2}_unpaired.fastq.gz`: Refering to the forward or reverse UNPAIRED reads.
  - `_assembledTrim.fastq.gz`: Refering to the ovelaping reads that were assembled by `PEAR`
  

All trimmed `FASTQs` are renamed as `{specimen_id}_{run_id}_{lane}_{R1,R2}_{paired,unpaired,assembledTrim}.fastq.gz`, thus removing the random string stemming from the sequencer or demultiplexing scripts. The run ID should be the same for all FASTQs and is provided in the `params.json`.

#### 5.3. Mapping

The mapping to the referent genome workflow consists of three integrated steps:

1. `bwa mem` aligns reads to the reference genome
2. `samtools` adds read group metadata tags
3. `picardtools` identifies and flags PCR duplicates

The pipeline will only output the BAM file with read groups (`<sample>_<seqrun>_<lane>_.bam`) and the marked duplicates (`<sample>_<seqrun>_<lane>_markdup.bam`).


  - `*{R1,R2}_PE_mapped.sam`: Refering to those forward and reverse paired 
  - *_UPF_mapped.sam"), 
  - path("*_UPR_mapped.sam"), 
  - path("*_SE_mapped.sam".

##### Adding read groups with samtools

Proper read group assignment is critical for downstream analysis, particularly for `GATK` workflows. The pipeline implements a flexible read group tagging system where all critical metadata is dynamically assigned (see [read group extraction](#51-read-group-extraction)):

|Tag | Description           | Source                   |
|----|-----------------------|--------------------------|
|`ID`|Unique read group ID   |FASTQ metadata extraction |
|`BC`|Barcode, if multiplexed|FASTQ metadata extraction |
|`CN`|Sequencing center      |Provided in `params.json` |
|`LB`|Library identifier     |Sample name               |
|`PL`|Sequencing platform    |Provided in `params.json` |
|`PM`|Platform model         |Provided in `params.json` |
|`SM`|Sample name            |FASTQ metadata extraction |

*Note:* The current implementation assumes all samples were sequenced using the same platform (e.g. Illumina NovaSeqX at our GTF center) with the same library preparation kit (Nextera). For studies **combining data from multiple sequencing platforms or centers**, process the samples per sequencing batch.

##### Mark duplicates

This pipeline uses Picard `MarkDuplicates` instead of `MarkDuplicatesSpark` (usually applied by the group). The Spark version caused frequent errors in Nextflow environments. This change affects only processing speed and resource usage, not the biological results—both tools use identical duplicate identification algorithms.

#### 5.4. Quality Control

The pipeline includes automated quality checks at key stages:

- Raw read counting
- `FastQC`: Raw and trimmed read quality assessment
- `Trimmomatic`: Trimming statistics and adapter removal metrics
- `GATK MarkDuplicates`: Duplicate rate and library complexity
- `Qualimap`: Coverage distribution and mapping quality on final BAMs
- `MultiQC`: Aggregated summary of all QC metrics

#### 5.5. Sex assignment based on genomic coverage stats

Sex determination is needed for genomic imputation with `glimpse`. More specifically, it informs ploidy settings for `bcftools call`. Since field-collected metadata may be missing or unreliable, we infer sex computationally using genomic coverage statistics. Currently, this process is implemented in the `assign_sex.py` script, which is specifically designed for the *Tyto alba* 2020 reference genome.

The script requires two inputs: **(1)** coverage statistics generated by `samtools coverage`, and **(2)** a reference genome identifier (e.g., "2020"). When the 2020 genome is specified, the script calculates the mean coverage depth of the sex chromosome (scaffold 13) relative to autosomes (excluding scaffolds 13 and 42). Females are identified by a lower sex-chromosome coverage ratio (expected ~0.5X due to ZZ/ZW heterogamety), while males show roughly equal coverage (homogametic, ~1X). The **current threshold is set to 0.65**, following analyses for sexing low coverage individuals in Eleonor's thesis (p. 274). Scaffold 42 is excluded from calculations due to its pseudoautosomal region (PAR), which could skew coverage-based sex determination.

The script outputs a tab-delimited file (SAMPLE\tSEX) for each individual. These files are used directly in the variant calling pipeline and later aggregated into a metadata table, which will also be used for the database.

💡 Future updates to the reference genome  will require modifications to the script. This includes adding new reference genome identifiers (e.g., "2025") and adjusting chromosome or scaffold names if they differ from the 2020 assembly. The logic for sex determination, however, will remain consistent—only the genomic coordinates and thresholds may need refinement.

#### 5.6. BAM List Output

The pipeline generates `bamlist_runIdGTF_2720.tsv`, a tab-separated file listing all BAMs produced:

| RingID     | Specimen ID   | Lane | BAM file                         |
|------------|---------------|------|----------------------------------|
| M026442    | M026442_3_3H  | L6   | M026442_3_3H_2720_L6_markdup.bam |
| **M038060**| M038060_5_12D | L6   | M038060_5_12D_2720_L6_markdup.bam|
| **M038060**| M038060_3_8B  | L6   | M038060_3_8B_2720_L6_markdup.bam |

This file enables merging of BAMs from the same individual across different lanes, runs, or libraries (e.g. *M038060*). See [Merging BAMs](#merging-bams) in the next section for database-compatible merge guidelines.

## Additional scripts & tools

### Merging BAMs

Example script from the 3000 owls survival project (/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/survival/) for merging BAMs from the same individual across different libraries. It is important to keep track of the files merged together for the databse. The code snippets here generate a table to conduct the merging, as well as a dataframe collecting the input and output of this merging process.

```python
import pandas as pd

# Path to the TSV file output by the fastq2bam pipeline containing the list of BAM files for low coverage individuals 
fastq2bam_out = '/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/survival/1_TrimmAlign/data/bamlist_runIdGTF_2720.tsv'

# Load the file
fastq2bam_df = pd.read_csv(fastq2bam_out, sep = '\t', header = None)
fastq2bam_df.rename(columns = {0:'source_organism_id', 
                               1:'lib_id', 
                               2: 'lane', 
                               3: 'file'
                               }, inplace = True)
duplicates = fastq2bam_df[fastq2bam_df['source_organism_id'].duplicated(keep=False)]
duplicates['file'] = '/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/survival/1_TrimmAlign/data/2.mapped/' + duplicates['file']

pairs = []

for source_id, group in duplicates.groupby('source_organism_id'):
    # Get all rows for this source_organism_id
    rows = group.reset_index(drop=True)
    
    # Create pairs (assuming you want all pairwise combinations)
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            pairs.append({
                'source_organism_id': source_id,
                'lib_id1': rows.loc[i, 'lib_id'],
                'lib_id2': rows.loc[j, 'lib_id'],
                'file1': rows.loc[i, 'file'],
                'file2': rows.loc[j, 'file']
            })

# Create the result DataFrame
paired = pd.DataFrame(pairs)

paired['bam_id'] = paired['source_organism_id'] + '_20250915'
paired['bam'] = '/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/survival/1_TrimmAlign/data/3.merged/' + paired['bam_id'] + '_markdup.bam'

# print a table to keep track of merged samples
paired[['source_organism_id', 'lib_id1', 'lib_id2', 'bam_id']].to_csv('/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/survival/1_TrimmAlign/data/DB_merged_samples.tsv', sep = '\t', index = False, header = True)

# print a table to help merging the BAMs
paired[[ 'bam', 'file1', 'file2']].to_csv('/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/survival/1_TrimmAlign/data/merging.tsv', sep = '\t', index = False, header = False)
```

Merging can then be executed with a simple loop either in an interactive node or a slurm job depending on the number of samples:

```bash
while read i; do eval "samtools merge $i"; done < merging.tsv
```

### Collect metadata for the database

- Genetic sex
- Read count (FASTQ)
- Coverage (from Qualimap)

### Gather QC for lab feedback

### Test dataset

- generating test datasets

```bash
cd data/

# Raw data directory
DIR="/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/survival/0_RawData"
# Grab names of files, removing the CTRL directory
ls $DIR | tail -n +2 > survival_fastq_list.txt

# get paths of selected IND
cat pheno_male_subsample.txt pheno_female_subsample.txt > subsample.txt
grep -f subsample.txt survival_fastq_list.txt > subsample_fastq.txt
rm subsample.txt

Sinteractive
module load seqtk
while read i; do
    eval "zcat "${DIR}${i}" | seqtk sample -s100 - 1000 | gzip > "FASTQ/${i}""
done < subsample_fastq.txt
```

### BAM2CRAM

For compressing BAMs that will not be used in the near future.
