# PughLab pipeline-suite (version 1.1)

## Introduction
This is a collection of pipelines to be used for NGS analyses, from alignment to variant calling.

Start by creating a clone of the repository:

<pre><code>cd /path/to/some/directory
git clone https://github.com/pughlab/pipeline-suite/
</code></pre>

## Set up config files
There are example config files located in the "examples" folder:
- data configs:
  - fastq_config.yaml, can be generated using create_fastq_yaml.pl
  - gatk_input_config.yaml, can be generated using shared/create_bam_yaml.pl

- tool configs:
  - all tool configs must specify:
    - desired versions of tools
    - reference (file or directory) and ref_type (hg19 or hg38)
    - path to desired output directory (initial run or resumed run)
    - optional steps (mark_dup; del_intermediate; create_output_yaml; dry_run), either Y or N
    - HPC driver (for job submission)
    - memory and run time parameters for each step

  - bwa_tool_config.yaml, specifies:
    - path to bwa-indexed reference
    - sequencing centre and platform
    - optional steps: mark_dup (either Y or N)

  - star_tool_config.yaml, specifies:
    - path to STAR reference directory
    - sequencing centre and platform
    - optional steps: mark_dup (either Y or N)

   - rsem_tool_config.yaml, specifies:
    - path/stem to RSEM reference directory

   - star_fusion_tool_config.yaml, specifies:
    - path/stem to STAR-Fusion reference directory

  - gatk_tool_config.yaml, works for both DNA- and RNA-Seq data, specifies:
    - path to reference genome (requires .fa, .dict and .fai files)
    - path to target intervals (exome capture kit, if defined)
    - path to dbSNP file (if undefined, a default will be used)

   - haplotype_caller_config.yaml, specifies:
    - path to reference genome (requires .fa, .dict and .fai files)
    - path to vcf2maf.pl (RNA mode)
    - path to VEP (tool/version, cache data) (RNA mode)
    - path to ExAC data (for filtering/annotating with population allele frequencies) (RNA mode)
    - path to target intervals (exome capture kit, if defined) (DNA mode)
    - path to dbSNP file (if undefined, a default will be used) (DNA mode)

## Running a pipeline
If you are running these pipelines on the cluster, be sure to first load perl!

### Prepare a config file containing paths to FASTQ files:
<pre><code>cd /path/to/some/directory/pipeline-suite/

module load perl
perl create_fastq_yaml.pl -i /path/to/sampleInfo.txt -d /path/to/fastq/directory/ -o /path/to/fastq_config.yaml -t dna
</code></pre>

Where sampleInfo.txt is a tab-separate table containing two columms, and each row represents a single sample:

| Patient.ID | Sample.ID  |
| ---------- | ---------- |
| Patient1   | Patient1-N |
| Patient1   | Patient1-T |
| Patient2   | Patient2-T |

The assumption is that each FASTQ file will follow a similar naming convention, with Sample.ID used to identify files
for this sample, and lane information is pulled from the file name 
(for example: Patient1-N_135546_D00355_0270_AATCCGTC_L001.R1.fastq.gz, is a normal, library Patient1-N, lane 135546_D00355_0270_AATCCGTC_L001, R1).

### DNA pipeline:
For your initial run, make sure **output_dir:** is defined, and **resume_dir:** is empty (undefined).

<pre><code>cd /path/to/some/directory/pipeline-suite/

module load perl

# run BWA to align to a reference genome
perl bwa.pl -t /path/to/bwa_tool_config.yaml -c /path/to/fastq_config.yaml > /path/to/output/bwa_submission_out.log

# run GATK indel realignment and base quality score recalibration
perl gatk.pl -t /path/to/gatk_tool_config.yaml -c /path/to/bam_config.yaml --dna > /path/to/output/gatk_submission_out.log
</code></pre>

# run GATK's HaplotypeCaller to produce gvcfs
perl haplotype_caller.pl -t /path/to/haplotype_caller_config.yaml -c /path/to/gatk_bam_config.yaml --dna > /path/to/output/haplotype_call_submission_out.log
</code></pre>

### RNA pipeline:
<pre><code>cd /path/to/some/directory/pipeline-suite/

module load perl

# run STAR to align to a reference genome
perl star.pl -t /path/to/star_tool_config.yaml -c /path/to/fastq_config.yaml > /path/to/output/star_submission_out.log

# run RSEM on STAR-aligned BAMs
perl rsem.pl -t /path/to/rsem_tool_config.yaml -c /path/to/bam_config.yaml > /path/to/output/rsem_submission_out.log

# run STAR-Fusion on STAR-aligned BAMs
perl star_fusion.pl -t /path/to/star_fusion_tool_config.yaml -c /path/to/bam_config.yaml > /path/to/output/star_fusion_submission_out.log

# run GATK split CIGAR, indel realignment and base quality score recalibration on MarkDup BAMs
perl gatk.pl -t /path/to/gatk_tool_config.yaml -c /path/to/bam_config.yaml --rna > /path/to/output/rna_gatk_submission_out.log

# run GATK HaplotypeCaller, variant filtration and annotataion
perl haplotype_caller.pl -t /path/to/haplotype_caller_config.yaml -c /path/to/gatk_bam_config.yaml --rna > /path/to/output/haplotype_caller_submission_out.log
</code></pre>

### Resuming a run:
If the initial run is unsuccessful or incomplete, check the logs to identify the problem - it is most likely due to insufficient memory or runtime allocation. 
In this case, update the necessary parameters for the affected stage in the tool_config.yaml. 
Next update the "resume_dir" option to point to the run directory, for example
**resume_dir:** /path/to/output_dir/DATE_BWA_VERSION

Now, rerun as above! In the event of another failure, increase the memory or time requirements and try again.

## Output
bwa.pl will produce the following directory structure and output files in output_dir, with gatk.pl, star.pl, etc. having similar structures:

```
.
└── DATE_TOOL_VERSION
    ├── bam_config.yaml
    ├── logs
    │   ├── stage1_Patient1-N
    │   │   ├── script.sh
    │   │   └── slurm
    │   │       └── slurm-jobid.out
    │   └── stage1_Patient1-T
    │       ├── script.sh
    │       └── slurm
    ├── PATIENT1
    │   ├── Patient1-N
    │   │   ├── fastq_links
    │   │   └── lane
    │   ├── Patient1-N_output.bam
    │   ├── Patient1-T
    │   │   ├── fastq_links
    │   │   └── lane
    │   └── Patient1-T_output.bam
    └── PATIENT2
        ├── Patient2-T
        │   ├── fastq_links
        │   └── lane
        └── Patient2-T_output.bam
```
Note, in addition to .bam, there will also be .bai and .bam.md5 files, as well as other, tool specific output.
bam_config.yaml is created by create_final_yaml.pl - it lists the final output BAMs and can be used for downstream steps.

