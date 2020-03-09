# PughLab pipeline-suite
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
  - bwa_tool_config.yaml, specifies:
    - desired versions of tools (bwa, samtools, picard)
    - path to bwa-indexed reference
    - path to desired output directory (initial run or resumed run)
    - sequencing centre and platform
    - optional steps (mark_dup; del_intermediate; dry_run), either Y or N
    - memory and run time parameters for each step (separate for tumour/normal)
  - gatk_tool_config.yaml, specifies:
    - desired versions of tools (gatk, samtools, picard)
    - path to reference genome
    - path to desired output directory (initial run or resumed run)
    - path to target intervals (exome capture kit, if defined)
    - sequencing centre and platform
    - optional steps (del_intermediate; dry_run), either Y or N
    - memory and run time parameters for each step
    
## Running a pipeline
If you are running these pipelines on the cluster, be sure to first load perl!

For your initial run, make sure **output_dir:** is defined, and **resume_dir:** is empty (undefined)
    
<pre><code>cd /path/to/some/directory/pipeline-suite/

module load perl
perl create_fastq_yaml.pl -d /path/to/fastq/directory/ -o /path/to/fastq_config.yaml

perl bwa.pl -t /path/to/bwa_tool_config.yaml -c /path/to/fastq_config.yaml > /path/to/output/bwa_submission_out.log

perl gatk.pl -t /path/to/gatk_tool_config.yaml -c /path/to/bam_config.yaml > /path/to/output/gatk_submission_out.log
</code></pre>

If the initial run is unsuccessful or incomplete, check the logs to identify the problem - it is most likely due to insufficient memory or runtime allocation. 
In this case, update the necessary parameters for the affected stage in the tool_config.yaml. 
Next update the "resume_dir" option to point to the run directory, for example
**resume_dir:** /path/to/output_dir/DATE_BWA_VERSION

Now, rerun as above! In the event of another failure, increase the memory or time requirements and try again.

## Output
bwa.pl and gatk.pl will produce the following directory structure and output files in output_dir:

```
.
└── DATE_TOOL_VERSION
    ├── bam_config.yaml
    ├── logs
    │   ├── stage1_patient1_sample1
    │   │   ├── script.sh
    │   │   └── slurm
    │   │       └── slurm-jobid.out
    │   └── stage1_patient1_sample2
    │       ├── script.sh
    │       └── slurm
    ├── PATIENT1
    │   ├── Patient1_sample1
    │   │   ├── fastq_links
    │   │   └── lane
    │   ├── patient1_sample1_output.bam
    │   ├── Patient1_sample2
    │   │   ├── fastq_links
    │   │   └── lane
    │   └── patient1_sample2_output.bam
    └── PATIENT2
        ├── Patient2_sample1
        │   ├── fastq_links
        │   └── lane
        └── patient2_sample1_output.bam
```
Note, in addition to .bam, there will also be .bai and .bam.md5 files, as well as other, tool specific output.
bam_config.yaml is created by create_bam_yaml.pl - it lists the final output BAMs and can be used for downstream steps.

