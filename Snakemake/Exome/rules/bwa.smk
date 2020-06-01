####################################################################################################
### PughLab Exome-Seq Pipeline #####################################################################
####################################################################################################

### ALIGNMENT ######################################################################################
# extract read1/read2 for input validation
def get_r1(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit),"read1"]

def get_r2(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit),"read2"]

# define function to generate readgroup for SAM/BAM header
def get_readgroup(wildcards):
    return "'@RG\tID:{sample}\tSM:{sample}\tPL:Illumina\tLB:{lib}\tPU:{lane}'".format(
        sample = wildcards.sample,
	lib = units.loc[(wildcards.sample, wildcards.unit), "library"],
	lane = units.loc[(wildcards.sample, wildcards.unit), "lane"]
	)

# define a function to list input files for merge
def get_bwa_lane_bams(wildcards):
    lane_bams = []
    for unit in units.loc[wildcards.sample].lane:
        lane_bams.append(
	    "{project_dir}/BWA/{sample}/{sample}_{unit}_bwamem_aligned_sorted.bam".format(
		project_dir = project_dir,
		sample = wildcards.sample,
                unit = unit
            )
        )
    return lane_bams

### run BWA
rule runBWA:
    input:
        ref = config["bwa_ref"],
        r1 = get_r1,
        r2 = get_r2
    params:
        runtime = "72:00:00",
        mem = config["bwa_mem"],
	extra_slurm_args = "-p all",
        tool = config["bwa_version"],
        readgroup = get_readgroup,
	jobname = "run_bwa_sort_index_{sample}_{unit}",
        out_dir = "{project_dir}/BWA/{sample}",
        stem = "{sample}_{unit}"
    output:
        bam = "{project_dir}/BWA/{sample}/{sample}_{unit}_bwamem_aligned_sorted.bam",
        md5 = "{project_dir}/BWA/{sample}/{sample}_{unit}_bwamem_aligned_sorted.bam.md5"
    threads: 4
    shell:
     """
     module load {params.tool}

     cd {params.out_dir}

     bwa mem -M -t4 \
     -R {params.readgroup} \
     {input.ref} {input.r1} {input.r2} \
     > {params.stem}_bwamem_aligned.sam

     ### run SORT SAM + BAM conversion
     samtools sort -@4 -O bam -o {params.stem}_bwamem_aligned_sorted.bam \
     -T {params.stem} \
     {params.stem}_bwamem_aligned.sam

     ### index BAM file
     samtools index {output.bam}

     md5sum {output.bam} > {output.md5}
     """

### run MERGE / Markduplicates
rule runMERGE:
    input: get_bwa_lane_bams
    params:
        runtime = "24:00:00",
        mem = config["merge_mem"],
        java_mem = config["merge_java_mem"],
 	extra_slurm_args = "-p all",
        tool = config["picard_version"],
	jobname = "run_merge_markduplicates_{sample}",
        out_dir = "{project_dir}/BWA/{sample}/",
        tmp_dir = "{project_dir}/BWA/{sample}/TEMP"
    output:
        bam = "{project_dir}/BWA/{sample}/{sample}_bwamem_aligned_sorted_merged_markdup.bam",
        metrics = "{project_dir}/BWA/{sample}/{sample}_bwamem_aligned_sorted_merged_markdup.bam.metrics"
    threads: 1
    shell:
     """
     module load {params.tool}

     cd {params.out_dir}

     java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} -jar $picard_dir/picard.jar MarkDuplicates \
     INPUT={input} \
     OUTPUT={output.bam} \
     METRICS_FILE={output.metrics} \
     TMP_DIR={params.tmp_dir} \
     ASSUME_SORTED=true CREATE_INDEX=true CREATE_MD5_FILE=true \
     MAX_RECORDS_IN_RAM=100000 VALIDATION_STRINGENCY=SILENT
    """
