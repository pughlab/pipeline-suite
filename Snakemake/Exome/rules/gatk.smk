####################################################################################################
### PughLab Exome-Seq Pipeline #####################################################################
####################################################################################################

### INDEL REALIGNMENT/RECALIBRATION ################################################################
# define function to generate output file names
#def get_lane_output(wildcards):
#    return "{sample}_{lane}".format(
#	sample = units.loc[(wildcards.sample, wildcards.unit), "sample"],
#	lane = units.loc[(wildcards.sample, wildcards.unit), "lane"]
#	)

patients = list(set(samples["patient"]))

# define a function to list input files for input to GATK
def get_bwa_patient_bams():
    patient_bams = []
    for patient in samples.loc[wildcards.sample].patient:
        patient_bams.append(
	    "{project_dir}/BWA/{sample}/{sample}_bwamem_aligned_sorted_merged_markdup.bam".format(
		project_dir = project_dir,
		sample = wildcards.sample
            )
        )
    return lane_bams

### run GATK indel realignment and base-quality score recalibration
rule runGATK:
    input:
        ref = config["gatk_ref"],
        bam = "{project_dir}/BWA/{sample}/{sample}_bwamem_aligned_sorted_merged_markdup.bam",
        known1 = known_1000G_indels,
        known2 = known_mills,
        known3 = known_1000G_snps,
        dbsnp = dbsnp
    params:
        runtime = "72:00:00",
	mem = config["gatk_mem"],
	java_mem = config["gatk_java_mem"],
	extra_slurm_args = "-p all",
        tool = config["gatk_version"],
        jobname = "run_GATK_{sample}",
	tmp_dir = "{project_dir}/GATK/{sample}/TEMP",
        n_samples = 
        intervals = config["intervals_bed"]
    output:
        targets = "{project_dir}/GATK/{sample}/{sample}_target.intervals",
        realign_bam = temp("{project_dir}/GATK/{sample}/{sample}_split_realigned.bam"),
        bqsr = "{project_dir}/GATK/{sample}/{sample}.recal_data.grp",
        recal_bam = "{project_dir}/GATK/{sample}/{sample}_realigned_recalibrated.bam"
    threads: 1
    shell:
     """
     module load {params.tool}

     java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} -jar $gatk_dir/GenomeAnalysisTK.jar \
     -T RealignerTargetCreator \
     -R {input.ref} \
     -I {output.bam} \
     -o {output.targets} \
     -known {input.known1} \
     -known {input.known2} \
     --disable_auto_index_creation_and_locking_when_reading_rods -nt {params.n_samples} -dt None \
     --intervals {params.intervals} \
     --interval_padding 100

     if [ -s {output.targets} ]; then
       md5sum {output.targets} > {output.targets}.md5
     fi

     java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} -jar $gatk_dir/GenomeAnalysisTK.jar -T IndelRealigner \
     -R {input.ref} \
     -I {output.split_bam} \
     -targetIntervals {output.targets} \
     -o {output.realign_bam} \
     --generate_md5

     java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} -jar $gatk_dir/GenomeAnalysisTK.jar -T BaseRecalibrator \
     -R {input.ref} \
     -I {output.realign_bam} \
     -knownSites {input.known3} \
     -knownSites {input.dbsnp} \
     -o {output.bqsr} 

     if [ -s {output.bqsr} ]; then
       md5sum {output.bqsr} > {output.bqsr}.md5
     fi

     java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} -jar $gatk_dir/GenomeAnalysisTK.jar -T PrintReads \
     -R {input.ref} \
     -I {output.realign_bam} \
     -BQSR {output.bqsr} \
     -o {output.recal_bam} \
     --generate_md5

     if [ -s {output.recal_bam}.md5 ]; then
       rm -rf {params.tmp_dir}
     fi
     """

