# PughLab pipeline-suite Report (version 0.4.1)

## Introduction
This section of the PughLab Pipeline-Suite will generate a pretty report summarizing the pipeline output. It can be run as part of the main pipeline (ie, perl pughlab_dnaseq_pipeline.pl {...} --create_report), or separately once everything has run successfully (see below).

### create final report
pughlab_pipeline_auto_report.pl requires as input:
- dna_pipeline_config.yaml OR rna_pipeline_config.yaml
- project start date (the date the pipeline was initially run)

Note: much of this pipeline requires that the input files (ie, final calls from each tool) already exist so running with the --dry-run option may fail if something is missing.

<pre><code>
perl pughlab_pipeline_auto_report.pl \
-t /path/to/dna_pipeline_config.yaml \
-d DATE \
-c slurm \
--dry-run { if this is a dry-run; NOTE that this will fail if the above pipeline has not completed } \
--no-wait { if not a dry-run and you don't want to wait around for it to finish }
</code></pre>

Alternatively, you may run each portion separately if you only want summary plots for a single section:

Reminder, the report generation portion of this tool requires installation of the BPG plotting package for R:
https://CRAN.R-project.org/package=BoutrosLab.plotting.general, which is currently installed on H4H R/3.6.1

<pre><code>
module load perl
module load R/3.6.1

Plot QC metrics for either DNA- or RNA-Seq:
Rscript /path/to/pipeline-suite/scripts/report/plot_qc_metrics.R \
--output_dir /path/to/output/directory \
--project PROJECT_NAME \
--seq_type exome { one of wgs, exome or rna } \
--correlations /path/to/germline/correlations { produced by collect_germline_genotypes.R or collect_rnaseqc_output.R } \
--coverage /path/to/coverage/summary { produced by collect_coverage_output.R or collect_rnaseqc_output.R } \

And, for DNA, at least one of (if both are provided, will use --contamination):
--contest /path/to/ContEst/output { produced by collect_contest_output.R } \
--contamination /path/to/contamination/estimates { produced by collect_sequencing_metrics.R }

Note: must be run from parent directory (such that BWA or STAR folders are visible) in order to match tumour/normal or multi-tumour samples for each patient.

</code></pre>

For RNA-Seq:
<pre><code>Rscript plot_expression_summary.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--rsem /path/to/rsem/expression/matrix { produced by collect_rsem_output.R; not z-scores }

Rscript write_rna_fusion_summary.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
And at least one of:
--starfus /path/to/star-fusion/output { produced by collect_star-fusion_output.R } \
--fusioncatcher /path/to/fusioncatcher/output { produced by collect_fusioncatcher_output.R }

Rscript plot_rna_snv_summary.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--mutations /path/to/snv/calls { produced by collect_snv_output.R }

Rscript plot_viral_counts.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--viral_counts /path/to/viral/counts/matrix { produced by collect_fusioncatcher_output.R } \

Output summary methods for each tool that was run ({tool}->{run} = Y):
perl write_rna_methods.pl \
-t rna_pipeline_config.yaml \
-d /path/to/output/directory
</code></pre>

For DNA-Seq:
<pre><code>Rscript plot_germline_snv_summary.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--maf /path/to/germline/snv/matrix { produced by collect_snv_output.R; ideally CPSR annotated } 

Rscript plot_seqz_cna_summary.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--cnas /path/to/cna/matrix { produced by collect_sequenza_output.R } \
--metrics /path/to/ploidy/metrics { produced by collect_sequenza_output.R } \
--scale ratio { or CN depending on file provided to --cnas }

Rscript plot_gatk_cna_summary.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--cnas /path/to/cna/matrix { produced by collect_gatk_cnv_output.R } \
--metrics /path/to/pga/metrics { produced by collect_gatk_cnv_output.R } \
--scale ratio { or CN depending on file provided to --cnas }

Rscript plot_sv_summary.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--mavis /path/to/mavis/output { produced by collect_mavis_output.R }

Rscript format_ensemble_mutations.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--mutect /path/to/mutect/output { produced by collect_snv_output.R } \
--mutect2 /path/to/mutect2/output { produced by collect_snv_output.R } \
--strelka /path/to/strelka/output { produced by collect_snv_output.R } \
--somaticsniper /path/to/somaticsniper/output { produced by collect_snv_output.R } \
--vardict /path/to/vardict/output { produced by collect_snv_output.R } \
--varscan /path/to/varscan/output { produced by collect_snv_output.R } \
--pindel /path/to/pindel/output { produced by collect_snv_output.R } \
--coverage minimum depth to be considered callable { default 20,15 for T/N }

Rscript plot_snv_tool_summary.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--input /path/to/combined/tool/output { RData produced by format_ensemble_mutations.R }

module load MutSigCV/1.4;
sh /cluster/tools/software/MutSigCV/1.4/run_MutSigCV.sh /cluster/tools/software/MCR/8.1/v81 \
/path/to/ensemble/matrix { produced by format_ensemble_mutations.R } \
/cluster/projects/pughlab/references/MutSigCV/exome_full192.coverage.txt \
/cluster/projects/pughlab/references/MutSigCV/gene.covariates.txt \
/path/to/output/directory/STEM_MutSigCV \
/cluster/projects/pughlab/references/MutSigCV/mutation_type_dictionary_file.txt \
/cluster/projects/pughlab/references/MutSigCV/chr_files_hg38

Rscript apply_cosmic_mutation_signatures.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--input /path/to/ensemble/matrix { produced by format_ensemble_mutations.R } \
--ref_type hg38 { or hg19 } \
--seq_type exome { or wgs } \
--vaf_threshold 0.1 \
--signatures /path/to/COSMIC/SBS/matrix

Rscript plot_snv_summary.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--input /path/to/ensemble/matrix { produced by format_ensemble_mutations.R } \
--callable /path/to/callable/bases { produced by count_callable_bases.R } \
--msi /path/to/msi/estimates { produed by collect_msi_output.R } \
--mutsig /path/to/output/directory/STEM_MutSigCV.sig_genes.txt \
--seq_type exome { or wgs }

Output summary methods for each tool that was run ({tool}->{run} = Y):
perl write_wgs_methods.pl | write_wxs_methods.pl \
-t dna_pipeline_config.yaml \
-d /path/to/output/directory
</code></pre>

It is also suggested to annotate your final ensemble SNV calls using OncoKB. Unfortunately, OncoKB requires internet access and is therefore not supported on our clusters.
To run OncoKB from Samwise/Hamfast/H4H:build partition:
<pre><code>module load python3
python3 /path/to/oncokb-annotator/MafAnnotator.py \
-i /path/to/snv/matirx { in MAF format } \
-o /path/to/output/file \
-t CANCER_TYPE { as ONCOTREE CODE; optional } \
-r GRCh38 { or GRCh37; optional if provided in input MAF } \
-b ACCESS TOKEN { see https://www.oncokb.org/apiAccess }
</code></pre>

You can then summarize the OncoKB output (this also filters any tumour-only samples, using OncoKB as a white-list):
<pre><code>Rscript apply_oncokb_filter_and_summarize.R \
--project PROJECT_NAME \
--output_dir /path/to/output/directory \
--input /path/to/oncokb/matrix { produced by format_ensemble_mutations.R and annotated by OncoKB } \
</code></pre>

Once all plotting/summary functions are complete, create the report pdf:
<pre><code>perl generate_report.pl \
--title PROJECT_NAME { or other title } \
--date DATE { project start date } \
--input_dir /path/to/plot/directory { location of above plots/summary } \
--output_dir /path/to/output/directory

module load texlive
pdflatex Report.tex { produced by generate_report.pl }
pdflatex Report.tex { run 2x to ensure proper indexing }
</code></pre>
