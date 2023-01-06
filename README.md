# PughLab pipeline-suite (version 0.8.1)

## Introduction
This is a collection of pipelines to be used for NGS (both DNA and RNA) analyses, from alignment to variant calling.

Start by creating a clone of the repository:

<pre><code>cd /path/to/some/directory
git clone git@github.com:pughlab/pipeline-suite.git
</code></pre>

## Dependencies
Pipeline-Suite was developed using perl and is currently compatible using the Slurm workload manager.

In addition, all of the desired tools (ie, BWA, GATK, etc) must be installed and be accessible using the 'module load \<tool\>' syntax.

Finally, a number of R packages are required:
- argparse, optparse, plyr, GenomicRanges
- org.Hs.eg.db, AnnotationDbi
- TxDb.Hsapiens.UCSC.hg19.knownGene, TxDb.Hsapiens.UCSC.hg38.knownGene
- for sequenza: sequenza
- for ichorCNA: ichorCNA, HMMcopy and GenomeInfoDb
- for ASCAT: ASCAT (https://github.com/VanLoo-lab/ascat) and maftools
- for panelCN.mops: panelcn.mops and CopyNumberPlots
- for summarization: xtable, deconstructSigs, BSgenome
- for visualizations: BPG (https://CRAN.R-project.org/package=BoutrosLab.plotting.general), UpSetR (https://cran.r-project.org/web/packages/UpSetR/index.html) and RCircos (https://github.com/hzhanghenry/RCircos)

## Workflow
See the [wiki](https://github.com/pughlab/pipeline-suite/wiki) for detailed instructions on how to set up a new project and run the pipeline.
