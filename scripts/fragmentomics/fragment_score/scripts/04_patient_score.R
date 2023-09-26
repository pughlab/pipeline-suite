# file: runFrag.R
# author: Derek Wong, Ph.D
# date: ??
# modified: June 28th, 2023 by sprokopec

library(optparse);

## Set script variables
option_list <- list(
	make_option(c("--id"), type = "character", help = "sample id. Required"),
	make_option(c("--bam"), type = "character", help = "Path to bam file. Required."),
	make_option(c("--outdir"), type = "character", help = "Path to output directory. Required."),
	make_option(c("--ref"), type = "character", help = "Path to reference set.",
		default = '/cluster/home/sprokope/git/fragmentomics/fragment_score/ref/'),
	make_option(c("--libdir"), type = "character",
		help = "Path to scripts (git/fragmentomics/fragment_score).",
		default = '/cluster/home/sprokope/git/fragmentomics/fragment_score')
	);

opt <- parse_args(OptionParser(option_list = option_list));
print(opt);
options(scipen = 0, stringsAsFactors = FALSE);

## Get variables from input script
id <- opt$id;
bam <- opt$bam;
ref <- opt$ref;
libdir <- opt$libdir;
outdir <- opt$outdir;

## Import functions
source(paste0(libdir,"/functions.R"));

## Read in files
reference <- read.delim(ref, header = FALSE);
reference <- as.vector(reference$V1);

## Run script
Patient_score <- GeneratePatientFS(reference, bam);
write.table(
	Patient_score,
	file.path(outdir, paste0(id, "_score.txt")),
	row.names = FALSE,
	col.names = TRUE,
	sep = "\t"
	);

q('no');
