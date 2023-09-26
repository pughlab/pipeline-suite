# file: runFrag.R
# author: Derek Wong, Ph.D
# date: October 5th, 2021
# modified: June 28th, 2023 by sprokopec

library(optparse);

## Set script variables
option_list <- list(
	make_option(c("--id"), type = "character", help = "sample id. Required"),
	make_option(c("--bam"), type = "character", help = "Path to bam file. Required."),
	make_option(c("--outdir"), type = "character", help = "Path to output directory. Required."),
	make_option(c("--libdir"), type = "character", help = "Path to scripts (git/fragmentomics/ratio).", 
		default = '/cluster/home/sprokope/git/fragmentomics/ratio'),
	make_option(c("--filters"), type = "character", help = "Path to genomic blacklist regions",
		default = "extdata/filters.hg38.rda"),
	make_option(c("--gaps"), type = "character", help = "Path to genome gaps.",
		default = "extdata/gaps.hg38.rda"),
	make_option(c("--tiles"), type = "character", help = "Path to 100kb tiled genome.",
		default = "extdata/hg38_tiles.bed"),
	make_option(c("--VNTRs"), type = "character", help = "Path to VNTRs.",
		default = "extdata/VNTRs.hg38.rda"),
	make_option(c("--healthy"), type = "character", help = "Path to panel of healthy controls.",
		default = "extdata/healthy.median.hg38.rda")
	);

opt <- parse_args(OptionParser(option_list = option_list));
print(opt);
options(scipen = 0, stringsAsFactors = FALSE);

## Load required packages
library(tidyverse);
library(Rsamtools);
library(GenomicAlignments);
library(biovizBase);
library(BSgenome.Hsapiens.UCSC.hg38);

## Get variables from input script
id	<- opt$id;
bam	<- opt$bam;
outdir	<- opt$outdir;
libdir	<- opt$libdir;

filters <- if (file.exists(opt$filters)) { opt$filters;
	} else if (file.exists(paste0(libdir, '/', opt$filters))) {
	paste0(libdir, '/', opt$filters);
	} else if (file.exists(paste0(libdir, '/extdata/filters.hg38.rda'))) {
	paste0(libdir, '/extdata/filters.hg38.rda');
	} else {
	stop("No 'filters' file found.");
	}

gaps <- if (file.exists(opt$gaps)) { opt$gaps;
	} else if (file.exists(paste0(libdir, '/', opt$gaps))) {
	paste0(libdir, '/', opt$gaps);
	} else if (file.exists(paste0(libdir, '/extdata/gaps.hg38.rda'))) {
	paste0(libdir, '/extdata/gaps.hg38.rda');
	} else {
	stop("No 'gaps' file found.");
	}

tiles <- if (file.exists(opt$tiles)) { opt$tiles;
	} else if (file.exists(paste0(libdir, '/', opt$tiles))) {
	paste0(libdir, '/', opt$tiles);
	} else if (file.exists(paste0(libdir, '/extdata/hg38_tiles.bed'))) {
	paste0(libdir, '/extdata/hg38_tiles.bed');
	} else {
	stop("No 'tiles' file found.");
	}

VNTRs <- if (file.exists(opt$VNTRs)) { opt$VNTRs;
	} else if (file.exists(paste0(libdir, '/', opt$VNTRs))) {
	paste0(libdir, '/', opt$VNTRs);
	} else if (file.exists(paste0(libdir, '/extdata/VNTRs.hg38.rda'))) {
	paste0(libdir, '/extdata/VNTRs.hg38.rda');
	} else {
	stop("No 'VNTRs' file found.");
	}

healthy <- if (file.exists(opt$healthy)) { opt$healthy;
	} else if (file.exists(paste0(libdir, '/', opt$healthy))) {
	paste0(libdir, '/', opt$healthy);
	} else if (file.exists(paste0(libdir, '/extdata/healthy.median.hg38.rda'))) {
	paste0(libdir, '/extdata/healthy.median.hg38.rda');
	} else {
	stop("No 'healthy' file found.");
	}

## Create output directory
if (!dir.exists(outdir)) {
	dir.create(outdir);
	}

## Run scripts
source(paste0(libdir,"/R/git_01-read_fragments.R"));
source(paste0(libdir,"/R/git_02-mito_frag.R"));
source(paste0(libdir,"/R/git_03-100kb_bins.R"));
source(paste0(libdir,"/R/git_04-5Mb_bins.R"));
source(paste0(libdir,"/R/git_05-summary.R"));
source(paste0(libdir,"/R/git_06-plotting.R"));

q('no');
