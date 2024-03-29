### runASCAT.R #####################################################################################
# R script to run ASCAT on a single tumour/normal pair
# Implemented and tested on HPC4H using R/4.1.0  module
# Original by JB; modified by SDP for use in pipeline-suite

### PREAMBLE #######################################################################################
library(optparse);

option_list <- list(
	make_option(c("-t", "--tumour_bam"), type="character", default=NULL, 
		help="path to tumour bam file [default= %default]", metavar="character"),
	make_option(c("-n", "--normal_bam"), type="character", default=NULL,
		help="path to matched normal bam file [default= %default]", metavar="character"),
	make_option(c("-s", "--sample_name"), type="character", default=NULL,
		help="Sample name [default= %default]", metavar="character"),
	make_option(c("-w", "--working_dir"), type="character", default=NULL,
		help="working directory (for temporary files) [default= %default]", metavar="character"),
	make_option(c("-o", "--out_dir"), type="character", default=NULL,
		help="output directory [default= %default]", metavar="character"),
	make_option(c("-b", "--genome_build"), type="character", default='hg38',
		help="reference type [default= %default]", metavar="character"),
	make_option(c("-r", "--ref_file"), type="character", default=NULL,
		help="reference file (SNP6 positions) [default= %default]", metavar="character"),
	make_option(c("-c", "--n_threads"), type="character", default=1,
		help="number of threads available for splitting [default= %default]", metavar="character"),
	make_option(c("-l", "--lib_paths"), type="character", default='/cluster/projects/pughlab/src/r_lib/library/4.1/',
		help="path to library if not installed in default [default= %default]", metavar="character")
	);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# load required libraries
if (!is.null(opt$lib_paths)) {
	.libPaths(c(opt$lib_paths,.libPaths()));
	}

library(maftools);
library(ASCAT);

### MAIN ###########################################################################################
# confirm sample ID
if (is.null(opt$sample_name)) {
	sample_name <- basename(opt$tumour_bam);
	} else {
	sample_name <- opt$sample_name;
	}

# move to output directory
setwd(opt$working_dir);

# read in marker loci
counts <- maftools::gtMarkers(
	t_bam = opt$tumour_bam,
	n_bam = opt$normal_bam,
	build = opt$genome_build,
	loci = opt$ref_file,
	nthreads = opt$n_threads
	);

# run prepAscat
ascat.bc <- maftools::prepAscat(
	t_counts = gsub(".bam","_nucleotide_counts.tsv",basename(opt$tumour_bam)),
	n_counts = gsub(".bam","_nucleotide_counts.tsv",basename(opt$normal_bam)),
	sample_name = sample_name
	);

# load prepared data
ascat.bc <- ASCAT::ascat.loadData(
	Tumor_LogR_file = paste0(sample_name, ".tumour.logR.txt"),
	Tumor_BAF_file = paste0(sample_name, ".tumour.BAF.txt"),
	Germline_LogR_file = paste0(sample_name, ".normal.logR.txt"),
	Germline_BAF_file = paste0(sample_name, ".normal.BAF.txt"),
	isTargetedSeq = FALSE
	);

# move to final output directory
setwd(opt$out_dir);

# plot raw data
ASCAT::ascat.plotRawData(ASCATobj = ascat.bc, img.prefix = "tumor");

# run ascat
ascat.bc <- ASCAT::ascat.aspcf(ascat.bc);
ASCAT::ascat.plotSegmentedData(ascat.bc);
ascat.output <- ASCAT::ascat.runAscat(ascat.bc, gamma = 1); # as per github documentation for HTS

# output purity and ploidy
ascat_tab <- data.frame(matrix(nrow = 1, ncol = 5));
colnames(ascat_tab) <- c('Sample','purity','ploidy','goodnessOfFit','ascatStatus');

ascat_tab$Sample <- sample_name;
if (!is.null(ascat.output$purity)) {
	ascat_tab$purity <- ascat.output$purity;
	ascat_tab$ploidy <- ascat.output$ploidy;
	ascat_tab$goodnessOfFit <- ascat.output$goodnessOfFit;
	ascat_tab$ascatStatus <- 'success';
	} else {
	ascat_tab$ascatStatus <- 'failed';
	}

write.table(
	ascat_tab,
	file = paste0(sample_name,"_ASCAT_purity_ploidy.txt"),
	row.names = FALSE,
	col.names = TRUE,
	quote = FALSE,
	sep = '\t'
	);

# Output absolute and raw CN seg files (formatted to look like sequenza output for sake of ease)
if (!is.null(ascat.output$segments)) {
	as_seg <- ascat.output$segments[,-1];
	colnames(as_seg) <- c("chromosome", "start.pos", "end.pos", "A", "B");
	as_seg$CNt <- as_seg$A + as_seg$B;

	raw_seg <- ascat.output$segments_raw[,-1];
	LogR_runs <- rle(ascat.bc$Tumor_LogR_segmented[,1]);

	seg <- raw_seg[,c("chr","chr","startpos","endpos")];
	colnames(seg) <- c("ID", "chrom", "loc.start", "loc.end");
	seg$ID <- sample_name;
	seg$num.mark <- LogR_runs$lengths;
	seg$seg.mean <- LogR_runs$values;

	# write outputs to file
	write.table(
		as_seg,
		file = paste0(sample_name, "_absolute_CN_segments.txt"),
		row.names = FALSE,
		col.names = TRUE,
		quote = FALSE,
		sep = '\t'	
		);

	write.table(
		seg,
		file = paste0(sample_name, "_Total_CN.seg"),
		row.names = FALSE,
		col.names = TRUE,
		quote = FALSE,
		sep = '\t'	
		);
	}
