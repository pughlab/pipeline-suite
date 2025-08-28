### summarize_global_methylation.R #################################################################
# Extract methylation metrics (percent methylated reads) for all sites and run comparisons

### FUNCTIONS ######################################################################################
# function to generate a standardized filename
generate.filename <- function(project.stem, file.core, extension, include.date = TRUE) {

	# build up the filename
	file.name <- paste(project.stem, file.core, sep = '_');
	file.name <- paste(file.name, extension, sep = '.');

	if (include.date) {
		file.name <- paste(Sys.Date(), file.name, sep = '_');
		}

	return(file.name);
	}

# function to write session profile to file
save.session.profile <- function(file.name) {

	# open the file
	sink(file = file.name, split = FALSE);

	# write memory usage to file
	cat('### MEMORY USAGE ###############################################################');
	print(proc.time());

	# write sessionInfo to file
	cat("\n### SESSION INFO ###############################################################");
	print(sessionInfo());

	# close the file
	sink();

	}

### PREAMBLE #######################################################################################
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-d', '--directory', type = 'character', help = 'path to data directory');
parser$add_argument('-s', '--sample_yaml', type = 'character', help = 'path to sample (BAM) yaml');
parser$add_argument('-t', '--target_bed', type = 'character', help = 'path to target regions (BED)');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'reference type', default = 'hg38');

arguments <- parser$parse_args();

setwd(arguments$directory);

# load libraries/functions
library(yaml);
library(GenomicRanges);
library(BiocParallel);
library(parallel);
library(bsseq);

### READ DATA ######################################################################################
# parse sample information from yaml file
project.yaml <- read_yaml(arguments$sample_yaml);

sample.info <- as.data.frame(matrix(nrow = 0, ncol = 3));
colnames(sample.info) <- c('Patient','Sample','Type');

patients <- names(project.yaml);

for (patient in patients) {
	normals <- names(project.yaml[[patient]]$normal);
	tumours <- names(project.yaml[[patient]]$tumour);
	sample.info <- rbind(sample.info, 
		data.frame(Patient = rep(patient, length(normals)), Sample = normals, Type = rep('normal', length(normals))),
		data.frame(Patient = rep(patient, length(tumours)), Sample = tumours, Type = rep('tumour', length(tumours)))
		);
	}

# find data files
data.files <- list.files(
	path = arguments$directory,
	pattern = 'CpG.bedGraph',
	recursive = TRUE,
	full.names = TRUE
	);

sample.order <- sapply(data.files, function(i) { gsub('_CpG.bedGraph','',basename(i)) } );

# determine sample covariates/order
sample.info$Type <- factor(sample.info$Type, levels = c('normal','tumour'));
sample.info$Sample <- factor(sample.info$Sample, levels = sample.order);
sample.info <- sample.info[order(sample.info$Sample),];
rownames(sample.info) <- sample.info$Sample;

# determine reference genome to use
if ('hg38' == arguments$ref_type) {
	library(BSgenome.Hsapiens.UCSC.hg38);
	genomedb <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38;
	} else if ('hg19' == arguments$ref_type) {
	library(BSgenome.Hsapiens.UCSC.hg19);
	genomedb <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19;
	} else {
	stop('Unrecognized ref_type; must be one of hg38 or hg19');
	}

# first find all candidate loci
cpg.loci <- findLoci(
	pattern = 'CG',
	subject = genomedb,
	include = paste0('chr',1:22)
	);

# read in target regions file
if (!is.null(arguments$target_bed)) {
	target.bed <- read.table(arguments$target_bed, header = F);
	colnames(target.bed) <- c('chrom','start','end','annotation');
	target.bed <- makeGRangesFromDataFrame(target.bed, starts.in.df.are.0based = TRUE,
		keep.extra.columns = TRUE);
	
	target.loci <- subsetByOverlaps(cpg.loci, target.bed);
	} else {
	target.loci <- cpg.loci;
	}

# read in files using bsseq (runs in parallel so quite a bit faster)
BSobj <- bsseq::read.bismark(
	files = data.files,
	loci = target.loci,
	BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
	colData = sample.info,
	strandCollapse = FALSE,
	#nThread = ifelse(detectCores() > 16,  4, 1),
	verbose = TRUE
	);

save(BSobj, target.loci, sample.info,
	file = generate.filename(arguments$project, 'collected_CpG_objects','RData')
	);

### FORMAT DATA ####################################################################################
# extract per-sample per-base methylation metrics
full.methylation.data <- getMeth(BSobj, type = 'raw', what = 'perBase');
colnames(full.methylation.data) <- rownames(pData(BSobj));
full.methylation.data <- cbind(
	as.data.frame(granges(BSobj)),
	full.methylation.data
	);

save(
	full.methylation.data,
	file = generate.filename(arguments$project, 'CpG_methylation_matrix','RData')
	);

# aggregate over chromosomes (WGS) or genes of interested (MultiMMR)
if (is.null(arguments$target_bed)) {

	# extract average methylation for each chromosome/sample
	per.chrom.methylation <- aggregate(
		full.methylation.data[,rownames(sample.info)],
		by = list(Chromosome = full.methylation.data$seqnames),
		mean,
		na.rm = TRUE
		);

	write.table(
		per.chrom.methylation,
		file = generate.filename(arguments$project, 'methylation_per_chromosome','tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);

	} else {

	# extract average methylation for each gene/sample
	average.methylation.per.region <- getMeth(
		BSobj, regions = target.bed, type = 'raw', what = 'perRegion');
	colnames(average.methylation.per.region) <- rownames(pData(BSobj));
	average.methylation.per.region <- cbind(
		as.data.frame(target.bed),
		average.methylation.per.region
		);

	write.table(
		average.methylation.per.region,
		file = generate.filename(arguments$project, 'methylation_per_target_region','tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectMethylDackelData','SessionProfile','txt'));

