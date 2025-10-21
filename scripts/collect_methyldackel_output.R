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
parser$add_argument('-m', '--mane', type = 'character', help = 'path to tsv containing MANE transcript IDs');

arguments <- parser$parse_args();

setwd(arguments$directory);

# load libraries/functions
library(yaml);
library(GenomicRanges);
library(BiocParallel);
library(parallel);
library(bsseq);
library(org.Hs.eg.db);

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
	target.bed <- read.table(arguments$target_bed, header = FALSE);
	colnames(target.bed)[1:3] <- c('chrom','start','end');
	target.bed <- makeGRangesFromDataFrame(target.bed, starts.in.df.are.0based = TRUE,
		keep.extra.columns = TRUE);
	
	target.loci <- subsetByOverlaps(cpg.loci, target.bed);
	} else {
	target.loci <- cpg.loci;
	}

# read in files using bsseq (runs in parallel so quite a bit faster)
# for some reason, this occasionally fails (but will subsequently succeed with no changes)
# typically for MultiMMR (targeted-panel), k = 1 (nThread = 1) is sufficient but for 
# WGS, parallelization is strongly recommended (k > 1 and nThread > 1)
BSobj <- NULL;
k <- c(16,8,4,2,1);
attempt <- 0;

while (is.null(BSobj) && attempt <= 5) {
	attempt <- attempt + 1;
	try(
		BSobj <- bsseq::read.bismark(
			files = data.files,
			loci = target.loci,
			BPPARAM = MulticoreParam(workers = k[attempt], progressbar = TRUE),
			colData = sample.info,
			strandCollapse = FALSE,
			nThread = ifelse(detectCores() > 16,  4, 1),
			verbose = TRUE
			)
		);
	}

save(BSobj, target.loci, sample.info,
	file = generate.filename(arguments$project, 'collected_CpG_objects','RData')
	);

### FORMAT ANNOTATION ##############################################################################
if (arguments$ref_type %in% c('hg38','GRCh38')) {
	library(TxDb.Hsapiens.UCSC.hg38.knownGene);
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;
	} else if (arguments$ref_type %in% c('hg19','GRCh37')) {
	library(TxDb.Hsapiens.UCSC.hg19.knownGene);
	txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;
	}
 
# extract gene positions and annotations
if (!is.null(arguments$target_bed)) {
	gene.positions <- data.frame(subsetByOverlaps(genes(txdb), target.bed));
	} else {
	gene.positions <- data.frame(genes(txdb));
	}

####
if (!is.null(arguments$mane)) {

	mane <- read.delim(arguments$mane, comment.char = '#');
	mane$TXNAME <- mane$Transcript.stable.ID.version;
	mane$REFSEQ <- sapply(mane$RefSeq.match.transcript..MANE.Select., function(i) {
		unlist(strsplit(i,'\\.'))[1] } );

	tx.positions <- data.frame(transcripts(txdb, columns = c('GENEID','TXNAME')));
	promoter.positions <- data.frame(promoters(txdb, columns = c('GENEID','TXNAME')));

	tx.positions$GENEID <- sapply(tx.positions$GENEID, function(i) { unlist(i)[1] } )
	promoter.positions$GENEID <- sapply(promoter.positions$GENEID, function(i) { unlist(i)[1] } )

	gene.annotations <- select(
		org.Hs.eg.db,
		keys = gene.positions$gene_id,
		keytype = 'ENTREZID',
		columns = c('SYMBOL','GENETYPE','REFSEQ')
		);

	gene.data <- merge(gene.annotations, tx.positions, by.x = 'ENTREZID', by.y = 'GENEID');
	gene.data <- merge(gene.data, mane[,c('TXNAME','REFSEQ')], by = c('TXNAME','REFSEQ'));

	gene.data <- gene.data[,c(3,4,1,2,5,6,7,8)];

	promoter.data <- merge(gene.annotations, promoter.positions, by.x = 'ENTREZID', by.y = 'GENEID');
	promoter.data <- merge(promoter.data, mane[,c('TXNAME','REFSEQ')], by = c('TXNAME','REFSEQ'));

	promoter.data <- promoter.data[,c(3,4,1,2,5,6,7,8)];

	} else {

	gene.annotations <- select(
		org.Hs.eg.db,
		keys = gene.positions$gene_id,
		keytype = 'ENTREZID',
		columns = c('SYMBOL','GENETYPE')
		);

	gene.data <- merge(gene.annotations,gene.positions,by.y = 'gene_id', by.x = 'ENTREZID');
	gene.data <- gene.data[,1:6];

	promoter.data <- merge(gene.annotations,gene.positions,by.y = 'gene_id', by.x = 'ENTREZID');

	promoter.data$Promoter.start <- promoter.data$start - 2000;
	promoter.data$Promoter.end <- promoter.data$start + 200;

	rev.idx <- which(promoter.data$strand == '-');
	promoter.data[rev.idx,]$Promoter.end <- promoter.data[rev.idx,]$end + 2000;
	promoter.data[rev.idx,]$Promoter.start <- promoter.data[rev.idx,]$end - 200;

	promoter.data <- promoter.data[,c(1:4,9:10)];
	colnames(promoter.data)[5:6] <- c('start','end');
	}

# do some additional formatting
gene.data$seqnames <- factor(gene.data$seqnames, levels = paste0('chr',c(1:22,'X','Y')));
gene.data <- gene.data[!is.na(gene.data$seqnames),];
gene.data <- gene.data[!is.na(gene.data$SYMBOL),];
gene.data <- gene.data[order(gene.data$seqnames, gene.data$start),];

gene.gr <- GRanges(gene.data);

colnames(gene.data)[which(colnames(gene.data) == 'seqnames')] <- 'Chromosome';
colnames(gene.data)[which(colnames(gene.data) == 'start')] <- 'Start';
colnames(gene.data)[which(colnames(gene.data) == 'end')] <- 'End';

promoter.data$seqnames <- factor(promoter.data$seqnames, levels = paste0('chr',c(1:22,'X','Y')));
promoter.data <- promoter.data[!is.na(promoter.data$seqnames),];
promoter.data <- promoter.data[!is.na(promoter.data$SYMBOL),];
promoter.data <- promoter.data[order(promoter.data$seqnames, promoter.data$start),];

promoter.gr <- GRanges(promoter.data);

colnames(promoter.data)[which(colnames(promoter.data) == 'seqnames')] <- 'Chromosome';
colnames(promoter.data)[which(colnames(promoter.data) == 'Promoter.start')] <- 'Start';
colnames(promoter.data)[which(colnames(promoter.data) == 'Promoter.end')] <- 'End';

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

# extract average methylation for each gene/sample
average.methylation.per.gene <- getMeth(
	BSobj, regions = gene.gr, type = 'raw', what = 'perRegion');
colnames(average.methylation.per.gene) <- rownames(pData(BSobj));
average.methylation.per.gene <- cbind(
	as.data.frame(gene.data),
	average.methylation.per.gene
	);

write.table(
	average.methylation.per.gene,
	file = generate.filename(arguments$project, 'methylation_per_gene','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# extract average methylation for each gene promoter region
average.methylation.per.promoter <- getMeth(
	BSobj, regions = promoter.gr, type = 'raw', what = 'perRegion');
colnames(average.methylation.per.promoter) <- rownames(pData(BSobj));
average.methylation.per.promoter <- cbind(
	as.data.frame(promoter.data),
	average.methylation.per.promoter
	);

write.table(
	average.methylation.per.promoter,
	file = generate.filename(arguments$project, 'methylation_per_promoter','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
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
