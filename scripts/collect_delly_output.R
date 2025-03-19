### collect_delly_output.R #########################################################################
# Collect SV and CNV calls from Delly

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
# load libraries
library(plyr);
library(org.Hs.eg.db);
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-d', '--directory', type = 'character', help = 'path to data directory',
	default = getwd());
parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'reference type', default = 'hg38');
parser$add_argument('-t', '--targets', type = 'character', help = 'path to target intervals');
parser$add_argument('-g', '--germline', type = 'logical', help = 'is this germline data?', default = FALSE);

arguments <- parser$parse_args();

setwd(arguments$directory);

### READ DATA ######################################################################################
# find data files
if (arguments$germline) {

	input.file <- list.files(pattern = 'filtered_germline_SVs.vcf');
	input.file <- rev(sort(input.file));

	# find VCF header
	tmp <- read.delim(input.file[1], nrow = 1000, header = FALSE);
	header <- max(grep('^##', tmp$V1));
	rm(tmp);

	# read callset
	delly.calls <- read.delim(input.file, skip = header);

	} else {

	input.files <- list.files(pattern = 'Delly_SVs_somatic_hc.vcf$', recursive = TRUE);

	# find VCF header
	tmp <- read.delim(input.files[1], nrow = 4000, header = FALSE);
	header <- max(grep('##', tmp$V1));
	rm(tmp);

	# read callset
	delly.calls <- join_all(lapply(input.files, read.delim, skip = header), type = 'full');
	}

# get common field ids
common.fields <- colnames(delly.calls)[1:9];

# get sample list
all.samples <- colnames(delly.calls)[10:ncol(delly.calls)];

# write function to extract from INFO field
extract.info <- function(x, metric) {
	tmp <- unlist(strsplit(x,';'));
	if (any(grepl(metric,tmp))) {
		tmp2 <- tmp[grepl(metric, tmp)];
		return(unlist(strsplit(tmp2,'='))[2]);
		} else {
		return(NA);
		}
	}

# write function to extract from FORMAT field
extract.format <- function(x, metric) {
	formats <- unlist(strsplit(x[1],':'));
	idx <- which(formats == metric);
	return(unlist(strsplit(x[2],':'))[idx]);
	}

# format SV info
sv.data <- data.frame(
	Break1_Chromosome = delly.calls$X.CHROM,
	Break1_Position = delly.calls$POS,
	Break1_Gene = NA,
	Break2_Chromosome = as.character(sapply(delly.calls$INFO, extract.info, metric = 'CHR2')),
	Break2_Position = as.numeric(sapply(delly.calls$INFO, extract.info, metric = 'END')), 
	Break2_Gene = NA,
	SV_Type = as.character(sapply(delly.calls$INFO, extract.info, metric = 'SVTYPE')),
	SV_Conf = as.character(sapply(delly.calls$INFO, function(i) { unlist(strsplit(i,';'))[1] } ))
	);

if (any(is.na(sv.data$Break2_Chromosome))) {
	sv.data[is.na(sv.data$Break2_Chromosome),]$Break2_Chromosome <- sv.data[is.na(sv.data$Break2_Chromosome),]$Break1_Chromosome;
	}

# annotate regions
if ('hg38' == arguments$ref_type) {
	library(TxDb.Hsapiens.UCSC.hg38.knownGene);
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;
	} else if ('hg19' == arguments$ref_type) {
	library(TxDb.Hsapiens.UCSC.hg19.knownGene);
	txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;
	} else {
	stop('Unrecognized ref_type; must be one of hg38 or hg19');
	}

gene.gr <- transcriptsBy(txdb, by = 'gene');
gene.info <- select(txdb,
	keys = names(gene.gr),
	keytype = 'GENEID',
	columns = c('TXNAME','TXCHROM','TXSTART','TXEND')
	);

gene.info <- gene.info[order(gene.info$TXCHROM, gene.info$TXSTART),];

gene.info$Symbol <- mapIds(org.Hs.eg.db,
	keys = gene.info$GENEID,
	keytype = "ENTREZID",
	column = "SYMBOL",
	multiVals = "first"
	);

gene.info$GeneType <- mapIds(org.Hs.eg.db,
	keys = gene.info$GENEID,
	keytype = "ENTREZID",
	column = "GENETYPE",
	multiVals = "first"
	);

gene.info$GeneType <- factor(gene.info$GeneType,
	levels = c('protein-coding','rRNA','scRNA','snRNA','snoRNA','ncRNA','pseudo')
	);

gene.gr <- makeGRangesFromDataFrame(gene.info,
	seqnames.field = 'TXCHROM',
	start.field = 'TXSTART',
	end.field = 'TXEND'
	);

break1.gr <- makeGRangesFromDataFrame(sv.data,
	seqnames.field = 'Break1_Chromosome',
	start.field = 'Break1_Position',
	end.field = 'Break1_Position'
	);
break2.gr <- makeGRangesFromDataFrame(sv.data,
	seqnames.field = 'Break2_Chromosome',
	start.field = 'Break2_Position',
	end.field = 'Break2_Position'
	);

break1.overlaps <- as.data.frame(findOverlaps(gene.gr, break1.gr));
break2.overlaps <- as.data.frame(findOverlaps(gene.gr, break2.gr));

for (i in unique(break1.overlaps$subjectHits)) {
	gene.idx <- break1.overlaps[which(break1.overlaps$subjectHits == i),]$queryHits;
	tmp <- gene.info[gene.idx,];
	gene <- as.character(tmp[order(tmp$GeneType, tmp$TXSTART),]$Symbol[1]);
	sv.data[i,]$Break1_Gene <- gene;
	}
for (i in unique(break2.overlaps$subjectHits)) {
	gene.idx <- break2.overlaps[which(break2.overlaps$subjectHits == i),]$queryHits;
	tmp <- gene.info[gene.idx,];
	gene <- as.character(tmp[order(tmp$GeneType, tmp$TXSTART),]$Symbol[1]);
	sv.data[i,]$Break2_Gene <- gene;
	}

# format per-sample
my.data <- list();

for (smp in all.samples) {

	tmp <- sv.data;
	tmp$Sample <- gsub('\\.','-',smp);
	tmp$Genotype <- apply(delly.calls[,c('FORMAT',smp)], 1, extract.format, metric = 'GT');
	tmp$Read_Counts <- as.numeric(apply(delly.calls[,c('FORMAT',smp)], 1, extract.format, metric = 'RC'));
	tmp$Ref_Read_Pairs <- as.numeric(apply(delly.calls[,c('FORMAT',smp)], 1, 
		extract.format, metric = 'DR'));
	tmp$Var_Read_Pairs <- as.numeric(apply(delly.calls[,c('FORMAT',smp)], 1, 
		extract.format, metric = 'DV'));
	tmp$Ref_Junction_Reads <- as.numeric(apply(delly.calls[,c('FORMAT',smp)], 1, 
		extract.format, metric = 'RR'));
	tmp$Var_Junction_Reads <- as.numeric(apply(delly.calls[,c('FORMAT',smp)], 1, 
		extract.format, metric = 'RV'));

	tmp$Filter <- apply(delly.calls[,c('FORMAT',smp)], 1, extract.format, metric = 'FT');
	tmp$CN <- if ('RDCN' %in% unlist(strsplit(delly.calls$FORMAT[1],':'))) {
		as.integer(apply(delly.calls[,c('FORMAT',smp)], 1, extract.format, metric = 'RDCN'));
		} else {
		as.integer(apply(delly.calls[,c('FORMAT',smp)], 1, extract.format, metric = 'CN'));
		}

	my.data[[smp]] <- tmp[!is.na(tmp$Genotype),];
	gc();
	}

combined.data <- do.call(rbind, my.data);

# remove non-variants
combined.data <- combined.data[which(! combined.data$Genotype %in% c('./.', '0/0')),];

# if target regions were provided
if (!is.null(arguments$targets)) {

	# get target regions
	target_bed <- read.delim(arguments$targets, header = FALSE, comment.char = '#');
	colnames(target_bed)[1:3] <- c('Chromosome','Start','End');
	target_bed$Start <- target_bed$Start - 100;
	target_bed$End <- target_bed$End + 100;

	# create genomic ranges object for target regions
	target.gr <- makeGRangesFromDataFrame(target_bed, starts.in.df.are.0based = TRUE);

	# create genomic ranges object for each breakpoint
	first_bp <- data.frame(
		Chromosome = sv.data$Break1_Chromosome,
		Start = sv.data$Break1_Position,
		End = sv.data$Break1_Position + 1
		);

	second_bp <- data.frame(
		Chromosome = sv.data$Break2_Chromosome,
		Start = sv.data$Break2_Position,
		End = sv.data$Break2_Position + 1
		);

	bp1.gr <- makeGRangesFromDataFrame(first_bp, starts.in.df.are.0based = FALSE);
	bp2.gr <- makeGRangesFromDataFrame(second_bp, starts.in.df.are.0based = FALSE);

	# find overlaps
	overlaps.p1 <- as.data.frame(findOverlaps(bp1.gr, target.gr));
	overlaps.p2 <- as.data.frame(findOverlaps(bp2.gr, target.gr));

	overlap.data <- merge(
		overlaps.p1,
		overlaps.p2,
		by = 'queryHits',
		suffixes = c('.1','.2'),
		all = TRUE
		);

	# only keep entries for which both breakpoints are within target regions
	to.remove <- which(is.na(overlap.data$subjectHits.1) | is.na(overlap.data$subjectHits.2));
	keep.idx <- unique(overlap.data[-to.remove,]$queryHits);

	# filter initial input
	sv.data.filtered <- sv.data[keep.idx,];

	# filter sample input
	sample.filtered <- merge(sv.data.filtered, combined.data);

	# filter by call quality
	sample.filtered <- sample.filtered[which(sample.filtered$Filter == 'PASS'),];

	} else {

	# filter by call quality
	sample.filtered <- combined.data[which(combined.data$Filter == 'PASS'),];
	}

# save filtered data
write.table(
	sample.filtered,
	file = generate.filename(arguments$project, 'Delly_output_filtered','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('FormatData','SessionProfile','txt'));
