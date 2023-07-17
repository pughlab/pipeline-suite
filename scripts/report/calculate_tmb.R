### calculate_tmb.R ################################################################################
# Calculate tumour mutation burden from total callable bases and mutation data

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

### PREPARE SESSION ################################################################################
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('-m', '--maf', type = 'character', help = 'path to input MAF file');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'hg38 or hg19', default = 'hg38');
parser$add_argument('-c', '--callable', type = 'character', help = 'path to callable bases file');
parser$add_argument('--method', type = 'character', default = 'txDb',
	help = "method to define genic/coding variants; one of 'VEP' (uses Variant_Classification field of maf), 'txDb' (uses UCSC transcript coordinates) or 'both'");

arguments <- parser$parse_args();

# import libraries
library(BoutrosLab.plotting.general);

### READ DATA ######################################################################################
# get mutation data
if (is.null(arguments$maf)) {
	stop('ERROR: No input MAF provided, please provide path to SNV/INDEL calls in MAF format.');
	} else {
	maf <- read.delim(arguments$maf, stringsAsFactors = FALSE, comment.char = '#');
	}

# get callable bases
if (is.null(arguments$callable)) {
	stop('ERROR: No callable bases provided, please provide path to total callable bases file.');
	} else {
	callable <- read.delim(arguments$callable);
	callable$Sample <- gsub('\\.','-', callable$Sample);
	}

# standardize methods
annotation.method <- tolower(arguments$method);

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# get sample set
all.samples <- intersect(callable$Sample, unique(maf$Tumor_Sample_Barcode));

# do some minor filtering
maf$t_vaf <- maf$t_alt_count / maf$t_depth;

maf <- maf[which(!maf$FLAG.low_coverage),];

# remove unnecessary fields from MAF
keep.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','HGVSp','HGVSp_Short','t_vaf');

maf <- unique(maf[,keep.fields]);

### ORGANIZE GENE DATA #############################################################################
if (arguments$method %in% c('txdb','both')) {

	library(org.Hs.eg.db);

	# find regions of interest (CDS/EXON)
	if (arguments$ref_type %in% c('hg38','GRCh38')) {
		library(TxDb.Hsapiens.UCSC.hg38.knownGene);
		txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;
		} else if (arguments$ref_type %in% c('hg19','GRCh37')) {
		library(TxDb.Hsapiens.UCSC.hg19.knownGene);
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;
		}

	# organize all transcripts by GeneID
	txdb_table <- transcriptsBy(txdb, by = 'gene');

	# extract GeneIDs
	tx_ids <- names(txdb_table);

	# extract 
	transcript_table <- select(
		TxDb.Hsapiens.UCSC.hg38.knownGene,
		keys = tx_ids,
		columns = c('GENEID','TXNAME','TXCHROM','TXSTART','TXEND'),
		keytype = 'GENEID'
		);
	colnames(transcript_table)[1] <- 'ENTREZID';

	# pull in extra annotations
	transcript_table$SYMBOL <- mapIds(
		org.Hs.eg.db,
		keys = transcript_table$TXNAME,
		keytype = 'UCSCKG',
		column = 'SYMBOL'
		);

	transcript_table$ENSEMBL <- mapIds(
		org.Hs.eg.db,
		keys = transcript_table$TXNAME,
		keytype = 'UCSCKG',
		column = 'ENSEMBL'
		);

	transcript_table$MAP <- mapIds(
		org.Hs.eg.db,
		keys = transcript_table$TXNAME,
		keytype = 'UCSCKG',
		column = 'MAP'
		);

	transcript_table$GENETYPE <- mapIds(
		org.Hs.eg.db,
		keys = transcript_table$TXNAME,
		keytype = 'UCSCKG',
		column = 'GENETYPE'
		);

	# remove unmatching cases
	transcript_table$MAPCHROM <- sapply(
		transcript_table$MAP,
		function(i) { paste0('chr', unlist(strsplit(i, 'p|q'))[1]) }
		);

	transcript_table <- transcript_table[which(transcript_table$TXCHROM == transcript_table$MAPCHROM),];

	# squish it to 1 entry per gene
	gene.annotation <- merge(
		aggregate(TXSTART ~ ENTREZID + SYMBOL + ENSEMBL + GENETYPE + TXCHROM, transcript_table, min),
		aggregate(TXEND ~ ENTREZID + SYMBOL + ENSEMBL + GENETYPE + TXCHROM, transcript_table, max),
		);
	colnames(gene.annotation)[5:7] <- c('Chromosome','Start','End');
	gene.annotation$Chromosome <- factor(gene.annotation$Chromosome, levels = paste0('chr',c(1:22,'X','Y')));
	gene.annotation <- gene.annotation[order(gene.annotation$Chromosome, gene.annotation$Start),];

	gene.gr <- GRanges(gene.annotation);

	# extract CDS info
	transcript_table <- select(
		TxDb.Hsapiens.UCSC.hg38.knownGene,
		keys = tx_ids,
		columns = c('GENEID','TXNAME','TXCHROM','TXSTART','TXEND','CDSID','CDSCHROM','CDSSTART','CDSEND'),
		keytype = 'GENEID'
		);

	transcript_table <- unique(transcript_table[!is.na(transcript_table$CDSID),2:5]);
	transcript_table$CDSCHROM <- factor(transcript_table$CDSCHROM, levels = paste0('chr',c(1:22,'X','Y')));
	transcript_table <- transcript_table[order(transcript_table$CDSCHROM, transcript_table$CDSSTART),];

	cds.gr <- GRanges(transcript_table[!is.na(transcript_table$CDSCHROM),]);
	}

### CALCULATE TMB ##################################################################################
# initiate list
tmb.data <- list();

for (smp in all.samples) {

	# get callable bases
	total <- callable[which(callable$Sample == smp),]$Patient.total / 10**6;
	gene <- callable[which(callable$Sample == smp),]$Patient.gene / 10**6;
	cds <- callable[which(callable$Sample == smp),]$Patient.cds / 10**6;

	# subset data
	for (vaf in c(0.05,0.1)) {

		smp.data <- maf[which(maf$Tumor_Sample_Barcode == smp & maf$t_vaf >= vaf),];

		snp.idx <- which(smp.data$Variant_Type == 'SNP');
		indel.idx <- which(smp.data$Variant_Type != 'SNP');

		# using UCSC txDb annotations
		if (arguments$method %in% c('txdb','both')) {

			anno.method <- 'UCSC txDb';

			# make a genomic ranges object
			smp.gr <- makeGRangesFromDataFrame(
				smp.data,
				keep.extra.columns = TRUE,
				seqnames.field = 'Chromosome',
				start.field = 'Start_Position',
				end.field = 'End_Position',
				ignore.strand = TRUE,
				starts.in.df.are.0based = FALSE
				);

			genic <- data.frame(subsetByOverlaps(smp.gr, gene.gr));
			coding <- data.frame(subsetByOverlaps(smp.gr, cds.gr));

			# classify variants
			genic.full <- which(
				genic$Variant_Classification != 'IGR' & 
				genic$Hugo_Symbol != 'Unknown'
				);
			genic.snp  <- intersect(genic.full, which(genic$Variant_Type == 'SNP'));
			genic.indel  <- intersect(genic.full, which(genic$Variant_Type != 'SNP'));

			coding.full <- which(!coding$Variant_Classification %in% c('Intron','Silent'));
			coding.snp  <- intersect(coding.full, which(coding$Variant_Type == 'SNP'));
			coding.indel  <- intersect(coding.full, which(coding$Variant_Type != 'SNP'));

			# organize results
			results <- data.frame(
				ID = smp,
				Method = anno.method,
				VAF = vaf,
				full = nrow(smp.data) / total,
				genic = length(genic.full) / gene,
				coding = length(coding.full) / cds,
				full.snp = nrow(smp.data[snp.idx,]) / total,
				genic.snp = length(genic.snp) / gene,
				coding.snp = length(coding.snp) / cds,
				full.indel = nrow(smp.data[indel.idx,]) / total,
				genic.indel = length(genic.indel) / gene,
				coding.indel = length(coding.indel) / cds
				);

			# push results to list
			tmb.data[[smp]] <- rbind(tmb.data[[smp]], results);
			gc();
			}

		# using VEP Variant_Classification field
		if (arguments$method %in% c('vep','both')) {

			anno.method <- 'VEP';

			# split by region
			genic.full <- smp.data[which(smp.data$Variant_Classification != 'IGR'),];
			genic.snp  <- which(genic.full$Variant_Type == 'SNP');
			genic.indel  <- which(genic.full$Variant_Type != 'SNP');

			coding.classes <- c('Missense_Mutation','Nonsense_Mutation','Splice_Site','Splice_Region','Translation_Start_Site','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins');

			coding.full <- genic.full[which(genic.full$Variant_Classification %in% coding.classes |
				(genic.full$HGVSp_Short != '' & !grepl('=$',genic.full$HGVSp_Short))),];

			coding.snp  <- which(coding.full$Variant_Type == 'SNP');
			coding.indel  <- which(coding.full$Variant_Type != 'SNP');

			# organize results
			results <- data.frame(
				ID = smp,
				Method = anno.method,
				VAF = vaf,
				full = nrow(smp.data) / total,
				genic = nrow(genic.full) / gene,
				coding = nrow(coding.full) / cds,
				full.snp = nrow(smp.data[snp.idx,]) / total,
				genic.snp = length(genic.snp) / gene,
				coding.snp = length(coding.snp) / cds,
				full.indel = nrow(smp.data[indel.idx,]) / total,
				genic.indel = length(genic.indel) / gene,
				coding.indel = length(coding.indel) / cds
				);

			# push results to list
			tmb.data[[smp]] <- rbind(tmb.data[[smp]], results);
			gc();
			}
		}
	}

tmb.results <- do.call(rbind, tmb.data);
tmb.results <- tmb.results[order(tmb.results$Method, tmb.results$VAF, tmb.results$ID),];

write.table(
	tmb.results,
	file = generate.filename(arguments$project, 'tumour_mutation_burden','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# make a symlink for use with final summary plots
unlink('tumour_mutation_burden.tsv');
file.symlink(
	generate.filename(arguments$project, 'tumour_mutation_burden', 'tsv'),
	'tumour_mutation_burden.tsv'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CalcTMB','SessionProfile','txt'));
