### count_callable_bases.R #########################################################################
# Finds and combines output from callable bases (get_coverage).

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
# import libraries
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-d', '--directory', type = 'character', help = 'path to data directory');
parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'reference build (hg19 or hg38)', default = 'hg38');

arguments <- parser$parse_args();

# what's the date?
date <- Sys.Date();

setwd(arguments$directory);

### MAIN ###########################################################################################
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

# extract EXON info
transcript_table <- select(
        TxDb.Hsapiens.UCSC.hg38.knownGene,
        keys = tx_ids,
        columns = c('GENEID','TXNAME','TXCHROM','TXSTART','TXEND','EXONID','EXONCHROM','EXONSTART','EXONEND'),
        keytype = 'GENEID'
        );

transcript_table <- unique(transcript_table[!is.na(transcript_table$EXONID),2:5]);
transcript_table$EXONCHROM <- factor(transcript_table$EXONCHROM, levels = paste0('chr',c(1:22,'X','Y')));
transcript_table <- transcript_table[order(transcript_table$EXONCHROM, transcript_table$EXONSTART),];

exon.gr <- GRanges(transcript_table[!is.na(transcript_table$EXONCHROM),]);


# find results files
cov.files <- list.files(pattern = '^CallableBases.tsv$', recursive = TRUE);

# read them in
cov.list <- list();
sample.list <- c();

for (file in cov.files) {
	# extract sample ID
	smp <- unlist(strsplit(file, '\\/'))[1];
	# store data in list
	cov.list[[smp]] <- read.delim(file);

	these.smps <- setdiff(colnames(cov.list[[smp]])[6:ncol(cov.list[[smp]])],'TargetRegions');
	sample.list <- c(sample.list, these.smps);
	}

# determine number of callable bases
callable.bases <- data.frame(
	Sample = sample.list,
	Sample.total = NA,
	Sample.exon = NA,
	Sample.cds = NA,
	Patient.total = NA,
	Patient.exon = NA,
	Patient.cds = NA
	);

if ('TargetRegions' %in% colnames(cov.list[[1]])) {
	callable.bases$TargetRegions <- NA;
	callable.bases$Sample.targeted <- NA;
	callable.bases$Sample.exon.targeted <- NA;
	callable.bases$Sample.cds.targeted <- NA;
	callable.bases$Patient.targeted <- NA;
	callable.bases$Patient.exon.targeted <- NA;
	callable.bases$Patient.cds.targeted <- NA;
	}

for (i in 1:length(cov.list)) {

	tmp <- cov.list[[i]];

	# find samples
	these.smps <- intersect(colnames(tmp), sample.list);

	# find per-sample metrics
	for (smp in these.smps) {

		smp.idx <- which(callable.bases$Sample == smp);

		# make a genomic ranges object
		smp.gr <- GRanges(tmp[which(tmp[,smp] == 1),c('chrom','start','end')]);

		# get total callable
		callable.bases[smp.idx,]$Sample.total <- sum(data.frame(smp.gr)$width);

		# get total callable in exons
		callable.bases[smp.idx,]$Sample.exon <- sum(data.frame(intersect(smp.gr, exon.gr))$width);

		# get total callable in CDS
		callable.bases[smp.idx,]$Sample.cds <- sum(data.frame(intersect(smp.gr, cds.gr))$width);

		# get overlap with target regions too if required
		if ('TargetRegions' %in% colnames(tmp)) {

			target.gr <- GRanges(tmp[which(tmp$TargetRegions == 1),c('chrom','start','end')]);

			# get total target
			callable.bases[smp.idx,]$TargetRegions <- sum(data.frame(target.gr)$width);

			target.smp.gr <- intersect(target.gr, smp.gr);

			# get total target + callable
			callable.bases[smp.idx,]$Sample.targeted <- sum(data.frame(target.smp.gr)$width);

			# get total target + callable in exons
			callable.bases[smp.idx,]$Sample.exon.targeted <- sum(data.frame(intersect(target.smp.gr, exon.gr))$width);

			# get total target + callable in CDS
			callable.bases[smp.idx,]$Sample.cds.targeted <- sum(data.frame(intersect(target.smp.gr, cds.gr))$width);
			}
		}

	# find per-patient metrics
	# this may be T/N pairs or T/T/N sets
	patient.idx <- which(callable.bases$Sample %in% these.smps);

	# make a genomic ranges object
	patient.gr <- GRanges(tmp[which(tmp$num == length(these.smps)),c('chrom','start','end')]);

	# get total callable
	callable.bases[patient.idx,]$Patient.total <- sum(data.frame(patient.gr)$width);

	# get total callable in exons
	callable.bases[patient.idx,]$Patient.exon <- sum(data.frame(intersect(patient.gr, exon.gr))$width);

	# get total callable in CDS
	callable.bases[patient.idx,]$Patient.cds <- sum(data.frame(intersect(patient.gr, cds.gr))$width);

	# get overlap with target regions too if required
	if ('TargetRegions' %in% colnames(tmp)) {

		target.patient.gr <- intersect(target.gr, patient.gr);

		# get total target + callable
		callable.bases[patient.idx,]$Patient.targeted <- sum(data.frame(target.patient.gr)$width);

		# get total target + callable in exons
		callable.bases[patient.idx,]$Patient.exon.targeted <- sum(data.frame(intersect(target.patient.gr, exon.gr))$width);

		# get total target + callable in CDS
		callable.bases[patient.idx,]$Patient.cds.targeted <- sum(data.frame(intersect(target.patient.gr, cds.gr))$width);
		}
	}

# save data to file
save(
	cov.list,
	file = generate.filename(arguments$project, 'CallableBases','RData')
	);

# save total bases to file
write.table(
	callable.bases,
	file = generate.filename(arguments$project, 'total_bases_covered','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectCallableBases','SessionProfile','txt'));
