### compareVariantFragmentSizes.R ##################################################################
# Wrapper script to compare variant fragment size on a single sample.
# This function loads a processed BAM file, subsets reads to target positions, annotated each read
# 	as REF or ALT, performs summary statisitcs for each variant and outputs visualizations.

# Example:
# step 1 (find_ms): find microsatellite sites within target regions
#	Rscript find_ts_msi.R --step find_sites -t /path/to/target.bed -o /path/to/annotated_target.bed 

# step 2 (get_readcounts): collect readcounts for each microsatellite site that contain an INS/DEL/REF allele
#	Rscript find_ts_msi.R --step get_readcounts -t /path/to/annotated_target.bed -b /path/to/sample.bam -o /path/to/sample_readcounts.tsv

# step 3 (make_pon): aggregate readcounts across normal samples to generate a panel of normals
#	Rscript find_ts_msi.R --step make_pon -t /path/to/annotated_target.bed --pon_dir /path/to/pon/directory

# step 4 (final): calculate difference between tumour and normal (matched or PoN) using a prop.test
#	Rscript find_ts_msi.R --step final -t /path/to/annotated_target.bed -f /path/to/tumour_readcounts.tsv -n /path/to/normal_readcounts.tsv

### PREPARE SESSION ################################################################################
# import libraries
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('--step', type = 'character', help = 'one of find_sites, get_readcounts, make_pon or final');
parser$add_argument('-t', '--targets', type = 'character', help = 'target bed file with annotations');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output file');
parser$add_argument('-b', '--bam', type = 'character', help = 'path to sample BAM (required if step is get_readcounts)');
parser$add_argument('-f', '--tumour_freq', type = 'character', help = 'path to tumour read counts (required if step is final)');
parser$add_argument('-n', '--normal_freq', type = 'character', help = 'path to normal read counts (required if step is final)');
parser$add_argument('-p', '--pon_dir', type = 'character', help = 'path to normal readcount files (required if step is make_pon)');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'reference genome to use (required if step is find_ms)', default = 'hg38');

arguments <- parser$parse_args();

# do some error checking
arguments$step <- tolower(arguments$step);
if (!arguments$step %in% c('find_sites','get_readcounts','make_pon','final')) {
	stop('Argument --step is required and must be one of find_ms, get_readcounts, make_pon or final.');
	}

if ( ('find_ms' == arguments$step) & (!arguments$ref_type %in% c('hg38','hg19')) ) {
	stop('Unrecognized ref_type; must be one of hg38 or hg19');
	}

if ( (is.null(arguments$targets)) & ('make_pon' != arguments$step) ) { 
	stop('Argument --targets is required; please provide path to target bed file.');
	}

if (is.null(arguments$output)) { 
	stop('Argument --output is required; please provide path to output file appropriate for --step.');
	}

if ( ('get_readcounts' == arguments$step) & (is.null(arguments$bam)) ) { 
	stop('Argument --bam is required for this step; please provide path to BAM file.');
	}

if ( ('make_pon' == arguments$step) & (is.null(arguments$pon_dir)) ) { 
	stop('Argument --pon_dir is required for this step; please provide path to directory containing normal readcount files.');
	}

if ('final' == arguments$step) {
	if (is.null(arguments$normal_freq)) {
		stop('Argument --normal_freq is required for this step; please provide path to directory containing normal (or PoN) readcount file.');
		}

	if (is.null(arguments$tumour_freq)) {
		stop('Argument --tumour_freq is required for this step; please provide path to directory containing tumour readcount file.');
		}
	}

### FUNCTIONS ######################################################################################
# function to find mononucleotide repeats with length > 13 within each target region
run_find_ms_sites <- function(ref_type, targets) {

	# determine reference genome to use
	if ('hg38' == ref_type) {
		library(BSgenome.Hsapiens.UCSC.hg38);
		genomedb <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens;
		} else if ('hg19' == ref_type) {
		library(BSgenome.Hsapiens.UCSC.hg19);
		genomedb <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens;
		}

	# read in target regions
	target.bed <- read.delim(targets, header = FALSE);
	colnames(target.bed)[1:3] <- c('chrom','start','end');

	# add some padding
	target.bed$start <- target.bed$start - 20;
	target.bed$end <- target.bed$end + 20;

	# add fields to fill in
	target.bed[,c('repeat_unit','repeat_times','repeat_pos1','repeat_pos2')] <- NA;

	# add repeat annotations
	for (idx in 1:nrow(target.bed)) {

		ref <- as.character(Biostrings::DNAStringSet(IRanges::Views(genomedb, GenomicRanges::GRanges(target.bed[idx,1:3]))));

		repeats <- gregexpr("([a-zA-Z])\\1+", ref, perl = TRUE);

		df <- data.frame(
			Seq = regmatches(ref, repeats)[[1]],
			Length = attr(repeats[[1]],'match.length'),
			Pos1 = as.numeric(target.bed[idx,]$start + attr(repeats[[1]],'capture.start'))
			);
		df$Pos2 <- (df$Pos1 + df$Length - 1);
		df <- df[order(df$Length, decreasing = TRUE),];

		if (max(df$Length) < 13) { next; }

		target.bed[idx,c('repeat_unit','repeat_times','repeat_pos1','repeat_pos2')] <- c(
			substr(df[1,]$Seq,0,1),
			df[1,c('Length','Pos1','Pos2')]
			);
		}

	# drop regions with NAs (likely di- or tri-nucleotide repeats)
	target.bed <- target.bed[!is.na(target.bed$repeat_unit),];

	# give each region a unique ID
	target.bed$ID <- paste0('region',1:nrow(target.bed));

	return(target.bed);
	}

# function to calculate readcounts with INS/DEL/REF allele for each target region MS site
run_get_readcounts <- function(targets, filePath) {

	# read in target regions
	target.bed <- read.delim(targets);

	# find the BAM index
	indexed.bam <- gsub(".bam$", ".bai", filePath);
	if (!file.exists(indexed.bam)) { 
		indexed.bam <- gsub("$", ".bai", filePath);
		Rsamtools::indexBam(filePath);
		}

	# loop over each variant
	read.counts <- list();

	for (idx in 1:nrow(target.bed)) {

		repeat.gr <- GenomicRanges::makeGRangesFromDataFrame(
			target.bed[idx,c('chrom','repeat_pos1','repeat_pos2')],
			start.field = 'repeat_pos1',
			end.field = 'repeat_pos2'
			);

		# indicate parameters to use
		param <- Rsamtools::ScanBamParam(
			mapqFilter = 30,
			which = repeat.gr,
			what = 'seq',
			flag = Rsamtools::scanBamFlag(
				isPaired = TRUE,
				isProperPair = TRUE,
				isDuplicate = FALSE,
				isUnmappedQuery = FALSE
				)
			);

		# create genomic alignments object
		ga <- IRanges::mergeByOverlaps(
			repeat.gr, 
			GenomicAlignments::readGAlignments(filePath, param = param, index = indexed.bam),
			type = 'within'
			);

		if (nrow(ga) == 0) { next; }

		ga$RelativeLength <- NA;

		for (j in 1:nrow(ga)) {

			# expand the dna sequence
			dna.seq <- as.character(GenomicAlignments::sequenceLayer(
				Biostrings::DNAStringSet(ga[j,]$seq),
				ga[j,]$cigar)[[1]]);

			# extract repeat information
			pos1 <- target.bed[idx,'repeat_pos1'] - ga[j,]$start + 1;
			pos2 <- target.bed[idx,'repeat_pos2'] - ga[j,]$start + 1;
			allele <- substr(dna.seq, pos1, pos2);

			# find repeat
			repeats <- gregexpr(paste0("([",target.bed[idx,]$repeat_unit,"])\\1+"),dna.seq);

			df <- data.frame(
				Seq = unlist(regmatches(dna.seq, repeats)),
				Length = attr(repeats[[1]],'match.length'),
				Pos1 = as.numeric(repeats[[1]]),
				Diff = abs(repeats[[1]] - pos1)
				);

			# find the (longest) sequence closest to the expected position
			this.seq <- df[order(df$Diff, -df$Length),]$Length[1];
			ga[j,]$RelativeLength <- this.seq - target.bed[idx,]$repeat_times;

			}

		if (!all(is.na(ga$RelativeLength))) {
			read.counts[[idx]] <- data.frame(table(ga$RelativeLength));
			colnames(read.counts[[idx]]) <- c('RelativeLength','Frequency');
			read.counts[[idx]]$Allele <- target.bed[idx,]$ID;
			}
		}

	frequency.data <- reshape(
		do.call(rbind, read.counts),
		direction = 'wide',
		idvar = 'RelativeLength',
		timevar = 'Allele'
		);

	colnames(frequency.data) <- gsub('Frequency.','',colnames(frequency.data));

	return(frequency.data);
	}

# function to find average normal readcounts for each target region MS site
run_make_panel_of_normals <- function(directory) {

	# read in data
	normal.counts <- do.call(rbind, lapply(list.files(path = directory, pattern = 'readcounts.tsv', full.names = TRUE), function(i) { 

		smp <- unlist(strsplit(basename(i),'_'))[1];
		normal.data <- read.delim(i);

		# format sample data
		tmp <- aggregate(normal.data[,-1],
			by = list(Status = sign(normal.data$RelativeLength)), FUN = sum, na.rm = TRUE);

		counts <- data.frame(cbind(
			apply(tmp[which(tmp$Status == -1),-1],2,sum, na.rm = TRUE),
			apply(tmp[which(tmp$Status == 0),-1],2,sum, na.rm = TRUE),
			apply(tmp[which(tmp$Status == 1),-1],2,sum, na.rm = TRUE)
			));
		colnames(counts) <- c('-1','0','1');
		counts$Region <- rownames(counts);
		counts$Sample <- smp;

		return(counts);
		}));

	tmp <- aggregate(normal.counts[,c('-1','0','1')], by = list(Region = normal.counts$Region), median);
	panelOfNormals <- data.frame(t(tmp));
	rownames(panelOfNormals)[1] <- 'RelativeLength';

	return(panelOfNormals);
	}

# function to compare tumour and normal readcounts for each target region MS site
run_compare_readcounts <- function(targets, tumour, normal) {

	# read in data
	results.data <- read.delim(targets);
	results.data[,c('tumour_prop','normal_prop','p.value')] <- NA;

	normal.data <- read.delim(normal);
	tumour.data <- read.delim(tumour);

	# format sample data
	normal_counts <- aggregate(normal.data[,-1],
		by = list(Status = sign(normal.data$RelativeLength)), FUN = sum, na.rm = TRUE);
	tumour_counts <- aggregate(tumour.data[,-1],
		by = list(Status = sign(tumour.data$RelativeLength)), FUN = sum, na.rm = TRUE);

	# compare proportion of MS-deletions between tumour and normal
	for (i in 1:nrow(results.data)) {

		id <- results.data[i,]$ID;

		if ( (!id %in% colnames(normal_counts)) | (!id %in% colnames(tumour_counts)) ) { next; }
		if (sum(normal_counts[which(normal_counts$Status < 1),id]) == 0) { next; }
		if (sum(tumour_counts[which(tumour_counts$Status < 1),id]) == 0) { next; }

		results.data[i,]$tumour_prop <- tumour_counts[which(tumour_counts$Status == -1),id] / 
			sum(tumour_counts[which(tumour_counts$Status < 1),id]);
		results.data[i,]$normal_prop <- normal_counts[which(normal_counts$Status == -1),id] / 
			sum(normal_counts[which(normal_counts$Status < 1),id]);

		results.data[i,]$p.value <- prop.test(
			c(
				tumour_counts[which(tumour_counts$Status == -1),id],
				normal_counts[which(normal_counts$Status == -1),id]
				),
			c(
				sum(tumour_counts[which(tumour_counts$Status < 1),id]),
				sum(normal_counts[which(normal_counts$Status < 1),id])
				),
			alternative = 'g'
			)$p.value;
		}

	return(results.data);
	}


### MAIN ###########################################################################################
if ('find_sites' == arguments$step) {

	# find mononucleotide repeats with length > 13 within each target region
	target.bed.annotated <-	run_find_ms_sites(
		ref_type = arguments$ref_type,
		targets = arguments$targets
		);

	# save to file
	write.table(
		target.bed.annotated,
		file = arguments$output,
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);

	} else if ('get_readcounts' == arguments$step) {

	# calculate readcounts with INS/DEL/REF allele for each target region MS site
	sample.data <- run_get_readcounts(
		targets = arguments$targets,
		filePath = arguments$bam
		);

	# save to file
	write.table(
		sample.data,
		file = arguments$output,
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);

	} else if ('make_pon' == arguments$step) {

	# calculate median readcounts across 'normal' samples
	sample.data <- run_make_panel_of_normals(
		directory = arguments$pon_dir
		);

	# save to file
	write.table(
		sample.data,
		file = arguments$output,
		row.names = TRUE,
		col.names = FALSE,
		quote = FALSE,
		sep = '\t'
		);

	} else if ('final' == arguments$step) {

	# compare tumour and normal readcounts for each target region MS site
	results.data <- run_compare_readcounts(
		targets = arguments$targets,
		tumour = arguments$tumour_freq,
		normal = arguments$normal_freq
		);	

	# save to file
	write.table(
		results.data,
		file = arguments$output,
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

