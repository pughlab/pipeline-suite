### collect_fragmentomics_outputs.R ################################################################
# Collect outputs from pughlab_fragmentomics_pipeline

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

parser$add_argument('-p', '--project', type = 'character', help = 'project name', default = 'PROJECT');
parser$add_argument('-d', '--directory', type = 'character', help = 'path to data directory',
	default = getwd());

arguments <- parser$parse_args();

setwd(arguments$directory);

### READ DATA ######################################################################################
# find data files
score.files		<- list.files(pattern = 'score.txt', recursive = TRUE);
insertsize.files	<- list.files(pattern = 'picard.txt', recursive = TRUE);
alignment.files		<- list.files(pattern = 'alignmentMetrics.txt', recursive = TRUE);
fs.ratio.files		<- list.files(pattern = '5Mb_bins.txt', recursive = TRUE);
nucleosome.files	<- list.files(pattern = 'peak_distance.txt', recursive = TRUE);
griffin.files		<- list.files(pattern = 'all_sites.coverage.txt', recursive = TRUE);
# motifs = *_motifs.txt, *_raw.txt
motif5.files		<- list.files(pattern = '5prime_motifs.txt', recursive = TRUE);
motif3.files		<- list.files(pattern = '3prime_motifs.txt', recursive = TRUE);
# bkpts = *_count.txt, *_frequency.txt, *_ratio.txt
breakpoint.files	<- list.files(pattern = 'frequency.txt', recursive = TRUE);
# dinus = *_len150_contexts.txt, *_len167_contexts.txt (dinuc x position)
dinucleotide.long	<- list.files(pattern = 'len167_contexts.txt', recursive = TRUE);
dinucleotide.short	<- list.files(pattern = 'len150_contexts.txt', recursive = TRUE);

# read in fragment scores
if (length(score.files) > 0) {
	score.data <- do.call(rbind, lapply(score.files, function(i) { 
		tmp <- read.delim(i);
		tmp$Sample <- sub('_score.txt','',basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- gsub('_downsample','',tmp$Sample);
			}
		return(tmp);
		}));

	colnames(score.data)[1] <- 'FS';
	}

# read in insert sizes
if (length(insertsize.files) > 0) {
	fragment.length.data <- do.call(rbind, lapply(insertsize.files, function(i) { 
		tmp <- read.delim(i, skip = 6, nrow = 1);
		tmp$Sample <- sub('_picard.txt','',basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- gsub('_downsample','',tmp$Sample);
			}
		return(tmp);
		}));

	insertsizes.data <- do.call(rbind, lapply(insertsize.files, function(i) { 
		tmp <- read.delim(i, skip = 12);
		tmp$Sample <- sub('_picard.txt','',basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- gsub('_downsample','',tmp$Sample);
			}
		return(tmp);
		}));
	}

# read in alignment metrics
if (length(alignment.files) > 0) {

	qc.data <- do.call(rbind, lapply(alignment.files, function(i) {
		tmp <- read.delim(i, skip = 6);
		tmp$Sample <- sub('_alignmentMetrics.txt','',basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- gsub('_downsample','',tmp$Sample);
			}
		idx <- which(colnames(tmp) == 'SAMPLE')-1;
		tmp <- tmp[,c('Sample',colnames(tmp)[1:idx])];
		return(tmp);
		}));
	}

# read in fragment ratios
if (length(fs.ratio.files) > 0) {
	ratio.data <- do.call(rbind, lapply(fs.ratio.files, function(i) { 
		tmp <- read.delim(i);
		if (grepl('downsample', basename(i))) {
			tmp$sample_id <- sub('_downsample','',tmp$sample_id);
			}
	#	return(tmp[,c('sample_id','bin','seqnames','arm','start','end','combined_centered')]);
		return(tmp[,c('sample_id','bin','seqnames','arm','start','end','ratio_centered')]);
		}));
	}

# read in nucleosome positioning data
if (length(nucleosome.files) > 0) {
	nucleosome.data <- do.call(rbind, lapply(nucleosome.files, function(i) { 
		tmp <- read.delim(i);
		tmp$Sample <- sub('_peak_distance.txt', '', basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- sub('_downsample','',tmp$Sample);
			}
		return(tmp);
		}));
	}

# read in nucleosome accessibility data
if (length(griffin.files) > 0) {
	griffin.data <- do.call(rbind, lapply(griffin.files, function(i) {
		tmp <- read.delim(i);
		type <- basename(gsub('/coverage/all_sites','',dirname(i)));
		tmp$site_type <- type;
		return(tmp);
		}));

	if (any(grepl('downsample', griffin.data$sample))) {
		griffin.data$sample <- gsub('_downsample','',griffin.data$sample);
		}
	}

# read in end motif scores
if (length(motif5.files) > 0) {
	motif5.data <- do.call(rbind, lapply(motif5.files, function(i) { 
		tmp <- read.delim(i);
		tmp$Sample <- sub('_motifs.txt','',basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- gsub('_downsample','',tmp$Sample);
			}
		return(tmp);
		}));
	}
if (length(motif3.files) > 0) {
	motif3.data <- do.call(rbind, lapply(motif3.files, function(i) { 
		tmp <- read.delim(i);
		tmp$Sample <- sub('_motifs.txt','',basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- gsub('_downsample','',tmp$Sample);
			}
		return(tmp);
		}));
	}

# read in breakpoint data
if (length(breakpoint.files) > 0) {
	breakpoint.data <- do.call(rbind, lapply(breakpoint.files, function(i) { 
		tmp <- read.delim(i);
		tmp$Sample <- sub('_frequency.txt', '', basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- sub('_downsample','',tmp$Sample);
			}
		return(tmp);
		}));
	}

# read in dinucleotide data
if (length(dinucleotide.short) > 0) {
	dinucleotide.150 <- do.call(rbind, lapply(dinucleotide.short, function(i) { 
		tmp <- read.delim(i);
		tmp$Sample <- sub('_len150_contexts.txt', '', basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- sub('_downsample','',tmp$Sample);
			}
		return(tmp);
		}));
	}

if (length(dinucleotide.long) > 0) {
	dinucleotide.167 <- do.call(rbind, lapply(dinucleotide.long, function(i) { 
		tmp <- read.delim(i);
		tmp$Sample <- sub('_len167_contexts.txt', '', basename(i));
		if (grepl('downsample', basename(i))) {
			tmp$Sample <- sub('_downsample','',tmp$Sample);
			}
		return(tmp);
		}));
	}

### FORMAT DATA ####################################################################################
## format insert size data
if (length(insertsize.files) > 0) {

	# reshape fragment counts
	histogram.data <- reshape(insertsizes.data[,c('insert_size','All_Reads.fr_count','Sample')],
		direction = 'wide',
		idvar = 'insert_size',
		timevar = 'Sample'
		);
	colnames(histogram.data) <- gsub('All_Reads.fr_count.','',colnames(histogram.data));

	# calculate some size metrics
	short.frags <- apply(
		histogram.data[which(histogram.data$insert_size > 10 & histogram.data$insert_size < 150),-1],
		2, sum, na.rm = TRUE);
	normal.frags <- apply(
		histogram.data[which(histogram.data$insert_size >= 150 & histogram.data$insert_size < 220),-1],
		2, sum, na.rm = TRUE);
	long.frags <- apply(
		histogram.data[which(histogram.data$insert_size >= 220 & histogram.data$insert_size < 600),-1],
		2, sum, na.rm = TRUE);
	total.frags <- apply(
		histogram.data[which(histogram.data$insert_size > 10 & histogram.data$insert_size < 600),-1],
		2, sum, na.rm = TRUE);

	tmp.metrics <- data.frame(
		Sample = names(short.frags),
		Total.Fragments = total.frags,
		Total.Short = short.frags,
		Total.Normal = normal.frags,
		Total.Long = long.frags,
		Proportion.Short = (short.frags / total.frags),
		Proportion.Long = (long.frags / total.frags),
		Ratio.Short.Normal = (short.frags / normal.frags)
		);

	# add metrics to summary data
	fragment.size.summary <- merge(
		fragment.length.data[,c('Sample','MEDIAN_INSERT_SIZE','MEDIAN_ABSOLUTE_DEVIATION','MEAN_INSERT_SIZE','STANDARD_DEVIATION','READ_PAIRS')],
		tmp.metrics,
		by = 'Sample'
		);

	histogram.data[is.na(histogram.data)] <- 0;
	}

## format fragment ratio data
if (length(fs.ratio.files) > 0) {

	# reshape per-bin ratios
	ratio.data.wide <- reshape(
		ratio.data,
		direction = 'wide',
		idvar = c('bin','seqnames','arm','start','end'),
		timevar = 'sample_id'
		);
	#colnames(ratio.data.wide) <- gsub('combined_centered.','',colnames(ratio.data.wide));
	colnames(ratio.data.wide) <- gsub('ratio_centered.','',colnames(ratio.data.wide));
	}

## format griffin data
if (length(griffin.files) > 0) {

	# select only 1 entry per site
	if (length(unique(griffin.data$GC_correction)) == 2) {
		griffin.data <- griffin.data[which(griffin.data$GC_correction == 'GC_corrected'),];
		}
	colnames(griffin.data) <- sub('sample', 'Sample', colnames(griffin.data));

	# find distance fields	
	distance.fields <- colnames(griffin.data)[grepl('X',colnames(griffin.data))];

	# reshape data
	griffin.long <- reshape(
		griffin.data[,c(distance.fields,'Sample','site_name','site_type')],
		direction = 'long',
		varying = list(1:length(distance.fields)),
		v.names = 'Coverage',
		timevar = 'Position',
		times = as.numeric(gsub('X','', gsub('\\.','-',distance.fields)))
		);

	griffin.long <- griffin.long[,1:5];
	griffin.long <- griffin.long[order(griffin.long$Sample, griffin.long$site_type,
		griffin.long$site_name, griffin.long$Position),];
	}

## format end motif frequencies
if (length(motif5.files) > 0) {

	# reshape frequency data
	motifs5.wide <- reshape(
		motif5.data,
		direction = 'wide',
		idvar = 'motif',
		timevar = 'Sample'
		);
	colnames(motifs5.wide) <- gsub('frequency.','',colnames(motifs5.wide));
	}

if (length(motif3.files) > 0) {
	# reshape frequency data
	motifs3.wide <- reshape(
		motif3.data,
		direction = 'wide',
		idvar = 'motif',
		timevar = 'Sample'
		);
	colnames(motifs3.wide) <- gsub('frequency.','',colnames(motifs3.wide));
	}

## format breakpoint nucleotide frequencies
if (length(breakpoint.files) > 0) {

	breakpoint.data <- breakpoint.data[,c('Sample','nucleotide',
		colnames(breakpoint.data)[grepl('X',colnames(breakpoint.data))])];

	breakpoint.long <- reshape(
		breakpoint.data,
		direction = 'long',
		varying = list(3:ncol(breakpoint.data)),
		v.names = 'Frequency',
		timevar = 'Position',
		times = gsub('X','',gsub('X\\.','-',colnames(breakpoint.data)[3:ncol(breakpoint.data)]))
		);
	}

## format dinucleotide frequencies
if (length(dinucleotide.short) > 0) {

	# reshape short data
	dinucleotide.data.short <- reshape(
		dinucleotide.150,
		direction = 'wide',
		idvar = 'context',
		timevar = 'Sample'
		);

	colnames(dinucleotide.data.short) <- gsub('freq.','',colnames(dinucleotide.data.short));
	colnames(dinucleotide.150) <- gsub('freq','Frequency',colnames(dinucleotide.150));
	}

# reshape long data
if (length(dinucleotide.long) > 0) {

	dinucleotide.data.normal <- reshape(
		dinucleotide.167,
		direction = 'long',
		varying = list(grep('X',colnames(dinucleotide.167))),
		v.names = 'Frequency',
		timevar = 'Position'
		);

	colnames(dinucleotide.data.normal)[5] <- 'Group';
	dinucleotide.data.normal$Group <- '167bp';

	dinucleotide.167.mean <- aggregate(
		Frequency ~ Sample + context, dinucleotide.data.normal, mean);

	}

# combine
if ( (length(dinucleotide.short) > 0) & (length(dinucleotide.long) > 0) ) {
	average.dinucleotide <- merge(dinucleotide.150,dinucleotide.167.mean,
		by = c('Sample','context'), suffixes = c('.len150','.len167'));
	}

### SAVE DATA ######################################################################################
# save fragment scores
if (exists('score.data')) {

	write.table(
		score.data[,c('Sample','FS')],
		file = generate.filename(arguments$project, 'fragment_scores', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

# save insert size summary
if (exists('fragment.size.summary')) {

	write.table(
		fragment.size.summary,
		file = generate.filename(arguments$project, 'insert_size_summary', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);

	write.table(
		histogram.data,
		file = generate.filename(arguments$project, 'insert_size_counts', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

if (exists('qc.data')) {

	write.table(
		qc.data,
		file = generate.filename(arguments$project, 'alignmentMetrics', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

# save per-bin fragment ratios
if (exists('ratio.data.wide')) {

	write.table(
		ratio.data.wide,
		file = generate.filename(arguments$project, 'per5Mb_fragment_ratios', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

# save per-chromosome nucleosome position summary
if (exists('nucleosome.data')) {

	write.table(
		nucleosome.data[,c('Sample','distance',paste0('chr',1:22))],
		file = generate.filename(arguments$project, 'nucleosome_peak_distances', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

# save nucleosome accessibility data
if (exists('griffin.long')) {

	write.table(
		griffin.long,
		file = generate.filename(arguments$project, 'chromatin_accessibility_distances', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

# save end-motif frequencies
if (exists('motifs5.wide')) {

	write.table(
		motifs5.wide,
		file = generate.filename(arguments$project, '5prime_endmotif_frequencies', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

if (exists('motifs3.wide')) {
	write.table(
		motifs3.wide,
		file = generate.filename(arguments$project, '3prime_endmotif_frequencies', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

# save breakpoint frequencies
if (exists('breakpoint.long')) {

	write.table(
		breakpoint.long[,1:4],
		file = generate.filename(arguments$project, 'breakpoint_frequencies', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

# save dinucleotide frequencies
if (exists('dinucleotide.data.short')) {

	write.table(
		dinucleotide.data.short,
		file = generate.filename(arguments$project, 'dinucleotide_frequencies__len150', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

if (exists('dinucleotide.data.normal')) {

	write.table(
		dinucleotide.data.normal[,c('Sample','Group','context','Position','Frequency')],
		file = generate.filename(arguments$project, 'dinucleotide_frequencies__len167', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

if (exists('average.dinucleotide')) {

	write.table(
		average.dinucleotide,
		file = generate.filename(arguments$project, 'average_dinucleotide_frequencies', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectFragmentomics','SessionProfile','txt'));
