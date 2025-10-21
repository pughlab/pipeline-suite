### summarize_mutation_data.R ######################################################################
# Summarize mutation characteristics (TMB, basechange frequency, functional consequences, etc).

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

# function to trim sample IDs
simplify.ids <- function(x) {
	match <- TRUE;
	if (length(x) == 1) {
		match <- FALSE;
		new.ids <- x;
		}
	index <- 1;
	while (match) {
		if (length(unique(sapply(x,function(i) { unlist(strsplit(i,''))[index] } ))) == 1) {
			index <- index + 1;
			} else {
			new.ids <- sapply(x,function(i) { substring(i,index,nchar(i)) } );
			match <- FALSE;
			}
		}
	return(new.ids);
	}

### VARIANT CODING
# 1 = missense, 2 = stop gain, 3 = stop loss, 4 = splicing, 5 = frameshift, 6 = in frame indel, 7 = tss
# 8 = RNA, 9 = other (up/downstream, UTR, intergenic, silent, intron), 10 = ITD
variant.codes <- data.frame(
	Classification = c("3'Flank", "5'Flank", "Intron", "RNA", "IGR", "3'UTR", "5'UTR", "Silent",
		"Missense_Mutation", "Splice_Region", "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
		"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation",
		"Translation_Start_Site", "ITD"),
	Group = c('other','other','other','RNA','other','other','other','other','missense',
		'splice_site','splice_site','in_frame_indel','in_frame_indel','frameshift_indel',
		'frameshirt_indel', 'nonsense', 'nonstop', 'tss', 'itd'),
	Code = c(9, 9, 9, 8, 9, 9, 9, 9, 1, 4, 4, 6, 6, 5, 5, 2, 3, 7, 10)
	);

### PREPARE SESSION ################################################################################
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('--maf', type = 'character', help = 'combined mutation calls in MAF format');
parser$add_argument('--msi', type = 'character', help = 'path to msi file');
parser$add_argument('--chord', type = 'character', help = 'path to CHORD output');
parser$add_argument('--hrdetect', type = 'character', help = 'path to HRDetect output');
parser$add_argument('-z', '--report', type = 'character', help = 'path to report directory',
	default = NULL);

arguments <- parser$parse_args();

# import libraries
library(BoutrosLab.plotting.general);
library(xtable);

# set some colour schemes
variant.colours <- c('darkseagreen4','darkorchid4','#9AA3F2','yellow','darkorange3','#F9B38E','turquoise1','plum','grey50')
names(variant.colours) <- c('missense','nonsense','nonstop','splicing','frameshift_indel','in_frame_indel','tss','RNA','noncoding');

# for these plots, we will ignore some variant types
variant.colours <- variant.colours[c(1:7,9)];
basechange.colours <- default.colours(8,'pastel')[-5];

### READ DATA ######################################################################################
# get data
if (is.null(arguments$maf)) {
	stop('ERROR: No input MAF provided, please provide path to SNV calls in MAF format.');
	} else {
	input.data <- read.delim(arguments$maf, stringsAsFactors = FALSE, comment.char = '#');
	}

if (is.null(arguments$msi)) {
	msi <- NULL; } else {
	msi <- read.delim(arguments$msi, row.names = 1);
	}

if (is.null(arguments$chord)) {
	chord <- NULL; } else {
	chord <- read.delim(arguments$chord, row.names = 1);
	}

if (is.null(arguments$hrdetect)) {
	hrdetect <- NULL; } else {
	hrdetect <- read.delim(arguments$hrdetect, row.names = 1);
	}

# collect list of all samples
all.samples <- sort(as.character(unique(input.data$Tumor_Sample_Barcode)));

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# calculate VAF
input.data$t_vaf <- input.data$t_alt_count / input.data$t_depth;

# indicate key fields
keep.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele2','t_vaf');

print(paste0("Total variants: ", nrow(input.data), " among ", length(all.samples), " samples."));

# add quick recurrence filter (likely false positives/artefacts)
recurrence.data <- aggregate(
	Tumor_Sample_Barcode ~ Chromosome + Start_Position + End_Position + Allele,
	input.data,
	length
	);
colnames(recurrence.data)[ncol(recurrence.data)] <- 'RecurrenceCount';

recurrence.threshold <- floor(length(all.samples)*0.9);
tmp <- merge(input.data, recurrence.data, all.x = TRUE);
tmp <- tmp[which(tmp$RecurrenceCount < recurrence.threshold),];

print(paste0("Removing ",
	nrow(recurrence.data[which(recurrence.data$RecurrenceCount > recurrence.threshold),]),
	" variants due to too high recurrence (>= 90% [n = ",
	recurrence.threshold,
	"] of samples)."
	));

mutation.data <- tmp[,keep.fields];

rm(tmp,recurrence.data);


## apply variant coding
mutation.data$Code <- variant.codes$Code[match(mutation.data$Variant_Classification, variant.codes$Classification)];
# drop variants classified as 'RNA'
mutation.data <- mutation.data[which(mutation.data$Code != 8),];
	
# classify each basechange
mutation.data$Basechange <- paste0(mutation.data$Reference_Allele, '>', mutation.data$Tumor_Seq_Allele2);
if (any(mutation.data$Variant_Type %in% c('INS','DEL'))) {
	mutation.data[which(mutation.data$Variant_Type %in% c('DEL', 'INS')),]$Basechange <- 'indel';
	}
if (any(mutation.data$Basechange == 'T>G')) {
	mutation.data[which(mutation.data$Basechange == 'T>G'),]$Basechange <- 'A>C';
	}
if (any(mutation.data$Basechange == 'T>C')) {
	mutation.data[which(mutation.data$Basechange == 'T>C'),]$Basechange <- 'A>G';
	}
if (any(mutation.data$Basechange == 'T>A')) {
	mutation.data[which(mutation.data$Basechange == 'T>A'),]$Basechange <- 'A>T';
	}
if (any(mutation.data$Basechange == 'C>T')) {
	mutation.data[which(mutation.data$Basechange == 'C>T'),]$Basechange <- 'G>A';
	}
if (any(mutation.data$Basechange == 'C>G')) {
	mutation.data[which(mutation.data$Basechange == 'C>G'),]$Basechange <- 'G>C';
	}
if (any(mutation.data$Basechange == 'C>A')) {
	mutation.data[which(mutation.data$Basechange == 'C>A'),]$Basechange <- 'G>T';
	}
mutation.data$Basechange <- factor(
	mutation.data$Basechange,
	levels = rev(c('A>C','A>G','A>T','G>A','G>C','G>T','indel'))
	);

# group data by VAF (for count ~ VAF | type barplot/histogram)
mutation.data$VAF_GROUP <- NA;
for (vaf in seq(0.05,1,0.01)) {
	idx <- which(mutation.data$t_vaf >= vaf-0.01 & mutation.data$t_vaf < vaf);
	if (length(idx) > 0) { mutation.data[idx,]$VAF_GROUP <- vaf; }
	}

# summarize mutation counts by Sample
mutation.summary <- list();

tmp <- aggregate(
	t_vaf ~ Tumor_Sample_Barcode,
	mutation.data,
	length
	);
colnames(tmp)[2] <- 'Count';

mutation.summary$overall <- tmp;

# summarize mutation counts by VAF
tmp <- aggregate(
	t_vaf ~ Tumor_Sample_Barcode + VAF_GROUP,
	mutation.data,
	length
	);
colnames(tmp)[3] <- 'Count';
tmp$VAF_GROUP <- factor(tmp$VAF_GROUP, levels = seq(0.05,1,0.01));

# fill in any zeros
for (smp in all.samples) {
	smp.idx <- which(tmp$Tumor_Sample_Barcode == smp);
	for (vaf in seq(0.05,1,0.01)) {
		vaf.idx <- which(tmp$VAF_GROUP == vaf);
		if (length(intersect(smp.idx, vaf.idx)) == 0) {
			tmp <- rbind(tmp,
				data.frame(
					Tumor_Sample_Barcode = smp, VAF_GROUP = vaf, Count = 0
					)
				);
			}
		}
	}

mutation.summary$per_vaf <- tmp;

rm(tmp);

# summarize mutations by functional consequence
functional.summary <- list();

tmp1 <- data.frame(table(mutation.data[,c('Tumor_Sample_Barcode','Code')]));
tmp2 <- aggregate(Freq ~ Tumor_Sample_Barcode,tmp1,sum);
tmp3 <- merge(tmp1, tmp2, by = 'Tumor_Sample_Barcode');
tmp3$Proportion <- tmp3$Freq.x / tmp3$Freq.y;
tmp3$Tumor_Sample_Barcode <- factor(tmp3$Tumor_Sample_Barcode,levels = all.samples);
tmp3$Code <- factor(tmp3$Code, levels = rev(c(1:7,9)));

functional.summary$overall <- tmp3;

rm(tmp1, tmp2, tmp3);

# summarize mutation functions by VAF
tmp1 <- data.frame(table(
	mutation.data[,c('Tumor_Sample_Barcode','VAF_GROUP','Code')]));
tmp2 <- aggregate(Freq ~ Tumor_Sample_Barcode + VAF_GROUP, tmp1, sum);
tmp3 <- merge(tmp1, tmp2, by = c('Tumor_Sample_Barcode','VAF_GROUP'));
tmp3$Proportion <- tmp3$Freq.x / tmp3$Freq.y;
tmp3$Code <- factor(tmp3$Code, levels = c(1:7,9));
tmp3$VAF_GROUP <- factor(tmp3$VAF_GROUP, levels = seq(0.05,1,0.01));
colnames(tmp3)[which(colnames(tmp3) == 'Freq.x')] <- 'Count';
colnames(tmp3)[which(colnames(tmp3) == 'Freq.y')] <- 'Total';

# fill in any zeros
for (smp in all.samples) {
	smp.idx <- which(tmp3$Tumor_Sample_Barcode == smp);
	for (vaf in seq(0.05,1,0.01)) {
		vaf.idx <- which(tmp3$VAF_GROUP == vaf);
		if (length(intersect(smp.idx, vaf.idx)) == 0) {
			tmp3 <- rbind(tmp3,
				data.frame(
					Tumor_Sample_Barcode = smp,
					VAF_GROUP = vaf, Code = 1, Count = 0, Total = 0, Proportion = 0
					)
				);
			}
		}
	}

functional.summary$per_vaf <- tmp3;

rm(tmp1, tmp2, tmp3, smp.idx, vaf.idx, smp, vaf);

# summarize mutations by basechanges
basechange.summary <- list();

tmp1 <- data.frame(table(mutation.data[,c('Tumor_Sample_Barcode','Basechange')]));
tmp2 <- aggregate(Freq ~ Tumor_Sample_Barcode, tmp1, sum);
tmp3 <- merge(tmp1, tmp2, by = 'Tumor_Sample_Barcode');
tmp3$Proportion <- tmp3$Freq.x / tmp3$Freq.y;
tmp3$Tumor_Sample_Barcode <- factor(tmp3$Tumor_Sample_Barcode, levels = all.samples);

basechange.summary$overall <- tmp3;

rm(tmp1, tmp2, tmp3);

# summarize basechange by VAF
tmp1 <- data.frame(table(
	mutation.data[,c('Tumor_Sample_Barcode','VAF_GROUP','Basechange')]));
tmp2 <- aggregate(Freq ~ Tumor_Sample_Barcode + VAF_GROUP, tmp1, sum);
tmp3 <- merge(tmp1, tmp2, by = c('Tumor_Sample_Barcode','VAF_GROUP'));
tmp3$Proportion <- tmp3$Freq.x / tmp3$Freq.y;
tmp3$VAF_GROUP <- factor(tmp3$VAF_GROUP, levels = seq(0.05,1,0.01));
colnames(tmp3)[which(colnames(tmp3) == 'Freq.x')] <- 'Count';
colnames(tmp3)[which(colnames(tmp3) == 'Freq.y')] <- 'Total';

# fill in any zeros
for (smp in all.samples) {
	smp.idx <- which(tmp3$Tumor_Sample_Barcode == smp);
	for (vaf in seq(0.05,1,0.01)) {
		vaf.idx <- which(tmp3$VAF_GROUP == vaf);
		if (length(intersect(smp.idx, vaf.idx)) == 0) {
			tmp3 <- rbind(tmp3,
				data.frame(
					Tumor_Sample_Barcode = smp,
					VAF_GROUP = vaf, Basechange = 'indel',
					Count = 0, Total = 0, Proportion = 0
					)
				);
			}
		}
	}

basechange.summary$per_vaf <- tmp3;

rm(tmp1, tmp2, tmp3, smp.idx, vaf.idx, smp, vaf);

# format mutation signature (msi, chord, hrdetect) calls
if (!is.null(msi)) {
	msi$MSI <- 0;
	if (any(msi$Percent > 10)) { msi[which(msi$Percent > 10),]$MSI <- 1; }
	if (any(msi$Percent > 20)) { msi[which(msi$Percent > 20),]$MSI <- 2; }
	msi.smps <- rownames(msi);
	} else {
	msi <- data.frame(
		matrix(data = 0, nrow = length(all.samples), ncol = 2,
		dimnames = list(all.samples, c('Percent','MSI')))
		);
	msi.smps <- NULL;
	}

if (!is.null(chord)) {
	chord$CHORD <- 0;
	if (any(grepl('deficient', chord$hr_status))) {
		chord[which(chord$hr_status == 'HR_deficient'),]$CHORD <- 2; }
	if (all(grepl('cannot_be_determined', chord$hr_status))) {
		chord.smps <- NULL;
		} else { 
		chord.smps <- rownames(chord[which(chord$hr_status != 'cannot_be_determined'),]);
		}
	} else {
	chord <- data.frame(
		matrix(data = 0, nrow = length(all.samples), ncol = 2,
			dimnames = list(all.samples, c('hr_status','CHORD')))
		);
	chord.smps <- NULL;
	}

if (!is.null(hrdetect)) {
	hrdetect$HRDetect <- 0;
	if (any(hrdetect$Probability >= 0.7)) { hrdetect[which(hrdetect$Probability >= 0.7),]$HRDetect <- 2; }
	hrdetect.smps <- rownames(hrdetect);
	} else {
	hrdetect <- data.frame(
		matrix(data = 0, nrow = length(all.samples), ncol = 2,
		dimnames = list(all.samples, c('Probability','HRDetect')))
		);
	hrdetect.smps <- NULL;
	}

sig.data <- list();
sig.data$bg <- merge(
	msi[,c('Percent','MSI')],
	merge(
		chord[,c('hr_status','CHORD')],
		hrdetect[,c('Probability','HRDetect')],
		by = 'row.names',
		all = TRUE
		),
	by.x = 'row.names',
	by.y = 'Row.names',
	all = TRUE
	);

rownames(sig.data$bg) <- sig.data$bg[,1];
sig.data$bg <- sig.data$bg[,c('MSI','CHORD','HRDetect')];

sig.data$call <- data.frame(matrix(data = 1, nrow = length(all.samples), ncol = 3,
	dimnames = list(all.samples, c('MSI','CHORD','HRDetect'))));

if (!is.null(msi.smps)) { sig.data$call[msi.smps,]$MSI <- 0; }
if (!is.null(chord.smps)) { sig.data$call[chord.smps,]$CHORD <- 0; }
if (!is.null(hrdetect.smps)) { sig.data$call[hrdetect.smps,]$HRDetect <- 0; }

sig.data$call[which(is.na(sig.data$call),arr.ind = TRUE)] <- 1;


### MAKE LEGENDS ###################################################################################
# make the plot legend (mutation type/consequence)
functional.legend <- legend.grob(
	legends = list(
		legend = list(
			colours = variant.colours,
			labels = names(variant.colours),
			title = 'Functional Consequence'
			)
		),
	title.cex = 0.7,
	title.fontface = 'plain',
	title.just = 'left',
	label.cex = 0.7,
	size = 1.5
	);

functional.legend.split <- list(
	rect = list(col = variant.colours, size = 1),
	text = list(lab = names(variant.colours), cex = 0.8),
	columns = length(variant.colours),
	between = 1,
	between.columns = 1,
	title = expression('Functional Consequence'),
	cex.title = 0.9
	);

basechange.legend <- legend.grob(
	legends = list(
		legend = list(
			colours = basechange.colours,
			labels = c(
				expression('A' %->% 'C,T' %->% 'G'),
				expression('A' %->% 'G,T' %->% 'C'),
				expression('A' %->% 'T,T' %->% 'A'),
				expression('G' %->% 'A,C' %->% 'T'),
				expression('G' %->% 'C,C' %->% 'G'),
				expression('G' %->% 'T,C' %->% 'A'),
				'indel'
				),
			title = 'Basechange'
			)
		),
	title.cex = 0.7,
	title.fontface = 'plain',
	title.just = 'left',
	label.cex = 0.7,
	size = 1.5
	);

basechange.legend.split <- list(
	rect = list(col = basechange.colours, size = 1),
	text = list(
		lab = c(
			expression('A' %->% 'C,T' %->% 'G'),
			expression('A' %->% 'G,T' %->% 'C'),
			expression('A' %->% 'T,T' %->% 'A'),
			expression('G' %->% 'A,C' %->% 'T'),
			expression('G' %->% 'C,C' %->% 'G'),
			expression('G' %->% 'T,C' %->% 'A'),
			'indel'
			),
		cex = 0.8
		),
	columns = length(basechange.colours),
	between = 1,
	between.columns = 1,
	title = expression('Basechange'),
	cex.title = 0.9
	);

sig.legend <- legend.grob(
	legends = list(
		legend = list(
			colours = c('white','grey50','black'),
			labels = c('MSS','MSI-L','MSI-H')
			),
		legend = list(
			colours = c('white','black','white'),
			labels = c('HR-proficient','HR-deficient','(X) undetermined')
			)
		),
	label.cex = 0.7,
	size = 1.5
	);

### PLOT RESULTS ###################################################################################
# create plots for each sample
for (smp in all.samples) {

	if (!dir.exists(smp)) {
		dir.create(smp);
		}

	setwd(smp);

	count.data <- mutation.summary$per_vaf;
	count.data <- count.data[which(count.data$Tumor_Sample_Barcode == smp),];

	functional.data <- functional.summary$per_vaf;
	functional.data <- functional.data[which(functional.data$Tumor_Sample_Barcode == smp),];
	if (any(is.na(functional.data$Proportion))) {
		functional.data[is.na(functional.data$Proportion),]$Proportion <- 0;
		}

	basechange.data <- basechange.summary$per_vaf;
	basechange.data <- basechange.data[which(basechange.data$Tumor_Sample_Barcode == smp),];
	if (any(is.na(basechange.data$Proportion))) {
		basechange.data[is.na(basechange.data$Proportion),]$Proportion <- 0;
		}

	# create plot for mutaiton counts
	mutation.count.plot <- create.barplot(
		Count ~ VAF_GROUP,
		count.data,
		xlimits = c(0.5,96.5),
		xat = c(1,seq(6,96,10)),
		xaxis.lab = c(0.05,seq(0.1,1,0.1)),
		yaxis.tck = c(0.5,0),
		xaxis.tck = c(0.5,0),
		ylab.label = 'Count',
		ylab.cex = 1.2,
		xlab.label = '',
		yaxis.cex = 1,
		xaxis.cex = 1,
		style = 'Nature'
		);

	# create plot for functional summary
	function.plot <- create.barplot(
		Proportion ~ VAF_GROUP,
		functional.data,
		groups = functional.data$Code,
		stack = TRUE,
		col = variant.colours,
		xlimits = c(0.5,96.5),
		xat = c(1,seq(6,96,10)),
		xaxis.lab = c(0.05,seq(0.1,1,0.1)),
		yaxis.tck = c(0.5,0),
		xaxis.tck = c(0.5,0),
		ylab.label = 'Proportion',
		ylab.cex = 1.2,
		xlab.label = '',
		ylimits = c(0,1),
		yat = seq(0,1,0.5),
		yaxis.cex = 1,
		xaxis.cex = 1,
		legend = list(
			inside = list(fun = draw.key,
				args = list(key = functional.legend.split), x = 0.1, y = 1.5)
			),
		style = 'Nature'
		);

	# create plot for basechange summary
	basechange.plot <- create.barplot(
		Proportion ~ VAF_GROUP,
		basechange.data,
		groups = basechange.data$Basechange,
		stack = TRUE,
		col = rev(basechange.colours),
		xlimits = c(0.5,96.5),
		xat = c(1,seq(6,96,10)),
		xaxis.lab = c(0.05,seq(0.1,1,0.1)),
		yaxis.tck = c(0.5,0),
		xaxis.tck = c(0.5,0),
		ylab.label = 'Proportion',
		ylab.cex = 1.2,
		xlab.label = expression('Variant Allele Frequency (x'['0']*', x'['1']*']'),
		xlab.cex = 1.2,
		ylimits = c(0,1),
		yat = seq(0,1,0.5),
		yaxis.cex = 1,
		xaxis.cex = 1,
		legend = list(
			inside = list(fun = draw.key,
				args = list(key = basechange.legend.split), x = 0.1, y = 1.5)
			),
		style = 'Nature'
		);

	# combine them
	create.multipanelplot(
		plot.objects = list(mutation.count.plot, function.plot, basechange.plot),
		plot.objects.heights = c(1,1,1.2),
		left.legend.padding = 0,
		right.legend.padding = 0,
		top.legend.padding = 0,
		bottom.legend.padding = 0,
		top.padding = -1,
		y.spacing = c(2,2),
		ylab.axis.padding = 3,
		xlab.axis.padding = 1,
		height = 6,
		width = 11, 
		resolution = 200,
		filename = generate.filename(smp, 'mutation_summary_by_VAF','png'),
		);

	setwd(arguments$output);
	}

# grab some parameters
axis.cex <- if (length(all.samples) <= 30) { 1
	} else if (length(all.samples) <= 50) { 0.75
	} else if (length(all.samples) <= 80) { 0.5
	} else { 0 };

# create plot for mutation signatures (msi, chord, hrdetect)
msi.plot <- create.dotmap(
	x = t(sig.data$call),
	bg.data = t(sig.data$bg),
	pch = 4,
	spot.size.function = function(x) { x },
	spot.colour.function = function(x) { c('black','transparent')[match(x, c(1,0))] },
	xaxis.lab = rep('',nrow(sig.data$call)),
	xaxis.tck = 0,
	yaxis.tck = c(0.5,0),
	yaxis.fontface = 'plain',
	yaxis.cex = 0.8,
	colour.scheme = c('white','grey50','black'),
	total.colours = 4,
	at = seq(0,2,1),
	bg.alpha = 1,
	lwd = 1, row.lwd = c(1.5,1), col.lwd = 1,
	row.colour = c('black','grey90'),
	col.colour = 'grey90', 
	legend = list(
		right = list(fun = sig.legend)
		)
	);

# create plot for functional summary
function.plot <- create.barplot(
	Proportion ~ Tumor_Sample_Barcode,
	functional.summary$overall,
	groups = functional.summary$overall$Code,
	stack = TRUE,
	col = rev(variant.colours),
	xaxis.lab = rep('',length(all.samples)),
	xlimits = c(0.5, length(all.samples)+0.5),
	xat = seq(1,length(all.samples)),
	yaxis.tck = c(0.5,0),
	xaxis.tck = 0,
	ylab.label = "\t\t\tProportion",
	ylab.cex = 1.5,
	xlab.label = NULL,
	ylimits = c(0,1),
	yat = seq(0,1,0.5),
	yaxis.cex = 1,
	style = 'Nature'
	);

# zoom in/exclude non-coding mutations
min.noncoding <- round(
	min(functional.summary$overall[which(functional.summary$overall$Code == 9),]$Proportion),3);

if (min.noncoding > 0.995) { # typical for WGS
	ylimits <- c(0.995, 1);
	yat <- seq(0.995, 1, 0.001);
	} else if (min.noncoding > 0.99) {
	ylimits <- c(0.99, 1);
	yat <- seq(0.99, 1, 0.002);
	} else if (min.noncoding > 0.98) {
	ylimits <- c(0.98, 1);
	yat <- seq(0.98, 1, 0.005);
	} else if (min.noncoding > 0.95) {
	ylimits <- c(0.95, 1);
	yat <- seq(0.95, 1, 0.01);
	} else if (min.noncoding > 0.9) {
	ylimits <- c(0.9, 1);
	yat <- seq(0.9, 1, 0.02);
	} else {
	ylimits <- c(floor(min.noncoding*10)/10, 1);
	yat <- seq(0,1,0.1);
	}

function.zoom.plot <- create.barplot(
	Proportion ~ Tumor_Sample_Barcode,
	functional.summary$overall,
	groups = functional.summary$overall$Code,
	stack = TRUE,
	col = rev(variant.colours),
	xaxis.lab = rep('',length(all.samples)),
	xlimits = c(0.5, length(all.samples)+0.5),
	xat = seq(1,length(all.samples)),
	yaxis.tck = c(0.5,0),
	xaxis.tck = 0,
	ylab.label = NULL,
	xlab.label = NULL,
	ylimits = ylimits,
	yat = yat,
	yaxis.cex = 1,
	style = 'Nature'
	);

# create plot for basechange summary
basechange.plot <- create.barplot(
	Proportion ~ Tumor_Sample_Barcode,
	basechange.summary$overall,
	groups = basechange.summary$overall$Basechange,
	stack = TRUE,
	col = rev(basechange.colours),
	xaxis.lab = simplify.ids(all.samples),
	xaxis.rot = if (length(all.samples) == 1) { 0 } else { 90 },
	xlimits = c(0.5, length(all.samples)+0.5),
	xat = seq(1,length(all.samples)),
	yaxis.tck = c(0.5,0),
	xaxis.tck = 0,
	ylab.label = 'Proportion',
	ylab.cex = 1.5,
	xlab.label = NULL,
	ylimits = c(0,1),
	yat = seq(0,1,0.5),
	yaxis.cex = 1,
	xaxis.cex = axis.cex,
	style = 'Nature',
	legend = list(
		right = list(fun = basechange.legend)
		)
	);

# combine them!
combined.plot <- create.multipanelplot(
	plot.objects = list(msi.plot, function.zoom.plot, function.plot, basechange.plot),
	plot.objects.heights = if (length(all.samples) <= 10) { c(0.6,1,1,2)
		 } else { c(0.6,1,1,2) },
	left.legend.padding = 0,
	right.legend.padding = 0,
	top.legend.padding = 0,
	bottom.legend.padding = 0,
	y.spacing =  if (length(all.samples) <= 10) { c(0.5,-1,0.5) } else { c(-1.3,-1.8,-1.8) },
	ylab.axis.padding = 0.8,
	right.padding = 3,
	legend = list(
		inside = list(
			fun = functional.legend,
			y = 0.65,
			x = if (length(all.samples) <= 10) { 0.815 } else { 0.863 }
			)
		),
	style = 'Nature'
	);

write.plot(
	combined.plot,
	height = 7,
	width = if (length(all.samples) <= 10) { 6 } else { 8 },
	resolution = 200,
	filename = generate.filename(arguments$project, 'mutation_summary','png'),
	);

# save key data to file
save(
	functional.summary,
	basechange.summary,
	sig.data,
	file = generate.filename(arguments$project, 'mutation_summary', 'RData')
	);

# run unlink, in case it exists from a previous run
if (!is.null(arguments$report)) {
	unlink(paste0(arguments$report, '/', 'mutation_summary.png'));
	file.symlink(
		paste0(arguments$output, '/',
			generate.filename(arguments$project, 'mutation_summary','png')),
		paste0(arguments$report, '/', 'mutation_summary.png')
		);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('MutationSummary','SessionProfile','txt'));
