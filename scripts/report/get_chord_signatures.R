### get_chord_signatures.R #########################################################################
# Extract HRD signatures using CHORD (Classifier of HOmologous Recombination Deficiency);
# uses somatic SNV/INDEL/SV calls

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

	# write key variables to file
	cat("\n### VARIABLES #################################################################\n");
	print(paste0('Input MAF: ', arguments$maf));
	print(paste0('Input SVs: ', arguments$sv));
	print(paste0('VAF threshold applied: ', vaf.threshold));

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

### PREPARE SESSION ################################################################################
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'hg38 or hg19', default = 'hg38');
parser$add_argument('-m', '--maf', type = 'character', help = 'mutation calls in MAF format');
parser$add_argument('-s', '--sv', type = 'character', help = 'combined SV calls from Mavis');
parser$add_argument('-t', '--vaf_threshold', type = 'character',
	help = 'threshold to filter variants', default = 0.1);
parser$add_argument('-l','--lib_paths', type = 'character',
	help = 'path to library if not installed in default',
	default = '/cluster/projects/pughlab/src/r_lib/library/4.1/');
parser$add_argument('-z', '--report', type = 'character', help = 'path to report directory',
	default = NULL);

arguments <- parser$parse_args();

# load required libraries
if (!is.null(arguments$lib_paths)) {
	.libPaths(c(arguments$lib_paths,.libPaths()));
	}

# import libraries
library(BoutrosLab.plotting.general);
library(xtable);

library(randomForest);
library(CHORD);
library(signature.tools.lib);

# get genome object
if (arguments$ref_type == 'hg38') {
	library("BSgenome.Hsapiens.UCSC.hg38");
	BSgenome.Hsapiens.UCSC <- 'BSgenome.Hsapiens.UCSC.hg38';
	} else if (arguments$ref_type == 'hg19') {
	library("BSgenome.Hsapiens.UCSC.hg19");
	BSgenome.Hsapiens.UCSC <- 'BSgenome.Hsapiens.UCSC.hg19';
	} else {
	stop('Unrecognized genome build requested.');
	}

# indicate filtering thresholds
vaf.threshold <- arguments$vaf_threshold;

### READ DATA ######################################################################################
# get mutation data
if (is.null(arguments$maf)) {
	stop('ERROR: No input MAF provided, please provide path to SNV/INDEL calls in MAF format.');
	} else {
	input.maf <- read.delim(arguments$maf, stringsAsFactors = FALSE, comment.char = '#');
	}

if (is.null(arguments$sv)) {
	stop('ERROR: No input SVs provided, please provide path to SV calls from MAVIS.');
	} else {
	sv.data <- read.delim(arguments$sv, stringsAsFactors = FALSE);
	}

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# get list of samples
all.samples <- sort(unique(input.maf$Tumor_Sample_Barcode));

# filter input MAF
input.maf$t_vaf <- input.maf$t_alt_count / input.maf$t_depth;
input.maf$n_vaf <- 1 - (input.maf$n_ref_count / input.maf$n_depth);
input.maf <- input.maf[which(input.maf$t_vaf >= vaf.threshold),];
input.maf <- input.maf[which(input.maf$t_depth >= 20 & input.maf$n_vaf < 0.05),];

# split SNVs
keep.fields <- c('Tumor_Sample_Barcode','Chromosome','Start_Position','Reference_Allele','Allele');
snv.data <- input.maf[which(input.maf$Variant_Type == 'SNP'),keep.fields];

# split and format INDELs
keep.fields <- c('Tumor_Sample_Barcode','Chromosome','Start_Position','flanking_bps');

del.data <- input.maf[which(input.maf$Variant_Type == 'DEL'),keep.fields];
del.data$Start_Position <- del.data$Start_Position - 1;
del.data$Reference_Allele <- sapply(del.data$flanking_bps, function(i) {
	substr(i, 2, nchar(i)-1) } );
del.data$Allele <- sapply(del.data$flanking_bps, function(i) {
	substr(i, 2, 2) } );

ins.data <- input.maf[which(input.maf$Variant_Type == 'INS'),c(keep.fields,'Allele')];
ins.data$Reference_Allele <- sapply(ins.data$flanking_bps, function(i) {
	substr(i, 2, 2) } );
ins.data$Allele <- paste0(ins.data$Reference_Allele, ins.data$Allele);

keep.fields <- c('Tumor_Sample_Barcode','Chromosome','Start_Position','Reference_Allele','Allele');
indel.data <- rbind(del.data[,keep.fields], ins.data[,keep.fields]);

# indicate key fields
colnames(snv.data) <- c('Sample','chr','pos','ref','alt');
colnames(indel.data) <- c('Sample','chr','pos','ref','alt');

# format SV data
sv.data <- sv.data[which(sv.data$library %in% all.samples),];

# subset only somatic variants
sv.states <- apply(
	sv.data[,grep('diseased_genome|normal_genome', colnames(sv.data))],
	1, function(i) {
		if (any(grepl('germline', i))) { return('germline') } else { return('somatic') } }
	);
sv.data <- sv.data[which(sv.states == 'somatic'),];

# add evidence filter
sv.data$Evidence <- apply(
	sv.data[,c('break1_split_reads','break2_split_reads','spanning_reads','linking_split_reads')],
	1, function(i) {
	if (any(i == 'None')) { i[which(i == 'None')] <- 0; }
	i <- sapply(i, function(j) { max(as.numeric(unlist(strsplit(as.character(j),';')))); });
	# if call method uses spanning reads
	if (i[3] >= 20) { return('PASS')
		# if call method uses split reads
		} else if (i[1] >= 10 & i[2] >= 10) { return('PASS')
		} else if (i[4] >= 10) { return('PASS')
		} else { return('FAIL');
		}
	});

sv.data <- sv.data[which(sv.data$Evidence == 'PASS' | grepl(';', sv.data$tools)),];

# subset/code variant types
sv.data$sv_type <- c('DEL','DUP','INV','TRA','TRA')[match(sv.data$event_type,
	c('deletion','duplication','inversion','translocation','inverted translocation'))];
sv.data <- sv.data[!is.na(sv.data$sv_type),];

# calculate variant lengths
sv.data$Length <- abs(sv.data$break2_position_start - sv.data$break1_position_start);
sv.data[which(sv.data$sv_type == 'TRA'),]$Length <- NA;
sv.data <- sv.data[which(sv.data$Length > 100 | is.na(sv.data$Length)),];

# get signature weights for each sample
mutation.context <- list();
sig.predictions <- list();

for (smp in all.samples) {

	mutation.context[[smp]] <- extractSigsChord(
		df.snv = snv.data[which(snv.data$Sample == smp),-1],
		df.indel = indel.data[which(indel.data$Sample == smp),-1],
		df.sv = sv.data[which(sv.data$Sample == smp),c('sv_type','Length')],
		ref.genome = BSgenome.Hsapiens.UCSC,
		sample.name = smp
		);

	sig.predictions[[smp]] <- chordPredict(mutation.context[[smp]],
		do.bootstrap = TRUE
		);
	}

contexts <- do.call(rbind, mutation.context);
predictions <- do.call(rbind, sig.predictions);

# save results
save(
	contexts,
	predictions,
	file = generate.filename(arguments$project, 'CHORD_signatures', 'RData')
	);

write.table(
	predictions,
	file = generate.filename(arguments$project, 'chord_predictions', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# make a symlink for use with final summary plots
unlink('CHORD_predictions.tsv');
file.symlink(
	generate.filename(arguments$project, 'chord_predictions', 'tsv'),
	'CHORD_predictions.tsv'
	);

### PLOTTING #######################################################################################
# format data
brca1 <- predictions[,c('sample','p_BRCA1.50%','p_BRCA1.5%','p_BRCA1.95%')];
brca2 <- predictions[,c('sample','p_BRCA2.50%','p_BRCA2.5%','p_BRCA2.95%')];
hrd <- predictions[,c('sample','p_hrd.50%','p_hrd.5%','p_hrd.95%')];

brca1$group <- 'BRCA1';
brca2$group <- 'BRCA2';
hrd$group <- 'HRD';

new.names <- c('sample','median','error.down','error.up','group');
colnames(brca1) <- new.names;
colnames(brca2) <- new.names;
colnames(hrd) <- new.names;

plot.data <- rbind(hrd, brca1, brca2);

plot.data$sample <- factor(plot.data$sample, levels = all.samples);
plot.data$x.pos <- as.numeric(plot.data$sample);
plot.data[which(plot.data$group == 'HRD'),]$x.pos <- plot.data[which(plot.data$group == 'HRD'),]$x.pos-0.2;
plot.data[which(plot.data$group == 'BRCA2'),]$x.pos <- plot.data[which(plot.data$group == 'BRCA2'),]$x.pos+0.2;

# grab some parameters
axis.cex <- if (length(all.samples) <= 30) { 1
	} else if (length(all.samples) <= 50) { 0.75
	} else if (length(all.samples) <= 80) { 0.5
	} else { 0 };

# set colour scheme
colour.scheme <- c('black',default.colours(3, 'spiral.dawn')[-2]);
names(colour.scheme) <- c('HRD','BRCA1','BRCA2');

# make a legend
type.key <- list(
	title = 'Score                      ',
	cex.title = 1.2,
	points = list(col = colour.scheme, pch = 19, cex = 1),
	text = list(lab = c('HRD','BRCA1-type HRD','BRCA2-type HRD'), cex = 1)
	);

ci.key <- list(
	points = list(col = 'black', pch = '|', cex = 1.5),
	text = list(lab = '90% CI', cex = 1)
	);

signif.key <- list(
	rect = list(col = 'grey80', size = 3, alpha = 0.5),
	text = list(lab = 'HR-deficient', cex = 1)
	);

# make the plot
create.scatterplot(
	median ~ x.pos,
	plot.data,
	col = colour.scheme[match(plot.data$group, names(colour.scheme))],
	y.error.up = plot.data$error.up - plot.data$median,
	y.error.down = plot.data$median - plot.data$error.down,
	y.error.bar.col = colour.scheme[match(plot.data$group, names(colour.scheme))],
	group.specific.colouring = FALSE,
	error.bar.length = 0,
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	ylab.label = expression('Signature Probability'),
	ylab.cex = 1.4,
	ylab.axis.padding = 3,
	xlimits = c(0.5,length(all.samples)+0.5),
	xat = seq(1,length(all.samples),1),
	xaxis.lab = simplify.ids(all.samples),
	xlab.label = NULL,
	xaxis.cex = axis.cex,
	yaxis.cex = 1,
	xaxis.rot = 90,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	add.rectangle = any(predictions$hr_status == 'HR_deficient'),
	xleft.rectangle = which(predictions$hr_status != 'HR_proficient')-0.5,
	xright.rectangle = which(predictions$hr_status != 'HR_proficient')+0.5,
	ytop.rectangle = 1.1,
	ybottom.rectangle = -0.1,
	col.rectangle = 'grey80',
	alpha.rectangle = 0.5,
	right.padding = 20,
	legend = list(
		inside = list(fun = draw.key, args = list(signif.key), x = 1, y = 0.95),
		inside = list(fun = draw.key, args = list(type.key), x = 1.02, y = 0.75),
		inside = list(fun = draw.key, args = list(ci.key), x = 1.01, y = 0.3)
		),
	height = 4,
	width = if (length(all.samples) < 10) { 6 } else { 8 },
	resolution = 600,
	filename = generate.filename(arguments$project, 'chord_predictions', 'png')
	);

predictions$hr_status <- factor(predictions$hr_status,
	levels = c('HR_deficient','HR_proficient','cannot_be_determined')
	);
predictions <- predictions[order(predictions$hr_status),];

successful.predictions <- which(predictions$hr_status != 'cannot_be_determined');

if (length(successful.predictions) >= 5) {
	to.write <- predictions[successful.predictions,];
	to.write <- to.write[1:min(c(nrow(to.write),12)),];
	} else {
	to.write <- predictions[1:min(c(nrow(predictions),12)),];
	}

to.write <- to.write[,c('sample','hr_status','p_BRCA1','p_BRCA2','p_hrd')];

# if making output for a report
if (!is.null(arguments$report)) {

	tex.file <- paste0(
		arguments$report,
		'/chord_signature_summary.tex'
		);

	# run unlink, in case it exists from a previous run
	unlink(paste0(arguments$report, '/', 'chord_predictions.png'));
	file.symlink(
		paste0(arguments$output, '/',
			generate.filename(arguments$project, 'chord_predictions', 'png')),
		paste0(arguments$report, '/', 'chord_predictions.png')
		);

	### LATEX ##################################################################################
	# write for latex
	write("\\section{Mutation Signatures}", file = tex.file);

	write(
		paste0('HRD signatures were applied to the current dataset using the CHORD (Classifier of HOmologous Recombination Deficiency) algorithm, using ENSEMBLE mutations with minimum a VAF threshold = ', vaf.threshold, '.'),
		file = tex.file,
		append = TRUE
		);

	# add mutation_signatures plot
	write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
	write("\\begin{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\includegraphics[height=0.3\\textheight]{",
		paste0(arguments$report, '/', 'chord_predictions.png'), '}'
		), file = tex.file, append = TRUE);
	write("\\end{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\caption{CHORD signature predictions; signature probabilities with 5-95\\% confidence intervals are shown for each signature for each sample.}"
		), file = tex.file, append = TRUE);
	write("\\end{figure}\n", file = tex.file, append = TRUE);

	# add table summary
	caption <- 'HRD status and signature probabilities, as generated by the CHORD package.';
	print(
		xtable(
			to.write,
			caption = caption 
			),
		file = tex.file,
		include.rownames = FALSE,
		append = TRUE
		);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CHORD','SessionProfile','txt'));
