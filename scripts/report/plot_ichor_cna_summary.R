### plot_cna_summary.R ############################################################################
# Plot CNA landscape (from ichorCNA).
# INPUT:
#	- copy number calls and ploidy/purity estimates (output by collect_ichorCNA_output.R)

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
	cat(paste0('Input: ', arguments$cnas));
	cat(paste0('CNA type: ', arguments$scale));

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
parser$add_argument('-c', '--cnas', type = 'character',
	help = 'path to gene by sample matrix of cna calls (generated by collect_ichorCNA_output.R)');
parser$add_argument('-m', '--metrics', type = 'character',
	help = 'path to ploidy and cellularity estimates (generated by collect_ichorCNA_output.R)');
parser$add_argument('-s', '--scale', type = 'character', default = 'CN',
	help = 'are we plotting CN or log2(ratio)? one of CN or ratio');
parser$add_argument('-z', '--report', type = 'character', help = 'path to report directory',
	default = NULL);

arguments <- parser$parse_args();

# import libraries
library(BoutrosLab.plotting.general);
library(xtable);

### READ DATA ######################################################################################
# get data
if (is.null(arguments$cnas)) {
	stop('ERROR: No CNA calls provided, please provide path to CNA matrix from IchorCNA.');
	} else {
	input.data <- read.delim(arguments$cnas, stringsAsFactors = FALSE);
	}

if (is.null(arguments$metrics)) {
	stop('ERROR: No CNA metrics provided, please provide path to metrics file from IchorCNA.');
	} else {
	ploidy.estimates <- read.delim(arguments$metrics, row.names = 1);
	}

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# collect list of all (ordered) samples
all.samples <- rev(rownames(ploidy.estimates));

ploidy.estimates$Order <- rev(1:nrow(ploidy.estimates));

# makes sure input.data is sorted
input.data$Chromosome <- factor(input.data$chr, levels = paste0('chr',c(1:22,'X','Y')));
input.data <- input.data[order(input.data$Chromosome, input.data$start),];
input.data <- input.data[which(input.data$Chromosome != 'chrY'),];

# reshape to bin x sample matrix
if ('CN' == arguments$scale) {
	tmp <- input.data[,c('Sample','Chromosome','start','end','Corrected_Copy_Number')];
	colnames(tmp)[5] <- 'CN';
	tmp$CN <- tmp$CN - 2;
	} else if ('ratio' == arguments$scale) {
	tmp <- input.data[,c('Sample','Chromosome','start','end','logR')];
	colnames(tmp)[5] <- 'Ratio';
	}

cna.wide <- reshape(
	tmp,
	direction = 'wide',
	idvar = c('Chromosome','start','end'),
	timevar = 'Sample'
	);

colnames(cna.wide) <- gsub('CN\\.|Ratio\\.','',colnames(cna.wide));
rownames(cna.wide) <- paste0('bin',1:nrow(cna.wide));

cna.data <- cna.wide[,all.samples];

if (length(all.samples) == 1) {
	cna.data <- data.frame(cbind(cna.data, cna.data));
	colnames(cna.data) <- rep(all.samples,2);
	}

# remove bins with 0 data
na.counts <- apply(cna.data,1,function(i) { length(i[!is.na(i)]); });
cna.data <- cna.data[which(na.counts > 0),]

### PLOT DATA ######################################################################################
# grab some parameters
axis.cex <- if (ncol(cna.data) <= 30) { 1
	} else if (ncol(cna.data) <= 50) { 0.75
	} else if (ncol(cna.data) <= 80) { 0.5
	} else if (ncol(cna.data) <= 100) { 0
	} else { 0 };

chromosomes <- gsub('chr','',tolower(cna.wide[which(na.counts > 0),]$Chromosome));
chr.breaks <- get.line.breaks(chromosomes);
covariate.colours <- force.colour.scheme(chromosomes, scheme = 'chromosome');

covariate.legends <- legend.grob(
	legends =  list(
		legend = list(
			colours = force.colour.scheme(c(1:22,'x'), scheme = 'chromosome'),
			labels = c(1:22,'X')
			)
		),
	size = 2
	);

# create heatmap
create.heatmap(
	cna.data,
	cluster.dimensions = 'none',
	covariates.top = list(
		rect = list(
			col = 'transparent',
			fill = covariate.colours,
			lwd = 0
			)
		),
	covariates.top.grid.border = list(col = 'black', lwd = 1),
	covariates.top.grid.col = list(col = 'black', lwd = 1),
	covariates.top.col.lines = chr.breaks,
	inside.legend = list(fun = covariate.legends, x = 1.02, y = 1),
	yaxis.lab = if (length(all.samples) == 1) { all.samples } else {
		gsub('\\.','-', simplify.ids(colnames(cna.data))) },
	yat = if (length(all.samples) == 1) { 1.5 } else { TRUE },
	xaxis.lab = rep('',nrow(cna.data)),
	yaxis.cex = axis.cex,
	xaxis.tck = 0,
	yaxis.tck = if (axis.cex == 0) { 0 } else { 0.2 },
	yaxis.fontface = 'plain',
	axes.lwd = 1,
	col.colour = 'black',
	grid.col = TRUE,
	force.grid.col = TRUE,
	col.lines = chr.breaks,
	print.colour.key = TRUE,
	fill.colour = 'grey50',
	at = if (arguments$scale == 'CN') { seq(-2.5,2.5,1) } else { seq(-2,2,0.1) },
#	colour.alpha = if (arguments$scale == 'CN') { 1 } else { 0.5 },
	colour.scheme = c('blue','white','red'),
	colourkey.labels.at = if (arguments$scale == 'CN') { seq(-2,2,1) } else { seq(-2,2,1) },
	colourkey.labels = if (arguments$scale == 'CN') { seq(-2,2,1) } else { seq(-2,2,1) },
	colourkey.cex = 1,
	height = if (length(all.samples) > 12) { 7 } else { 5 },
	width = 12,
	resolution = 1200,
	right.padding = 12,
	filename = generate.filename(arguments$project, 'ichorCNA_landscape','png')
	);

# create plot for cellulariy (purity) estimates
purity.plot <- create.barplot(
	Order ~ Tumour.Fraction,
	ploidy.estimates,
	yaxis.lab = rev(simplify.ids(rownames(ploidy.estimates))),
	ylimits = c(0.5, length(all.samples)+0.5),
	yat = seq(1,length(all.samples)),
	xaxis.tck = c(0.5,0),
	yaxis.tck = if (axis.cex == 0) { 0 } else { c(0.2,0) },
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xlab.label = "Tumour\nFraction",
	xlab.cex = 1,
	ylab.label = NULL,
	xlimits = c(0,1),
	xat = seq(0,1,0.5),
	xaxis.cex = 1,
	yaxis.cex = axis.cex,
	axes.lwd = 1,
	top.padding = 0,
	plot.horizontal = TRUE
	);

# create plot for ploidy estimates
ploidy.plot <- create.barplot(
	Order ~ Ploidy,
	ploidy.estimates,
	yaxis.lab = rep('', length(all.samples)),
	ylimits = c(0.5, length(all.samples)+0.5),
	yat = seq(1, length(all.samples)),
	xaxis.tck = c(0.5,0),
	yaxis.tck = 0,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xlab.label = "Ploidy\nEstimate",
	xlab.cex = 1,
	ylab.label = NULL,
	xlimits = c(0,max(ceiling(ploidy.estimates$Ploidy))+0.1),
	xat = seq(0,max(ceiling(ploidy.estimates$Ploidy)), 1),
	xaxis.cex = 1,
	yaxis.cex = axis.cex,
	axes.lwd = 1,
	top.padding = 0,
	plot.horizontal = TRUE
	);

# create plot for subclone fraction
subclone.plot <- create.barplot(
	Order ~ Estimated.Subclone.Fraction,
	ploidy.estimates,
	yaxis.lab = rep('', length(all.samples)),
	ylimits = c(0.5, length(all.samples)+0.5),
	yat = seq(1,length(all.samples)),
	xaxis.tck = c(0.5,0),
	yaxis.tck = 0,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xlab.label = "Subclone\nFraction",
	xlab.cex = 1,
	ylab.label = NULL,
	xlimits = c(0,1),
	xat = seq(0,1,0.5),
	xaxis.cex = 1,
	yaxis.cex = if (max(nchar(all.samples))>10) { 0.75*axis.cex } else { axis.cex },
	axes.lwd = 1,
	top.padding = 0,
	plot.horizontal = TRUE
	);

# combine them!
create.multipanelplot(
	plot.objects = list(purity.plot, ploidy.plot, subclone.plot),
	height = if (nrow(ploidy.estimates) <= 30) { 6 } else { 8 },
	width = 7,
	resolution = 200,
	filename = generate.filename(arguments$project, 'ichorCNA_metrics','png'),
	layout.width = 3,
	layout.height = 1,
	plot.objects.widths = c(1.2,1,1), #if (axis.cex > 0) { c(1.5,1,1) } else {  c(1.2,1,1) },
	left.legend.padding = 0,
	right.legend.padding = 0,
	top.legend.padding = 0,
	bottom.legend.padding = 0,
	x.spacing = 1,
	top.padding = -1
	);

save(
	cna.wide,
	ploidy.estimates,
	file = generate.filename(arguments$project, 'ichorCNA_data', 'RData')
	);

### LATEX ##########################################################################################
# if making output for a report
if (!is.null(arguments$report)) {

	tex.file <- paste0(
		arguments$report,
		'/cna_summary_ichor.tex'
		);

	# create symlinks for plots
	unlink(paste0(arguments$report,'/ichorCNA_landscape.png'));
	file.symlink(
		paste0(arguments$output,'/',
			generate.filename(arguments$project, 'ichorCNA_landscape','png')),
		paste0(arguments$report,'/ichorCNA_landscape.png')
		);

	unlink(paste0(arguments$report,'/ichorCNA_metrics.png'));
	file.symlink(
		paste0(arguments$output,'/',
			generate.filename(arguments$project, 'ichorCNA_metrics','png')),
		paste0(arguments$report,'/ichorCNA_metrics.png')
		);

	# write some captions
	summary.caption <- "Estimated tumour fraction (left), ploidy (middle) and subclonal fraction (right) estimates for each sample; only available for samples with a matched normal."; 
	if (arguments$scale == 'CN') {
		cna.caption <- "Ploidy-adjusted CNA profile, ordered by genomic position (chr1 [left] to chrX [right]) and sample (alphabetical [top to bottom]); only available for tumour samples with a matched normal.";
		} else {
		cna.caption <- "Heatmap shows log$_2$(ratio) per bin from ichorCNA; ordered by genomic position (chr1 [left] to chrX [right]) and sample (alphabetical [top to bottom]); only available for tumour samples with a matched normal.";
		}

	# write for latex
	write("\\section{SCNA Summary}", file = tex.file);
	write("Copy number calls from ichorCNA:", file = tex.file, append = TRUE);

	# CNA landscape plot
	write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
	write("\\begin{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\includegraphics[width=0.9\\textwidth]{",
		paste0(arguments$report, '/', 'ichorCNA_landscape.png'), '}'
		), file = tex.file, append = TRUE);
	write("\\end{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\caption{", cna.caption, "}"
		), file = tex.file, append = TRUE);
	write("\\end{figure}\n", file = tex.file, append = TRUE);

	# Purity/Ploidy plot
	if ((nrow(ploidy.estimates) <= 20)) {
		write("\n\\begin{figure}[h!]", file = tex.file, append = TRUE);
		write("\\begin{center}", file = tex.file, append = TRUE);
		write(paste0(
			"\\includegraphics[width=0.5\\textwidth]{",
			paste0(arguments$report, '/', 'ichorCNA_metrics.png'), '}'
			), file = tex.file, append = TRUE);
		} else {
		write("\\pagebreak\n\\begin{figure}[h!]", file = tex.file, append = TRUE);
		write("\\begin{center}", file = tex.file, append = TRUE);
		write(paste0(
			"\\includegraphics[width=0.9\\textwidth]{",
			paste0(arguments$report, '/', 'ichorCNA_metrics.png'), '}'
			), file = tex.file, append = TRUE);
		}

	write("\\end{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\caption{", summary.caption, "}"
		), file = tex.file, append = TRUE);
	write("\\end{figure}\n", file = tex.file, append = TRUE);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('IchorCNASummary','SessionProfile','txt'));
