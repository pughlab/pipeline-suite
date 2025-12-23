### plot_qc_metrics.R #############################################################################
# Plot various DNA-Seq or EM-Seq QC metrics for the MultiMMR panel

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

# function to format data for output table
format.mean <- function(x) {
	paste0(round(median(x),2), ' (', round(quantile(x,0.25),2), ' - ', round(quantile(x,0.75),2), ')')
	}

### PREPARE SESSION ################################################################################
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-o', '--output_dir', type = 'character', help = 'path to output directory');
parser$add_argument('-r', '--report_dir', type = 'character', help = 'path to report directory');
parser$add_argument('-s', '--sample_yaml', type = 'character', help = 'path to sample (BAM) yaml');
parser$add_argument('-g', '--correlations', type = 'character', 
	help = 'path to germline correlation file');
parser$add_argument('-c', '--coverage', type = 'character', 
	help = 'path to coverage metrics (HSMetrics or WGSMetrics)');
parser$add_argument('-t', '--seq_type', type = 'character', help = 'either wgs, dnaseq or emseq', 
	default = 'dnaseq');

arguments <- parser$parse_args();

# load libraries
library(BoutrosLab.plotting.general);
library(yaml);
library(xtable);

# what's the date?
date <- Sys.Date();

# get input files
cov.file <- arguments$coverage;
cor.file <- arguments$correlations;
output.dir <- arguments$output_dir;
reports.dir <- arguments$report_dir;

# create (if necessary) and move to output directory
if (!dir.exists(output.dir)) {
	dir.create(output.dir);
	}

if (!dir.exists(reports.dir)) {
	dir.create(reports.dir);
	}

### MAIN ##########################################################################################
## read in data files
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

sample.info$Patient <- as.character(sample.info$Patient);
sample.info$Sample <- as.character(sample.info$Sample);

sample.info <- sample.info[order(sample.info$Patient, sample.info$Sample),];
smp.names <- sample.info$Sample;


# read in genotype concordance file
# account for differences in header depending on version of bcftools used
concordance.data <- if ('dnaseq' == arguments$seq_type) { 
	read.delim(cor.file, skip = 19);
	} else { NULL; }

if (!is.null(concordance.data) && (ncol(concordance.data) == 1)) {
	concordance.data <- read.delim(cor.file, skip = 32);
	}

if (!is.null(concordance.data)) {
	colnames(concordance.data) <- c('Metric','Query.Sample','Genotyped.Sample','Discordance',
		'pvalue','Total.Sites','N.matched')[1:ncol(concordance.data)];

	# convert N mismatches (discordance) to a proportion of total sites
	concordance.data$Proportion <- 1 - (concordance.data$Discordance / concordance.data$Total.Sites);

	# reshape concordance data
	cor.data <- reshape(
		concordance.data[,c('Query.Sample','Genotyped.Sample','Proportion')],
		direction = 'wide',
		idvar = 'Genotyped.Sample',
		timevar = 'Query.Sample'
		);

	rownames(cor.data) <- cor.data[,1];
	colnames(cor.data) <- gsub('Proportion\\.','',colnames(cor.data));
	cor.data <- cor.data[,-1];

	# current data is one-sided, so mirror this
	for (i in 1:nrow(concordance.data)) {
		query.smp <- concordance.data[i,]$Query.Sample;
		genotyped.smp <- concordance.data[i,]$Genotyped.Sample;
		cor.data[query.smp,genotyped.smp] <- cor.data[genotyped.smp,query.smp];
		}

	for (smp in smp.names) { cor.data[smp,smp] <- 1; }

	} else {
	cor.data <- data.frame();
	}

# add in missing samples
missing.samples <- setdiff(smp.names, rownames(cor.data));
cor.data[missing.samples,] <- NA;

missing.samples <- setdiff(smp.names, colnames(cor.data));
cor.data[,missing.samples] <- NA;

cor.data <- cor.data[smp.names,smp.names];


# read in coverage metrics
metric.data <- read.delim(cov.file, row.names = 1);

if ('wgs' == arguments$seq_type) {

	target_coverage <- '1';
	metric.data$percent_bases_above_targetX <- metric.data[,paste0('PCT_', target_coverage, 'X')];

	colnames(metric.data)[which(colnames(metric.data) == 'MEAN_COVERAGE')] <- 'mean';
	metric.data$lower_cov <- metric.data$mean - metric.data$SD_COVERAGE; 
	metric.data$upper_cov <- metric.data$mean + metric.data$SD_COVERAGE; 
	metric.data <- metric.data[,c('mean','lower_cov','upper_cov','percent_bases_above_targetX')];

	} else {

	target_coverage <- '100';
	metric.data$percent_on_bait <- metric.data$ON_BAIT_BASES / metric.data$PF_BASES_ALIGNED;
	metric.data$percent_off_bait <- metric.data$OFF_BAIT_BASES / metric.data$PF_BASES_ALIGNED;
	metric.data$percent_near_bait <- metric.data$NEAR_BAIT_BASES / metric.data$PF_BASES_ALIGNED;

	metric.data <- metric.data[,c('MEAN_TARGET_COVERAGE','PCT_EXC_DUPE','PCT_TARGET_BASES_100X','percent_on_bait','percent_off_bait','percent_near_bait')];
	colnames(metric.data)[1:3] <- c('mean','duplication_rate','percent_bases_above_targetX');
	metric.data[,c('lower_cov','upper_cov')] <- NA;

	}

# format data for plotting
metric.data$Sample <- rownames(metric.data);
metric.data$Order <- 1:nrow(metric.data);


# find line breaks (separate patients with multiple samples)
line.breaks <- get.line.breaks(sample.info$Patient);

# determine some parameters
add.rectangle <- TRUE;
if (length(smp.names) > 50) { add.rectangle <- FALSE; }
if (length(smp.names) < 6) { add.rectangle <- FALSE; }

axis.cex <- if (length(smp.names) <= 30) { 1
	} else if (length(smp.names) <= 80) { 0.8
	} else if (length(smp.names) <= 100) { 0.6
	} else if (length(smp.names) <= 140) { 0.5
	} else if (length(smp.names) <= 180) { 0.4
	} else { 0 };

# find plot limits	
max.cov <- max(c(metric.data$upper_cov,metric.data$mean), na.rm = TRUE);
if ( (max.cov == 1) || (max.cov == 500) ) { max.cov <- max(metric.data$mean, na.rm = TRUE); }

if ('wgs' == arguments$seq_type) {
	total.cov.limits <- c(0,10);
	total.cov.at <- seq(0,10,2);
	} else if (max.cov <= 100) {
	total.cov.limits <- c(0, 100);
	total.cov.at <- seq(0, 100, 25);
	} else if (max.cov <= 200) {
	total.cov.limits <- c(0, 200);
	total.cov.at <- seq(0, 200, 50);
	} else if (max.cov <= 400) {
	total.cov.limits <- c(0, 400);
	total.cov.at <- seq(0, 400, 100);
	} else if (max.cov > 400 & max.cov < 850) {
	total.cov.limits <- c(0, ceiling(max.cov/100)*100);
	total.cov.at <- seq(0, total.cov.limits[2], length.out = 5);
	} else if (max.cov <= 1000) {
	total.cov.limits <- c(0, 1000);
	total.cov.at <- seq(0, 1000, 250);
	} else if (max.cov <= 2000) {
	total.cov.limits <- c(0, 2000);
	total.cov.at <- seq(0, 2000, 500);
	} else {
	total.cov.limits <- c(0, ceiling(max.cov/1000)*1000);
	total.cov.at <- seq(0, total.cov.limits[2], length.out = 5);
	}

# move to output directory
setwd(output.dir);

output.stem <- if ('wgs' == arguments$seq_type) {
	paste0(arguments$project, '_sWGS');
	} else if ('emseq' == arguments$seq_type) {
	paste0(arguments$project, '_EMSeq');
	} else {
	paste0(arguments$project, '_DNASeq');
	}

# plot coverage metrics
metric.data$percent_target <- metric.data$percent_bases_above_targetX * total.cov.limits[2];
key.lab <- if ('wgs' == arguments$seq_type) { 'mean read depth ( \u00B1 SD )';
	} else { 'mean target coverage' };

cov.key <- list(
	fun = draw.key,
	args = list(
		key = list(
			lines = list(cex = 1, col = 'black', pch = 19, 
				type = if ('wgs' == arguments$seq_type) { 'b' } else { 'p' } ),
			text = list(lab = key.lab, cex = 1),
			divide = 1
			)
		),
	x = 0.5, y = 1
	);

total.coverage.plot <- create.scatterplot(
	mean ~ Order,
	metric.data,
	xaxis.lab = rep('',nrow(metric.data)),
	xlimits = c(0.5, length(smp.names)+0.5),
	xat = seq(1,length(smp.names)),
	yaxis.tck = c(0.5,0.5),
	xaxis.tck = 0,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	ylab.label = 'Read depth',
	ylab.cex = 1.5,
	xlab.label = '',
	ylimits = total.cov.limits,
	yat = total.cov.at,
	yaxis.cex = 1.3,
	cex = max(c(axis.cex,0.2)),
	axes.lwd = 1,
	y.error.up = metric.data$upper_cov - metric.data$mean,
	y.error.down = metric.data$mean - metric.data$lower_cov,
	error.whisker.angle = 0,
	abline.v = line.breaks,
	abline.col = 'grey80',
	add.rectangle = TRUE,
	xleft.rectangle = seq(0.5, length(smp.names)-0.5, 1),
	xright.rectangle = seq(1.5, length(smp.names)+0.5,1),
	ybottom.rectangle = -1,
	ytop.rectangle = metric.data$percent_target,
	col.rectangle = 'red',
	alpha.rectangle = 0.5,
	top.padding = 1,
	key = list(
		space = 'right',
		text = list(
			lab = rev(c('0 %','25 %','50 %','75 %','100 %')),
			cex = 1.3,
			col = 'black'
			),
		padding.text = 8.2
		),
	legend = list(top = cov.key)
	);

if ('dnaseq' != arguments$seq_type) {

	total.coverage.plot$x.scales$labels <- simplify.ids(smp.names);
	total.coverage.plot$x.scales$rot <- c(90,0);
	total.coverage.plot$x.scales$cex <- axis.cex;
	total.coverage.plot$ylab$fontface <- 'plain';
	total.coverage.plot$legend$right$args$key$padding.text <- 12;
	total.coverage.plot$legend$right$args$key$text$col <- 'red';
	total.coverage.plot$par.settings$layout.widths$ylab.axis.padding <- 3;
	total.coverage.plot$par.settings$layout.widths$key.right <- 0.8;
	total.coverage.plot$par.settings$layout.widths$axis.key.padding <- 0;	

	write.plot(
		total.coverage.plot,
		filename = generate.filename(output.stem, 'coverage_summary','png'),
		height = 4,
		width = 11,
		resolution = 200
		);

	} else {

	# plot the results (correlations of germline genotypes)
	cor.heatmap <- create.heatmap(
		cor.data,
		cluster.dimensions = 'none',
		ylab.label = 'GT concordance',
		ylab.cex = 1.5,
		xaxis.lab = NA,
		yaxis.lab = NA,
		xaxis.cex = axis.cex,
		yaxis.cex = axis.cex,
		xaxis.tck = if (nrow(cor.data) <= 50) { 0.2 } else { 0 },
		yaxis.tck = if (nrow(cor.data) <= 50) { 0.2 } else { 0 },
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		colourkey.cex = 1,
		colour.scheme = c('white','black'),
		at = seq(0.8,1,0.001),
		colourkey.labels.at = c(0.8,0.9,0.95,1),
		colourkey.labels = c('< 0.8','0.90', '0.95', '1'),
		grid.row = TRUE,
		force.grid.row = TRUE,
		row.lines = line.breaks,
		row.colour = 'grey80',
		col.colour = 'grey80',
		grid.col = TRUE,
		force.grid.col = TRUE,
		col.lines = line.breaks,
		axes.lwd = 1,
		fill.colour = 'grey80'
		);

	# combine them!
	create.multipanelplot(
		plot.objects = list(total.coverage.plot, cor.heatmap),
		filename = generate.filename(output.stem, 'coverage_summary','png'),
		height = 12,
		width = 11,
		resolution = 200,
		layout.width = 1,
		layout.height = 2,
		plot.objects.heights = c(1,4),
		left.legend.padding = 0,
		right.legend.padding = 0,
		top.legend.padding = 1,
		bottom.legend.padding = 0,
		xlab.axis.padding = 1,
		ylab.axis.padding = 1,
		y.spacing = 0,
		legend = list(
			inside = list(
				fun = draw.key,
				args = list(
					key = list(
						rect = list(col = 'red', alpha = 0.5, size = 1),
						text = list(
							lab = paste0('Percent bases above ', 
								target_coverage, 'x'),
							cex = 1
							)
						)
					),
				x = 0.78,
				y = 0.96
				)
			)
		);
	}

# if this is target sequencing, plot the on/off/near target bases
if ('wgs' != arguments$seq_type) {

	stacked.data <- reshape(
		metric.data[,c('Sample','percent_on_bait','percent_off_bait','percent_near_bait')],
		direction = 'long',
		varying = list(2:4),
		v.names = 'Percent',
		times = c('percent_on_bait','percent_off_bait','percent_near_bait'),
		timevar = 'Group',
		idvar = 'Sample'
		);

	stacked.data$Group <- factor(
		stacked.data$Group,
		levels = rev(c('percent_on_bait','percent_near_bait','percent_off_bait'))
		);

	stack.legend <- list(fun = draw.key, args = list(key = list(
		space = 'top',
		rect = list(col = rev(rev(default.colours(12))[1:3]), size = 1),
		text = list(lab = c('on target', 'near target', 'off target'), cex = 1),
		columns = 3,
		between = 0.7
		)), x = 0.02, y = 1);

	# plot bait efficiency
	efficiency.plot <- create.barplot(
		Percent ~ Sample,
		stacked.data,
		groups = stacked.data$Group,
		stack = TRUE,
		col = rev(default.colours(12))[1:3],
		xlab.label = NULL,
		ylab.label = '% Aligned Bases',
		ylab.cex = 1.5,
		xaxis.lab = simplify.ids(stacked.data$Sample),
		xaxis.rot = 90,
		xaxis.tck = c(0.5,0),
		yaxis.tck = c(0.5,0),
		xaxis.cex = axis.cex,
		yaxis.cex = 1,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		ylimits = c(0,1.1),
		yat = seq(0,1,0.2),
		ylab.axis.padding = 2,
		right.padding = 0,
		legend = list(inside = stack.legend),
		top.padding = 1,
		style = 'Nature'
		);

	write.plot(
		efficiency.plot,
		filename = generate.filename(output.stem, 'bait_efficiency_metrics','png'),
		height = 4,
		width = 11,
		resolution = 200
		);

	# plot duplication rates
	duplication.plot <- create.barplot(
		duplication_rate ~ Sample,
		metric.data,
		xlab.label = NULL,
		ylab.label = 'Duplication Rate',
		ylab.cex = 1.5,
		ylab.axis.padding = 2,
		right.padding = 0,
		xaxis.lab = simplify.ids(metric.data$Sample),
		xaxis.rot = 90,
		xaxis.tck = c(0.5,0),
		yaxis.tck = c(0.5,0),
		xaxis.cex = axis.cex,
		yaxis.cex = 1,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		ylimits = c(0,1),
		yat = seq(0,1,0.2),
		top.padding = 1,
		style = 'Nature'
		);

	write.plot(
		duplication.plot,
		filename = generate.filename(output.stem, 'duplication_summary','png'),
		height = 4,
		width = 11,
		resolution = 200
		);
	}

save(
	sample.info,
	metric.data,
	cor.data,	
	file = generate.filename(output.stem, 'qc_plot_data','RData')
	);

### PER-SAMPLE SUMMARIES ###########################################################################
# find all sWGS files
all.files <- if ('wgs' == arguments$seq_type) {
	list.files(pattern = 'sWGS_qc_plot_data.RData');
	} else if ('emseq' == arguments$seq_type) {
	list.files(pattern = 'EMSeq_qc_plot_data.RData');
	} else {
	list.files(pattern = 'DNASeq_qc_plot_data.RData');
	}

orig.file <- generate.filename(output.stem, 'qc_plot_data','RData');
all.files <- setdiff(all.files, orig.file);

orig.sample.info <- sample.info;
orig.metric.data <- metric.data;

combined.metrics <- metric.data;

for (file in all.files) {
	load(file);
	combined.metrics <- unique(rbind(combined.metrics, metric.data));
	}

setwd(reports.dir);

# initiate report object
tex.file <- if ('wgs' == arguments$seq_type) {
	'qc_summary_sWGS.tex';
	} else if ('emseq' == arguments$seq_type) {
	'qc_summary_EMSeq.tex';
	} else {
	'qc_summary_DNASeq.tex';
	}

for (patient in unique(orig.sample.info$Patient)) {

	if (!dir.exists(patient)) {
		dir.create(patient);
		}

	setwd(patient);

	these.smps <- orig.sample.info[which(orig.sample.info$Patient == patient),]$Sample;
	smp.data <- combined.metrics[which(combined.metrics$Sample %in% these.smps),];

	if ('wgs' == arguments$seq_type) {
		to.write <- smp.data[,c('Sample','mean')];
		colnames(to.write)[2] <- 'Mean Coverage';

		to.write[,'Historical Average'] <- format.mean(combined.metrics$mean);
		caption <- "Mean coverage per sample compared to historical average (all previously processed samples; reported as median (IQR) ).";

		} else {

		to.write <- smp.data[,c('Sample','duplication_rate','mean','percent_bases_above_targetX')];
		colnames(to.write)[2:4] <- c('Duplication Rate','Mean Coverage','Fraction Bases Above 100X');
		to.write[,2:4] <- round(to.write[,2:4],2);

		to.write[nrow(to.write)+1,] <- c('Historical Average',0,0,0,0);
		to.write[nrow(to.write),2] <- format.mean(combined.metrics$duplication_rate);
		to.write[nrow(to.write),3] <- format.mean(combined.metrics$mean);
		to.write[nrow(to.write),4] <- format.mean(combined.metrics$percent_bases_above_targetX);

		caption <- "HS coverage metrics per sample compared to historical average (all previously processed samples; reported as median + IQR).";
		}

	if ('wgs' == arguments$seq_type) {
		write("\\subsection{sWGS}", file = tex.file);
		} else if ('emseq' == arguments$seq_type) {
		write("\\subsection{EM-Seq}", file = tex.file);
		} else {
		write("\\subsection{DNA-Seq}", file = tex.file);
		}

	# add data to tex file
	print(
		xtable(
			to.write,
			align = c("r","l",rep("c", ncol(to.write)-1)),
			caption = caption
			),
		file = tex.file,
		append = TRUE,
		include.rownames = FALSE,
		latex.environments = ""
		);

	setwd(reports.dir);
	}

### SAVE SESSION INFO ##############################################################################
setwd(output.dir);
save.session.profile(generate.filename(output.stem, 'PlotQC_SessionProfile','txt'));
