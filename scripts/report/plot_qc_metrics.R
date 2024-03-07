### plot_qc_metrics.R #############################################################################
# Plot various DNA-Seq QC metrics: summarize covereage, contamination estimates and cross-sample
#	similarity (germline correlations)
# OR
# Plot various RNA-Seq QC metrics: summarize genome coverage (regions) + depth metrics and cross-
# 	sample similary (correlation of genes.rpkm.gct)

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

parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-o', '--output_dir', type = 'character', help = 'path to output directory');
parser$add_argument('-s', '--sample_yaml', type = 'character', help = 'path to sample (BAM) yaml');
parser$add_argument('-t', '--seq_type', type = 'character', help = 'either wgs, exome, targeted or rna', 
	default = 'exome');
parser$add_argument('-g', '--correlations', type = 'character', help = 'path to germline correlation file');
parser$add_argument('-c', '--contest', type = 'character', help = 'path to contamination file (ContEst)', default = NULL);
parser$add_argument('-m', '--contamination', type = 'character', help = 'path to contamination file (CalculateContamination)',
	default = NULL);
parser$add_argument('-v', '--coverage', type = 'character', help = 'path to coverage metrics (depthOfCoverage)');
parser$add_argument('-w', '--wgscoverage', type = 'character', help = 'path to WGS coverage metrics (WGSMetrics)');
parser$add_argument('-e', '--hscoverage', type = 'character', help = 'path to HS coverage metrics (HSMetrics)');
parser$add_argument('-z', '--report', type = 'character', help = 'path to report directory',
	default = NULL);

arguments <- parser$parse_args();

# load libraries
library(BoutrosLab.plotting.general);
library(yaml);
library(xtable);

# what's the date?
date <- Sys.Date();

# clarify input files
is.dna <- !(arguments$seq_type == 'rna');
summary.file <- ifelse(!is.null(arguments$wgscoverage), arguments$wgscoverage, arguments$coverage);
contest.file <- arguments$contest;
other.contamination <- arguments$contamination;
cor.file <- arguments$correlations;
output.dir <- arguments$output_dir;

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


# read in correlation heatmap data
cor.data <- read.delim(cor.file, row.names = 1);


# read in coverage metrics
metric.data  <- read.delim(summary.file, row.names = 1);

if (!is.null(arguments$wgscoverage)) {

	target_coverage <- '15';
	metric.data$percent_bases_above_targetX <- metric.data[,paste0('PCT_', target_coverage, 'X')];

	colnames(metric.data)[which(colnames(metric.data) == 'MEAN_COVERAGE')] <- 'mean';
	metric.data$lower_cov <- metric.data$mean - metric.data$SD_COVERAGE; 
	metric.data$upper_cov <- metric.data$mean + metric.data$SD_COVERAGE; 
	metric.data <- metric.data[,c('mean','lower_cov','upper_cov','percent_bases_above_targetX')];

	cov.method <- "Picard's WGSMetrics was used to estimate coverage.";

	} else if (is.dna) {

	target_coverage <- as.numeric(unlist(
		strsplit(colnames(metric.data)[grepl('percent_bases_above_', colnames(metric.data))],'_')
		)[4]);

	metric.data$percent_bases_above_targetX <- metric.data[,paste0('percent_bases_above_', target_coverage)]/100;

	colnames(metric.data)[which(colnames(metric.data) == 'granular_third_quartile')] <- 'upper_cov';
	colnames(metric.data)[which(colnames(metric.data) == 'granular_first_quartile')] <- 'lower_cov';
	metric.data <- metric.data[,c('mean','lower_cov','upper_cov','percent_bases_above_targetX')];

	cov.method <- "GATK's depthOfCoverage was used to estimate coverage.";
	} else {
	cov.method <- NULL;
	}

if (!is.null(arguments$hscoverage)) {

	hs.data <- read.delim(arguments$hscoverage, row.names = 1);
	hs.data$percent_on_bait <- hs.data$ON_BAIT_BASES / hs.data$PF_BASES_ALIGNED;
	hs.data$percent_off_bait <- hs.data$OFF_BAIT_BASES / hs.data$PF_BASES_ALIGNED;
	hs.data$percent_near_bait <- hs.data$NEAR_BAIT_BASES / hs.data$PF_BASES_ALIGNED;

	metric.data <- merge(metric.data, hs.data[,c('percent_on_bait','percent_off_bait','percent_near_bait')],
		by = 'row.names', all.x = TRUE);
	rownames(metric.data) <- metric.data$Row.names;
	metric.data <- metric.data[,-1];

	}

# read in contamination metrics
if (is.dna) {
	if (!is.null(contest.file)) {
		contest.data <- read.delim(contest.file);
		contest.data <- contest.data[which(contest.data$name == 'META'),];
		rownames(contest.data) <- contest.data$Sample;
		contest.data <- contest.data[,c('contamination','confidence_interval_95_low','confidence_interval_95_high')];
		contest.data$Method <- 'ContEst';
		}

	if (!is.null(other.contamination)) {
		pileup.contam <- read.delim(other.contamination);
		pileup.contam$contamination <- pileup.contam$contamination * 100;
		pileup.contam$confidence_interval_95_low <- pileup.contam$contamination - (pileup.contam$error*100);
		pileup.contam$confidence_interval_95_high <- pileup.contam$contamination + (pileup.contam$error*100);
		rownames(pileup.contam) <- pileup.contam$ID;
		pileup.contam <- pileup.contam[,c('contamination','confidence_interval_95_low','confidence_interval_95_high')];
		pileup.contam$Method <- 'pileup';
		}

	if (!is.null(other.contamination)) {
		contam.data <- pileup.contam; } else if (!is.null(contest.file)) {
		contam.data <- contest.data; } else {
		contam.data <- NULL;
		}

	contam.method <- if (!is.null(other.contamination)) {
		"GATK's GetPileupSummaries and CalculateContamination were used to estimate contamination."; } else if (!is.null(contest.file)) {
		"GATK's ContEst was used to estimate contamination for T/N pairs only."; } else {
		NULL;
		}
	}

# format data for plotting
if (is.dna) {
	qc.metrics <- merge(metric.data, contam.data, by = 'row.names', all.x = TRUE);
	qc.metrics <- qc.metrics[which(qc.metrics$Row.names %in% smp.names),];
	qc.metrics$Sample <- factor(qc.metrics$Row.names, levels = smp.names);
	} else {
	keep.metrics <- c('Total.Purity.Filtered.Reads.Sequenced','Mapping.Rate','Exonic.Rate','Intronic.Rate','Intergenic.Rate','rRNA.rate','Mean.Per.Base.Cov','Mean.CV','Duplication.Rate.of.Mapped','Base.Mismatch.Rate');
	qc.metrics <- metric.data[intersect(rownames(metric.data), smp.names),keep.metrics];
	qc.metrics$Sample <- factor(rownames(qc.metrics), levels = smp.names);
	qc.metrics$min <- qc.metrics$Mean.Per.Base.Cov - qc.metrics$Mean.CV/2;
	qc.metrics$max <- qc.metrics$Mean.Per.Base.Cov + qc.metrics$Mean.CV/2;
	}

qc.metrics <- qc.metrics[order(qc.metrics$Sample),];
qc.metrics$Order <- 1:nrow(qc.metrics);

if (!any(smp.names %in% rownames(cor.data))) {
        if ( (any(grepl('-', smp.names))) && (!any(grepl('-', rownames(cor.data)))) ) {
		new.smps <- gsub('-', '.', smp.names);
		}
	if ( (any(grepl("^[[:digit:]]", smp.names))) && (any(grepl('^X',rownames(cor.data)))) ) {
		new.smps <- paste0('X', new.smps);
		}
	if (!any(new.smps %in% rownames(cor.data))) {
		stop('Correlation sample names do not match names from other tables');
		}
        cor.data <- cor.data[new.smps,new.smps];
	colnames(cor.data) <- smp.names;
	rownames(cor.data) <- smp.names;
        } else {
        cor.data <- cor.data[smp.names,smp.names];
        }

# find line breaks (separate patients with multiple samples)
line.breaks <- get.line.breaks(sample.info$Patient);

# determine some parameters
add.rectangle <- TRUE;
if (length(smp.names) > 50) { add.rectangle <- FALSE; }
if (length(smp.names) < 6) { add.rectangle <- FALSE; }

axis.cex <- if (nrow(cor.data) <= 30) { 1
	} else if (nrow(cor.data) <= 80) { 0.8
	} else if (nrow(cor.data) <= 100) { 0.6
	} else if (nrow(cor.data) <= 140) { 0.5
	} else if (nrow(cor.data) <= 180) { 0.4
	} else { 0 };

if (is.dna) {

	max.cov <- max(qc.metrics$upper_cov, na.rm = TRUE);
	if ( (max.cov == 1) || (max.cov == 500) ) { max.cov <- max(qc.metrics$mean, na.rm = TRUE); }

	if (max.cov <= 100) {
		total.cov.limits <- c(0, 100);
		total.cov.at <- seq(0, 100, 25);
		} else if (max.cov <= 200) {
		total.cov.limits <- c(0, 200);
		total.cov.at <- seq(0, 200, 50);
		} else if (max.cov <= 400) {
		total.cov.limits <- c(0, 400);
		total.cov.at <- seq(0, 400, 100);
		} else {
		total.cov.limits <- c(0, ceiling(max.cov/10)*10);
		total.cov.at <- seq(0, total.cov.limits[2], length.out = 5);
		}

	max.contamination <- max(qc.metrics$contamination + qc.metrics$confidence_interval_95_high, na.rm = TRUE);
	if (max.contamination <= 3) {
		contest.limits <- c(0, 3);
		contest.at <- seq(0, 3, 1);
		} else if (max.contamination <= 8)  {
		contest.limits <- c(0, 8);
		contest.at <- seq(0, 8, 2);
		} else if (max.contamination <= 12)  {
		contest.limits <- c(0, 12);
		contest.at <- seq(0, 12, 4);
		} else {
		contest.limits <- NULL;
		contest.at <- TRUE;
		}

	} else {

	total.reads <- max(qc.metrics$Total.Purity.Filtered.Reads.Sequenced, na.rm = TRUE)/10**6;

	if (total.reads <= 100) {
		total.cov.limits <- c(0, 100);
		total.cov.at <- seq(0, 100, 20);
		} else if (total.reads <= 200) {
		total.cov.limits <- c(0, 200);
		total.cov.at <- seq(0, 200, 50);
		} else if (total.reads <= 300) {
		total.cov.limits <- c(0, 300);
		total.cov.at <- seq(0, 300, 100);
		} else if (total.reads <= 400) {
		total.cov.limits <- c(0, 400);
		total.cov.at <- seq(0, 400, 100);
		} else if (total.reads <= 600) {
		total.cov.limits <- c(0, 600);
		total.cov.at <- seq(0, 600, 200);
		} else if (total.reads <= 1000) {
		total.cov.limits <- c(0, 1000);
		total.cov.at <- seq(0, 1000, 200);
		} else {
		total.cov.limits <- NULL;
		total.cov.at <- TRUE;
		}

	max.cov <- max(qc.metrics$max , na.rm = TRUE);
	if (max.cov <= 1) {
		cov.limits <- c(0, 1);
		cov.at <- seq(0, 1, 0.2);
		} else if (max.cov <= 4) {
		cov.limits <- c(0, 4);
		cov.at <- seq(0, 4, 1);
		} else if (max.cov <= 10) {
		cov.limits <- c(0, 10);
		cov.at <- seq(0, 10, 2);
		} else if (max.cov <= 20) {
		cov.limits <- c(0, 20);
		cov.at <- seq(0, 20, 5);
		} else if (max.cov <= 50) {
		cov.limits <- c(0, 50);
		cov.at <- seq(0, 50, 10);
		} else {
		cov.limits <- NULL;
		cov.at <- TRUE;
		}
	}

# move to output directory
setwd(output.dir);

# plot results
if (is.dna) {

	# plot the results (correlations of germline genotypes)
	cor.heatmap <- create.heatmap(
		cor.data,
		cluster.dimensions = 'none',
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
		at = seq(0.5,1,0.001),
		colourkey.labels.at = c(0.5,0.75,0.9,0.95,1),
		colourkey.labels = c('< 0.5', '0.75','0.90', '0.95', '1'),
		grid.row = TRUE,
		force.grid.row = TRUE,
		row.lines = line.breaks,
		row.colour = 'grey80',
		col.colour = 'grey80',
		grid.col = TRUE,
		force.grid.col = TRUE,
		col.lines = line.breaks,
		axes.lwd = 1
		);

	# identify any suspect cases
	suspect.cases <- data.frame(matrix(ncol = 2, nrow = 0));
	colnames(suspect.cases) <- c('Sample','Suspect.Relation');
	tmp <- cor.data;
	tmp[which(tmp < 0.8, arr.ind = TRUE)] <- NA;
	colnames(tmp) <- gsub('\\.','-',colnames(tmp));
	rownames(tmp) <- gsub('\\.','-',rownames(tmp));
	for (patient in unique(sample.info$Patient)) {
		smps <- sample.info[which(sample.info$Patient == patient),]$Sample;
		for (smp in smps) {
			smp.data <- tmp[,smp];
			similar.smps <- setdiff(rownames(tmp)[!is.na(smp.data)], smps);
			if (length(similar.smps) > 0) {
				suspect.cases[nrow(suspect.cases)+1,] <- c(
					smp,
					paste(similar.smps, collapse = ', ')
					);
				}
			}
		}

	# plot coverage metrics
	qc.metrics$percent_target <- qc.metrics$percent_bases_above_targetX * total.cov.limits[2];
	key.lab <- if (!is.null(arguments$wgscoverage)) { 'mean read depth ( \u00B1 SD )';
		} else if (!is.null(arguments$hscoverage)) { 'mean target coverage';
		} else { 'mean read depth (1st - 3rd quartile)'; }

	cov.key <- list(fun = draw.key, args = list(key = list(
		lines = list(cex = 1, col = 'black', pch = 19, type = if (!is.null(arguments$hscoverage)) { 'p' } else { 'b' } ),
		text = list(lab = key.lab, cex = 1),
		divide = 1)),
		x = 0.5, y = 1
		);

	total.coverage.plot <- create.scatterplot(
		mean ~ Order,
		qc.metrics,
		xaxis.lab = rep('',nrow(qc.metrics)),
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
		y.error.up = qc.metrics$upper_cov - qc.metrics$mean,
		y.error.down = qc.metrics$mean - qc.metrics$lower_cov,
		error.whisker.angle = 0,
		abline.v = line.breaks,
		abline.col = 'grey80',
		add.rectangle = TRUE,
		xleft.rectangle = seq(0.5, length(smp.names)-0.5, 1),
		xright.rectangle = seq(1.5, length(smp.names)+0.5,1),
		ybottom.rectangle = -1,
		ytop.rectangle = qc.metrics$percent_target,
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

	# identify cases with poor coverage
	poor.coverage <- qc.metrics[which(qc.metrics$percent_bases_above_targetX < 0.8),c(1,2,5)];
	coverage.caption <- paste0(
		"Samples with low coverage ($<$80\\% bases with at least ", target_coverage,
		"x coverage)."
		);
	if (nrow(poor.coverage) > 20) {
		poor.coverage <- qc.metrics[which(qc.metrics$percent_bases_above_targetX < 0.5),c(1,2,5)];
		coverage.caption <- paste0(
			"Samples with low coverage ($<$50\\% bases with at least ",
			target_coverage, "x coverage)."
			);
		}
	colnames(poor.coverage) <- c('Sample','Mean Coverage','fraction_above_targetX');
	
	# if this is exome/target sequencing, plot the on/off/near target bases
	if (!is.null(arguments$hscoverage)) {

		stacked.data <- reshape(
			qc.metrics[,c('Sample','percent_on_bait','percent_off_bait','percent_near_bait')],
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

		efficiency.plot <- create.barplot(
			Percent ~ Sample,
			stacked.data,
			groups = stacked.data$Group,
			stack = TRUE,
			col = rev(default.colours(12))[1:3],
			xlab.label = NULL,
			ylab.label = '% Aligned Bases',
			ylab.cex = 1.5,
			xaxis.lab = NA,
			xaxis.rot = 90,
			xaxis.tck = c(0.5,0),
			yaxis.tck = c(0.5,0),
			xaxis.cex = axis.cex,
			yaxis.cex = 1,
			xaxis.fontface = 'plain',
			yaxis.fontface = 'plain',
			ylimits = c(0,1.1),
			yat = seq(0,1,0.2),
			legend = list(inside = stack.legend),
			top.padding = 1,
			style = 'Nature'
			);

		write.plot(
			efficiency.plot,
			filename = generate.filename(arguments$project, 'hs_efficiency_metrics','png'),
			height = 4,
			width = 11,
			resolution = 200
			);
		}

	# plot contamination estimates
	contest.plot <- create.segplot(
		Order ~ confidence_interval_95_low + confidence_interval_95_high,
		qc.metrics,
		yaxis.lab = rep('',nrow(qc.metrics)),
		yaxis.tck = 0,
		ylimits = c(0.5, length(smp.names)+0.5),
		yat = seq(1,length(smp.names)),
		xaxis.tck = c(0.5,0),
		xaxis.cex = 1.3,
		xlab.label = if (axis.cex > 0) {
			'Contamination\nEstimate (%)' } else {
			'\nContamination\nEstimate (%)' },
		xlab.cex = 1.5,
		ylab.label = NULL,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		xlimits = contest.limits,
		xat = contest.at,
		centers = qc.metrics$contamination,
		symbol.cex = max(c(axis.cex,0.2)),
		abline.h = line.breaks,
		abline.col = 'grey80',
		top.padding = 0.1,
		add.rectangle = add.rectangle,
		xleft.rectangle = -1,
		ybottom.rectangle = c(0.5, line.breaks),
		xright.rectangle = 100,
		ytop.rectangle = c(line.breaks, length(smp.names)+0.5),
		col.rectangle = c('white', 'grey90'),
		alpha.rectangle = 0.8,
		axes.lwd = 1
		);

	# identify cases with high contamination estimates
	high.contamination <- qc.metrics[which(qc.metrics$contamination > 3),c(1,6)];
	colnames(high.contamination) <- c('Sample','Contamination_Estimate');

	# combine them!
	create.multipanelplot(
		plot.objects = list(total.coverage.plot, cor.heatmap, contest.plot),
		filename = generate.filename(arguments$project, 'dnaseqc_metrics','png'),
		height = 12,
		width = 12,
		resolution = 200,
		layout.width = 2,
		layout.height = 2,
		layout.skip = c(FALSE,TRUE,FALSE,FALSE),
		plot.objects.widths = c(5,1),
		plot.objects.heights = c(1,4.5),
		left.legend.padding = 0,
		right.legend.padding = 1.5,
		top.legend.padding = 0,
		bottom.legend.padding = 0,
		xlab.axis.padding = -4,
		y.spacing = 1,
		x.spacing = -4,
		legend = list(
			inside = list(
				fun = draw.key,
				args = list(
					key = list(
						rect = list(col = 'red', alpha = 0.5, size = 1),
						text = list(
							lab = paste0('Percent bases\nabove ', 
								target_coverage, 'x'),
							cex = 1
							)
						)
					),
				x = 0.93,
				y = 0.9
				)
			)
		);

	} else {

	# plot the results (correlations of genes.rpkm.gct)
	cor.heatmap <- create.heatmap(
		cor.data,
		cluster.dimensions = 'none',
		xaxis.lab = smp.names,
		yaxis.lab = smp.names,
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
		colourkey.labels = c('< 0.8', '0.9', '0.95', '1'),
		grid.row = TRUE,
		force.grid.row = TRUE,
		row.lines = line.breaks,
		row.colour = 'grey80',
		col.colour = 'grey80',
		grid.col = TRUE,
		force.grid.col = TRUE,
		col.lines = line.breaks,
		axes.lwd = 1
		);

	# identify any suspect cases
	suspect.cases <- data.frame(matrix(ncol = 2, nrow = 0));
	colnames(suspect.cases) <- c('Sample','Suspect.Relation');
	tmp <- cor.data;
	tmp[which(tmp < 0.9, arr.ind = TRUE)] <- NA;
	colnames(tmp) <- gsub('^X','',gsub('\\.','-',colnames(tmp)));
	rownames(tmp) <- gsub('^X','',gsub('\\.','-',rownames(tmp)));
	for (patient in unique(sample.info$Patient)) {
		smps <- sample.info[which(sample.info$Patient == patient),]$Sample;
		for (smp in smps) {
			smp.data <- tmp[,smp];
			similar.smps <- setdiff(rownames(tmp)[!is.na(smp.data)], smps);
			if (length(similar.smps) > 0) {
				suspect.cases[nrow(suspect.cases)+1,] <- c(
					smp,
					paste(similar.smps, collapse = ', ')
					);
				}
			}
		}

	# plot region summary
	stacked.data <- data.frame(
		Order = rep(qc.metrics$Order, times = 5),
		ind   = rep(
			c('Unmapped','rRNA','Exonic','Intronic','Intergenic'),
			each = nrow(qc.metrics)
			),
		counts = NA
		);

	for (i in 1:nrow(qc.metrics)) {
		total <- qc.metrics[i,]$Total.Purity.Filtered.Reads.Sequenced/10**6;
		mapped <- qc.metrics[i,]$Mapping.Rate*total;

		unmapped <- total - mapped;
		rrna <- qc.metrics[i,]$rRNA.rate * total;
		exonic <- qc.metrics[i,]$Exonic.Rate * mapped;
		intronic <- qc.metrics[i,]$Intronic.Rate * mapped;
		intergenic <- qc.metrics[i,]$Intergenic.Rate * mapped;
		stacked.data[which(stacked.data$Order == i & stacked.data$ind == 'Unmapped'),]$counts <- unmapped;
		stacked.data[which(stacked.data$Order == i & stacked.data$ind == 'rRNA'),]$counts <- rrna;
		stacked.data[which(stacked.data$Order == i & stacked.data$ind == 'Exonic'),]$counts <- exonic;
		stacked.data[which(stacked.data$Order == i & stacked.data$ind == 'Intronic'),]$counts <- intronic;
		stacked.data[which(stacked.data$Order == i & stacked.data$ind == 'Intergenic'),]$counts <- intergenic;

		}

	stacked.data$ind <- factor(
		stacked.data$ind,
		levels = rev(c('Exonic','Intronic','Intergenic','rRNA','Unmapped'))
		);

	my.legend.grob <- legend.grob(
		legends = list(
			legend = list(
				colours = default.colours(6,'pastel')[-5],
				labels = c('Exonic','Intronic','Intergenic','rRNA','Unmapped/duplicate'),
				title = 'Region'
				)
			),
		title.just = 'left',
		size = 2
		);

	library.plot <- create.barplot(
		counts ~ Order,
		stacked.data,
		groups = stacked.data$ind,
		stack = TRUE,
		xaxis.lab = rep('',nrow(qc.metrics)),
		xlimits = c(0.5, length(smp.names)+0.5),
		xat = seq(1,length(smp.names)),
		yaxis.tck = c(0.5,0.5),
		xaxis.tck = 0,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		ylab.label = 'Total Reads',
		ylab.cex = 1.5,
		xlab.label = '',
		ylimits = total.cov.limits,
		yat = total.cov.at,
		yaxis.cex = 1.3,
		axes.lwd = 1,
		abline.v = line.breaks,
		abline.col = 'grey80',
		col = rev(default.colours(6, 'pastel')[-5]),
		top.padding = 1
		);

	# extract rate data
	rate.data <- qc.metrics[,grepl('Rate|rate',colnames(qc.metrics))];
	rate.data <- rate.data[which(
		rate.data$Mapping.Rate < quantile(rate.data$Mapping.Rate, probs = 0.1) |
		rate.data$Duplication.Rate.of.Mapped > quantile(rate.data$Duplication.Rate.of.Mapped, probs = 0.9) |
		rate.data$Base.Mismatch.Rate > quantile(rate.data$Base.Mismatch.Rate, probs = 0.9) |
		rate.data$rRNA.rate > quantile(rate.data$rRNA.rate, probs = 0.9)
		),];
	colnames(rate.data)[which(colnames(rate.data) == 'Duplication.Rate.of.Mapped')] <- 'Duplication.Rate';

	# plot coverage data
	cov.plot <- create.segplot(
		Order ~ min + max,
		qc.metrics,
		yaxis.lab = rep('',nrow(qc.metrics)),
		yaxis.tck = 0,
		ylimits = c(0.5, length(smp.names)+0.5),
		yat = seq(1,length(smp.names)),
		xaxis.tck = c(0.5,0),
		xaxis.cex = 1.3,
		xlab.label = 'Mean per-base\nCoverage',
		xlab.cex = 1.5,
		ylab.label = NULL,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		xlimits = cov.limits,
		xat = cov.at,
		centers = qc.metrics$Mean.Per.Base.Cov,
		symbol.cex = axis.cex,
		abline.h = line.breaks,
		abline.col = 'grey80',
		top.padding = 0.1,
		add.rectangle = add.rectangle,
		xleft.rectangle = -1,
		ybottom.rectangle = c(0.5, line.breaks),
		xright.rectangle = max.cov+10,
		ytop.rectangle = c(line.breaks, length(smp.names)+0.5),
		col.rectangle = c('white', 'grey90'),
		alpha.rectangle = 0.8,
		axes.lwd = 1
		);

	# identify cases with low coverage
	poor.coverage <- qc.metrics[which(qc.metrics$Mean.Per.Base.Cov < 10),c('Mean.Per.Base.Cov','Mean.CV')];

	# combine them!
	create.multipanelplot(
		plot.objects = list(library.plot, cor.heatmap, cov.plot),
		filename = generate.filename(arguments$project, 'rnaseqc_metrics','png'),
		height = 12,
		width = 12,
		resolution = 200,
		layout.width = 2,
		layout.height = 2,
		layout.skip = c(FALSE,TRUE,FALSE,FALSE),
		plot.objects.widths = c(5,1),
		plot.objects.heights = c(1,4.5),
		left.legend.padding = 0,
		right.legend.padding = 1.5,
		top.legend.padding = 0,
		bottom.legend.padding = 0,
		xlab.axis.padding = -4,
		y.spacing = 1,
		x.spacing = 0,
		legend = list(
			inside = list(
				fun = my.legend.grob,
				x = 0.9, y = 0.9
				)
			)
		);
	}

save(
	sample.info,
	qc.metrics,
	cor.data,	
	file = generate.filename(arguments$project, 'qc_plot_data','RData')
	);

write.table(
	sample.info,
	file = 'sample_info.txt',
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### LATEX ##########################################################################################
# if making output for a report
if (!is.null(arguments$report)) {

	# write a caption for this plot
	caption <- if (arguments$seq_type == 'rna') {
		"Top plot: Breakdown of total reads across the genome (including total mapped to exonic, intronic and intergenic regions, rRNA sequences and total unmapped and/or duplicate reads). Right plot: Mean per-base coverage across the middle 1000 expressed transcripts. Central heatmap: Pearson's correlation of RPKM values for all genes across all samples.";
		} else if (!is.null(arguments$wgscoverage)) {
		"Top plot: Summary of coverage across the genome. Points indicate mean read depth while red bars show percent of bases with $>$15x coverage. Right plot: Contamination estimate (as a percent $\\pm$ 95\\% CI if ContEst [available for T/N pairs only] or $\\pm$ error if GATK:Pileup was used). Central heatmap: Spearman's correlation ($\\rho$) of germline variant genotypes across all samples.";
		} else if (arguments$seq_type == 'targeted') {
		"Top plot: Summary of coverage across target bases. Points indicate mean read depth while red bars show percent of target bases with $>$200x coverage. Right plot: Contamination estimate (as a percent $\\pm$ 95\\% CI if ContEst [available for T/N pairs only] or $\\pm$ error if GATK:Pileup was used). Central heatmap: Spearman's correlation ($\\rho$) of germline variant genotypes across all samples.";
		} else if (arguments$seq_type == 'exome') {
		"Top plot: Summary of coverage across target bases. Points indicate median read depth (range from first to third quartiles) while red bars show percent of target bases with $>$15x coverage. Right plot: Contamination estimate (as a percent $\\pm$ 95\\% CI if ContEst [available for T/N pairs only] or $\\pm$ error if GATK:Pileup was used). Central heatmap: Spearman's correlation ($\\rho$) of germline variant genotypes across all samples.";
		} 

	# initiate report object
	tex.file <- paste0(arguments$report, '/qc_summary.tex');

	# create symlinks for plots
	unlink(paste0(arguments$report, '/', 'qc_metrics.png'));
	file.symlink(
		paste0(output.dir, '/', ifelse (is.dna,
			generate.filename(arguments$project, 'dnaseqc_metrics','png'),
			generate.filename(arguments$project, 'rnaseqc_metrics','png')
			)),
		paste0(arguments$report, '/', 'qc_metrics.png')
		);

	if (!is.null(arguments$hscoverage)) {
		unlink(paste0(arguments$report, '/', 'hs_efficiency_plot.png'));
		file.symlink(
			paste0(output.dir, '/', generate.filename(arguments$project, 'hs_efficiency_metrics','png')),
			paste0(arguments$report, '/', 'hs_efficiency_plot.png')
			);
		}

	# write for latex
	write("\\section{QC}", file = tex.file);
	write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
	write("\\begin{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\includegraphics[width=0.9\\textwidth]{",
		paste0(arguments$report, '/', 'qc_metrics.png'), '}'
		), file = tex.file, append = TRUE);
	write("\\end{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\caption{", caption, "}"
		), file = tex.file, append = TRUE);
	write("\\end{figure}\n", file = tex.file, append = TRUE);
	write("\\pagebreak\n", file = tex.file, append = TRUE);

	# write concerning metrics to file
	if (is.dna) {

		# record correlations
		write("\\subsection{Germline Correlations}", file = tex.file);
		if (nrow(suspect.cases) > 20) {
			write("Many cases of high germline correlations detected - this is likely an artefact of using a targeted panel.", file = tex.file, append = TRUE);
			} else if (nrow(suspect.cases) > 0) {
			suspect.cases <- xtable(
				suspect.cases,
				caption = "Samples with high levels of germline genotype concordance ($\\rho >$0.8); may suggest relatedness, cross-sample contamination or label mix-ups.",
				align = rep('c',ncol(suspect.cases)+1)
				);
			print(suspect.cases, file = tex.file, append = TRUE, include.rownames = FALSE, latex.environments = "");
			} else {
			write("No unexpected cases of germline concordance detected.", file = tex.file, append = TRUE);
			}

		# record coverage
		write("\\subsection{Coverage}", file = tex.file, append = TRUE);
		if (nrow(poor.coverage) > 0) {
			poor.coverage <- xtable(
				poor.coverage,
				caption = paste0(coverage.caption, " Low coverage of either a tumour or normal sample may affect variant callablility."),
				align = rep('c',ncol(poor.coverage)+1),
				digits = c(0,0,1,1)
				);
			print(poor.coverage, file = tex.file, append = TRUE, include.rownames = FALSE, latex.environments = "");
			} else {
			write("No samples had unusually low coverage.", file = tex.file, append = TRUE);
			}

		# record contamination
		write("\\subsection{Contamination}", file = tex.file, append = TRUE);
		if (!is.null(contam.method)) {
			write(contam.method, file = tex.file, append = TRUE);
			}
		if (nrow(high.contamination) > 0) {
			high.contamination <- xtable(
				high.contamination,
				caption = "Samples with estimated high levels of cross-sample contamination ($>$3\\%). It is suggested that samples with a cross-sample contamination beyond this threshold be removed from downstream analyses.",
				align = rep('c', ncol(high.contamination)+1)
				);
			print(high.contamination, file = tex.file, append = TRUE, include.rownames = FALSE, latex.environments = "");
			} else {
			write("No samples had unusually high cross-sample contamination estimates.", file = tex.file, append = TRUE);
			}

		# add HS efficiency if available
		if (!is.null(arguments$hscoverage)) {
			write("\\pagebreak", file = tex.file, append = TRUE);
			write("\\subsection{HS Efficiency}", file = tex.file, append = TRUE);
			write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
			write("\\begin{center}", file = tex.file, append = TRUE);
			write(paste0(
				"\\includegraphics[width=0.9\\textwidth]{",
				paste0(arguments$report, '/', 'hs_efficiency_plot.png'), '}'
				), file = tex.file, append = TRUE);
			write("\\end{center}", file = tex.file, append = TRUE);
			write(paste0(
				"\\caption{Fraction of bases aligned to on-target, near-target or off-target regions as determined using Picard's CollectHsMetrics tool.}"
				), file = tex.file, append = TRUE);
			write("\\end{figure}\n", file = tex.file, append = TRUE);
			write("\\pagebreak\n", file = tex.file, append = TRUE);
			}

		# for RNA
		} else {

		# record correlations
		write("\\subsection{Sample Correlations}", file = tex.file);
		if (nrow(suspect.cases) > 0 & nrow(suspect.cases) <= 10) {
			suspect.cases <- xtable(
				suspect.cases,
				caption = "Samples with high levels of concordance across the transcriptome (Pearson's correlation $>$0.9); may suggest relatedness, cross-sample contamination or label mix-ups.",
				align = rep('c',ncol(suspect.cases)+1)
				);
			print(suspect.cases, file = tex.file, append = TRUE, include.rownames = FALSE, latex.environments = "");
			} else if (nrow(suspect.cases) > 10) {
			write("Cohort is highly correlated; not listing highly concordant samples.", file = tex.file, append = TRUE);
			} else {
			write("No unexpected cases of transcriptome concordance detected.", file = tex.file, append = TRUE);
			}

		# record coverage
		write("\\subsection{Coverage by Base}", file = tex.file, append = TRUE);
		if (nrow(poor.coverage) > 0 & nrow(poor.coverage) <= 10) {
			poor.coverage <- xtable(
				poor.coverage,
				caption = "Samples with low coverage (mean per-base coverage $<$10x); table shows mean per-base coverage and mean coefficient of variation of the middle 1000 expressed genes.",
				align = rep('c',ncol(poor.coverage)+1),
				digits = c(0,2,2)
				);
			print(poor.coverage, file = tex.file, append = TRUE, latex.environments = "");
			} else if (nrow(poor.coverage) > 10) {
			poor.coverage <- xtable(
				poor.coverage[1:10,],
				caption = "Subset of samples with low coverage (mean per-base coverage $<$10x); table shows mean per-base coverage and mean coefficient of variation of the middle 1000 expressed genes.",
				align = rep('c',ncol(poor.coverage)+1),
				digits = c(0,2,2)
				);
			print(poor.coverage, file = tex.file, append = TRUE, latex.environments = "");
			} else {
			write("No samples had unusually low coverage.", file = tex.file, append = TRUE);
			}

		# record region metrics
		write("\\pagebreak\n\\subsection{Coverage by Region}", file = tex.file, append = TRUE);
		if (nrow(rate.data) > 0) {
			if (nrow(rate.data) > 15) { rate.data <- rate.data[1:15,]; }
			rate.data1 <- xtable(
				rate.data[,c('Mapping.Rate','Duplication.Rate','Base.Mismatch.Rate','rRNA.rate')],
				align = rep('c',5),
				digits = c(0,2,2,2,2)
				);

			rate.data2 <- xtable(
				rate.data[,!colnames(rate.data) %in% c('Mapping.Rate','Duplication.Rate','Base.Mismatch.Rate','rRNA.rate')],
				caption = "Samples detected as possible outliers due to one or more of the following criteria: low rate of mapped reads; high duplication rate, high base mismatch rate or high rRNA rate (from total reads). Per-region rates for the above samples (per mapped unique reads).",
				align = rep('c',4),
				digits = c(0,2,2,2)
				);

			print(rate.data1, file = tex.file, append = TRUE, latex.environments = "");
			print(rate.data2, file = tex.file, append = TRUE, latex.environments = "", table.placement = 'h!');
			} else {
			write("No samples had unusual coverage metrics.", file = tex.file, append = TRUE);
			}
		}
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('PlotQCdata','SessionProfile','txt'));
