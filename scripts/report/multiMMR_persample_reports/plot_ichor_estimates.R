### plot_ichor_metrics.R ###########################################################################
# Plot sWGS ichorCNA tumour fraction estimates

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
parser$add_argument('-i', '--ichor_dir', type = 'character', help = 'path to ichorCNA directory');
parser$add_argument('-r', '--report_dir', type = 'character', help = 'path to report directory');
parser$add_argument('-s', '--sample_yaml', type = 'character', help = 'path to sample (BAM) yaml');
parser$add_argument('-t', '--tumour_fraction', type = 'character', 
	help = 'path to ichorCNA tumour fraction estimates');

arguments <- parser$parse_args();

# load libraries
library(BoutrosLab.plotting.general);
library(yaml);
library(xtable);

# what's the date?
date <- Sys.Date();

# get input files
tf.file <- arguments$tumour_fraction;
ichor.dir <- arguments$ichor_dir;
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


# read in TF estimates
tf.data <- read.delim(tf.file);

# determine some parameters
axis.cex <- if (nrow(tf.data) <= 30) { 1
	} else if (nrow(tf.data) <= 80) { 0.8
	} else if (nrow(tf.data) <= 100) { 0.6
	} else if (nrow(tf.data) <= 140) { 0.5
	} else if (nrow(tf.data) <= 180) { 0.4
	} else { 0 };

# move to output directory
setwd(output.dir);

sex.colours <- c('tomato1','steelblue2');
sex.key <- list(
	legend = list(
		colours = 'white',
		border = 'transparent',
		labels = 'Inferred Sex:',
		size = 0.1
		),
	legend = list(
		colours = sex.colours[1],
		labels = 'female',
		size = 1.5
		),
	legend = list(
		colours = sex.colours[2],
		labels = 'male',
		size = 1.5
		)
	);

my.legend <- legend.grob(
	legends = sex.key,
	layout = c(3,1)
	);

# plot TF estimates
tf.plot <- create.barplot(
	Tumour.Fraction ~ ID,
	tf.data,
	col = sex.colours[match(tf.data$Sex, c('female','male'))],
	xlab.label = NULL,
	ylab.label = 'Tumour Fraction',
	ylab.cex = 1.5,
	ylab.axis.padding = 2,
	right.padding = 0,
	xaxis.lab = simplify.ids(tf.data$ID),
	xaxis.rot = 90,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xaxis.cex = axis.cex,
	yaxis.cex = 1,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	top.padding = 3,
	legend = list(inside = list(fun = my.legend, x = 0, y = 1.15)),
	abline.h = 0.1,
	abline.lty = 2,
	style = 'Nature'
	);

write.plot(
	tf.plot,
	filename = generate.filename(arguments$project, 'tumour_fraction_estimates','png'),
	height = 4,
	width = 11,
	resolution = 200
	);

save(
	sample.info,
	tf.data,
	file = generate.filename(arguments$project, 'ichorCNA_tf_data','RData')
	);

### PER-SAMPLE SUMMARIES ###########################################################################
# find all sWGS files
all.files <- list.files(pattern = 'ichorCNA_tf_data.RData');

orig.file <- generate.filename(arguments$project, 'ichorCNA_tf_data','RData');
all.files <- setdiff(all.files, orig.file);

orig.sample.info <- sample.info;
orig.metric.data <- tf.data;

combined.metrics <- tf.data;

for (file in all.files) {
	load(file);
	combined.metrics <- unique(rbind(combined.metrics, metric.data));
	}

# calculate kernal density for smoothed histogram
tf.density <- data.frame(
	x = density(combined.metrics$Tumour.Fraction)$x,
	y = density(combined.metrics$Tumour.Fraction)$y
	);

label.positions <- seq(0.95,0.7,-0.03);
get.label.pos <- function(x) { 
	y <- x + 0.04;
	if (any(x >= 1)) { y[which(x >= 1)] <- x[which(x >= 1)] }
	return(y);
	}

setwd(reports.dir);

# initiate report object
tex.file <- 'ichor_estimates_sWGS.tex';

for (patient in unique(orig.sample.info$Patient)) {

	if (!dir.exists(patient)) {
		dir.create(patient);
		}

	setwd(patient);

	these.smps <- orig.sample.info[which(orig.sample.info$Patient == patient),]$Sample;
	smp.data <- combined.metrics[which(combined.metrics$ID %in% these.smps),];
	smp.data <- smp.data[order(smp.data$Tumour.Fraction),];

	to.write <- smp.data[,c('ID','Sex','Tumour.Fraction')];
	colnames(to.write) <- c('Sample','Inferred Sex','Tumour Fraction');

	to.write[,3] <- round(to.write[,3],3);

	# add in historical data
	to.write[nrow(to.write)+1,] <- c('Historical Average','',0);
	to.write[nrow(to.write),3] <- format.mean(combined.metrics$Tumour.Fraction);

	caption <- "Tumour fraction estimates by sWGS:ichorCNA for this sample, compared to historical average (all previously processed samples; reported as median + IQR).";

	# plot per-sample TF estimates
	smp.plot <- create.polygonplot(
		formula = NA ~ x,
		data = tf.density,
		max = tf.density$y,
		min = rep(0,nrow(tf.density)),
		col = 'grey80',
		xlab.label = 'Tumour Fraction',
		ylab.label = 'Population Frequency',
		xlab.cex = 1.5,
		ylab.cex = 1.5,
		xlimits = c(-0.05,1.05),
		xat = seq(0,1,0.2),
		ylimits = c(0,ceiling(max(tf.density$y))),
		xaxis.tck = c(0.5,0),
		yaxis.tck = c(0.5,0),
		xaxis.cex = 1,
		yaxis.cex = 1,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		top.padding = 3,
		bottom.padding = 0,
		right.padding = 0,
		left.padding = 0,
		ylab.axis.padding = 2,
		add.text = TRUE,
		text.labels = simplify.ids(smp.data$ID),
		text.x = get.label.pos(smp.data$Tumour.Fraction), 
		text.y = ceiling(max(tf.density$y))*label.positions[1:nrow(smp.data)],
		text.cex = 0.8,
		text.fontface = 'plain',
		legend = list(inside = list(fun = my.legend, x = 0, y = 1.1)),
		abline.v = smp.data$Tumour.Fraction,
		abline.lty = 1,
		abline.col = sex.colours[match(smp.data$Sex, c('female','male'))],
		xgrid.at = NULL,
		ygrid.at = NULL,
		style = 'Nature'
		);

	write.plot(
		smp.plot,
		filename = generate.filename(arguments$project, 'tumour_fraction_estimates','png'),
		height = 4,
		width = 9,
		resolution = 200
		);	

	# add data to tex file
	write("\\subsubsection{ichorCNA}\n", file = tex.file);

	print(
		xtable(
			to.write,
			caption = caption
			),
		file = tex.file,
		append = TRUE,
		include.rownames = FALSE,
		latex.environments = ""
		);

	write("", file = tex.file, append = TRUE);
	write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
	write("\\begin{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\includegraphics[width=0.9\\textwidth]{",
		generate.filename(arguments$project, 'tumour_fraction_estimates','png'), '}'
		), file = tex.file, append = TRUE);
	write("\\end{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\caption{", caption, "}"
		), file = tex.file, append = TRUE);
	write("\\end{figure}\n", file = tex.file, append = TRUE);
	write("\\pagebreak\n", file = tex.file, append = TRUE);

	ichor.plots <- list.files(
		path = paste0(ichor.dir,'/',patient),
		pattern = 'genomeWide.pdf', recursive = TRUE, full.names = TRUE);

	for (plot in ichor.plots) {

		write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
		write("\\begin{center}", file = tex.file, append = TRUE);
		write(paste0(
			"\\includegraphics[width=1.1\\textwidth]{",
			plot, '}'
			), file = tex.file, append = TRUE);
		write("\\end{center}", file = tex.file, append = TRUE);
		write("\\end{figure}\n", file = tex.file, append = TRUE);

		}

	setwd(reports.dir);
	}

### SAVE SESSION INFO ##############################################################################
setwd(output.dir);
save.session.profile(generate.filename(arguments$project, 'PlotTF_SessionProfile','txt'));
