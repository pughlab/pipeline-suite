### apply_cosmic_mutation_signatures.R #############################################################
# apply COSMIC mutation signatures to ensemble SNVs

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
	cat(paste0('Input MAF: ', arguments$maf));
	cat(paste0('VAF threshold applied: ', vaf.threshold));
	cat(paste0('Minimum number of SNPs required: ', min.snps));

	# write sessionInfo to file
	cat("\n### SESSION INFO ###############################################################\n");
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
parser$add_argument('-t', '--seq_type', type = 'character', help = 'exome or wgs', default = 'exome');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'hg38 or hg19', default = 'hg38');
parser$add_argument('-m', '--maf', type = 'character', help = 'mutation calls in MAF format');
parser$add_argument('-v', '--vaf_threshold', type = 'character',
	help = 'threshold to filter variants', default = 0.1);
parser$add_argument('-s', '--signatures', type = 'character', default = NULL,
	help = 'path to signatures matrix (96 trinucleotide x N sigs)');
parser$add_argument('-z', '--report', type = 'character', help = 'path to report directory',
	default = NULL);

arguments <- parser$parse_args();

# load required libraries
library(BoutrosLab.plotting.general);
library(xtable);
library(deconstructSigs);
library(BSgenome);

# get genome object
BSgenome.Hsapiens.UCSC <- if (arguments$ref_type == 'hg38') {
	getBSgenome("BSgenome.Hsapiens.UCSC.hg38");
	} else if (arguments$ref_type == 'hg19') {
	getBSgenome("BSgenome.Hsapiens.UCSC.hg19");
	}

# indicate normalization method (assuming signatures matrix was developed from WGS)
normalization.method <- 'default';
if (arguments$seq_type %in% c('exome','targeted')) {
	normalization.method <- 'exome2genome';
	}

# indicate filtering thresholds
vaf.threshold <- arguments$vaf_threshold;
min.snps <- 50;

### READ DATA ######################################################################################
# get data
if (is.null(arguments$maf)) {
	stop('ERROR: No input MAF provided, please provide path to SNV calls in MAF format.');
	} else {
	input.maf <- read.delim(arguments$maf, stringsAsFactors = FALSE, comment.char = '#');
	}

# read in mutation signatures
signatures.to.apply <- if (is.null(arguments$signatures)) { signatures.cosmic } else {
	as.data.frame(t(read.delim(arguments$signatures, row.names = 1)));
	}

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

# run unlink, in case it exists from a previous run
if (!is.null(arguments$signatures)) {
	unlink( basename(arguments$signatures) );
	file.symlink(arguments$signatures, basename(arguments$signatures));
	}

### FORMAT DATA ####################################################################################
# remove INDELs
input.maf <- input.maf[which(input.maf$Variant_Type == 'SNP'),];

# calculate VAF
input.maf$VAF <- 1 - (input.maf$t_ref_count / input.maf$t_depth);

# remove low VAF (if threshold > 0)
input.maf <- input.maf[which(input.maf$VAF > vaf.threshold),];

# indicate key fields
snv.data <- input.maf[,c('Tumor_Sample_Barcode','Chromosome','Start_Position','Reference_Allele','Allele')];
colnames(snv.data) <- c('Sample','chr','pos','ref','alt');

# get list of samples
all.samples <- sort(unique(snv.data$Sample));

# convert to signatures object
sigs.input <- mut.to.sigs.input(
	mut.ref = snv.data,
	bsg = BSgenome.Hsapiens.UCSC
	);

# get signature weights for each sample
assigned.sigs <- list();

for (smp in all.samples) {

	if (!dir.exists(smp)) {
		dir.create(smp);
		}

	setwd(smp);

	if (nrow(snv.data[which(snv.data$Sample == smp),]) < 50) {
		assigned.sigs[[smp]]$weights <- rep(NA, ncol(signatures.to.apply));
		next;
		}

	assigned.sigs[[smp]] <- whichSignatures(
		tumor.ref = sigs.input,
		signatures.ref = signatures.to.apply,
		sample.id = smp,
		contexts.needed = TRUE,
		tri.counts.method = normalization.method
		);

	sig.weights.ordered <- sort(assigned.sigs[[smp]]$weights, decreasing = TRUE);

	# organize signature weights (scores) for plotting
	plot.data <- data.frame(
		sig = colnames(sig.weights.ordered),
		score = as.numeric(sig.weights.ordered)
		);

	plot.data$sig <- factor(plot.data$sig, levels = c(colnames(sig.weights.ordered),'Unknown'));
	plot.data <- rbind(plot.data, c('Unknown', as.numeric(assigned.sigs[[smp]]$unknown)));
	plot.data$score <- as.numeric(plot.data$score);

	# organize mutation data (trinucleotide counts) for plotting
	tnc.data <- data.frame(
		context = colnames(assigned.sigs[[smp]]$tumor),
		tumour = as.numeric(assigned.sigs[[smp]]$tumor),
		reconstructed = as.numeric(assigned.sigs[[smp]]$product)
		);

	tnc.data$Group <- as.character(
		sapply(as.character(tnc.data$context), function(i) { substr(i,3,5) } ));
	tnc.data$Group <- factor(
		tnc.data$Group, levels = rev(c('T>G','T>C','T>A','C>T','C>G','C>A')));
	tnc.data <- tnc.data[order(tnc.data$Group, tnc.data$context),];
	tnc.data$context <- factor(tnc.data$context, levels = as.character(tnc.data$context));

	sse <- sum(assigned.sigs[[smp]]$diff)^2;

	stats <- list(text = list(lab = c(paste0('SSE = ', round(sse,3))),cex = 1));

	sbs.plot <- create.barplot(
		score ~ sig,
		plot.data,
		xlab.label = NULL,
		ylab.label = 'Weight',
		ylab.cex = 1,
		xaxis.cex = 0.6,
		yaxis.cex = 1,
		axes.lwd = 1,
		xaxis.rot = 90,
		yaxis.fontface = 'plain',
		xaxis.fontface = 'plain',
		yaxis.tck = c(0.5,0),
		xaxis.tck = c(0.5,0)
		);

	tnc.plot <- create.barplot(
		reconstructed ~ context,
		tnc.data,
		xlab.label = NULL,
		ylab.label = 'Fraction',
		ylab.cex = 1,
		xaxis.cex = 0.6,
		yaxis.cex = 1,
		axes.lwd = 1,
		xaxis.rot = 90,
		yaxis.fontface = 'plain',
		xaxis.fontface = 'plain',
		yaxis.tck = c(0.5,0),
		xaxis.tck = c(0.5,0),
		col = rep(rev(default.colours(8,'pastel')[-5])[-1],each = 16),
		add.rectangle = TRUE,
		xleft.rectangle = c(0.5,get.line.breaks(tnc.data$Group)),
		xright.rectangle = c(get.line.breaks(tnc.data$Group),nrow(tnc.data)+0.5),
		ytop.rectangle = 1,
		ybottom.rectangle = 0,
		col.rectangle = rev(default.colours(8,'pastel')[-5])[-1],
		alpha.rectangle = 0.2
		);

	create.multipanelplot(
		plot.objects = list(sbs.plot, tnc.plot),
		layout.width = 1,
		layout.height = 2,
		plot.objects.heights = c(1.8,3),
		y.spacing = -2.8,
		left.legend.padding = 0,
		right.legend.padding = 0,
		bottom.legend.padding = 0,
		top.legend.padding = -2,
		legend = list(
			inside = list(fun = draw.key, args = list(key = stats),x = 0.13, y = 0.58)
			),
		height = 4,
		width = 13,
		resolution = 600,
		filename = generate.filename(smp, 'mutation_counts_per_signature', 'png')
		);

	setwd(arguments$output);
	}

# extract weights
sig.weights.list <- list();
for (smp in all.samples) {
	sig.weights.list[[smp]] <- assigned.sigs[[smp]]$weights;
	}

sig.weights <- do.call(rbind, sig.weights.list);

# do some error checking
if (all(is.na(sig.weights))) {
	warning('All signature weights are NA - sample(s) probably have too few mutations.');
	should.plot <- FALSE;
	} else {

	save(
		signatures.to.apply,
		assigned.sigs,
		file = generate.filename(arguments$project, 'assignedSignatures', 'RData')
		);

	write.table(
		sig.weights,
		file = generate.filename(arguments$project, 'mutation_signature_weights', 'tsv'),
		row.names = TRUE,
		col.names = NA,
		sep = '\t'
	 	);

	# make a symlink for use with final summary plots
	unlink('mutation_signature_weights.tsv');
	file.symlink(
		generate.filename(arguments$project, 'mutation_signature_weights', 'tsv'),
		'mutation_signature_weights.tsv'
		);

	# order samples by recurrence for prettier heatmap
	heatmap.data <- t(sig.weights);
	heatmap.data[which(heatmap.data == 0)] <- NA;
	heatmap.data[!is.na(heatmap.data)] <- 1;

	heatmap.data <- heatmap.data[names(sort(apply(heatmap.data,1,sum,na.rm = T), decreasing = T)),];

	heatmap.data <- t(heatmap.data);
	heatmap.data <- heatmap.data[do.call(order, transform(heatmap.data)),];

	sample.order <- rownames(heatmap.data);
	sig.order <- colnames(heatmap.data);

	# extract top signatures
	top.signatures <- sig.order[1:10];
	sig.summary <- data.frame(
		Signature = top.signatures,
		N.Samples = NA
		);	

	for (i in 1:10) {
		sig <- top.signatures[i];
		sig.summary[i,2] <- nrow(sig.weights[which(sig.weights[,sig] > 0.05),]);
		}

	should.plot <- TRUE;
	}

### PLOTTING #######################################################################################
if (should.plot) {

	# grab some parameters
	axis.cex <- if (length(all.samples) <= 30) { 1
		} else if (length(all.samples) <= 50) { 0.75
		} else if (length(all.samples) <= 80) { 0.5
		} else { 0 };

	# make the heatmap
	create.heatmap(
		sig.weights[sample.order,sig.order],
		same.as.matrix = TRUE,
		cluster.dimensions = 'none',
		colour.scheme = c('white','red'),
		xaxis.cex = 0.5,
		xaxis.lab.top = colnames(heatmap.data),
		x.alternating = 2,
		yaxis.cex = axis.cex,
		yaxis.lab = simplify.ids(rownames(heatmap.data)),
		axes.lwd = 1,
		yaxis.fontface = 'plain',
		xaxis.fontface = 'plain',
		yaxis.tck = c(0.2,0),
		right.padding = 1,
		grid.row = TRUE,
		force.grid.row = TRUE,
		row.colour = 'grey80',
		col.colour = 'grey80',
		row.lwd = if (length(all.samples) < 30) { 3 } else { 1 },
		col.lwd = if (length(all.samples) < 30) { 3 } else { 1 },
		grid.col = TRUE,
		force.grid.col = TRUE,
		fill.colour = 'white',
		colourkey.cex = 1,
		at = seq(0,1,0.001),
		height = if (length(all.samples) < 10) { 5 } else { 7 },
		width = 12,
		resolution = 600,
		filename = generate.filename(arguments$project, 'mutation_signatures', 'png')
		);
	}

### LATEX ##########################################################################################
# if making output for a report
if (!is.null(arguments$report)) {

	tex.file <- paste0(
		arguments$report,
		'/mutation_signature_summary.tex'
		);

	# write for latex
	write("\\section{Mutation Signatures}", file = tex.file);

	if (!should.plot) {

		write(
			'All signature weights are NA - sample(s) probably have too few mutations.',
			file = tex.file, append = TRUE
			);

		} else {
		sig.version <- if (is.null(arguments$signatures)) {
			'v2; \\url{https://cancer.sanger.ac.uk/signatures/signatures\\_v2/}'
			} else if (grepl('3.2', basename(arguments$signatures))) {
			'v3.2; \\url{https://cancer.sanger.ac.uk/signatures/sbs/}';
			} else {
			paste0(gsub('_','\\_',basename(arguments$signatures)));
			}

		write(
			paste0('COSMIC SBS mutation signatures (', sig.version, ') were applied to the current dataset using ENSEMBLE mutations, with minimum a VAF threshold = ', vaf.threshold, '.'),
			file = tex.file,
			append = TRUE
			);

		# run unlink, in case it exists from a previous run
		unlink('mutation_signatures.png');
		file.symlink(
			paste0(arguments$output, '/', 
				generate.filename(arguments$project, 'mutation_signatures', 'png')),
			paste0(arguments$report, '/', 'mutation_signatures.png')
			);

		# add mutation_signatures plot
		write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
		write("\\begin{center}", file = tex.file, append = TRUE);
		write(paste0(
			"\\includegraphics[height=0.3\\textheight]{",
			paste0(arguments$report, '/', 'mutation_signatures.png'), '}'
			), file = tex.file, append = TRUE);
		write("\\end{center}", file = tex.file, append = TRUE);
		write(paste0(
			"\\caption{COSMIC mutation signatures were applied to the current dataset using deconstructSig. Heatmap shows weights attributed to each signature for each sample.}"
			), file = tex.file, append = TRUE);
		write("\\end{figure}\n", file = tex.file, append = TRUE);

		# add table summary
		caption <- 'Summary of top 10 recurrent signatures across the cohort. Table shows number of samples with signature weight $>$ 0.05.';
		print(
			xtable(
				sig.summary,
				caption = caption 
				),
			file = tex.file,
			include.rownames = FALSE,
			append = TRUE
			);
		}
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('ApplyMutSigs','SessionProfile','txt'));
