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
	help = 'path to signatures matrix (96 trinucleotide x N sigs); takes priority over --cosmic_version');
parser$add_argument('-c', '--cosmic_version', type = 'character', default = 3.2,
	help = 'COSMIC SBS signature version to use (see R package: cosmicsig)');
parser$add_argument('-z', '--report', type = 'character', help = 'path to report directory',
	default = NULL);

arguments <- parser$parse_args();

# load required libraries
library(BoutrosLab.plotting.general);
library(cosmicsig);
library(xtable);
library(deconstructSigs);
library(BSgenome);

# get genome object
ref_type <- if (arguments$ref_type == 'hg38') { 'GRCh38'
	} else if (arguments$ref_type == 'hg19') { 'GRCh37'
	}

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
if (!is.null(arguments$signatures)) {
	signatures.to.apply <- as.data.frame(t(read.delim(arguments$signatures, row.names = 1)));
	signature.etiologies <- NULL;
	} else if (arguments$cosmic_version %in% c('3.0','3.1','3.2','3.3')) {
	signatures.to.apply <- get(paste0('COSMIC_v',arguments$cosmic_version))$signature[[ref_type]]$SBS96;
	signatures.to.apply <- as.data.frame(t(signatures.to.apply));
	colnames(signatures.to.apply) <- sapply(colnames(signatures.to.apply), function(i) {
		i <- unlist(strsplit(i,''));
		paste0(i[1],'[',i[2],'>',i[4],']',i[3])
		} );
	signature.etiologies <- data.frame(get(paste0('COSMIC_v',arguments$cosmic_version))$etiology$SBS96);
	} else {
	signatures.to.apply <- deconstructSigs::signatures.cosmic;
	signature.etiologies <- NULL;
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
		signature.etiologies,
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

	plot.data <- merge(t(sig.weights), signature.etiologies, by = 'row.names');
	plot.data$Count <- apply(plot.data[,sample.order],1,function(i) { length(i[which(i > 0)]) } );

	plot.data$Group <- 'other';
	plot.data[grepl('exposure',tolower(plot.data$proposed_etiology)),]$Group <- 'other exposure';
	plot.data[grepl('ultraviolet',tolower(plot.data$proposed_etiology)),]$Group <- 'UV';
	plot.data[grepl('tobacco',tolower(plot.data$proposed_etiology)),]$Group <- 'tobacco';
	plot.data[grepl('AID|APOBEC',plot.data$proposed_etiology),]$Group <- 'APOBEC';
	plot.data[grepl('treatment|chemotherapy',tolower(plot.data$proposed_etiology)),]$Group <- 'treatment';
	plot.data[grepl('deamination',tolower(plot.data$proposed_etiology)),]$Group <- 'age';
	plot.data[grepl('excision',tolower(plot.data$proposed_etiology)),]$Group <- 'BER';
	plot.data[grepl('mismatch repair',tolower(plot.data$proposed_etiology)),]$Group <- 'MMR';
	plot.data[grepl('HR|BRCA',plot.data$proposed_etiology),]$Group <- 'HRD';
	plot.data[grepl('POLD|POLE|polymerase',plot.data$proposed_etiology),]$Group <- 'polymerase';
	plot.data[grepl('reactive',tolower(plot.data$proposed_etiology)),]$Group <- 'ROS';
	plot.data[grepl('artefact',tolower(plot.data$proposed_etiology)),]$Group <- 'junk';
	plot.data[which(plot.data$proposed_etiology == ''),]$Group <- 'unknown';
	plot.data$Group <- factor(plot.data$Group, levels = rev(c('age','APOBEC','BER','HRD','MMR','polymerase','ROS','tobacco','treatment','other exposure','UV','other','unknown','junk')));

	plot.data <- plot.data[order(plot.data$Group, plot.data$Count),];

	# extract top signatures
	top.signatures <- sig.order[1:10];
	sig.summary <- data.frame(
		Signature = top.signatures,
		N.Samples = apply(sig.weights[,top.signatures],2,function(i) { length(i[which(i > 0.05)]) } ),
		Etiology = if (!is.null(signature.etiologies)) { signature.etiologies[top.signatures,1]
			} else { rep('', 10) }
		);

	should.plot <- TRUE;
	}

### PLOTTING #######################################################################################
if (should.plot) {

	# grab some parameters
	axis.cex <- if (length(all.samples) <= 30) { 1
		} else if (length(all.samples) <= 50) { 0.75
		} else if (length(all.samples) <= 80) { 0.5
		} else { 0 };

	# make some covariates
	sig.colours <- force.colour.scheme(scheme = 'chromosome', return.scheme = TRUE)$scheme$colours;
	sig.colours <- rev(sig.colours[1:nlevels(plot.data$Group)]);

	sig.covariates <- list(
		rect = list(col = 'white', lwd = 0, 
			fill = rev(sig.colours[match(plot.data$Group, levels(plot.data$Group))])
			)
		);

	sig.legends <- legend.grob(
		legends = list(
			legend = list(colours = rev(sig.colours), labels = rev(levels(plot.data$Group)))
			),
		label.cex = 0.8,
		size = 1.5
		);

	# make the heatmap
	create.heatmap(
		sig.weights[sample.order,rev(plot.data$Row.names)],
		same.as.matrix = TRUE,
		cluster.dimensions = 'none',
		colour.scheme = c('white','red'),
		covariates.top = sig.covariates,
		covariates.top.grid.border = list(col = 'black', lwd = 1),
		covariates.top.col.lines = get.line.breaks(rev(plot.data$Group))-0.5,
		covariates.top.grid.col = list(col = 'grey80', lwd = 2),
		inside.legend = list(fun = sig.legends, x = 1.01, y = 1),
		right.padding = 15,
		xaxis.cex = 0.5,
		xaxis.lab = rev(plot.data$Row.names),
		yaxis.cex = axis.cex,
		yaxis.lab = simplify.ids(rownames(heatmap.data)),
		axes.lwd = 1,
		yaxis.fontface = 'plain',
		xaxis.fontface = 'plain',
		yaxis.tck = c(0.2,0),
		grid.row = TRUE,
		force.grid.row = TRUE,
		row.colour = 'grey80',
		col.colour = 'grey80',
		row.lwd = if (length(all.samples) < 30) { 3 } else { 1 },
		col.lwd = if (length(all.samples) < 30) { 3 } else { 1 },
		col.lines = get.line.breaks(rev(plot.data$Group)),
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
		unlink(paste0(arguments$report, '/', 'mutation_signatures.png'));
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
