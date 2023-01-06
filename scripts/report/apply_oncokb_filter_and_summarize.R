### filter_ensemble_mutations.R ####################################################################
# FILTERs ensemble mutations (using OncoKB annotations as a white list) and plots results
# replaces output from plot_snv_summary.R for use in report

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

	# write iemory usage to file
	cat('### MEMORY USAGE ###############################################################');
	print(proc.time());

	# write key variables to file
	cat("\n### VARIABLES #################################################################\n");
	cat(paste0('Input: ', arguments$input));
	cat(paste0('Treat low coverage as tumour-only: ', arguments$set_tumour_only));
	cat(paste0('Target coverage: ', arguments$coverage));

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

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('-i', '--input', type = 'character', help = 'mutation calls in MAF format');
parser$add_argument('-c', '--coverage', type = 'character', default = '20,15',
	help = 'minimum depth for tumour and normal to be considered callable; length 1 or 2 (default: 20,15)');
parser$add_argument('-t', '--set_tumour_only', type = 'logical', default = FALSE,
	help = 'should somatic mutations with low coverage in matched normal be processed as tumour-only?');
parser$add_argument('-z', '--report', type = 'character', help = 'path to report directory',
	default = NULL);

arguments <- parser$parse_args();

# import libraries
library(BoutrosLab.plotting.general);
library(xtable);

### READ DATA #####################################################################################
# get data
if (is.null(arguments$input)) {
	stop('ERROR: No input file provided, please provide path to SNV/INDELs calls in MAF format.');
	} else {
	input.data <- read.delim(arguments$input, stringsAsFactors = FALSE, comment.char = '#');
	}

# stop if no OncoKB annotations available
if (! 'GENE_IN_ONCOKB' %in% colnames(input.data)) {
	stop('No OncoKB annotations available in input maf; please run OncoKB first.');
	}

# collect list of all samples
all.samples <- sort(as.character(unique(input.data$Tumor_Sample_Barcode)));

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

### FILTER DATA ####################################################################################
# oncokb seems to select genes *purely* on gene symbol (including alternate names)
# for example, TEP1 is selected because of PTEN (aka TEP1), even though they are
# on different chromosomes
oncokb.columns <- c(grep('GENE_IN_ONCOKB', colnames(input.data)):ncol(input.data));
known.bugs <- which(input.data$Hugo_Symbol %in% c('PRMT9','RAD1','TEP1','FAH'));
if (length(known.bugs) > 0) {
	input.data[known.bugs,oncokb.columns] <- NA;
	}

# should somatic mutations with low coverage in matched normal be processed as tumour-only?
if (arguments$set_tumour_only) {

	target.cov.normal <- if (grepl(',', arguments$coverage)) {
		as.numeric(unlist(strsplit(arguments$coverage,','))[2]);
		} else if (!is.null(arguments$coverage)) {
		as.numeric(arguments$coverage);
		} else { 20; }

	low.n <- which(input.data$n_depth < target.cov.normal);
	input.data[low.n,]$FLAG.tumour_only <- TRUE;
	}

# find probable false positives to remove
remove.these <- which(input.data$FLAG.tumour_only & input.data$FLAG.high_pop);
also.remove.these <- which(
	input.data$FLAG.tumour_only & 
	input.data$FLAG.low_vaf
	);

# but whitelist these
keep.these <- grep('Oncogenic', input.data$ONCOGENIC);

# apply filter
if (length(setdiff(unique(c(remove.these,also.remove.these)), keep.these)) > 0) {
	remove.idx <- setdiff(unique(c(remove.these,also.remove.these)), keep.these);
	print(paste('Removing', length(remove.idx), 'variants as they are low confidence calls and not deemed oncogenic by oncoKB.'));

	input.data <- input.data[-remove.idx,];
	}

write.table(
	input.data,
	file = sub('.tsv','_filtered.tsv', basename(arguments$input)),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# extract key hits
significant.hits <- which(
	input.data$ONCOGENIC == 'Oncogenic' |
	input.data$MUTATION_EFFECT %in% c('Loss-of-function','Likely Loss-of-function') |
	input.data$MUTATION_EFFECT %in% c('Switch-of-function') |
	input.data$MUTATION_EFFECT %in% c('Gain-of-function','Likely Gain-of-function')
	);

if (length(significant.hits) > 0) {

	# format for latex
	oncokb.hits <- input.data[significant.hits,];

	oncokb.to.print <- data.frame(
		ID = oncokb.hits$Tumor_Sample_Barcode,
		Gene = oncokb.hits$Hugo_Symbol,
#		Position = oncokb.hits$Start_Position,
		HGVSp = oncokb.hits$HGVSp_Short,
		VAF = round(oncokb.hits$t_alt_count / oncokb.hits$t_depth,4),
		VARIANT_IN_ONCOKB = as.logical(oncokb.hits$VARIANT_IN_ONCOKB),
		Effect = oncokb.hits$MUTATION_EFFECT
		);

	# order by count and median VAF
	gene.counts <- aggregate(VAF ~ Gene, oncokb.to.print, length);
	colnames(gene.counts)[2] <- 'Count';
	gene.counts$VAF <- aggregate(VAF ~ Gene, oncokb.to.print, median)$VAF;
	gene.counts <- gene.counts[order(-gene.counts$Count, -gene.counts$VAF),];
	gene.counts <- gene.counts[1:min(40,nrow(gene.counts)),];

	oncokb.to.print$Gene <- factor(oncokb.to.print$Gene, levels = gene.counts$Gene);
	oncokb.to.print <- oncokb.to.print[!is.na(oncokb.to.print$Gene),];
	oncokb.to.print <- oncokb.to.print[order(
		oncokb.to.print$Gene,
		oncokb.to.print$VAF
		),];

	# save results
	write.table(
		oncokb.to.print,
		file = generate.filename(arguments$project, 'OncoKB_hits', 'tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);

	### MAKE PLOTS #############################################################################
	# plot oncoKB hits
	effect.colours <- c(default.colours(8,'div')[c(2,4,6,8)],'purple4');
	names(effect.colours) <- c('Loss-of-function','Likely Loss-of-function','Likely Gain-of-function','Gain-of-function','Switch-of-function');

	create.boxplot(
		VAF ~ Gene,
		oncokb.to.print,
		add.stripplot = TRUE,
		points.col = effect.colours[match(oncokb.to.print$Effect, names(effect.colours))],
		points.alpha = 0.8,
		points.cex = 0.8,
		xaxis.rot = 90,
		yaxis.tck = c(0.5,0),
		xaxis.tck = c(0.5,0),
		ylab.label = 'Variant Allele Frequency',
		ylab.cex = 1.2,
		ylab.axis.padding = 1,
		xlab.label = NULL,
		ylimits = c(0,1),
		yat = seq(0,1,0.2),
		yaxis.cex = 1,
		xaxis.cex = if (nrow(gene.counts) < 30) { 1 } else { 0.75 },
		style = 'Nature',
		top.padding = 2,
		right.padding =if (nrow(gene.counts) >= 30) { 3 } else { 15 },
		key = list(
			points = list(col = effect.colours, pch = 19, cex = 1),
			text = list(lab = names(effect.colours), cex = 0.8, col = 'black'),
			title = 'Mutation Effect',
			cex.title = 1,
			between = 0.3,
			x = if (nrow(gene.counts) >= 30) { 0.82 } else { 0.95 },
			y = if (nrow(gene.counts) >= 30) { 1.1 } else { 1 }
			),
		filename = generate.filename(arguments$project, 'oncoKB_hits', 'png'),
		height = 4,
		width = 8
		);
	}

### LATEX ##########################################################################################
# if making output for a report
if (!is.null(arguments$report)) {

	tex.file <- paste0(
		arguments$report,
		'/oncokb_summary.tex'
		);

	# write for latex
	write("\\section{OncoKB Results}", file = tex.file);

	# summarize key OncoKB findings
	oncokb.caption <- paste0(
		'OncoKB was used to annotate known/suspected oncogenic variants; this identified ',
		length(unique(as.character(oncokb.hits$Hugo_Symbol))),
		' genes (', nrow(oncokb.hits), ' total mutations) across the cohort.'
		);

	if (length(significant.hits) ==  0) {

		write(
			'No significant mutations found by OncoKB.',
			file = tex.file, append = TRUE
			);

		} else {
		# create symlinks for plots
		unlink(paste0(arguments$report,'/oncoKB_hits.png'));
		file.symlink(
			paste0(arguments$output,'/',
				generate.filename(arguments$project, 'oncoKB_hits','png')),
			paste0(arguments$report,'/oncoKB_hits.png')
			);

		if (nrow(gene.counts) > 5) {
			write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
			write("\\begin{center}", file = tex.file, append = TRUE);
			write(paste0(
				"\\includegraphics[width=0.85\\textwidth]{",
				paste0(arguments$report, '/', 'oncoKB_hits.png'), '}'
				), file = tex.file, append = TRUE);
			write("\\end{center}", file = tex.file, append = TRUE);
			write(paste0(
				"\\caption{", paste0(oncokb.caption, " Top hits are shown, sorted by recurrence and median VAF."), "}"
				), file = tex.file, append = TRUE);
			write("\\end{figure}\n", file = tex.file, append = TRUE);
			write("\\vspace{1.0cm}\n", file = tex.file, append = TRUE);
			}

		oncokb.to.print <- oncokb.to.print[order(oncokb.to.print$Gene, -oncokb.to.print$VAF),];
		idx <- if (nrow(oncokb.to.print) <= 15) { idx <- 1:nrow(oncokb.to.print); } else {
			idx <- which(oncokb.to.print$VARIANT_IN_ONCOKB == TRUE);
			idx <- idx[1:min(length(idx),15)];
			}

		print(
			xtable(
				oncokb.to.print[idx,],
				caption = oncokb.caption
				),
			file = tex.file,
			include.rownames = FALSE,
			size = 'scriptsize',
			append = TRUE
			);
		}
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SummarizeOncoKB','SessionProfile','txt'));
