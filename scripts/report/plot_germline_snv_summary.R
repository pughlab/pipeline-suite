### plot_snv_summary.R ############################################################################
# Identify and plot recurrent SNVs.
# INPUT:
#	- germline mutation calls (final list from annotate_germline.pl)

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
parser$add_argument('-m', '--maf', type = 'character', help = 'mutation calls in MAF format');
parser$add_argument('-z', '--report', type = 'character', help = 'path to report directory',
	default = NULL);

arguments <- parser$parse_args();

# import libraries
library(BoutrosLab.plotting.general);
library(xtable);

### READ DATA ######################################################################################
# if there happens to only be 'T' within any allele field (ie, Reference_Allele, Tumor_Seq_Allele),
# R will interpret this as 'TRUE' - set the field classes to character to avoid this

# edit: unfortunately, this sometimes causes errors when trying to read in the table...
# since we don't actually use allele information here, just ignore the problem for now

#maf.classes <- rep('character',132);
#maf.classes[c(2,6,7,40:45,58,60,77:85,100:107,112:122,124:132)] <- 'numeric';

# get data
if (is.null(arguments$maf)) {
	stop('ERROR: No input MAF provided, please provide path to SNV calls in MAF format.');
	} else {
	input.data <- read.delim(arguments$maf, stringsAsFactors = FALSE, comment.char = '#');
	}

# collect list of all samples
all.samples <- sort(as.character(unique(input.data$Tumor_Sample_Barcode)));

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
variant.colours <- c('darkseagreen4','darkorchid4','#9AA3F2','yellow','darkorange3','#F9B38E','turquoise1','plum','grey50')
names(variant.colours) <- c('missense','nonsense','nonstop','splicing','frameshift_indel','in_frame_indel','tss','RNA','other');

# for these plots, we will ignore some variant types
variant.colours <- variant.colours[c(1:6,9)];

# make the plot legend (mutation type/consequence)
functional.legend <- legend.grob(
	legends = list(
		legend = list(
			colours = variant.colours,
			labels = names(variant.colours)
			)
		),
	title.just = 'left',
	label.cex = 0.7,
	size = 1
	);

### manually update annotation ###
# vcf2maf mis-annotates a common downstream/regulartory_region TERC mutation (rs2293607) 
# as a downstream variant in ACTRT3; therefore we will manually update this for plotting
idx <- which(input.data$dbSNP_RS == 'rs2293607' & input.data$SYMBOL == 'ACTRT3');
if (length(idx) > 0) {
	input.data[idx,]$Hugo_Symbol <- 'TERC';
	input.data[idx,]$Entrez_Gene_Id <- '7012';
	}
##################################

# indicate key fields
keep.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2', 't_depth','t_ref_count','n_depth','n_ref_count');

mutation.data <- input.data[,keep.fields];

# calculate tumour vaf
mutation.data$t_vaf <- 1-(mutation.data$t_ref_count/mutation.data$t_depth);
mutation.data$n_vaf <- if (all(is.na(mutation.data$n_depth))) { rep(NA, nrow(mutation.data)); } else {
	1-(mutation.data$n_ref_count/mutation.data$n_depth);
	}

# apply variant coding
mutation.data$Code <- variant.codes$Code[match(mutation.data$Variant_Classification, variant.codes$Classification)];
if (any(mutation.data$Code > 6)) {
	mutation.data[which(mutation.data$Code > 6),]$Code <- 9;
	}

# add some basic filters to remove low coverage variants
mutation.data <- mutation.data[which(mutation.data$n_depth >= 15 | is.na(mutation.data$n_vaf)),];
mutation.data <- mutation.data[which(
	(mutation.data$t_depth >= 15 & is.na(mutation.data$n_vaf)) | 
	!is.na(mutation.data$n_vaf) ),];

# true germline variants should have VAF of 0.5 or 1, so add a VAF filter
mutation.data <- mutation.data[which(
	(is.na(mutation.data$n_vaf) & mutation.data$t_vaf > 0.4) | 
	(mutation.data$n_vaf > 0.4)),];

# if there are fewer than 5 germline variants here, write a table to output
show.germline <- mutation.data[,c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','t_vaf','n_vaf')];
colnames(show.germline) <- c('Sample','Symbol','Class','TumourVAF','NormalVAF');

if (nrow(mutation.data) == 0) {
	plot.data <- data.frame(matrix(nrow = 0, ncol = 1));
	} else {
	# reduce to 1 mutation per gene per sample [taking the higher priority code]
	mutation.data.trimmed <- aggregate(
		Code ~ Tumor_Sample_Barcode + Hugo_Symbol + Chromosome,
		mutation.data,
		min
		);

	mutation.data.trimmed <- merge(
		mutation.data.trimmed,
		mutation.data[,c('Tumor_Sample_Barcode','Hugo_Symbol','Code','t_vaf','n_vaf')],
		all.x = TRUE
		);

	# reshape data
	plot.data <- reshape(
		mutation.data.trimmed[,c('Tumor_Sample_Barcode','Chromosome','Hugo_Symbol','Code')],
		direction = 'wide',
		timevar = 'Tumor_Sample_Barcode',
		idvar = c('Chromosome','Hugo_Symbol')
		);
	colnames(plot.data) <- gsub('Code.','',colnames(plot.data));
	}

### MAKE PLOTS #####################################################################################
if (nrow(plot.data) > 1) {

	vaf.data <- reshape(
		mutation.data.trimmed[,c('Tumor_Sample_Barcode','Chromosome','Hugo_Symbol','t_vaf')],
		direction = 'wide',
		timevar = 'Tumor_Sample_Barcode',
		idvar = c('Chromosome','Hugo_Symbol')
		);
	colnames(vaf.data) <- gsub('t_vaf.','',colnames(vaf.data));

	# add in any missing samples (those with no variants)
	missing.samples <- setdiff(all.samples, colnames(plot.data));
	plot.data[,missing.samples] <- NA;
	vaf.data[,missing.samples] <- NA;

	# get per-gene sample counts
	plot.data <- plot.data[,c('Hugo_Symbol', 'Chromosome', all.samples)];
	plot.data$Count <- apply(plot.data[,all.samples],1,function(i) { length(i[!is.na(i)]) } );
	plot.data <- plot.data[order(-plot.data$Count),];

	# trim down to top recurrently mutated genes
	top.count <- min(nrow(plot.data),20);
	plot.data <- plot.data[seq(1,top.count,1),];

	# clean up for plotting
	plot.data$Label <- as.character(plot.data$Hugo_Symbol);
	if (any(duplicated(plot.data$Hugo_Symbol))) {
		dup.symbols <- unique(plot.data$Hugo_Symbol[duplicated(plot.data$Hugo_Symbol)]);
		for (gene in dup.symbols) {
			tmp <- plot.data[which(plot.data$Hugo_Symbol == gene),1:2];
			new.symbols <- paste0(plot.data$Hugo_Symbol, '_', plot.data$Chromosome);
			plot.data[which(plot.data$Hugo_Symbol == gene),]$Label <- new.symbols;
			}
		}

	# grab some parameters
	axis.cex <- if (length(all.samples) <= 30) { 1
		} else if (length(all.samples) <= 50) { 0.75
		} else if (length(all.samples) <= 80) { 0.5
		} else if (length(all.samples) <= 100) { 0.4
		} else { 0 };

	# create heatmap for recurrent genes (ordered by recurrence)
	heatmap.data <- t(plot.data[,all.samples]);
	heatmap.data[!is.na(heatmap.data)] <- 1;
	heatmap.data <- heatmap.data[do.call(order, transform(heatmap.data)),];

	vaf.data <- vaf.data[colnames(heatmap.data),rownames(heatmap.data)];
	vaf.data[is.na(vaf.data)] <- 0;

	# create function to determine spot size
	modifier <- if (length(all.samples) < 12) { 1.5 } else { 1 }
	spot.size.vaf <- function(x) { abs(x)* modifier; }
	spot.colour.vaf <- function(x) {
		sapply(x, function(i) { if (i >= 0.5) { 'black' } else { 'grey70' } } )
		}

	dot.key <- list(
		points = list(
			col = c('black','black','black','grey70','grey70','grey70'),
			cex = spot.size.vaf(c(1,0.8,0.5,0.2,0.1,0)),
			pch = 19
			),
		text = list(
			lab = c('1.0','0.8','0.5','0.2','0.1','0'),
			cex = 0.7
			),
		padding.text = 1.05,
		title = 'Tumour VAF',
		cex.title = 0.8
		);

	create.dotmap(
		x = vaf.data,
		bg.data = plot.data[rownames(vaf.data),colnames(vaf.data)],
		spot.size.function = spot.size.vaf,
		spot.colour.function = spot.colour.vaf,
		pch = 19,
		colour.scheme = variant.colours,
		at = seq(0,length(variant.colours),1),
		total.colours = length(variant.colours)+1,
		colourkey = FALSE,
		legend = list(
			right = list(fun = draw.key, args = list(key = dot.key)),
			inside = list(fun = functional.legend, x = 1.02, y = 1)
			),
		right.padding = 5,
		xaxis.lab = rownames(heatmap.data),
		yaxis.lab = plot.data$Label,
		xaxis.cex = axis.cex,
		yaxis.cex = 1,
		xaxis.tck = 0,
		yaxis.tck = 0,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'italic',
		xaxis.rot = 90,
		bg.alpha = 1,
		lwd = 1,
		col.lwd = 2,
		row.lwd = 2,
		col.colour = 'grey80',
		row.colour = 'grey80',
		height = if (nrow(plot.data) > 20) { 8 } else { 6 },
		width = if (length(all.samples) > 30) { 11 } else { 8 },
		resolution = 200,
		filename = generate.filename(arguments$project, 'germline_snvs','png')
		);
	}

save(
	show.germline,
	mutation.data,
	plot.data,
	variant.colours,
	file = generate.filename(arguments$project, 'germline_mutation_summary', 'RData')
	);

### LATEX ##########################################################################################
# if making output for a report
if (!is.null(arguments$report)) {

	tex.file <- paste0(
		arguments$report,
		'/germline_snv_summary.tex'
		);

	# write some captions
	recurrence.caption <- "Summary of pathogenic germline variants (SNVs and INDELs). Figure shows most frequently mutated genes; background colours indicate predicted functional consequence (white = no mutation detected) while circles indicate the variant was detected in the tumour (such that absence suggests the germline mutation may have reverted); black indicates VAF $>=$ 0.5. For T/N pairs, germline variants are detected in the normal sample. For tumour-only samples, variants are those detected in any tumour for a given patient (ie, for patients with multiple tumours, the variant may be in any them).";

	# write for latex
	write("\\section{Germline Mutation Summary}", file = tex.file);

	# first, check for mutation_overlap plot
	if (nrow(mutation.data) == 0) {
		write("No potentially pathogenic germline mutations were detected in the cohort.", file = tex.file, append = TRUE);
		} else if (nrow(plot.data) == 0) {
		write("No genes were recurrently mutated across the cohort.", file = tex.file, append = TRUE);
		} else if (nrow(input.data) < 5) {

		print(
			xtable(
				show.germline,
				caption = 'List of high-confidence germline mutations across the cohort.'
				),
			file = tex.file,
			include.rownames = FALSE,
			append = TRUE
			);
		
		} else {

		# create symlinks for plots
		unlink(paste0(arguments$report,'/germline_snvs.png'));
		file.symlink(
			paste0(arguments$output,'/',
				generate.filename(arguments$project, 'germline_snvs','png')),
			paste0(arguments$report,'/germline_snvs.png')
			);		

		write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
		write("\\begin{center}", file = tex.file, append = TRUE);
		write(paste0(
			"\\includegraphics[width=0.9\\textwidth]{",
			paste0(arguments$report, '/', 'germline_snvs.png'), '}'
			), file = tex.file, append = TRUE);
		write("\\end{center}", file = tex.file, append = TRUE);
		write(paste0(
			"\\caption{", recurrence.caption, "}"
			), file = tex.file, append = TRUE);
		write("\\end{figure}\n", file = tex.file, append = TRUE);
		}
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('GermlineSummary','SessionProfile','txt'));
