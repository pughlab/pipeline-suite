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

### PREPARE SESSION ################################################################################
# import libraries
library(xtable);
library(BoutrosLab.plotting.general);
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('-m', '--maf', type = 'character', help = 'mutation calls in MAF format');

arguments <- parser$parse_args();

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

#variant.colours <- c('#673AB7','#2196F3','#F44336','#00BCD4','#E91E63','#8BC34A','#FFC107','#03A9F4','grey50');
variant.colours <- c('darkseagreen4','darkorchid4','#9AA3F2','yellow','darkorange3','#F9B38E','turquoise1','plum','grey50')
names(variant.colours) <- c('missense','nonsense','nonstop','splicing','frameshift_indel','in_frame_indel','tss','RNA','other');

# for these plots, we will ignore some variant types
variant.colours <- variant.colours[c(1:6,9)];

### READ DATA ######################################################################################
# get data
input.data <- read.delim(arguments$maf);

# collect list of all samples
all.samples <- sort(as.character(unique(input.data$Tumor_Sample_Barcode)));

# move to output directory
setwd(arguments$output);

### FORMAT DATA ####################################################################################
# indicate key fields
keep.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2', 't_depth','t_ref_count');
mutation.data <- input.data[,keep.fields];

# calculate tumour vaf
mutation.data$t_vaf <- 1-(mutation.data$t_ref_count/mutation.data$t_depth);

# apply variant coding
mutation.data$Code <- variant.codes$Code[match(mutation.data$Variant_Classification, variant.codes$Classification)];
mutation.data[which(mutation.data$Code > 6),]$Code <- 9;

# get per-sample mutation counts
sample.counts <- data.frame(cbind(
	table(mutation.data$Tumor_Sample_Barcode),
	table(mutation.data[which(mutation.data$Code != 9),]$Tumor_Sample_Barcode)
	));
colnames(sample.counts) <- c('Total','Non.Silent');

sample.counts <- sample.counts[order(sample.counts$Total, decreasing = TRUE),];
sample.counts$Order <- 1:nrow(sample.counts);

# reduce to 1 mutation per gene per sample [taking the higher priority code]
mutation.data <- aggregate(
	Code ~ Tumor_Sample_Barcode + Hugo_Symbol + Chromosome + t_vaf,
	mutation.data,
	min
	);

# reshape data
plot.data <- reshape(
	mutation.data[,c('Tumor_Sample_Barcode','Chromosome','Hugo_Symbol','Code')],
	direction = 'wide',
	timevar = 'Tumor_Sample_Barcode',
	idvar = c('Chromosome','Hugo_Symbol')
	);
colnames(plot.data) <- gsub('Code.','',colnames(plot.data));

vaf.data <- reshape(
	mutation.data[,c('Tumor_Sample_Barcode','Chromosome','Hugo_Symbol','t_vaf')],
	direction = 'wide',
	timevar = 'Tumor_Sample_Barcode',
	idvar = c('Chromosome','Hugo_Symbol')
	);
colnames(vaf.data) <- gsub('t_vaf.','',colnames(vaf.data));

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

vaf.data <- vaf.data[rev(colnames(heatmap.data)),rownames(heatmap.data)];
vaf.data[is.na(vaf.data)] <- 0;

create.heatmap(
	plot.data[,rownames(heatmap.data)],
	cluster.dimensions = 'none',
	same.as.matrix = TRUE,
	xaxis.lab = rownames(heatmap.data),
	yaxis.lab = plot.data$Label,
	xaxis.cex = axis.cex,
	yaxis.cex = 1,
	xaxis.tck = if (axis.cex == 0) { 0 } else { 0.2 },
	yaxis.tck = 0.2,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	axes.lwd = 1,
	grid.row = TRUE,
	force.grid.row = TRUE,
	row.colour = 'grey80',
	col.colour = 'grey80',
	row.lwd = if (length(all.samples) < 30) { 3 } else { 1 },
	col.lwd = if (length(all.samples) < 30) { 3 } else { 1 },
	grid.col = TRUE,
	force.grid.col = TRUE,
	print.colour.key = FALSE,
	fill.colour = 'white',
	at = seq(0,length(variant.colours),1),
	total.colours = length(variant.colours)+1,
	colour.scheme = variant.colours,
	row.pos = which(vaf.data > 0, arr.ind = TRUE)[,1],
	col.pos = which(vaf.data > 0, arr.ind = TRUE)[,2],
	cell.text = rep(expression("\u25CF"), times = sum(vaf.data > 0)),
	text.col = 'black',
	text.cex = 0.4,
	inside.legend = list(fun = functional.legend, x = 1.02, y = 1),
	right.padding = 15,
	height = 6,
	width = 8,
	resolution = 200,
	filename = generate.filename(arguments$project, 'germline_snvs','png')
	);

save(
	sample.counts,
	mutation.data,
	plot.data,
	variant.colours,
	file = generate.filename(arguments$project, 'germline_mutation_summary', 'RData')
	);

# write some captions
recurrence.caption <- "Summary of short germline variants (SNVs and INDELs). Figure shows most frequently mutated genes; background colours indicate predicted functional consequence (white = no mutation detected) while black circles indicate the variant was also detected in the tumour (such that absence suggests the germline mutation may have been lost).";

# write for latex
write("\\section{Germline SNV Summary}", file = 'germline_snv_summary.tex');

# first, check for mutation_overlap plot
if (nrow(plot.data) == 0) {
	write("\\pagebreak\nNo genes were recurrently mutated across the cohort.", file = 'germline_snv_summary.tex', append = TRUE);
	} else {
	write("\\pagebreak\n\\begin{figure}[h!]", file = 'germline_snv_summary.tex', append = TRUE);
	write("\\begin{center}", file = 'germline_snv_summary.tex', append = TRUE);
	write(paste0(
		"\\includegraphics[width=0.9\\textwidth]{",
		getwd(), '/',
		generate.filename(arguments$project, 'germline_snvs','png'), '}'
		), file = 'germline_snv_summary.tex', append = TRUE);
	write("\\end{center}", file = 'germline_snv_summary.tex', append = TRUE);
	write(paste0(
		"\\caption{", recurrence.caption, "}"
		), file = 'germline_snv_summary.tex', append = TRUE);
	write("\\end{figure}\n", file = 'germline_snv_summary.tex', append = TRUE);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('GermlineSummary','SessionProfile','txt'));
