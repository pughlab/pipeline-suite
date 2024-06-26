### plot_rna_snv_summary.R #########################################################################
# Identify and plot recurrent SNVs from RNA-Seq data. This is a trimmed down version of 
# plot_snv_summary.R used for DNA because of the increased SNV count from RNA.
# INPUT:
#	- mutation calls (variant x variant matrix output by collect_snv_output.R)

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
parser$add_argument('-m', '--mutations', type = 'character',
	help = 'mutation calls (combined maf output by collect_snv_output.R)');

arguments <- parser$parse_args();

# import libraries
library(BoutrosLab.plotting.general);
library(xtable);

variant.colours <- c('darkseagreen4','darkorchid4','#9AA3F2','yellow','darkorange3','#F9B38E','turquoise1','plum','grey50')
names(variant.colours) <- c('missense','nonsense','nonstop','splicing','frameshift_indel','in_frame_indel','tss','RNA','other');

# for these plots, we will ignore some variant types
variant.colours <- variant.colours[c(1:6,9)];
basechange.colours <- default.colours(8,'pastel')[-5];

### READ DATA ######################################################################################
# get data
if (is.null(arguments$mutations)) {
	stop('ERROR: No input MAF provided, please provide path to SNV calls in MAF format.');
	} else {
	input.data <- read.delim(arguments$mutations, stringsAsFactors = FALSE, comment.char = '#');
	}

# collect list of all samples
all.samples <- unique(input.data$Tumor_Sample_Barcode);

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# get per-sample mutation counts
sample.counts <- data.frame(
	Sample = gsub('\\.','-',all.samples),
	Total = as.numeric(table(input.data$Tumor_Sample_Barcode)[all.samples]),
	SNVs = as.numeric(table(input.data[which(input.data$Variant_Type == 'SNP'),]$Tumor_Sample_Barcode)[all.samples]),
	INDELs = as.numeric(table(input.data[which(!input.data$Variant_Type == 'SNP'),]$Tumor_Sample_Barcode)[all.samples]) 
	);

sample.counts <- sample.counts[order(sample.counts$Total, decreasing = F),];

plot.data <- reshape(
	sample.counts,
	direction = 'long',
	v.names = 'Count',
	varying = list(3:4),
	timevar = 'Type',
	times = c('SNVs','INDELs')
	);

if (any(is.na(plot.data$Count))) {
	plot.data[is.na(plot.data$Count),]$Count <- 0;
	}

plot.data$Sample <- factor(plot.data$Sample, levels = sample.counts$Sample);

### PLOT DATA ######################################################################################
# grab some parameters
axis.cex <- if (length(all.samples) <= 30) { 1
	} else if (length(all.samples) <= 50) { 0.75
	} else if (length(all.samples) <= 80) { 0.5
	} else if (length(all.samples) <= 100) { 0
	} else { 0 };

# create plot for mutation rates
create.barplot(
	Sample ~ Count,
	plot.data,
	plot.horizontal = TRUE,
	stack = TRUE,
	groups = plot.data$Type,
	col = default.colours(2),
	yaxis.cex = axis.cex,
	xaxis.cex = 1,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.2,0),
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	axes.lwd = 1,
	xlab.label = 'Count',
	xlab.cex = 1.2,
	ylab.label = NULL,
	height = 10,
	width = 8,
	resolution = 200,
	filename = generate.filename(arguments$project, 'mutation_summary','png'),
	);

save(
	sample.counts,
	file = generate.filename(arguments$project, 'mutation_summary', 'RData')
	);

### LATEX ##########################################################################################
# initiate file
tex.file <- 'snv_summary.tex';
write("\\section{SNV Summary}", file = tex.file);

# write some captions
summary.caption <- "Summary of short somatic variants (SNVs, indels) for each sample: orange = INDELs; green = SNVs.";

write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
write("\\begin{center}", file = tex.file, append = TRUE);
write(paste0(
	"\\includegraphics[width=0.9\\textwidth]{",
	getwd(), '/',
	generate.filename(arguments$project, 'mutation_summary','png'), '}'
	), file = tex.file, append = TRUE);
write("\\end{center}", file = tex.file, append = TRUE);
write(paste0(
	"\\caption{", summary.caption, "}"
	), file = tex.file, append = TRUE);
write("\\end{figure}\n", file = tex.file, append = TRUE);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('MutationSummary','SessionProfile','txt'));
