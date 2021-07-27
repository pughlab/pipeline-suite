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

	# write sessionInfo to file
	cat("\n### SESSION INFO ###############################################################");
	print(sessionInfo());

	# close the file
	sink();

	}

# function to trim sample IDs
simplify.ids <- function(x) {
	match <- TRUE;
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
# import libraries
library(BoutrosLab.plotting.general);
library(deconstructSigs);
library(BSgenome);
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('-t', '--seq_type', type = 'character', help = 'exome or wgs', default = 'exome');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'hg38 or hg19', default = 'hg38');
parser$add_argument('-i', '--input', type = 'character', help = 'mutation calls in MAF format');
parser$add_argument('-s', '--signatures', type = 'character', help = 'path to signatures matrix (96 trinucleotide x N sigs)', default = NULL);

arguments <- parser$parse_args();

# get genome object
BSgenome.Hsapiens.UCSC <- if (arguments$ref_type == 'hg38') {
	getBSgenome("BSgenome.Hsapiens.UCSC.hg38");
	} else if (arguments$ref_type == 'hg19') {
	getBSgenome("BSgenome.Hsapiens.UCSC.hg19");
	}

# indicate normalization method (assuming signatures matrix was developed from WGS)
normalization.method <- 'default';
if (arguments$seq_type == 'exome') {
	normalization.method <- 'exome2genome';
	}

### READ DATA ######################################################################################
# get data
input.maf <- read.delim(arguments$input, stringsAsFactors = FALSE);

# read in mutation signatures
signatures.to.apply <- if (is.null(arguments$signatures)) { signatures.cosmic } else {
	as.data.frame(t(read.delim(arguments$signatures, row.names = 1)));
	}

### FORMAT DATA ####################################################################################
# remove INDELs
input.maf <- input.maf[which(input.maf$Variant_Type == 'SNP'),];

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
	assigned.sigs[[smp]] <- whichSignatures(
		tumor.ref = sigs.input,
		signatures.ref = signatures.to.apply,
		sample.id = smp,
		contexts.needed = TRUE,
		tri.counts.method = normalization.method
		);
	}

save(
	signatures.to.apply,
	assigned.sigs,
	file = generate.filename(arguments$project, 'assignedSignatures', 'RData')
	);

# extract weights
sig.weights.list <- list();
for (smp in all.samples) {
	sig.weights.list[[smp]] <- assigned.sigs[[smp]]$weights;
	}

sig.weights <- do.call(rbind, sig.weights.list);

write.table(
	sig.weights,
	file = generate.filename(arguments$project, 'mutation_signature_weights', 'tsv'),
	row.names = TRUE,
	col.names = NA,
	sep = '\t'
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

### PLOTTING #######################################################################################
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
	height = if (length(all.samples) < 10) { 6 } else { 8 },
	width = 12,
	resolution = 600,
	filename = generate.filename(arguments$project, 'mutation_signatures', 'png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('ApplyMutSigs','SessionProfile','txt'));


