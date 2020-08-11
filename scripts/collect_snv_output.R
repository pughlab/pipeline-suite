### collect_snv_output.R ###########################################################################
# find, collate and format output from SNV caller

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

# function to identify symlinks
is.symlink <- function(file) {
	isTRUE(nzchar(Sys.readlink(file), keepNA = TRUE));
	}

### PREPARE SESSION ################################################################################
# import libraries
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-d', '--directory', type = 'character', help = 'path to data directory');
parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-g', '--gtf', type = 'character', help = 'annotation gtf (refGene)',
	default = '/cluster/projects/pughlab/references/gencode/GRCh38/gencode.v31.GRCh38.genes.gtf');

arguments <- parser$parse_args();

# what's the date?
date <- Sys.Date();

setwd(arguments$directory);


### VARIANT CODING
# 1 = missense, 2 = stop gain/loss, 3 = splicing, 4 = frameshift, 5 = in frame indel, 6 = tss
# 7 = RNA, 8 = other (up/downstream, UTR, intergenic, silent, intron)
variant.codes <- data.frame(
	Classification = c("3'Flank", "5'Flank", "Intron", "RNA", "IGR", "3'UTR", "5'UTR", "Silent", "Missense_Mutation",
		"Splice_Region", "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins",
		"Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "ITD"),
	Group = c('other','other','other','RNA','other','other','other','other','missense',
		'splice_site','splice_site','in_frame_indel','in_frame_indel','frameshift_indel','frameshirt_indel',
		'nonsense', 'nonstop', 'tss', 'itd'),
	Code = c(9, 9, 9, 8, 9, 9, 9, 9, 1,
		4, 4, 6, 6, 5, 5,
		2, 3, 7, 10)
	);

variant.colours <- c("#F44336", "#E91E63", "grey70", "#673AB7", "#2196F3", "#03A9F4", 
	"#00BCD4", "#8BC34A", "#CDDC39", "#FFC107"
	);
names(variant.colours) <- c("nonstop","frameshift_indel","other", "missense", "nonsense", "RNA",
	"splicing", "in_frame_indel", "itd", "tss"
	);

# if this is RNA-Seq, only keep variants in coding regions
is.rnaseq <- grepl('RNASeq', getwd());
is.exome <- grepl('Exome', getwd());
classes.to.keep <- c('RNA','Missense_Mutation','Splice_Region','Splice_Site','In_Frame_Del','In_Frame_Ins',
	'Frame_Shift_Del','Frame_Shift_Ins','Nonsense_Mutation','Nonstop_Mutation','Translation_Start_Site');

###

### FORMAT ANNOTATION
# using refGene, indicate genes to keep (coding genes) for gene x patient matrix
refGene <- read.delim(arguments$gtf, header = F, comment.char = '#');

refGene <- droplevels(refGene[which(refGene$V3 == 'gene'),c(1,4,5,9)]);
colnames(refGene) <- c('Chromosome','Start','End','INFO');

refGene$GeneID <- sapply(
	refGene$INFO,
	function(i) {
		parts <- unlist(strsplit(as.character(i), ';'));
		gene_name <- unlist(strsplit(parts[grepl('gene_id', parts)], ' '));
		gene_id  <- gene_name[length(gene_name)];
		return(gene_id);
		}
	);

refGene$Symbol <- sapply(
	refGene$INFO,
	function(i) {
		parts <- unlist(strsplit(as.character(i), ';'));
		gene_name <- unlist(strsplit(parts[grepl('gene_name', parts)], ' '));
		gene_symbol <- gene_name[length(gene_name)];
		return(gene_symbol);
		}
	);

refGene$Chromosome <- factor(refGene$Chromosome, levels = paste0('chr',c(1:22,'X','Y','M')));

### MAIN ###########################################################################################
# find results files
maf.files <- list.files(pattern = '*.maf$', recursive = TRUE);
maf.files <- maf.files[!sapply(maf.files, is.symlink)]; 

# read them in and store them
samples <- c();

variant.data <- data.frame();
maf.fields <- c(
	'Chromosome',
	'Start_Position',
	'End_Position',
	'Reference_Allele',
	'Allele',
	'dbSNP_RS',
	'Hugo_Symbol'
	);

for (i in 1:length(maf.files)) {

	file <- maf.files[i];

	# extract sample ID
	if (is.rnaseq) { smp <- unlist(strsplit(file,'\\/'))[3]; } else { smp <- unlist(strsplit(file,'\\/'))[2]; }
	samples <- unique(c(samples, smp));

	# read in data
	tmp <- read.delim(file, comment.char = "#");

	## do some filtering
	# if RNA-Seq or Exome-Seq, focus on genic regions only
	if (is.rnaseq || is.exome) {
		tmp <- tmp[which(tmp$Variant_Classification %in% classes.to.keep),];
		}
	tmp <- tmp[which(tmp$FILTER == 'PASS'),];
	tmp <- tmp[which(tmp$t_depth > 20),];

	if (!all(is.na(tmp$n_depth))) {
		tmp <- tmp[which(tmp$n_depth > 15),];
		}

	# filter out known (dbSNP) variants
	#tmp <- tmp[!grepl('rs', tmp$dbSNP_RS),];

	if (nrow(tmp) == 0) {
		next;
		}

	# store it
	if (nrow(variant.data) == 0) {

		variant.data <- tmp[,maf.fields];
		colnames(variant.data) <- c('Chromosome','Start','End','Ref','Alt','dbSNP','Symbol');

		variant.data[,smp] <- variant.codes$Code[match(tmp$Variant_Classification, variant.codes$Classification)];

		} else {

		tmp2 <- tmp[,maf.fields];
		colnames(tmp2) <- c('Chromosome','Start','End','Ref','Alt','dbSNP','Symbol');
		tmp2[,smp] <- variant.codes$Code[match(tmp$Variant_Classification, variant.codes$Classification)];

		variant.data <- merge(variant.data, tmp2, all = TRUE);

		rm(tmp2);
		}

	rm(tmp);
	print("");
	print(file);
	print(gc());
	print("");
	}

# add in empty samples
for (smp in samples) {
	if (smp %in% colnames(variant.data)) { next; }
	variant.data[,smp] <- NA;
	}

# save data to file
write.table(
	variant.data,
	file = generate.filename(arguments$project, 'variant_by_patient', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# collapse to per-gene
gene.by.patient <- aggregate(
	variant.data[,samples],
	by = list(
		Symbol = variant.data$Symbol,
		Chromosome = variant.data$Chromosome
		),
	FUN = min,
	na.rm = TRUE
	);

# filter data
gene.by.patient <- merge(
	refGene[,c('Symbol','Chromosome','Start','End')],
	gene.by.patient,
	all.x = TRUE
	);

for (smp in samples) {
	smp.codes <- gene.by.patient[,smp];
	smp.codes[is.na(smp.codes)] <- 0;
	smp.codes[which(smp.codes == 'Inf')] <- 0;
	gene.by.patient[,smp] <- smp.codes;
	}

gene.by.patient <- gene.by.patient[order(gene.by.patient$Chromosome, gene.by.patient$Start),];

write.table(
	gene.by.patient,
	file = generate.filename(arguments$project, 'gene_by_patient', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectVariantData','SessionProfile','txt'));
