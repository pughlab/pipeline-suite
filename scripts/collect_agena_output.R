### collect_agena_output.R #########################################################################
# Calculate VAF for Agena SNPs and get sample correlations to estimate sample similarity 
# (cohort contamination).

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
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-i', '--input', type = 'character', help = 'path to multi-sample VCF of Agena SNPs');
parser$add_argument('-a', '--agena', type = 'character', help = 'path to Agena SNP information (/path/to/pipeline-suite/data/agena_snps_annotated.txt)');

arguments <- parser$parse_args();

# what's the date?
date <- Sys.Date();

data.directory <- dirname(arguments$input);
setwd(data.directory);

### MAIN ###########################################################################################
agena <- read.delim(arguments$agena);

header <- read.delim(arguments$input, nrow = 5000, header = FALSE);
header <- header[grepl('##', header[,1]),];
smp.data <- read.delim(arguments$input, skip = length(header), stringsAsFactors = FALSE);
if (! colnames(smp.data)[1] == 'X.CHROM') {
	stop('Could not read input VCF');
	} else {
	colnames(smp.data)[1] <- 'CHROM';
	}

standard.fields <- c('CHROM','POS','REF','ALT','ID','FORMAT','QUAL','INFO','FILTER');
all.samples <- setdiff(colnames(smp.data),standard.fields);

# format output table
agena[,all.samples] <- NA;

# for each SNP
for (snp.idx in 1:nrow(agena)) {
	chr <- agena[snp.idx,]$Chromosome;
	pos <- agena[snp.idx,]$Position;
	vcf.idx <- which(smp.data$CHROM == chr & smp.data$POS == pos);

	if (length(vcf.idx) == 0) { next; }

	# identify index for genotype and allele depth (ref,alt)
	tmp <- unlist(strsplit(smp.data[vcf.idx,]$FORMAT,':'));
	gt.idx <- which(tmp == 'GT');
	ad.idx <- which(tmp == 'AD');

	# for each sample
	for (smp in all.samples) {
		genotype <- unlist(strsplit(smp.data[vcf.idx,smp],':'))[gt.idx];
		if (! genotype %in% c('0/0','0/1','1/1')) { next; }
		tmp <- unlist(strsplit(smp.data[vcf.idx,smp],':'))[ad.idx];
		ref <- as.numeric(unlist(strsplit(tmp,','))[1]);
		alt <- as.numeric(unlist(strsplit(tmp,','))[2]);
		agena[snp.idx,smp] <- alt/(ref + alt);
		}
	}

# get sample correlations
cor.data <- cor(agena[,all.samples], use = 'pairwise');

# save combined/formatted data to file
write.table(
	cor.data,
	file = generate.filename(arguments$project, '_agena_snps__correlations','tsv'),
	row.names = TRUE,
	col.names = NA,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectAgena','SessionProfile','txt'));
