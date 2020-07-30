### collect_rnaseqc_output.R #######################################################################
# find, collate and format output from RNASeqC

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

parser$add_argument('-d', '--directory', type = 'character', help = 'path to data directory');
parser$add_argument('-p', '--project', type = 'character', help = 'project name');

arguments <- parser$parse_args();

# what's the date?
date <- Sys.Date();

setwd(arguments$directory);


### RNASeQC fields
# where End.1...Sense and End.2...Sense = End 1 % Sense and End 2 % Sense
# where Num..Gaps and Gap.. = Number of Gaps and Gap % (cumulative gap length/cumulative transcript length)
rnaseqc.fields <- list(
	overall.mapping.qc = c('Total.Purity.Filtered.Reads.Sequenced','Alternative.Aligments','Failed.Vendor.QC.Check','Read.Length','Estimated.Library.Size'),
	aligned.mapping.qc = c('Mapped','Mapping.Rate','Mapped.Unique','Mapped.Unique.Rate.of.Total','Unique.Rate.of.Mapped','Duplication.Rate.of.Mapped','Base.Mismatch.Rate','rRNA','rRNA.rate'),
	mate.mapping.qc = c('Mapped.Pairs','Unpaired.Reads','End.1.Mapping.Rate','End.2.Mapping.Rate','End.1.Mismatch.Rate','End.2.Mismatch.Rate','Fragment.Length.Mean','Fragment.Length.StdDev','Chimeric.Pairs'),
	transcript.mapping.qc = c('Intragenic.Rate','Intronic.Rate','Exonic.Rate','Intergenic.Rate','Split.Reads','Expression.Profiling.Efficiency','Genes.Detected','Transcripts.Detected'),
	strand.specific.qc = c('End.1.Sense','End.1.Antisense','End.2.Sense','End.2.Antisense','End.1...Sense','End.2...Sense'),
	coverage.qc = c('Mean.CV','Mean.Per.Base.Cov.','Num..Gaps','Gap..')
	);

### MAIN ###########################################################################################
# find results files 
qc.files <- list.files(pattern = 'metrics.tsv', recursive = TRUE);
gct.files <- list.files(pattern = 'genes.rpkm.gct', recursive = TRUE);

# read them in and store them

### IF all samples were run together
if (length(qc.files) == 1) {

	combined.qc <- read.delim(qc.files[1]);

	smp.columns <- gsub('-', '.', unique(combined.qc$Sample));

	# get correlations (these are Pearson's Correlation)
	corr.data <- combined.qc[,c('Sample', smp.columns)];
	corr.data$Sample <- gsub('-', '.', corr.data$Sample);

	}

if (length(qc.files) > 1) {

	# read them in and store them
	qc.data <- list();
	gct.data <- list();

	# qc metrics
	for (i in 1:length(qc.files)) {

		file <- qc.files[i];

		tmp <- read.delim(file);
		tmp <- tmp[,c('Sample',unlist(rnaseqc.fields))];

		# read in data
		qc.data[[i]] <- tmp;
		}

	combined.qc <- do.call(rbind, qc.data);

	# genes.rpkm.gct data for correlation
	for (i in 1:length(gct.files)) {

		file <- gct.files[i];

		tmp <- read.delim(file, skip = 2, row.names = 1);

		# read in data
		gct.data[[i]] <- t(tmp[,-1]);
		}

	combined.gct <- t(do.call(rbind, gct.data));
	corr.data <- cor(combined.gct, method = 'pearson');
	}

# format data
formatted.qc <- combined.qc[,c('Sample',unlist(rnaseqc.fields))];
colnames(formatted.qc)[which(colnames(formatted.qc) == 'End.1...Sense')] <- 'Percent.Intragenic.End.1.Sense';
colnames(formatted.qc)[which(colnames(formatted.qc) == 'End.2...Sense')] <- 'Percent.Intragenic.End.2.Sense';
colnames(formatted.qc)[which(colnames(formatted.qc) == 'Mean.Per.Base.Cov.')] <- 'Mean.Per.Base.Cov';
colnames(formatted.qc)[which(colnames(formatted.qc) == 'Num..Gaps')] <- 'Num.Gaps';
colnames(formatted.qc)[which(colnames(formatted.qc) == 'Gap..')] <- 'Gap.Percentage';

# save data to file
write.table(
	corr.data,
	file = generate.filename(arguments$project, 'rnaseqc_Pearson_correlations', 'tsv'),
	row.names = TRUE,
	col.names = NA,
	sep = '\t'
	);

write.table(
	formatted.qc,
	file = generate.filename(arguments$project, 'rnaseqc_output', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectRNASeQC','SessionProfile','txt'));
