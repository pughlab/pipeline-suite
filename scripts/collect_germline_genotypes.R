### collect_germline_genotypes.R ###################################################################
# combine and compare germline variants (genotypes) across a dataset

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
library(plyr);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-d', '--directory', type = 'character', help = 'path to data directory');
parser$add_argument('-p', '--project', type = 'character', help = 'project name');

arguments <- parser$parse_args();

# what's the date?
date <- Sys.Date();

setwd(arguments$directory);

### MAIN ###########################################################################################
# find results files
vcf.files <- list.files(pattern = 'vcf.gz$');

# create a list to store the raw data
germline.data <- list();

# read in each file
for (file in vcf.files) {

	# extract sample id
	smp <- unlist(strsplit(file,'_'))[1];

	# do some quick filtering
	sys.command <- paste0("zcat ", file, " | grep -v '##' | cut -f3,6-9 --complement > ", smp, ".vcf");
	system(sys.command);

	# read in the file and update header
	tmp <- read.delim(paste0(smp, '.vcf'));
	colnames(tmp)[1:4] <- c('Chromosome','Position','Ref','Alt');

	# for each sample, extract the genotypes
	for (field in colnames(tmp)[5:ncol(tmp)]) {
		tmp[,field] <- sapply(
			tmp[,field],
			function(i) {
				unlist(strsplit(as.character(i),':'))[1];
				}
			);
		}

	# save it in the list
	germline.data[[smp]] <- tmp;

	# remove junk
	system(paste0("rm ", smp, ".vcf"));
	rm(smp,tmp,sys.command,field);

	}

# merge per-sample results
combined <- join_all(
	germline.data,
	by = c('Chromosome','Position','Ref','Alt'),
	type = 'full',
	match = 'first'
	);

# and write to file for future use (ie, plotting)
write.table(
	combined,
	file = generate.filename(arguments$project, 'germline_genotypes', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# do a quick summarization
plot.data <- combined[,5:ncol(combined)];

plot.data[which(plot.data == '0/0', arr.ind = TRUE)] <- 0;
plot.data[which(plot.data == '0/1', arr.ind = TRUE)] <- 1;
plot.data[which(plot.data == '1/1', arr.ind = TRUE)] <- 2;

for (i in 1:ncol(plot.data)) {
	plot.data[,i] <- as.numeric(plot.data[,i]);
	}

cor.data <- cor(plot.data, method = 'spearman', use = 'pairwise');

for (i in rownames(cor.data)) {

	cor.data[i,i] <- NA;
	similar.smps <- colnames(cor.data)[which(cor.data[i,] > 0.8 & cor.data[i,] < 1)];

	print(paste("Sample", i, "shows high similarity to sample(s):", paste(similar.smps, collapse = ', ')));
	}

# and write to file for future use (ie, plotting)
write.table(
	cor.data,
	file = generate.filename(arguments$project, 'germline_correlation', 'tsv'),
	row.names = TRUE,
	col.names = NA,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectVariantData','SessionProfile','txt'));
