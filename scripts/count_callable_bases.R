### count_callable_bases.R #########################################################################
# Finds and combines output from callable bases (get_coverage).

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

### MAIN ###########################################################################################
# find results files
cov.files <- list.files(pattern = '^CallableBases.tsv$', recursive = TRUE);

# read them in
cov.list <- list();
sample.list <- c();

for (file in cov.files) {
	# extract sample ID
	smp <- unlist(strsplit(file, '\\/'))[1];
	# store data in list
	cov.list[[smp]] <- read.delim(file);

	these.smps <- setdiff(colnames(cov.list[[smp]])[6:ncol(cov.list[[smp]])],'TargetRegions');
	sample.list <- c(sample.list, these.smps);
	}

# determine number of callable bases
callable.bases <- data.frame(
	Sample = sample.list,
	Total = NA,
	Patient.total = NA
	);

if ('TargetRegions' %in% colnames(cov.list[[1]])) {
	callable.bases$TargetRegions <- NA;
	callable.bases$Total.targeted <- NA;
	callable.bases$Patient.targeted <- NA;
	}

for (i in 1:length(cov.list)) {

	# count interval length
	cov.list[[i]]$Length <- cov.list[[i]]$end - cov.list[[i]]$start;

	tmp <- cov.list[[i]];

	# find samples
	these.smps <- intersect(colnames(tmp), sample.list);
	smp.idx <- which(callable.bases$Sample %in% these.smps);

	# find total bases covered for each sample and patient
	if (length(these.smps) == 1) {
		total.bases <- sum(tmp[,these.smps]*tmp$Length);
		total.patient <- total.bases;
		} else {
		total.bases <- apply(
			tmp[,these.smps],
			2,
			function(i) { sum(i*tmp$Length) }
			);

		total.patient <- sum(
			apply(tmp[,these.smps], 1, prod) * tmp$Length
			);
		}

	callable.bases[smp.idx,]$Total <- total.bases;
	callable.bases[smp.idx,]$Patient.total <- total.patient;

	# find total targeted bases covered for each sample and patient (if available)
	if ('TargetRegions' %in% colnames(tmp)) {

		if (length(these.smps)  == 1) {
			total.targeted <- sum(tmp[,these.smps]*tmp$Length*tmp$TargetRegions);
			patient.targeted <- total.targeted;
			} else {
			total.patient <- sum(
				apply(tmp[,these.smps], 1, prod) * tmp$Length
				);

			total.targeted <- apply(
				tmp[,these.smps],
				2,
				function(i) { sum(i*tmp$Length*tmp$TargetRegions) }
				);

			patient.targeted <- sum(
				apply(tmp[,these.smps], 1, prod) * tmp$Length * tmp$TargetRegions
				);
			}

		callable.bases[smp.idx,]$TargetRegions <- sum(tmp$TargetRegions * tmp$Length);
		callable.bases[smp.idx,]$Total.targeted <- total.targeted;
		callable.bases[smp.idx,]$Patient.targeted <- patient.targeted;
		}
	}

# save data to file
save(
	cov.list,
	file = generate.filename(arguments$project, 'CallableBases','RData')
	);

# save total bases to file
write.table(
	callable.bases,
	file = generate.filename(arguments$project, 'total_bases_covered','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectCallableBases','SessionProfile','txt'));
