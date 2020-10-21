### format_ensemble_mutations.R ####################################################################
# Collect and filter mutation calls from multiple tools using an ensemble approach.
# INPUT (must be in MAF format):
#	tool-specific mutation calls (output by collect_snv_output.R [usually DATE_PROJECT_mutations_for_cbioportal.tsv]) OR
#	2-column, tab-delimited file listing Tool and Path [ordered by priority]

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

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('--mutect', type = 'character', help = 'path to combined mutect output');
parser$add_argument('--mutect2', type = 'character', help = 'path to combined mutect2 output');
parser$add_argument('--strelka', type = 'character', help = 'path to combined strelka output');
parser$add_argument('--varscan', type = 'character', help = 'path to combined varscan output');
parser$add_argument('-n', '--n_tools', type = 'character', help = 'minimum number of tools required to call a variant', default = 3);
parser$add_argument('-i', '--input', type = 'character', help = 'tab-delimited file listing Tool and Path to tool output; overrides tool-specific input (--mutect, --mutect2, --strelka, --varscan)');

arguments <- parser$parse_args();

# do some quick error checks
run.mutect <- !is.null(arguments$mutect);
run.mutect2 <- !is.null(arguments$mutect2);
run.strelka <- !is.null(arguments$strelka);
run.varscan <- !is.null(arguments$varscan);

tool.count <- sum(run.mutect,run.mutect2,run.strelka,run.varscan);

if (is.null(arguments$input) & (tool.count < 2)) {
	stop('Must provide path to input file or paths to 2+ tool-specific outputs');
	}

### READ DATA ######################################################################################
# get data
mutation.data <- list();

if (!is.null(arguments$input)) {

	input <- read.delim(arguments$input, as.is = TRUE);

	for (i in 1:nrow(input)) {
		tool <- as.character(input[i,1]);
		mutation.data[[tool]] <- read.delim(input[i,2]);
		}

	} else {

	if (!is.null(arguments$mutect2)) {
		mutation.data[['Mutect2']] <- read.delim(arguments$mutect2);
		}
	if (!is.null(arguments$mutect)) {
		mutation.data[['MuTect']] <- read.delim(arguments$mutect);
		}
	if (!is.null(arguments$strelka)) {
		mutation.data[['Strelka']] <- read.delim(arguments$strelka);
		}
	if (!is.null(arguments$varscan)) {
		mutation.data[['VarScan']] <- read.delim(arguments$varscan);
		}
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# create minimal tables for overlap
keep.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Entrez_Gene_Id','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2');

combined.data <- data.frame();

for (tool in names(mutation.data)) {

	tmp <- mutation.data[[tool]][,keep.fields];
	tmp[,tool] <- 1;

	if (nrow(combined.data) == 0) {
		combined.data <- tmp;
		} else {

		# merge datasets
		combined.data <- unique(merge(
			combined.data,
			tmp,
			all = TRUE
			));
		}
	}

combined.data$Count <- apply(
	combined.data[,names(mutation.data)],
	1,
	sum,
	na.rm = TRUE
	);

# filter to min callset
filtered.calls <- unique(combined.data[which(combined.data$Count >= arguments$n_tools),]);

filtered.calls$Chromosome <- factor(filtered.calls$Chromosome, levels = paste0('chr',c(1:22,'X','Y')));
filtered.calls <- filtered.calls[order(filtered.calls$Chromosome, filtered.calls$Start_Position),];

# extract unique entries for each of these filtered calls
keep.data <- list();

for (i in 1:nrow(filtered.calls)) {

	smp <- as.character(filtered.calls[i,]$Tumor_Sample_Barcode);
	chr <- as.character(filtered.calls[i,]$Chromosome);
	start <- filtered.calls[i,]$Start_Position;

	for (tool in names(mutation.data)) {

		tool.idx <- which(
			mutation.data[[tool]]$Tumor_Sample_Barcode == smp & 
			mutation.data[[tool]]$Chromosome == chr & 
			mutation.data[[tool]]$Start_Position == start
			);

		if (length(tool.idx) == 1) {
			keep.data[[i]] <- mutation.data[[tool]][tool.idx,];
			break;
			} else { next; }
		}

	# indicate which methods called each variant
	tmp <- filtered.calls[i,names(mutation.data)];
	keep.data[[i]]$Called_By <- paste(colnames(tmp)[which(tmp == 1)], collapse = ',');

	}

keep.calls <- do.call(rbind, keep.data);

# remove unwanted fields
keep.calls <- keep.calls[,!grepl('variant',colnames(keep.calls))];

# write data
write.table(
	keep.calls,
	file = generate.filename(arguments$project, 'ensemble_mutation_data', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

file.symlink(
	generate.filename(arguments$project, 'ensemble_mutation_data', 'tsv'),
	'ensemble_mutation_data.tsv'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('format_ensemble_mutations','SessionProfile','txt'));
