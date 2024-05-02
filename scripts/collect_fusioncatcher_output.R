### collect_fusioncatcher_output.R #################################################################
# find, collate and format output from Fusioncatcher

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
fusion.files <- list.files(pattern = 'final-list_candidate-fusion-genes.txt', recursive = TRUE);
virus.files <- list.files(pattern = 'viruses_bacteria_phages.txt', recursive = TRUE);

# read them in and store them
fusion.data <- list();
viral.data <- list();

fusions.empty.samples <- c();

for (file in fusion.files) {

	# extract sample ID
	smp <- unlist(strsplit(file,'\\/'))[2];
	# read in data
	tmp <- read.delim(file);

	if (nrow(tmp) == 0) {
		fusions.empty.samples <- c(fusions.empty.samples, smp);
		next;
		}

	# store it
	tmp$Sample <- smp;
	fusion.data[[length(fusion.data)+1]] <- tmp;
	}

viral.empty.samples <- c();

for (file in virus.files) {

	# extract sample ID
	smp <- unlist(strsplit(file,'\\/'))[2];
	# read in data
	tmp <- read.table(file, header = TRUE);

	if (nrow(tmp) == 0) {
		viral.empty.samples <- c(viral.empty.samples, smp);
		next;
		}

	# store it
	tmp$Sample <- smp;
	viral.data[[length(viral.data)+1]] <- tmp;
	}

# format data
combined.fusions <- do.call(rbind, fusion.data);
combined.viral 	 <- do.call(rbind, viral.data);

formatted.fusions <- unique(combined.fusions[,c(17,1:16)]);

tmp <- unique(formatted.fusions[,1:3]);
tmp$Call <- 1;

formatted.fusions.by.patient <- reshape(
	tmp,
	direction = 'wide',
	timevar = 'Sample',
	idvar = c('Gene_1_symbol.5end_fusion_partner.','Gene_2_symbol.3end_fusion_partner.')
	);
colnames(formatted.fusions.by.patient) <- gsub('Call.','',colnames(formatted.fusions.by.patient));

# add in empty samples
for (smp in fusions.empty.samples) {
	formatted.fusions.by.patient[,smp] <- NA;
	}

# now for viral count data
formatted.viral <- reshape(
	combined.viral,
	direction = 'wide',
	timevar = 'Sample',
	idvar = 'Virus.Bacteria.Phage'
	);
colnames(formatted.viral) <- gsub('Counts_of_mapping_reads.','',colnames(formatted.viral));

# add in empty samples
for (smp in viral.empty.samples) {
	formatted.viral[,smp] <- NA;
	}

# save data to file
write.table(
	formatted.fusions,
	file = generate.filename(arguments$project, 'fusioncatcher_output_long', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

write.table(
	formatted.fusions.by.patient,
	file = generate.filename(arguments$project, 'fusioncatcher_output_wide', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

write.table(
	formatted.viral,
	file = generate.filename(arguments$project, 'fusioncatcher_viral_counts', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# format for cbioportal
tmp <- formatted.fusions[which(formatted.fusions$Spanning_unique_reads > 0),];
tmp <- tmp[which(tmp$Predicted_effect %in% c('in-frame','out-of-frame')),];

cbio.data <- data.frame(
	Hugo_Symbol = c(
		as.character(tmp$Gene_1_symbol.5end_fusion_partner.),
		as.character(tmp$Gene_2_symbol.3end_fusion_partner.)
		),
	Entrez_Gene_Id = NA,
	Center = NA,
	Tumor_Sample_Barcode = rep(tmp$Sample, times = 2),
	Fusion = NA,
	DNA_support = 'no',
	RNA_support = 'yes',
	Method = 'Fusioncatcher',
	Frame = rep(tmp$Predicted_effect, times = 2),
	Fusion_Status = '',
	stringsAsFactors = FALSE
	);

cbio.data$Fusion <- rep(
	paste0(tmp$Gene_1_symbol.5end_fusion_partner., '--', tmp$Gene_2_symbol.3end_fusion_partner.),
	times = 2
	);

cbio.data$Frame <- factor(cbio.data$Frame, levels = c('in-frame','out-of-frame'), labels = c('inframe','frameshift'));

write.table(
	unique(cbio.data),
	file = generate.filename(arguments$project, 'fusioncatcher_for_cbioportal', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	quote = FALSE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectFusionCatcher','SessionProfile','txt'));
