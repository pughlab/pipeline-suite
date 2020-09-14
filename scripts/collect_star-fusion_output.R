### collect_star-fusion_output.R ###################################################################
# find, collate and format output from STAR-Fusion

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
fusion.files <- list.files(pattern = 'star-fusion.fusion_predictions.abridged.tsv', recursive = TRUE);

# read them in and store them
fusion.data <- list();

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

# format data
combined.fusions <- do.call(rbind, fusion.data);

formatted.fusions <- unique(combined.fusions[,c(16,1:15)]);

tmp <- unique(formatted.fusions[,1:2]);

tmp$Gene1 <- sapply(
	tmp$X.FusionName,
	function(i) {
		unlist(strsplit(as.character(i),'--'))[1];
		}
	);

tmp$Gene2 <- sapply(
	tmp$X.FusionName,
	function(i) {
		unlist(strsplit(as.character(i),'--'))[2];
		}
	);

tmp$Call <- 1;

formatted.fusions.by.patient <- reshape(
	tmp[,-2],
	direction = 'wide',
	timevar = 'Sample',
	idvar = c('Gene1','Gene2')
	);
colnames(formatted.fusions.by.patient) <- gsub('Call.','',colnames(formatted.fusions.by.patient));

# add in empty samples
for (smp in fusions.empty.samples) {
	formatted.fusions.by.patient[,smp] <- NA;
	}

# save data to file
write.table(
	formatted.fusions,
	file = generate.filename(arguments$project, 'star-fusion_output_long', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

write.table(
	formatted.fusions.by.patient,
	file = generate.filename(arguments$project, 'star-fusion_output_wide', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# format data for cbioportal
tmp <- formatted.fusions[which(formatted.fusions$SpanningFragCount > 0),];

cbio.data <- data.frame(
	Hugo_Symbol = c(
		sapply(tmp$X.FusionName, function(i) { unlist(strsplit(as.character(i),'--'))[1] } ),
		sapply(tmp$X.FusionName, function(i) { unlist(strsplit(as.character(i),'--'))[2] } )
		),
	Entrez_Gene_Id = NA,
	Center = NA,
	Tumor_Sample_Barcode = rep(tmp$Sample, times = 2),
	Fusion = rep(tmp$X.FusionName, times = 2),
	DNA_support = 'no',
	RNA_support = 'yes',
	Method = 'STAR-Fusion',
	Frame = 'in-frame',
	Fusion_Status = '',
	stringsAsFactors = FALSE
	);

tmp <- tmp[which(tmp$SpliceType == 'INCL_NON_REF_SPLICE'),];
for (i in 1:nrow(tmp)) {
	idx <- which(cbio.data$Tumor_Sample_Barcode == tmp[i,]$Sample & cbio.data$Fusion == tmp[i,]$X.FusionName);
	cbio.data[idx,]$Frame <- 'frameshift';
	}

write.table(
	unique(cbio.data),
	file = generate.filename(arguments$project, 'star-fusion_for_cbioportal', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectSTARFusion','SessionProfile','txt'));
