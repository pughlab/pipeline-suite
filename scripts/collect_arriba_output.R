### collect_arriba_output.R ########################################################################
# find, collate and format output from Arriba

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
fusion.files <- list.files(pattern = 'fusions.tsv', recursive = TRUE);
virus.files  <- list.files(pattern = 'virus_expression.tsv', recursive = TRUE);

# read them in and store them
fusion.data <- list();

fusions.empty.samples <- c();
viral.empty.samples <- c();

for (file in fusion.files) {

	# extract sample ID
	smp <- unlist(strsplit(file,'\\/'))[2];
	# read in data
	tmp <- read.delim(file, stringsAsFactors = FALSE);

	if (nrow(tmp) == 0) {
		fusions.empty.samples <- c(fusions.empty.samples, smp);
		next;
		}

	# store it
	tmp$Sample <- smp;
	idx <- which(colnames(tmp) == 'read_identifiers'); # remove read_identifiers field
	fusion.data[[length(fusion.data)+1]] <- tmp[,-idx];
	}

virus.data <- list();

for (file in virus.files) {

	# extract sample ID
	smp <- unlist(strsplit(file,'\\/'))[2];

	# read in data
	tmp <- read.delim(file, stringsAsFactors = FALSE);
	if (nrow(tmp) == 0) {
		viral.empty.samples <- c(viral.empty.samples, smp);
		next;
		}

	# store it
	tmp$Sample <- smp;
	virus.data[[length(virus.data)+1]] <- tmp;
	}

# format data
combined.fusions <- do.call(rbind, fusion.data);
combined.virus <- do.call(rbind, virus.data);

colnames(combined.fusions)[1] <- 'gene1';

keep.fields <- c('Sample','gene1','gene2','gene_id1','gene_id2','breakpoint1','breakpoint2','site1','site2','type','split_reads1','split_reads2','confidence','reading_frame');
virus.fields <- setdiff(colnames(combined.virus),'Sample');

formatted.fusions <- unique(combined.fusions[,keep.fields]);

# organize data
formatted.fusions[which(formatted.fusions$site1 == 'intergenic'),]$gene1 <- 'None';
formatted.fusions[which(formatted.fusions$site2 == 'intergenic'),]$gene2 <- 'None';

formatted.fusions$Fusion <- apply(
	formatted.fusions[,c('gene1','gene2')],
	1,
	function(i) {
		paste(sort(i), collapse = '--');
		}
	);
formatted.fusions$SVType <- sapply(formatted.fusions$type, function(i) { unlist(strsplit(i,'\\/'))[1] } );
formatted.fusions$confidence <- factor(formatted.fusions$confidence, levels = c('high','medium','low'));

formatted.fusions <- formatted.fusions[order(formatted.fusions$confidence),];

# create fusion x patient matrix (for recurrence)
tmp <- unique(formatted.fusions[!grepl('None', formatted.fusions$Fusion),c('Sample','Fusion')]);
tmp$Gene1 <- sapply(tmp$Fusion, function(i) { unlist(strsplit(i,'--'))[1] } );
tmp$Gene2 <- sapply(tmp$Fusion, function(i) { unlist(strsplit(i,'--'))[2] } );
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

# now for viral count data
formatted.viral <- reshape(
	combined.virus[,c('Sample','VIRUS','RPKM')],
	direction = 'wide',
	timevar = 'Sample',
	idvar = 'VIRUS'
	);

colnames(formatted.viral) <- gsub('RPKM.','',colnames(formatted.viral));
colnames(formatted.viral)[1] <- 'Species';

# add in empty samples
for (smp in viral.empty.samples) {
	formatted.viral[,smp] <- NA;
	}

# save data to file
write.table(
	formatted.fusions,
	file = generate.filename(arguments$project, 'arriba_output_long', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

write.table(
	formatted.fusions.by.patient,
	file = generate.filename(arguments$project, 'arriba_output_wide', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

write.table(
	formatted.viral, #combined.virus[,c('Sample',virus.fields)],
	file = generate.filename(arguments$project, 'arriba_viral_expression', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# format data for cbioportal
tmp <- formatted.fusions[which(formatted.fusions$confidence == 'high'),];

cbio.data <- data.frame(
	Hugo_Symbol = c(
		sapply(tmp$Fusion, function(i) { unlist(strsplit(as.character(i),'--'))[1] } ),
		sapply(tmp$Fusion, function(i) { unlist(strsplit(as.character(i),'--'))[2] } )
		),
	Entrez_Gene_Id = NA,
	Center = NA,
	Tumor_Sample_Barcode = rep(tmp$Sample, times = 2),
	Fusion = rep(tmp$Fusion, times = 2), # gsub('--','-',rep(tmp$Fusion, times = 2)),
	DNA_support = 'no',
	RNA_support = 'yes',
	Method = 'arriba',
	Frame = rep(tmp$reading_frame, times = 2),
	Fusion_Status = NA,
	stringsAsFactors = FALSE
	);

cbio.data$Frame <- factor(
	cbio.data$Frame,
	levels = c('out-of-frame','in-frame'),
	labels = c('frameshift','inframe')
	);

write.table(
	unique(cbio.data),
	file = generate.filename(arguments$project, 'arriba_for_cbioportal', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	quote = FALSE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectArriba','SessionProfile','txt'));
