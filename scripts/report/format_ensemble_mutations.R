### format_ensemble_mutations.R ####################################################################
# Collect and filter mutation calls from multiple tools using an ensemble approach.
# INPUT (must be in MAF format):
#	tool-specific mutation calls (output by collect_snv_output.R [usually DATE_PROJECT_mutations_for_cbioportal.tsv from each tool directory])

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

	# write key variables to file
	cat("\n### VARIABLES #################################################################\n");
	cat(paste0('Experiment type: ', arguments$seq_type));
	cat(paste0('Target coverage: ', arguments$coverage));
	cat(paste0('VAF threshold for final filter: ', vaf.threshold));

	# write sessionInfo to file
	cat("\n### SESSION INFO ###############################################################");
	print(sessionInfo());

	# close the file
	sink();

	}

### PREPARE SESSION ################################################################################
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('--mutect', type = 'character', help = 'path to combined mutect output');
parser$add_argument('--mutect2', type = 'character', help = 'path to combined mutect2 output');
parser$add_argument('--pindel', type = 'character', help = 'path to combined pindel output');
parser$add_argument('--strelka', type = 'character', help = 'path to combined strelka output');
parser$add_argument('--somaticsniper', type = 'character', help = 'path to combined somaticsniper output');
parser$add_argument('--varscan', type = 'character', help = 'path to combined varscan output');
parser$add_argument('--vardict', type = 'character', help = 'path to combined vardict output');
parser$add_argument('-t', '--seq_type', type = 'character', help = 'type of sequencing experiment');
parser$add_argument('-c', '--coverage', type = 'character', help = 'minimum depth for tumour and normal to be considered callable; length 1 or 2 (comma separated)', default = '20,15');

arguments <- parser$parse_args();

# import libraries
library(BoutrosLab.plotting.general);

# do some quick error checks
run.mutect <- !is.null(arguments$mutect);
run.mutect2 <- !is.null(arguments$mutect2);
run.pindel <- !is.null(arguments$pindel);
run.strelka <- !is.null(arguments$strelka);
run.somaticsniper <- !is.null(arguments$somaticsniper);
run.varscan <- !is.null(arguments$varscan);
run.vardict <- !is.null(arguments$vardict);

tool.count <- sum(run.mutect,run.mutect2,run.pindel,run.strelka,run.varscan,run.somaticsniper,run.vardict);
snp.tool.count <- sum(run.mutect,run.mutect2,run.strelka,run.varscan,run.somaticsniper,run.vardict);
indel.tool.count <- sum(run.mutect2,run.pindel,run.strelka,run.varscan,run.vardict);

if (tool.count < 2) {
	stop('Must provide path to input file or paths to 2+ tool-specific outputs');
	}

# determine minimum depth for callable flag
if (is.null(arguments$coverage) || is.na(arguments$coverage)) {
	t_depth <- 0;
	n_depth <- 0;
	} else if (length(unlist(strsplit(arguments$coverage,','))) == 2) {
	t_depth <- as.numeric(unlist(strsplit(arguments$coverage,','))[1]);
	n_depth <- as.numeric(unlist(strsplit(arguments$coverage,','))[2]);
	} else if (length(unlist(strsplit(arguments$coverage,','))) == 1) {
	t_depth <- as.numeric(unlist(strsplit(arguments$coverage,','))[1]);
	n_depth <- t_depth;
	}

if (is.na(t_depth)) { t_depth <- 0; }
if (is.na(n_depth)) { n_depth <- 0; }

### READ DATA ######################################################################################
# get data
mutation.data <- list();

if (!is.null(arguments$strelka)) {
	mutation.data[['Strelka']] <- read.delim(arguments$strelka, comment.char = '#');
	}
if (!is.null(arguments$somaticsniper)) {
	mutation.data[['SomaticSniper']] <- read.delim(arguments$somaticsniper, comment.char = '#');
	}
if (!is.null(arguments$varscan)) {
	mutation.data[['VarScan']] <- read.delim(arguments$varscan, comment.char = '#');
	}
if (!is.null(arguments$vardict)) {
	mutation.data[['VarDict']] <- read.delim(arguments$vardict, comment.char = '#');
	}
if (!is.null(arguments$mutect2)) {
	mutation.data[['MuTect2']] <- read.delim(arguments$mutect2, comment.char = '#');
	}
if (!is.null(arguments$mutect)) {
	mutation.data[['MuTect']] <- read.delim(arguments$mutect, comment.char = '#');
	}
if (!is.null(arguments$pindel)) {
	mutation.data[['Pindel']] <- read.delim(arguments$pindel, comment.char = '#');
	}

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# create minimal tables for overlap
keep.fields <- c('Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode','Hugo_Symbol','Entrez_Gene_Id','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type');

combined.data <- data.frame();

for (tool in names(mutation.data)) {

	tool.data <- mutation.data[[tool]];

	tmp <- tool.data[,c(keep.fields,'t_depth','t_alt_count')];
	tmp$VAF <- tmp$t_alt_count / tmp$t_depth;
	tmp <- tmp[,c(keep.fields,'VAF')];
	colnames(tmp)[which(colnames(tmp) == 'VAF')] <- tool;

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

	rm(tool.data);
	}

# how many times was this variant/sample called?
combined.data$Count <- apply(
	combined.data[,names(mutation.data)],
	1,
	function(i) { length(i[which(i > 0)]) }
	);

### FILTER VARIANTS ################################################################################
# apply the following criteria to filter variants:
#	1) is snv and called by a minimum 4+ tools (or 50% of tools)
#		EXCEPTION is if this is a tumour-only sample because somaticsniper doesn't run these
#	2) is indel and called by 3+ tools (or 50% of indel-callers) (because mutect does not call indels)
#	3) is called by MuTect2 and VAF is < 0.1
#	4) keep position if it passes any of the above in any sample AND called by MuTect2

combined.data[,c('FILTER.n_tools','FILTER.m2_vaf')] <- FALSE;

min.snp.count <- ceiling(snp.tool.count*0.55);
min.indel.count <- ceiling(indel.tool.count*0.55);

# 1) is snv and called by a minimum n_tools
is.snp <- combined.data$Variant_Type == 'SNP';
callers.min <- combined.data$Count >= min.snp.count;

if (length(which(is.snp & callers.min)) > 0) {
	combined.data[which(is.snp & callers.min),]$FILTER.n_tools <- TRUE;
	}

# reduce n_tools if tumour_only (due to SomaticSniper) [only required if cohort is mixed T/N and T-only]
if (run.somaticsniper) {
	callers.min.mod <- combined.data$Count >= ceiling((snp.tool.count-1)*0.55);
	t.only.idx <- is.na(combined.data$Matched_Norm_Sample_Barcode);
	if (length(which(is.snp & callers.min.mod & t.only.idx)) > 0) {
		combined.data[which(is.snp & callers.min.mod & is.na(combined.data$Matched_Norm_Sample_Barcode)),]$FILTER.n_tools <- TRUE;
		}
	}

# 2) is indel and called by 2+ tools (because mutect and somaticsniper do not call indels)
is.indel <- combined.data$Variant_Type != 'SNP';
callers.min <- combined.data$Count >= min.indel.count;

if (length(which(is.indel & callers.min)) > 0) {
	combined.data[which(is.indel & callers.min),]$FILTER.n_tools <- TRUE;
	}

# 3) is called by MuTect2 and VAF is < 0.1 (because many tools don't call low vaf/coverage variants)
if ('MuTect2' %in% colnames(combined.data)) {
	low.vaf.idx <- which(combined.data$MuTect2 < 0.1 & combined.data$MuTect2 > 0);
	if (length(low.vaf.idx) > 0) { combined.data[low.vaf.idx,]$FILTER.m2_vaf <- TRUE; }
	}

# 4) keep position if it passes any of the above in any sample from this patient
# this should catch most instances of a variant with low VAF/coverage in 1 part of a multi-region tumour
passed.variants <- unique(
	combined.data[which(combined.data$FILTER.n_tools | combined.data$FILTER.m2_vaf),
	c('Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type','Matched_Norm_Sample_Barcode')]
	);

passed.variants <- passed.variants[!is.na(passed.variants$Matched_Norm_Sample_Barcode),];

if (nrow(passed.variants) > 0) {

	passed.variants$FILTER.patient_evidence <- TRUE;

	combined.data <- merge(
		combined.data,
		passed.variants,
		all.x = TRUE
		);
	}

combined.data$FILTER <- NA;
combined.data[which(combined.data$FILTER.n_tools |
	combined.data$FILTER.m2_vaf |
	combined.data$FILTER.patient_evidence),]$FILTER <- 'PASS';

save(combined.data, file = generate.filename(arguments$project, 'CombinedMutationData','RData'));

# run unlink, in case it exists from a previous run
unlink('CombinedMutationData.RData');

file.symlink(
	generate.filename(arguments$project, 'CombinedMutationData', 'RData'),
	'CombinedMutationData.RData'
	);

# filter data
filtered.calls <- unique(combined.data[which(combined.data$FILTER == 'PASS'),]);

# order data for easier formatting
filtered.calls$Chromosome <- factor(filtered.calls$Chromosome, levels = paste0('chr',c(1:22,'X','Y')));
filtered.calls <- filtered.calls[order(filtered.calls$Chromosome, filtered.calls$Start_Position),];

# add column indicating which tools called it
filtered.calls$Called_By <- apply(
	filtered.calls[,names(mutation.data)],
	1,
	function(i) { paste(names(mutation.data)[!is.na(i)], collapse = ',') }
	);

# pull out entries from base mafs 
tmp <- filtered.calls;
keep.data <- list();

# this will prioritize tools in the order they were read in (ie, Strelka first)
# for each tool
for (tool in names(mutation.data)) {

	# extract calls from the tool
	keep.data[[tool]] <- merge(
		tmp[!is.na(tmp[,tool]),c('Chromosome','Start_Position','Variant_Type','Tumor_Sample_Barcode','Called_By')],
		mutation.data[[tool]],
		all.x = TRUE
		);

	# then remove them from the list
	tmp <- tmp[is.na(tmp[,tool]),];

	# repeat until tmp is empty
	if (nrow(tmp) == 0) { break; }
	}

keep.calls <- do.call(rbind, keep.data);

# order columns
keep.calls <- keep.calls[,c(colnames(mutation.data[[1]]), 'Called_By')];


# add some quick QC flags
annotated.data <- keep.calls;
annotated.data[,c('FLAG.low_vaf','FLAG.low_coverage','FLAG.high_pop','FLAG.tumour_only')] <- FALSE;

# low VAF flag
vaf <- annotated.data$t_alt_count / annotated.data$t_depth;
if (any(na.omit(vaf) < 0.05)) {
	annotated.data[which(vaf < 0.05),]$FLAG.low_vaf <- TRUE;
	}

# low COVERAGE (non-callable base)
if (any(as.numeric(annotated.data$t_depth) < t_depth)) {
	annotated.data[which(as.numeric(annotated.data$t_depth) < t_depth),]$FLAG.low_coverage <- TRUE;
	}
if (!all(is.na(annotated.data$n_depth)) & (any(as.numeric(annotated.data$n_depth) < n_depth))) {
	annotated.data[which(as.numeric(annotated.data$n_depth) < n_depth),]$FLAG.low_coverage <- TRUE;
	}

# high population frequency
af.fields <- colnames(annotated.data)[grepl('AF', colnames(annotated.data))];
pop.freq <- apply(annotated.data[,af.fields],1,function(i) { max(i,na.rm = TRUE) > 0.001 } );
if (any(pop.freq)) {
	annotated.data[which(pop.freq),]$FLAG.high_pop <- TRUE;
	}

# and did this sample have a matched normal
if (any(is.na(annotated.data$Matched_Norm_Sample_Barcode))) {
	annotated.data[is.na(annotated.data$Matched_Norm_Sample_Barcode),]$FLAG.tumour_only <- TRUE;
	}

### FINAL FORMATTING ###############################################################################
# do some minor formatting for cbioportal (might already be done)
for (field in c('t_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count')) {
	if (any(is.na(annotated.data[,field]))) {
		annotated.data[is.na(annotated.data[,field]),field] <- '';
		}
	}

annotated.data <- unique(annotated.data);

# save to file
write.table(
	annotated.data,
	file = generate.filename(arguments$project, 'ensemble_mutation_data', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	quote = FALSE,
	sep = '\t'
	);

# apply some basic filters
# higher depth = lower VAF threshold
vaf.threshold <- if ('targeted' == arguments$seq_type) { 0.005;
	} else if ('exome' == arguments$seq_type) { 0.01
	} else { 0.05; }

for (field in c('t_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count')) {
	annotated.data[,field] <- as.numeric(annotated.data[,field]);
	}

tumour.keep <- which(
	annotated.data$t_depth > 10 & 
	(annotated.data$t_alt_count / annotated.data$t_depth) > vaf.threshold
	);

normal.keep <- which(
	annotated.data$FLAG.tumour_only |
	(annotated.data$n_depth > 10 & (1 - ( annotated.data$n_ref_count / annotated.data$n_depth ) < 0.1) )
	);

annotated.filtered <- annotated.data[intersect(tumour.keep, normal.keep),];

# do some minor formatting for cbioportal (might already be done)
for (field in c('t_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count')) {
	if (any(is.na(annotated.filtered[,field]))) {
		annotated.filtered[is.na(annotated.filtered[,field]),field] <- '';
		}
	}

# save to file
write('# ENSEMBLE MAF', file = generate.filename(arguments$project, 'ensemble_mutation_data_filtered', 'tsv'));
write.table(
	annotated.filtered,
	file = generate.filename(arguments$project, 'ensemble_mutation_data_filtered', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	quote = FALSE,
	append = TRUE,
	sep = '\t'
	);

# run unlink, in case it exists from a previous run
unlink('ensemble_mutation_data.tsv');

file.symlink(
	generate.filename(arguments$project, 'ensemble_mutation_data_filtered', 'tsv'),
	'ensemble_mutation_data.tsv'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('format_ensemble_mutations','SessionProfile','txt'));
