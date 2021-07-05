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

        # write sessionInfo to file
        cat("\n### SESSION INFO ###############################################################");
        print(sessionInfo());

        # close the file
        sink();

        }

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('--mutect', type = 'character', help = 'path to combined mutect output');
parser$add_argument('--mutect2', type = 'character', help = 'path to combined mutect2 output');
parser$add_argument('--strelka', type = 'character', help = 'path to combined strelka output');
parser$add_argument('--somaticsniper', type = 'character', help = 'path to combined somaticsniper output');
parser$add_argument('--varscan', type = 'character', help = 'path to combined varscan output');
parser$add_argument('--vardict', type = 'character', help = 'path to combined vardict output');
parser$add_argument('-c', '--coverage', type = 'character', help = 'minimum depth for tumour and normal to be considered callable; length 1 or 2 (comma separated)', default = '20,15');

arguments <- parser$parse_args();

# do some quick error checks
run.mutect <- !is.null(arguments$mutect);
run.mutect2 <- !is.null(arguments$mutect2);
run.strelka <- !is.null(arguments$strelka);
run.somaticsniper <- !is.null(arguments$somaticsniper);
run.varscan <- !is.null(arguments$varscan);
run.vardict <- !is.null(arguments$vardict);

tool.count <- sum(run.mutect,run.mutect2,run.strelka,run.varscan,run.somaticsniper,run.vardict);

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
	mutation.data[['Strelka']] <- read.delim(arguments$strelka);
	}
if (!is.null(arguments$somaticsniper)) {
	mutation.data[['SomaticSniper']] <- read.delim(arguments$somaticsniper);
	}
if (!is.null(arguments$varscan)) {
	mutation.data[['VarScan']] <- read.delim(arguments$varscan);
	}
if (!is.null(arguments$vardict)) {
	mutation.data[['VarDict']] <- read.delim(arguments$vardict);
	}
if (!is.null(arguments$mutect2)) {
	mutation.data[['MuTect2']] <- read.delim(arguments$mutect2);
	}
if (!is.null(arguments$mutect)) {
	mutation.data[['MuTect']] <- read.delim(arguments$mutect);
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# create minimal tables for overlap
keep.fields <- c('Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode','Hugo_Symbol','Entrez_Gene_Id','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type');

combined.data <- data.frame();

for (tool in names(mutation.data)) {

	tool.data <- mutation.data[[tool]];

	if (tool == 'MuTect2') {
		tmp <- tool.data[,c(keep.fields,'t_depth','t_alt_count')];
		tmp$VAF <- tmp$t_alt_count / tmp$t_depth;
		tmp <- tmp[,c(keep.fields,'VAF')];
		tmp[,tool] <- 1;
		} else {
		tmp <- tool.data[,keep.fields];
		tmp[,tool] <- 1;
		}

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
	sum,
	na.rm = TRUE
	);

### FILTER VARIANTS ################################################################################
# apply the following criteria to filter variants:
#	1) is snv and called by a minimum n_tools (or 50% of tools)
#		EXCEPTION is if this is a tumour-only sample because somaticsniper doesn't run these
#	2) is indel and called by 2+ tools (or 50% of indel-callers) (because mutect does not call indels)
#	3) is called by MuTect2 and VAF is < 0.1
#	4) keep position if it passes any of the above in any sample AND called by MuTect2

combined.data$FILTER <- NA;

min.snp.count <- if (tool.count == 6) { 4
	} else { ceiling(tool.count*0.5)
	}

min.indel.count <- if (sum(run.mutect2, run.varscan, run.strelka, run.vardict) == 4) { 3
	} else { ceiling(sum(run.mutect2, run.varscan, run.strelka, run.vardict)*0.5)
	}

# 1) is snv and called by a minimum n_tools
is.snp <- combined.data$Variant_Type == 'SNP';
callers.min <- combined.data$Count >= min.snp.count;

if (length(which(is.snp & callers.min)) > 0) {
	combined.data[which(is.snp & callers.min),]$FILTER <- 'PASS';
	}

# reduce n_tools if tumour_only (due to SomaticSniper)
if (run.somaticsniper & tool.count <= 4) {
	callers.min.mod <- combined.data$Count >= (min.snp.count-1);
	if (length(which(is.snp & callers.min)) > 0) {
		combined.data[which(is.snp & callers.min.mod & is.na(combined.data$Matched_Norm_Sample_Barcode)),]$FILTER <- 'PASS';
		}
	}

# 2) is indel and called by 2+ tools (because mutect and somaticsniper do not call indels)
is.indel <- combined.data$Variant_Type != 'SNP';
callers.min <- combined.data$Count >= min.indel.count;

if (length(which(is.indel & callers.min)) > 0) {
	combined.data[which(is.indel & callers.min),]$FILTER <- 'PASS';
	}

# 3) is called by MuTect2 and VAF is < 0.1 (because many tools don't call low vaf/coverage variants)
if ('MuTect2' %in% colnames(combined.data)) {
	combined.data[which(combined.data$MuTect2 == 1 & combined.data$VAF < 0.1),]$FILTER <- 'PASS';
	}

# 4) keep position if it passes any of the above in any sample AND called by MuTect2
# this should catch most instances of a variant with low VAF/coverage in 1 part of a multi-region tumour
passed.variants <- unique(combined.data[which(combined.data$FILTER == 'PASS'),c('Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type','Matched_Norm_Sample_Barcode')]);
passed.variants <- passed.variants[!is.na(passed.variants$Matched_Norm_Sample_Barcode),];

if (nrow(passed.variants) > 0) {

	passed.variants$C4 <- 1;

	combined.data <- merge(
		combined.data,
		passed.variants,
		all.x = TRUE
		);

	combined.data[which(combined.data$C4 == 1),]$FILTER <- 'PASS';
	}

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
if (any(vaf < 0.05)) {
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

# save to file
write.table(
	annotated.data,
	file = generate.filename(arguments$project, 'ensemble_mutation_data', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# run unlink, in case it exists from a previous run
unlink('ensemble_mutation_data.tsv');

file.symlink(
	generate.filename(arguments$project, 'ensemble_mutation_data', 'tsv'),
	'ensemble_mutation_data.tsv'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('format_ensemble_mutations','SessionProfile','txt'));
