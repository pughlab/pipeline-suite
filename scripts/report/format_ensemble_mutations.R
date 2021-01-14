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

# function to trim sample IDs
simplify.ids <- function(x) {
	match <- TRUE;
	index <- 1;
	while (match) {
		if (length(unique(sapply(x,function(i) { unlist(strsplit(i,''))[index] } ))) == 1) {
			index <- index + 1;
			} else {
			new.ids <- sapply(x,function(i) { substring(i,index,nchar(i)) } );
			match <- FALSE;
			}
		}
	return(new.ids);
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
parser$add_argument('-n', '--n_tools', type = 'character', help = 'minimum number of tools required to call a variant');

parser$add_argument('-i', '--input', type = 'character', help = 'tab-delimited file listing Tool and Path to tool output; overrides tool-specific input (--mutect, --mutect2, --strelka, --varscan --vardict --somaticsniper)');

arguments <- parser$parse_args();

# do some quick error checks
run.mutect <- !is.null(arguments$mutect);
run.mutect2 <- !is.null(arguments$mutect2);
run.strelka <- !is.null(arguments$strelka);
run.somaticsniper <- !is.null(arguments$somaticsniper);
run.varscan <- !is.null(arguments$varscan);
run.vardict <- !is.null(arguments$vardict);

tool.count <- sum(run.mutect,run.mutect2,run.strelka,run.varscan,run.somaticsniper,run.vardict);

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
		mutation.data[['MuTect2']] <- read.delim(arguments$mutect2);
		}
	if (!is.null(arguments$mutect)) {
		mutation.data[['MuTect']] <- read.delim(arguments$mutect);
		}
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
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# create minimal tables for overlap
keep.fields <- c('Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode','Hugo_Symbol','Entrez_Gene_Id','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type');

combined.data <- data.frame();

for (tool in names(mutation.data)) {

	tool.data <- mutation.data[[tool]];
	if (any(is.na(tool.data$n_depth))) {
		tool.data[is.na(tool.data$n_depth),]$Matched_Norm_Sample_Barcode <- NA;
		}

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

min.snp.count <- if (!is.null(arguments$n_tools)) { arguments$n_tools
	} else { ceiling(length(mutation.data)*0.5)
	}

min.indel.count <- if (sum(run.mutect2, run.varscan, run.strelka, run.vardict) < 4) { 2
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
passed.variants <- unique(combined.data[which(combined.data$FILTER == 'PASS'),c('Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type')]);

passed.variants$C4 <- 1;

combined.data <- merge(
	combined.data,
	passed.variants,
	all.x = TRUE
	);

combined.data[which(combined.data$C4 == 1 & combined.data$MuTect2 == 1),]$FILTER <- 'PASS';


# filter data
filtered.calls <- unique(combined.data[which(combined.data$FILTER == 'PASS'),]);

# order data for easier formatting
filtered.calls$Chromosome <- factor(filtered.calls$Chromosome, levels = paste0('chr',c(1:22,'X','Y')));
filtered.calls <- filtered.calls[order(filtered.calls$Chromosome, filtered.calls$Start_Position),];

# extract unique entries for each of these filtered calls
keep.data <- list();

for (i in 1:nrow(filtered.calls)) {

	smp <- as.character(filtered.calls[i,]$Tumor_Sample_Barcode);
	chr <- as.character(filtered.calls[i,]$Chromosome);
	start <- filtered.calls[i,]$Start_Position;
	type <- as.character(filtered.calls[i,]$Variant_Type);

	for (tool in names(mutation.data)) {
		tool.idx <- which(
			mutation.data[[tool]]$Tumor_Sample_Barcode == smp & 
			mutation.data[[tool]]$Chromosome == chr & 
			mutation.data[[tool]]$Start_Position == start &
			mutation.data[[tool]]$Variant_Type == type
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

### PLOT DATA ######################################################################################
# if all six tools were used, plot overlap summary
if (6 == tool.count) {

	# initiate plot data
	plot.data <- list();
	plot.data$per_tool <- merge(
		aggregate(MuTect2 ~ Tumor_Sample_Barcode, combined.data, sum),
		merge(
			aggregate(MuTect ~ Tumor_Sample_Barcode, combined.data, sum),
			merge(
				aggregate(Strelka ~ Tumor_Sample_Barcode, combined.data, sum),
				merge(
					aggregate(VarScan ~ Tumor_Sample_Barcode, combined.data, sum),
					merge(
						aggregate(VarDict ~ Tumor_Sample_Barcode, combined.data, sum),
						aggregate(SomaticSniper ~ Tumor_Sample_Barcode, combined.data, sum),
						all = TRUE
						),
					all = TRUE
					),
				all = TRUE
				),
			all = TRUE
			),
		all = TRUE
		);

	plot.data$overlap <- merge(
		aggregate(Hugo_Symbol ~ Tumor_Sample_Barcode + Count, combined.data, length),
		aggregate(Hugo_Symbol ~ Tumor_Sample_Barcode, combined.data, length),
		by = 'Tumor_Sample_Barcode'
		);
	colnames(plot.data$overlap) <- c('ID','Overlap','Count','Total');
	plot.data$overlap$Overlap <- factor(plot.data$overlap$Overlap, levels = c(1,2,3,4,5,6));
	plot.data$overlap$Percent <- plot.data$overlap$Count / plot.data$overlap$Total * 100;

	save(
		plot.data,
		file = generate.filename(arguments$project, 'mutation_overlap', 'RData')
		);

	# grab some parameters
	axis.cex <- if (nrow(plot.data$per_tool) <= 30) { 1
		} else if (nrow(plot.data$per_tool) <= 50) { 0.75
		} else if (nrow(plot.data$per_tool) <= 80) { 0.5
		} else { 0 };

	# create some plots to summarize overlap
	plot.objects <- list();

	# plot per-tool counts
	for (tool in c('MuTect2','MuTect','Strelka','VarScan','SomaticSniper','VarDict')) {

		plot.objects[[tool]] <- create.barplot(
			get(tool) ~ Tumor_Sample_Barcode | tool,
			plot.data$per_tool,
			ylab.label = NULL,
			xlab.label = NULL,
			yaxis.cex = 0.8,
			xaxis.lab = rep('',nrow(plot.data$per_tool)),
			yaxis.tck = c(0.5,0),
			xaxis.tck = 0,
			yaxis.fontface = 'plain'
			);
		}

	overlap.legend <- legend.grob(
		legends = list(
			legend = list(
				colours = rev(default.colours(6,'seq.greenblue')),
				labels = rev(c('1','2','3','4','5','6')),
				title = 'Tool Count'
				)
			),
		title.just = 'left',
		size = 2
		);

	# plot overlap
	plot.objects[['overlap']] <- create.barplot(
		Percent ~ ID,
		plot.data$overlap,
		groups = plot.data$overlap$Overlap,
		stack = TRUE,
		col = default.colours(6,'seq.greenblue'),
		ylab.label = '% of Total',
		xlab.label = NULL,
		ylimits = c(0,100),
		yat = seq(0,100,20),
		ylab.cex = 1,
		yaxis.cex = 1,
		xaxis.lab = simplify.ids(levels(plot.data$overlap$ID)),
		xaxis.cex = axis.cex,
		yaxis.tck = c(0.5,0),
		xaxis.tck = if (axis.cex == 0) { 0 } else { c(0.5, 0) },
		xaxis.rot = 90,
		yaxis.fontface = 'plain',
		xaxis.fontface = 'plain',
		legend = list(
			right = list(fun = overlap.legend)
			)			
		);

	# combine them!
	create.multipanelplot(
		plot.objects = plot.objects,
		height = 10,
		width = 8,
		resolution = 200,
		filename = generate.filename(arguments$project, 'mutation_overlap','png'),
		plot.objects.heights = c(rep(1,6),2.5),
		left.legend.padding = 0,
		right.legend.padding = 0,
		top.legend.padding = 0,
		bottom.legend.padding = 0,
		y.spacing = c(rep(-2.5,5),-1),
		ylab.label = paste(
			paste(rep(' ', 50), collapse = ''),
			'Somatic Variant Count (SNVs and INDELs)'
			),
		ylab.cex = 1,
		left.padding = 0,
		ylab.axis.padding = 0,
		bottom.padding = 0,
		top.padding = -1,
		right.padding = 0
		);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('format_ensemble_mutations','SessionProfile','txt'));
