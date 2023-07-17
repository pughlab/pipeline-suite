### get_hrdetect_signatures.R ######################################################################
# Extract HRD signatures using HRDetect (uses SNV/INDEL/SV/CNV calls)

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
	print(paste0('Input VCFs: ', arguments$vcf_dir));
	print(paste0('Input SVs: ', arguments$sv));
	print(paste0('Input CNVs: ', arguments$cn));

	# write sessionInfo to file
	cat("\n### SESSION INFO ###############################################################");
	print(sessionInfo());

	# close the file
	sink();

	}

# function to trim sample IDs
simplify.ids <- function(x) {
	match <- TRUE;
	if (length(x) == 1) {
		match <- FALSE;
		new.ids <- x;
		}
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
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'hg38 or hg19', default = 'hg38');
parser$add_argument('-v', '--vcf_dir', type = 'character', help = 'path to directory containing VCFs');
parser$add_argument('-s', '--sv', type = 'character', help = 'combined SV calls from Mavis');
parser$add_argument('-c', '--cnv', type = 'character', help = 'cna segments file');

parser$add_argument('-l','--lib_paths', type = 'character',
	help = 'path to library if not installed in default',
	default = '/cluster/projects/pughlab/src/r_lib/library/4.1/');
parser$add_argument('-z', '--report', type = 'character', help = 'path to report directory',
	default = NULL);

arguments <- parser$parse_args();

# load required libraries
if (!is.null(arguments$lib_paths)) {
	.libPaths(c(arguments$lib_paths,.libPaths()));
	}

# import libraries
library(BoutrosLab.plotting.general);
library(xtable);

library(signature.tools.lib);

### READ DATA ######################################################################################
# find VCF files
vcf.files <- list.files(path = arguments$vcf_dir, pattern = '.vcf.gz$', full.names = TRUE);
vcf.files <- vcf.files[!grepl('ensemble_mutation_data',vcf.files)];

if (length(vcf.files) == 0) {
	stop('ERROR: No VCFs found in input directory, please provide path to VCF files.');
	}

all.samples <- as.character(sapply(vcf.files,function(i) { unlist(strsplit(basename(i),'_'))[1] } ));

# get mutation data
if (is.null(arguments$sv)) {
	stop('ERROR: No input SVs provided, please provide path to SV calls from MAVIS.');
	} else {
	sv.data <- read.delim(arguments$sv, stringsAsFactors = FALSE);
	}

if (is.null(arguments$cnv)) {
	stop('ERROR: No CNA segments provided, please provide path to CNA segment file from Sequenza.');
	} else {
	seg.data <- read.delim(arguments$cnv, stringsAsFactors = FALSE);
	}

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output)) {
	dir.create(arguments$output);
	}

setwd(arguments$output);

### FORMAT DATA ####################################################################################
# format CNV data
seg.data <- seg.data[which(seg.data$Sample %in% all.samples),];

cnv.data <- data.frame(
	Sample = seg.data$Sample,
	seg_no = NA,
	Chromosome = gsub('chr','', seg.data$chromosome),
	chromStart = seg.data$start.pos,
	chromEnd = seg.data$end.pos,
	total.copy.number.inNormal = rep(2, nrow(seg.data)),
	minor.copy.number.inNormal = rep(1, nrow(seg.data)),
	total.copy.number.inTumour = seg.data$CNt,
	minor.copy.number.inTumour = seg.data$B
	);

for (smp in all.samples) {
	idx <- which(cnv.data$Sample == smp);
	if (length(idx) == 0) { next; }

	cnv.data[idx,]$seg_no <- 1:nrow(cnv.data[idx,]);
	}

# format SV data
sv.data <- sv.data[which(sv.data$library %in% all.samples),];

# subset only somatic variants
sv.states <- apply(
	sv.data[,grep('diseased_genome|normal_genome', colnames(sv.data))],
	1, function(i) {
		if (any(grepl('germline', i))) { return('germline') } else { return('somatic') } }
	);
sv.data <- sv.data[which(sv.states == 'somatic'),];

# remove short SVs (INDELs)
sv.data$Length <- abs(sv.data$break2_position_end - sv.data$break1_position_start);
sv.data[grepl('translocation', sv.data$event_type),]$Length <- NA;
sv.data <- sv.data[which(is.na(sv.data$Length) | sv.data$Length > 200),];

# add evidence filter
sv.data$Evidence <- apply(
	sv.data[,c('break1_split_reads','break2_split_reads','spanning_reads','linking_split_reads')],
	1, function(i) {
	if (any(i == 'None')) { i[which(i == 'None')] <- 0; }
	i <- sapply(i, function(j) { max(as.numeric(unlist(strsplit(as.character(j),';')))); });
	# if call method uses spanning reads
	if (i[3] >= 20) { return('PASS')
		# if call method uses split reads
		} else if (i[1] >= 10 & i[2] >= 10) { return('PASS')
		} else if (i[4] >= 10) { return('PASS')
		} else { return('FAIL');
		}
	});

sv.data <- sv.data[which(sv.data$Evidence == 'PASS' | grepl(';', sv.data$tools)),];

# format for HRDetect
somatic.svs <- data.frame(
	chrom1 = gsub('chr','',sv.data$break1_chromosome),
	start1 = sv.data$break1_position_start,
	end1 = sv.data$break1_position_end,
	chrom2 = gsub('chr','',sv.data$break2_chromosome),
	start2 = sv.data$break2_position_start,
	end2 = sv.data$break2_position_end,
	sample = sv.data$library,
	svclass = gsub('inverted translocation','translocation', sv.data$event_type)
	);

# indicate signatures to use (fit the 12 breast cancer signatures using the bootstrap signature fit approach)
sigsToUse <- c(1,2,3,5,6,8,13,17,18,20,26,30);
hrdetect.features <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8");

hrdetect.outputs <- list();

# loop over each sample
for (smp in all.samples) {

	# make sample-specific directory
	if (!dir.exists(smp)) { dir.create(smp); }
	setwd(smp);

	# get CNA data
	cna.idx <- which(cnv.data$Sample == smp);
	cnv.filename <- paste0(smp, "_seg_data.tsv");
	names(cnv.filename) <- smp;
	write.table(cnv.data[cna.idx,-1], file = cnv.filename, row.names = F, col.names = T, quote = F, sep = '\t');

	# get SV data
	sv.idx <- which(somatic.svs$sample == smp);
	sv.filename <- paste0(smp, "_sv.bedpe");
	names(sv.filename) <- smp;
	write.table(somatic.svs[sv.idx,], file = sv.filename, row.names = F, col.names = T, quote = F, sep = '\t');

	# get SNV/INDEL data
	snv.filename <- vcf.files[grep(smp, vcf.files)];
	names(snv.filename) <- smp;

	# extract substitution signatures
	res <- vcfToSNVcatalogue(snv.filename, genome.v = arguments$ref_type);
	colnames(res$catalogue) <- smp;
	SNV_catalogues <- res$catalogue;

	plotSubsSignatures(
		signature_data_matrix = SNV_catalogues,
		plot_sum = TRUE,
		output_file = "SNV_catalogues.pdf"
		);

	subs_fit_res <- Fit(
       		catalogues = SNV_catalogues,
        	signatures = COSMIC30_subs_signatures[,sigsToUse],
	        useBootstrap = TRUE,
	        nboot = 100,
	        nparallel = 4
        	);

	# the signature exposures can be found here and correspond to the median
	# of the bootstrapped runs followed by false positive filters. See
	# ?Fit for details
	snv_exp <- subs_fit_res$exposures;

	# calculate indel classification
	indel_class <- vcfToIndelsClassification(snv.filename, smp, genome.v = arguments$ref_type);

	# calculate HRD LOH count
	hrdloh <- ascatToHRDLOH(cnv.data[cna.idx,-1], smp);

	# initialise feature matrix
	input_matrix <- matrix(nrow = 1, ncol = length(hrdetect.features), dimnames = list(smp,hrdetect.features));

	# supply SNV signatures in the sample via the input matrix
	input_matrix[rownames(snv_exp),"SNV3"] <- snv_exp[smp,"Signature3"];
	input_matrix[rownames(snv_exp),"SNV8"] <- snv_exp[smp,"Signature8"];

	# supply del.mh.prop via input matrix
	input_matrix[rownames(snv_exp),"del.mh.prop"] <- indel_class$count_proportion$del.mh.prop;

	# supply hrd loh count via input matrix
	input_matrix[rownames(snv_exp),"hrd"] <- hrdloh[[1]];

	# run the HRDetect pipeline
	hrd.results <- HRDetect_pipeline(
		input_matrix,
		genome.v = arguments$ref_type,
		SV_bedpe_files = sv.filename
		);

	# output SV catalogue and signatures from (Degasperi et al. 2020, Nature Cancer;
	# https://www.nature.com/articles/s43018-020-0027-5)
	save(
		hrd.results,
		file = generate.filename(smp, 'HRDetect_SV_catalogue_and_exposures', 'RData')
		);

	# output HRDetect scores
	output <- as.data.frame(hrd.results$hrdetect_output);
	output$Sample <- smp;
	output <- output[,c("Sample", "intercept", hrdetect.features, "Probability")];

	hrdetect.outputs[[smp]] <- output;

	rm(hrd.results,output,input_matrix,cna.idx,cnv.filename,sv.idx,sv.filename,snv.filename,res,SNV_catalogues,subs_fit_res,snv_exp,indel_class,hrdloh);

	# return to top-level output directory
	setwd(arguments$output);

	}

scores <- do.call(rbind, hrdetect.outputs);

# save to file
write.table(
	scores,
	file = generate.filename(arguments$project, 'HRDetect_scores', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# make a symlink for use with final summary plots
unlink('HRDetect_scores.tsv');
file.symlink(
	generate.filename(arguments$project, 'HRDetect_scores', 'tsv'),
	'HRDetect_scores.tsv'
	);

### PLOTTING #######################################################################################
# grab some parameters
axis.cex <- if (length(all.samples) <= 30) { 1
	} else if (length(all.samples) <= 50) { 0.75
	} else if (length(all.samples) <= 80) { 0.5
	} else { 0 };

# make the plot
create.barplot(
	Probability ~ Sample,
	scores,
	col = 'grey60',
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	ylab.label = expression('Signature Probability'),
	ylab.cex = 1.4,
	ylab.axis.padding = 3,
	xaxis.lab = simplify.ids(scores$Sample),
	xlab.label = NULL,
	xaxis.cex = axis.cex,
	yaxis.cex = 1,
	xaxis.rot = 90,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	abline.h = 0.7,
	abline.lty = 2,
	abline.col = 'red',
	height = 4,
	width = if (length(all.samples) < 10) { 6 } else { 8 },
	resolution = 600,
	filename = generate.filename(arguments$project, 'HRDetect_scores', 'png')
	);

scores <- scores[order(scores$Probability, decreasing = TRUE),];

# if making output for a report
if (!is.null(arguments$report)) {

	tex.file <- paste0(
		arguments$report,
		'/hrdetect_signature_summary.tex'
		);

	# run unlink, in case it exists from a previous run
	unlink(paste0(arguments$report, '/', 'hrdetect_scores.png'));
	file.symlink(
		paste0(arguments$output, '/',
			generate.filename(arguments$project, 'HRDetect_scores', 'png')),
		paste0(arguments$report, '/', 'hrdetect_scores.png')
		);

	### LATEX ##################################################################################
	# write for latex
	write("\\section{Mutation Signatures}", file = tex.file);

	write(
		'HRD signatures were applied to the current dataset by HRDetect, using ENSEMBLE mutations, high confidence SV calls and copy-number profiles.',
		file = tex.file,
		append = TRUE
		);

	# add mutation_signatures plot
	write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
	write("\\begin{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\includegraphics[height=0.3\\textheight]{",
		paste0(arguments$report, '/', 'hrdetect_scores.png'), '}'
		), file = tex.file, append = TRUE);
	write("\\end{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\caption{HRDetect was applied to the current dataset. Samples with signature probabilities $\\geq$ 0.7 are considered HR-deficient.}"
		), file = tex.file, append = TRUE);
	write("\\end{figure}\n", file = tex.file, append = TRUE);

	# add table summary
	caption <- 'Probabilities of HR-deficiency generated by HRDetect.';
	print(
		xtable(
			scores[1:(min(c(nrow(scores),12))),c('Sample','Probability')],
			caption = caption 
			),
		file = tex.file,
		include.rownames = FALSE,
		append = TRUE
		);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('HRDetect','SessionProfile','txt'));
