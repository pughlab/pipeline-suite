# file:   ichorCNA.R
# authors: Gavin Ha, Ph.D.
#          Fred Hutch
# contact: <gha@fredhutch.org>
#
#         Justin Rhoades
#          Broad Institute
# contact: <rhoades@broadinstitute.org>

# ichorCNA: https://github.com/broadinstitute/ichorCNA
# date:   July 24, 2019
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

### Modified by sprokopec for use with pipeline-suite

### PREAMBLE ######################################################################################
library(optparse);

option_list <- list(
	make_option("--id", type = "character", help = "Patient ID. Required."),
	make_option("--WIG", type = "character", help = "Path to tumor WIG file. Required."),
	make_option("--NORMWIG", type = "character", default = NULL, help = "Path to normal WIG file. Default: [%default]"),
	make_option("--normal", type = "character", default = "c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)", help = "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [%default]"),
	make_option("--scStates", type = "character", default = "c(1,3)", help = "Subclonal states to consider. Default: [%default]"),
	make_option("--ploidy", type = "character", default = "c(2,3)", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
	make_option("--maxCN", type = "numeric", default = 5, help = "Total clonal CN states. Default: [%default]"),
	make_option("--estimateNormal", type = "logical", default = TRUE, help = "Estimate normal. Default: [%default]"),
	make_option("--estimateScPrevalence", type = "logical", default = TRUE, help = "Estimate subclonal prevalence. Default: [%default]"),
	make_option("--estimatePloidy", type = "logical", default = TRUE, help = "Estimate tumour ploidy. Default: [%default]"),
	make_option("--chrTrain", type = "character", default = "c(1:22)", help = "Specify chromosomes to estimate params. Default: [%default]"),
	make_option("--chrNormalize", type = "character", default = "c(1:22,'X')", help = "Specify chromosomes to normalize GC/mappability biases. Default: [%default]"),
	make_option("--chrs", type = "character", default = "c(1:22,'X')", help = "Specify chromosomes to analyze. Default: [%default]"),
	make_option("--genomeBuild", type = "character", default = "hg38", help="Geome build. Default: [%default]"),
	make_option("--genomeStyle", type = "character", default = "UCSC", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option("--includeHOMD", type = "logical", default = FALSE, help = "If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [%default]"),
	make_option("--txnE", type = "numeric", default = 0.9999, help = "Self-transition probability. Increase to decrease number of segments. Default: [%default]"),
	make_option("--txnStrength", type = "numeric", default = 10000, help = "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [%default]"),
	make_option("--plotYLim", type = "character", default = "c(-2,2)", help = "ylim to use for chromosome plots. Default: [%default]"),
	make_option("--outDir", type = "character", default = "./", help = "Output Directory. Default: [%default]"),
	make_option("--target_bed", type = "character", default = NULL, help = "Path to target bed file (ie, for whole-exome sequencing"),
	make_option("--make_pon", type = "logical", default = FALSE, help = "Should a panel of normals be created? Requires --normal_list (input) and --normalPanel (output)"),
	make_option("--normal_list", type = "character", default = NULL, help = "list of normal WIG files for PoN"),
	make_option("--normalPanel", type="character", default = NULL, help = "Median corrected depth from panel of normals. Default: [%default]")
	);

parseobj <- OptionParser(option_list = option_list);
opt <- parse_args(parseobj);
options(scipen = 0, stringsAsFactors = FALSE);

# input files
patientID	<- opt$id;
tumour_file	<- opt$WIG;
normal_file	<- opt$NORMWIG;

# is this run for PoN creation?
create_pon	<- opt$make_pon;
normal_list	<- opt$normal_list;

# genome arguments
genomeBuild	<- opt$genomeBuild;
genomeStyle	<- opt$genomeStyle;
chrs		<- as.character(eval(parse(text = opt$chrs)));
chrTrain	<- as.character(eval(parse(text = opt$chrTrain)));
chrNormalize	<- as.character(eval(parse(text = opt$chrNormalize)));

# goals
normal		<- eval(parse(text = opt$normal));
scStates	<- eval(parse(text = opt$scStates));
ploidy		<- eval(parse(text = opt$ploidy));
estimateNormal	<- opt$estimateNormal;
estimatePloidy	<- opt$estimatePloidy;
estimateScPrevalence	<- opt$estimateScPrevalence;

# parameters
maxCN		<- opt$maxCN;
txnE		<- opt$txnE;
txnStrength	<- opt$txnStrength;
includeHOMD	<- as.logical(opt$includeHOMD);

# output parameters
outDir		<- sub('/$','',opt$outDir);
plotYLim	<- eval(parse(text = opt$plotYLim));
outImage	<- paste0(outDir,"/", patientID,".RData");
plotFileType	<- 'pdf';

# other/ignored/default only
minMapScore	<- 0.9;
flankLength	<- 1e5;
coverage	<- NULL;
lambda		<- NULL;
lambdaScaleHyperParam	<- 3;
maxFracCNASubclone	<- 0.7;
maxFracGenomeSubclone	<- 0.5;
minSegmentBins		<- 50;
altFracThreshold	<- 0.05;
normalizeMaleX		<- TRUE;
minTumFracToCorrect	<- 0.1;
fracReadsInChrYForMale	<- 0.001;
chrXMedianForMale	<- -0.1;
maleChrXLogRThres	<- -0.80;
ponMethod		<- 'median';
gender		<- NULL;

# load required libraries
library(HMMcopy);
library(GenomicRanges);
library(GenomeInfoDb);
library(ichorCNA);
options(stringsAsFactors = FALSE);
options(bitmapType = 'cairo');

# find built-in reference files
if ('hg38' == opt$genomeBuild) {
	gcWig <- system.file('extdata', 'gc_hg38_1000kb.wig', package = 'ichorCNA');
	mapWig <- system.file('extdata', 'map_hg38_1000kb.wig', package = 'ichorCNA');
	normal_panel <- system.file('extdata',
		'HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds', package = 'ichorCNA');
	centromere <- system.file('extdata',
		'GRCh38.GCA_000001405.2_centromere_acen.txt', package = 'ichorCNA');
	} else if ('hg19' == opt$genomeBuild) {
	gcWig <- system.file('extdata', 'gc_hg19_1000kb.wig', package = 'ichorCNA');
	mapWig <- system.file('extdata', 'map_hg19_1000kb.wig', package = 'ichorCNA');
	normal_panel <- system.file('extdata',
		'HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds', package = 'ichorCNA');
	centromere <- system.file('extdata',
		'GRCh37.p13_centromere_UCSC-gapTable.txt', package = 'ichorCNA');
	} else {
	stop('Unrecognized genomeBuild provided. Must be one of hg38 or hg19.');
	}

if (!is.null(opt$normalPanel)) {
	normal_panel <- opt$normalPanel;
	}

# based purely on built-in gc/map wigs, hg19 = NCBI and hg39 = UCSC
seqlevelsStyle(chrs)		<- genomeStyle;
seqlevelsStyle(chrNormalize)	<- genomeStyle;
seqlevelsStyle(chrTrain)	<- genomeStyle;

# load seqinfo 
seqinfo <- getSeqInfo(genomeBuild, genomeStyle);

# get centromere data
centromere <- read.delim(centromere);
seqlevelsStyle(centromere$Chr) <- genomeStyle;

# if target bed is provided
targetSequences <- NULL;
if (!is.null(opt$target_bed)) {
	targetSequences <- read.delim(opt$target_bed, header = FALSE);
	}

### CREATE PANEL OF NORMALS #######################################################################
if (create_pon) {

	# indicate output stem
	outfile <- gsub('_median.rds','',normal_panel);

	message("Reading GC and mappability files");
	gc <- wigToGRanges(gcWig);
	map <- wigToGRanges(mapWig);

	files <- read.delim(normal_list, header = FALSE, stringsAsFactors = FALSE)[,1];
	normalGR <- NULL;

	for (normalFile in files) {

		norm.id <- gsub('.wig','',basename(normalFile));
		message("Loading normal file:", normalFile);

		normal_reads <- wigToGRanges(normalFile);

		counts <- loadReadCountsFromWig(normal_reads, chrs = chrs, gc = gc, map = map, genomeStyle = genomeStyle,
			centromere = centromere, targetedSequences = targetSequences, chrNormalize = chrNormalize);

		gender.normal <- counts$gender;

		message("Correcting ", norm.id);
		if (is.null(normalGR)) {
			normalGR <- counts$counts;
			values(normalGR) <- values(normalGR)$copy;
			colnames(values(normalGR)) <- norm.id;
			} else {
			values(normalGR)[[norm.id]] <- counts$counts$copy;
			}

		chrXMedian <- gender.normal$chrXMedian;
		chrXStr	<- grep("X", chrs, value = TRUE);
		chrXInd <- as.character(seqnames(normalGR)) == chrXStr;

		values(normalGR)[[norm.id]][chrXInd] <- values(normalGR)[[norm.id]][chrXInd] - chrXMedian;

		}

	mat <- values(normalGR);
	medianVal <- apply(mat, 1, median, na.rm = TRUE);
	values(normalGR)[['Median']] <- medianVal;

	write.table(
		as.data.frame(normalGR[,'Median']),
		file = paste0(outfile, '_median.txt'),
		row.names = FALSE,
		col.names = TRUE,
		quote = FALSE,
		sep = '\t'
		);

	saveRDS(normalGR, file = paste0(outfile, '_median.rds'));
	q();
	}

### MAIN ##########################################################################################
## LOAD IN WIG FILES ##
tumour_copy <- list();
id <- patientID;

## create output directories for each sample ##
dir.create(paste0(outDir, "/", id, "/"), recursive = TRUE);

### LOAD TUMOUR AND NORMAL FILES ###
message("Loading tumour file:", patientID);
tumour_reads <- wigToGRanges(tumour_file);

## LOAD GC/MAP WIG FILES ###
# find the bin size and load corresponding wig files #
binSize <- as.data.frame(tumour_reads[1,])$width;

message("Reading GC and mappability files");
gc <- wigToGRanges(gcWig);
map <- wigToGRanges(mapWig);

message("Correcting Tumour");
counts <- loadReadCountsFromWig(tumour_reads, chrs = chrs, gc = gc, map = map, 
		centromere = centromere, flankLength = flankLength, 
		targetedSequences = targetSequences, chrXMedianForMale = chrXMedianForMale,
		genomeStyle = genomeStyle, fracReadsInChrYForMale = fracReadsInChrYForMale,
		chrNormalize = chrNormalize, mapScoreThres = minMapScore);

tumour_copy[[id]] <- counts$counts;
gender <- counts$gender;

## load in normal file if provided 
if (!is.null(normal_file) && normal_file != "None" && normal_file != "NULL") {
	message("Loading normal file:", normal_file);
	normal_reads <- wigToGRanges(normal_file);

	message("Correcting Normal");
	counts <- loadReadCountsFromWig(normal_reads, chrs = chrs, gc = gc, map = map, 
		centromere = centromere, flankLength = flankLength, targetedSequences = targetSequences,
		genomeStyle = genomeStyle, chrNormalize = chrNormalize, mapScoreThres = minMapScore);

	normal_copy <- counts$counts;
	gender.normal <- counts$gender;
	} else{
	normal_copy <- NULL;
	}

### DETERMINE GENDER ###
## if normal file not given, use chrY, else use chrX
message("Determining gender...", appendLF = FALSE);
gender.mismatch <- FALSE;

if (!is.null(normal_copy)) {
	if (gender$gender != gender.normal$gender) {
		# check if normal is same gender as tumour
		gender.mismatch <- TRUE;
		}
	}

message("Gender: ", gender$gender);

### there is a bug in normalizeByPanelOrMatchedNormal for male samples
# if tumour_copy has different intervals than the normal_panel, both are filtered to
# overlapping regions and then chrX normalization is applied (but the chrX indices are created
# on the unfiltered tumour_copy. To fix this, we will filter first.
chrs <- setdiff(chrs, 'chrY');
panel <- readRDS(normal_panel);
seqlevelsStyle(panel) <- genomeStyle;

hits <- findOverlaps(tumour_copy[[id]], panel, type="equal");
tumour_copy[[id]] <- tumour_copy[[id]][queryHits(hits),];

## NORMALIZE GENOME-WIDE BY MATCHED NORMAL OR NORMAL PANEL (MEDIAN) ##
tumour_copy[[id]] <- normalizeByPanelOrMatchedNormal(tumour_copy[[id]], chrs = chrs, 
	normal_panel = normal_panel, normal_copy = normal_copy, 
	gender = gender$gender, normalizeMaleX = normalizeMaleX);

### OUTPUT FILE ###
### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
outMat <- as.data.frame(tumour_copy[[id]]);
outMat <- outMat[,c("seqnames","start","end","copy")];
colnames(outMat) <- c("chr","start","end","log2_TNratio_corrected");
outFile <- paste0(outDir,"/",id,".correctedDepth.txt");
message(paste("Outputting to:", outFile));

write.table(
	outMat,
	file = outFile,
	row.names = FALSE,
	col.names = TRUE,
	quote = FALSE,
	sep = "\t"
	);

# get indices
chrInd <- as.character(seqnames(tumour_copy[[id]])) %in% chrTrain;

## get positions that are valid
valid <- tumour_copy[[1]]$valid;
if (length(tumour_copy) >= 2) {
	for (i in 2:length(tumour_copy)) {
		valid <- valid & tumour_copy[[i]]$valid;
		}
	}

### RUN HMM ###
## store the results for different normal and ploidy solutions ##
ptmTotalSolutions <- proc.time(); # start total timer
results <- list();
loglik <- as.data.frame(matrix(NA, nrow = length(normal) * length(ploidy), ncol = 7, 
                 dimnames = list(c(), c("init", "n_est", "phi_est", "BIC", 
		"Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))));
counter <- 1;
compNames <- rep(NA, nrow(loglik));
mainName <- rep(NA, length(normal) * length(ploidy));

#### restart for purity and ploidy values ####
for (n in normal) {
	for (p in ploidy) {
		if (n == 0.95 & p != 2) { next; }

		# NEED TO EXCLUDE CHR X #
		logR <- as.data.frame(lapply(tumour_copy, function(x) { x$copy }));
		param <- getDefaultParameters(logR[valid & chrInd, , drop=F], maxCN = maxCN,
				includeHOMD = includeHOMD, ct.sc = scStates, ploidy = floor(p),
				e = txnE, e.same = 50, strength = txnStrength);

		param$phi_0 <- p;
		param$n_0 <- n;

		############################################
		######## CUSTOM PARAMETER SETTINGS #########
		############################################
		# 0.1x cfDNA #
		if (is.null(lambda)) {
			logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2);
			param$lambda <- rep(logR.var, length(param$ct));
			param$lambda[param$ct %in% c(2)] <- logR.var;
			param$lambda[param$ct %in% c(1,3)] <- logR.var; 
			param$lambda[param$ct >= 4] <- logR.var / 5;
			param$lambda[param$ct == max(param$ct)] <- logR.var / 15;
			param$lambda[param$ct.sc.status] <- logR.var / 10;
			} else {
			param$lambda[param$ct %in% c(2)] <- lambda[2];
			param$lambda[param$ct %in% c(1)] <- lambda[1];
			param$lambda[param$ct %in% c(3)] <- lambda[3];
			param$lambda[param$ct >= 4] <- lambda[4];
			param$lambda[param$ct == max(param$ct)] <- lambda[2] / 15;
			param$lambda[param$ct.sc.status] <- lambda[2] / 10;
			}

		param$alphaLambda <- rep(lambdaScaleHyperParam, length(param$ct));

		#############################################
		################ RUN HMM ####################
		#############################################
		hmmResults.cor <- HMMsegment(tumour_copy, valid, dataType = "copy", 
				param = param, chrTrain = chrTrain, maxiter = 50,
				estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
				estimateSubclone = estimateScPrevalence, verbose = TRUE);

		iter <- hmmResults.cor$results$iter;
		id <- names(hmmResults.cor$cna)[1];

		## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
		## check if there is an altered segment that has at least a minimum # of bins
		segsS <- hmmResults.cor$results$segs[[1]];
		segsS <- segsS[segsS$chr %in% chrTrain, ];
		segAltInd <- which(segsS$event != "NEUT");
		maxBinLength <- -Inf;

		if (sum(segAltInd) > 0) {
			maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1);
			maxSegRD <- GRanges(
				seqnames = segsS$chr[segAltInd[maxInd]],
				ranges = IRanges(
					start = segsS$start[segAltInd[maxInd]],
					end = segsS$end[segAltInd[maxInd]]
					)
				);

			hits <- findOverlaps(query = maxSegRD, subject = tumour_copy[[id]][valid, ]);
			maxBinLength <- length(subjectHits(hits));
			}

		## check if there are proportion of total bins altered 
		# if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
		cnaS <- hmmResults.cor$cna[[1]];
		altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT";
		altFrac <- sum(!altInd, na.rm = TRUE) / length(altInd);

		if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)) {
			hmmResults.cor$results$n[1, iter] <- 1.0;
			}

		# correct integer copy number based on estimated purity and ploidy
		correctedResults <- correctIntegerCN(
			cn = hmmResults.cor$cna[[1]],
			segs = hmmResults.cor$results$segs[[1]], 
			purity = 1 - hmmResults.cor$results$n[1, iter],
			ploidy = hmmResults.cor$results$phi[1, iter],
			cellPrev = 1 - hmmResults.cor$results$sp[1, iter], 
			maxCNtoCorrect.autosomes = maxCN,
			maxCNtoCorrect.X = maxCN,
			minPurityToCorrect = minTumFracToCorrect, 
			gender = gender$gender,
			chrs = chrs,
			correctHOMD = includeHOMD
			);

		hmmResults.cor$results$segs[[1]] <- correctedResults$segs;
		hmmResults.cor$cna[[1]] <- correctedResults$cn;

		## plot solution ##
		outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_", "n", n, "-p", p);
		mainName[counter] <- paste0(id, ", n: ", n, ", p: ", p, ", log likelihood: ",
			signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4));

		plotGWSolution(hmmResults.cor, s = 1, outPlotFile = outPlotFile, plotFileType = plotFileType, 
			logR.column = "logR", call.column = "Corrected_Call", plotYLim = plotYLim,
			estimateScPrevalence = estimateScPrevalence, seqinfo = seqinfo, main = mainName[counter]);

		iter <- hmmResults.cor$results$iter;
		results[[counter]] <- hmmResults.cor;
		loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4);
		subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }));
		fracGenomeSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }));
		fracAltSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna,
			function(x){ sum(x$copy.number != 2) }));
		fracAltSub <- lapply(fracAltSub, function(x){if (is.na(x)) { 0 } else { x }} );
		loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub, digits = 2),
			collapse=",");
		loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSub), digits = 2),
			collapse=",");
		loglik[counter, "init"] <- paste0("n", n, "-p", p);
		loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2),
			collapse = ",");
		loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4),
			collapse = ",");

		counter <- counter + 1;
		}
	}

## get total time for all solutions ##
elapsedTimeSolutions <- proc.time() - ptmTotalSolutions;
message("Total ULP-WGS HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.");

### SAVE R IMAGE ###
save.image(outImage);

### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###
loglik <- loglik[!is.na(loglik$init), ];

if (estimateScPrevalence) {
	## sort but excluding solutions with too large % subclonal 
	fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone & 
			loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone);

	if (length(fracInd) > 0) { ## if there is a solution satisfying % subclonal
		ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing = TRUE)];
		} else { # otherwise just take largest likelihood
		ind <- order(as.numeric(loglik[, "loglik"]), decreasing = TRUE);
		}

	} else { #sort by likelihood only
	ind <- order(as.numeric(loglik[, "loglik"]), decreasing = TRUE);
	}

# new loop by order of solutions (ind)
outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_all_sols");

for (i in 1:length(ind)) {
	hmmResults.cor <- results[[ind[i]]];
	turnDevOff <- FALSE;
	turnDevOn <- FALSE;
	if (i == 1) { turnDevOn <- TRUE; }
	if (i == length(ind)) { turnDevOff <- TRUE; }

	plotGWSolution(hmmResults.cor, s = 1, outPlotFile = outPlotFile, plotFileType = "pdf", 
		logR.column = "logR", call.column = "Corrected_Call", seqinfo = seqinfo,
		plotYLim = plotYLim, estimateScPrevalence = estimateScPrevalence, 
		turnDevOn = turnDevOn, turnDevOff = turnDevOff, main = mainName[ind[i]]);
	}

hmmResults.cor <- results[[ind[1]]];
hmmResults.cor$results$loglik <- as.data.frame(loglik);
hmmResults.cor$results$gender <- gender$gender;
hmmResults.cor$results$chrYCov <- gender$chrYCovRatio;
hmmResults.cor$results$chrXMedian <- gender$chrXMedian;
hmmResults.cor$results$coverage <- coverage;

outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
	results = hmmResults.cor$results, patientID = patientID, outDir = outDir);

outFile <- paste0(outDir, "/", patientID, ".params.txt");
outputParametersToFile(hmmResults.cor, file = outFile);

## plot solutions for all samples 
plotSolutions(hmmResults.cor, tumour_copy, chrs, outDir, numSamples = 1,
              logR.column = "logR", call.column = "Corrected_Call",
              plotFileType = plotFileType, plotYLim = plotYLim, seqinfo = seqinfo,
              estimateScPrevalence = estimateScPrevalence, maxCN = maxCN);

# write for pughlab-pipeline-suite
final.iteration <- hmmResults.cor$results$iter;
estimated.tumour.fraction <- signif(1 - hmmResults.cor$results$n[,final.iteration], digits = 4);
estimated.tumour.ploidy <- signif(hmmResults.cor$results$phi[,final.iteration], digits = 4);
estimated.sex <- hmmResults.cor$results$gender;
fraction.genome.subclonal <- sum(hmmResults.cor$cna[[1]]$subclone.status)/nrow(hmmResults.cor$cna[[1]]);
fraction.cna.subclonal <- sum(hmmResults.cor$cna[[1]]$subclone.status)/sum(
	hmmResults.cor$cna[[1]]$copy.num != 2
	);
estimated.subclone.fraction <- if (fraction.genome.subclonal == 0) { NA; } else {
	signif(1 - hmmResults.cor$results$sp[,final.iteration], digits = 4);
	}

final.output <- data.frame(
	ID = patientID,
	Sex = estimated.sex,
	Tumour.Fraction = estimated.tumour.fraction,
	Ploidy = estimated.tumour.ploidy,
	Estimated.Subclone.Fraction = estimated.subclone.fraction,
	Fraction.Genome.Subclonal = fraction.genome.subclonal,
	Fraction.CNA.Subclonal = fraction.cna.subclonal
	);

write.table(
	final.output,
	file = paste0(outDir, '/', patientID, '_final_metrics.txt'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);
