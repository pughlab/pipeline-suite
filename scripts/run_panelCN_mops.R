## Copy-Number calling of targeted panel DNA-seq data using panelcn.mops
# Jeff Bruce - jeffpbruce@gmail.com
# Script to take in a set of test bams and reference ('normal') bams along
# with a target bed file and output figures and copy-number calls

.libPaths(.libPaths()[!grepl('home', .libPaths())]);

### PREAMBLE #######################################################################################
# load required libraries
library(optparse);
library(panelcn.mops);

####################################
### Get inputs from command line argumets
####################################
option_list <- list(
	make_option(c("-s", "--sample_list"), default = NULL,
		help = "Required; path to tab seperated file with 3 columns (no header): 
		1) Sample ID 
		2) Sample Type (values must be either 'control' or 'tumour')
		3) Path to bam file
		"),
	make_option(c("-b", "--bed_file"), default = NULL,
		help = "Required; standard format bed file (chromsome<tab>start<tab>end) with 
		no header row and probe/gene name in 4th. In order to correctly 
		parse gene and exon names from probe ID, must be in this format: 
		'GeneName.ExonName'"
		),
	make_option(c("-l", "--read_length"), type = "integer", default = 100,
		help = "Read length of input data. [default %default]"
		),
	make_option(c("-o", "--output_directory"), type = "character", default = ".",
		help = "Path to output directory; defaults to current working directory."
		),
	make_option(c("-g", "--genome_build"), type = "character", default = "hg38",
		help = "Genome build; one of hg19 or hg38."
		),
	make_option(c("-p", "--pon"), type = "character", default = NULL,
		help = "Path to panel of normals file."
		),
	make_option(c("-m", "--make_pon"), type = "logical", default = FALSE,
		help = "Does the provided sample_list ONLY contain control samples? [default %default]"
		)
	);

# set arguments
opt <- parse_args(OptionParser(option_list=option_list));

sample_file		<- opt$sample_list;
bed_file		<- opt$bed_file;
output_directory	<- opt$output_directory;
read_length		<- opt$read_length;
ref_type		<- opt$genome_build;
pon			<- opt$pon;
make_pon		<- opt$make_pon;

# move to output directory
setwd(output_directory);

### MAIN ###########################################################################################
# format target regions into count windows
if (make_pon) {

	# add annotations
	library(GenomicRanges);

	tmp_bed <- read.delim(bed_file, header = FALSE);
        colnames(tmp_bed)[1:3] <- c('Chromosome','Start','End');
	has_chr <- grepl('chr', tmp_bed[1,1]);
        target.gr <- reduce(GRanges(tmp_bed));

	if ('hg19' == ref_type) {
		library('TxDb.Hsapiens.UCSC.hg19.knownGene');
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;

		} else if ('hg38' == ref_type) {
		library('TxDb.Hsapiens.UCSC.hg38.knownGene');
		txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;

		} else {
		stop('Unknown reference provided; must be one of hg19 or hg38.');
		}

	gene.gr <- transcriptsBy(txdb, by = 'gene');

	gene.ids <- names(gene.gr);
	gene.info <- select(txdb,
		keys = gene.ids,
		keytype = 'GENEID',
		columns = c('TXNAME','TXCHROM','TXSTART','TXEND','EXONRANK','EXONSTART','EXONEND')
		);

	library(org.Hs.eg.db);
	gene.info$Symbol <- mapIds(org.Hs.eg.db,
		keys = gene.info$GENEID,
		keytype = "ENTREZID",
		column = "SYMBOL",
		multiVals = "first"
		);

	# reduce to gene coordinates (from transcript)
	gene.info <- gene.info[order(gene.info$GENEID, gene.info$TXCHROM, gene.info$TXNAME),];

	gene.data <- gene.info[!duplicated(gene.info$GENEID),c(1,9,6,5,7,8)];
	gene.data <- gene.data[order(gene.data$TXCHROM, gene.data$TXSTART, gene.data$TXEND),];

	exon.data <- merge(gene.data[,1:4],gene.info);
	exon.data <- exon.data[order(exon.data$TXCHROM, exon.data$EXONSTART, exon.data$EXONEND),];

	# add gene annotations
	gene.gr <- makeGRangesFromDataFrame(gene.data,
		seqnames.field = 'TXCHROM',
		start.field = 'TXSTART',
		end.field = 'TXEND',
		keep.extra.columns = TRUE
		);
	exon.gr <- makeGRangesFromDataFrame(exon.data,
		seqnames.field = 'TXCHROM',
		start.field = 'EXONSTART',
		end.field = 'EXONEND'
		);

	gene.overlaps <- as.data.frame(findOverlaps(gene.gr, target.gr));
	exon.overlaps <- as.data.frame(findOverlaps(exon.gr, target.gr));

	# panelcn.mops expects gene 'names' to be in the format
	# SYMBOL.exon.chr.start.end
	# so fill in as much as possible
	annotated_bed <- as.data.frame(target.gr);
	annotated_bed$seqnames <- as.character(annotated_bed$seqnames);
	annotated_bed$V4 <- NA;
	for (i in 1:nrow(annotated_bed)) {
		if (i %in% gene.overlaps$subjectHits) {
			gene.idx <- gene.overlaps[which(gene.overlaps$subjectHits == i),]$queryHits;
			gene <- sort(unique(gene.data[gene.idx,]$Symbol))[1];
			} else { gene <- 'GENE';
			}
		if (i %in% exon.overlaps$subjectHits) {
			exon.idx <- exon.overlaps[which(exon.overlaps$subjectHits == i),]$queryHits;
			exon <- paste0('E',unique(exon.data[exon.idx,]$EXONRANK)[1]);
			} else { exon <- 'other';
			}

		annotated_bed[i,]$V4 <- paste(c(gene, exon, annotated_bed[i,1:3]), collapse = '.');
		}

	new_bed_file <- sub('.bed','_annotated.bed',basename(bed_file));
	write.table(
		annotated_bed[,c('seqnames', 'start', 'end', 'V4')],
		file = new_bed_file,
		row.names = FALSE,
		col.names = FALSE,
		quote = FALSE,
		sep = '\t'
		);

	countWindows <- getWindows(new_bed_file, chr = has_chr);
	if (any(countWindows$gene == 'GENE')) {
		countWindows[which(countWindows$gene == 'GENE'),]$gene <- NA;
		}
	if (any(countWindows$exon == 'other')) {
		countWindows[which(countWindows$exon == 'other'),]$exon <- NA;
		}
	countWindows$name <- gsub('GENE.other','region',countWindows$name);

	write.table(
		countWindows,
		file = 'formatted_countWindows.bed',
		row.names = FALSE,
		col.names = TRUE,
		quote = FALSE,
		sep = '\t'
		);
	
	} else {

	countWindows <- read.delim(bed_file);

	}

## PROCESS SAMPLES ##
# get sample list (ID, type, path)
sample_list <- read.delim(sample_file, header = FALSE);

sampleBams <- sample_list[which(sample_list$V2 == 'tumour'),]$V3;
sampleNames <- sample_list[which(sample_list$V2 == 'tumour'),]$V1;
refBams <- sample_list[which(sample_list$V2 == 'control'),]$V3;

# IF processing control samples ONLY
if (make_pon) {

	# Load control files
	control <- countBamListInGRanges(
		countWindows = countWindows,
		bam.files = refBams,
		read.width = read_length
		);
 
	# save reference for future use
	save(control, file = 'merged_GRanges_count_obj_for_panelcn.RData');

	# ELSE if processing tumour samples ONLY
	} else {

	# load control files
	if (file.exists(pon)) {
		message(paste0("Loading previously processed reference count-set: ", pon));
  		load(pon);
		} else {
		stop('Panel of normals file does not exist; please confirm and try again.');
		}

	# process sample Bam
	if (length(sampleBams) == 1) {
		print(paste("Sample:",sampleNames));
		output.file <- paste0(sampleNames, "_panelcn.mops_results.tsv");
		} else {
		print(paste("Processing cohort:", length(sampleBams), "samples."));
		output.file <- "cohort_panelcn.mops_results.tsv";
		}

	# get sample read counts
	XandCB <- countBamListInGRanges(
		countWindows = countWindows,
		bam.files = sampleBams,
		read.width = read_length
		);

	colnames(elementMetadata(XandCB)) <- sampleNames;
	elementMetadata(XandCB) <- cbind(
		elementMetadata(XandCB),
		elementMetadata(control)
		);

	# run panelCN.mops
	resultlist <- runPanelcnMops(
		XandCB,
		testiv = 1:length(sampleNames),
		countWindows = countWindows,
		selectedGenes = NULL
		);

	# save raw results object
	save(resultlist, file = 'panelcn_mops_resultlist_object.RData');

	# extract/format results table
	results.table <- createResultTable(
		resultlist = resultlist,
		XandCB = XandCB,
		countWindows = countWindows,
		sampleNames = sampleNames
		);

	output.results <- do.call(rbind, results.table);

	# write results to file
	write.table(
		output.results,
		file = output.file,
		row.names = FALSE,
		col.names = TRUE,
		quote = FALSE,
		sep = '\t'
		);

	# output some plots
	library(CopyNumberPlots);
	plot.params <- getDefaultPlotParams(plot.type = 3);		  
	plot.params$leftmargin <- 0.1;
	plot.params$ideogramheight <- 30;
	plot.params$data2height <- 100;

	for (i in 1:length(sampleBams)) {

		gr_ob <- resultlist[[i]]@gr;

		gr_ob$Chromosome <- seqnames(gr_ob);
		gr_ob$baf <- 1;
		gr_ob$lrr <- log2(results.table[[i]]$RC.norm / results.table[[i]]$medRC.norm);
		gr_ob$cn <- as.numeric(gsub("CN","",resultlist[[i]]@integerCopyNumber[,1]));
		gr_ob$loh <- NA;
		gr_ob$segment.value <- gr_ob$cn;
		gr_ob$CopyNumberInteger <- gr_ob$cn;

		cnplot_dat <- loadSNPData(gr_ob);

		# output plots for whole genome and per-chromosome
		chrs <- c('all', as.character(unique(gr_ob$Chromosome)));

		output.file <- paste0(sampleNames[i], "_panelcn.mops_per-chrom_plots.pdf");

		# initiate pdf file
		pdf(output.file, height = 6, width = 12);

		for (chr in chrs) {

			kp <- plotKaryotype(
				ref_type,
				chromosomes = if (chr == 'all') { 'auto' } else { chr },
				plot.type = 3,
				labels.plotter = NULL, 
				main = if (chr == 'all') { '' } else { chr },
				cex = 1,
				plot.params = plot.params
				);

			if (chr == 'all') {
				kpAddChromosomeNames(kp, srt = 90, cex = 1);
				kpAddChromosomeSeparators(kp, lwd = 2, col = "#666666");
				}

			plotLRR(
				kp,
				cnplot_dat,
				ymin = -1.5,
				ymax = 1.5,
				labels = NA,
				line.at.0 = TRUE,
				line.at.0.col = "black",
				points.col = "dodgerblue2",
				points.cex = 0.8,
				add.axis = FALSE,
				data.panel = 1
				);

			kpAxis(kp, ymin = -1.5, ymax = 1.5, col = "gray50", cex = 1, numticks = 5);
			kpAddLabels(
				kp,
				labels = expression('log'['2']*'(Depth Ratio)'),
				srt = 90, 
				cex = 1,
				pos = 3,
				label.margin = 0.07
				);

			plotCopyNumberCallsAsLines(
				kp,
				cn.calls = cnplot_dat,
				col = "red",
				lwd = if (chr == 'all') { 3 } else { 1 },
				data.panel = 2,
				style = if (chr == 'all') { 'segments' } else { 'line' },
				labels = NA,
				ymin = 0,
				ymax = 4,
				add.axis = FALSE,
				);

			kpAxis(kp, ymin = 0, ymax = 4, col = "gray50", cex = 1, numticks = 5, data.panel = 2);
			kpAddLabels(
				kp,
				labels = "Copy-Number Estimate",
				srt = 90, 
				cex = 1,
				pos = 3,
				label.margin = 0.07,
				data.panel = 2
				);
			}

		dev.off();

		}
	}
