#############################
# script to run             #
# sequenza on snps and cnv  #
# files from Varscan        #
# including a step to run   #
# many possible gamma values#
# and choose the value that #
# optimizes signal/noise    #
############################

## SequnzaSingleSample_TuneB.R
## This script reads in a previously generated sequenza input
## and runs sequenza_extract using user-defined gamma values.
## This script is meant to be run in parrallel using many gamma
## values in order to find the optimal value

library(sequenza);
library(optparse);
library(GenomicRanges);

#############
# FUNCTIONS #
#############
plotRawGenome <- function(sequenza.extract, cellularity, ploidy, CNt.max = 7, main = "",
	seg_col = adjustcolor("red", alpha.f = 0.75),
	dot_col = adjustcolor("black", alpha.f = 0.5), add = FALSE,
	error_col = adjustcolor("lightblue", alpha.f = 1), ylim = c(0,2.5),
	...) {

	max.end <- sapply(sequenza.extract$ratio, FUN = function(x) max(x$end, na.rm = T));
	max.end <- c(0, cumsum(as.numeric(max.end)));
	chrs <- names(sequenza.extract$ratio);
	coords.names <- (max.end + c(diff(max.end)/2,0))[1:length(chrs)];
	new.coords <- function(win.list, max.end) {
		lapply(1:length(win.list), FUN = function(x) {
			y <- win.list[[x]]
			y$start <- y$start + max.end[x]
			y$end <- y$end + max.end[x]
			y
		})}

	new.coords.segs <- function(segs, max.end) {
		lapply(1:length(segs), FUN = function(x) {
			y <- segs[[x]]
			y$start.pos <- y$start.pos + max.end[x]
			y$end.pos <- y$end.pos + max.end[x]
			y
		})}

	ratio.new <- new.coords(sequenza.extract$ratio,max.end);
	BAF.new   <- new.coords(sequenza.extract$BAF,max.end);
	segs.new  <- do.call(rbind, new.coords.segs(sequenza.extract$segments,max.end));
	avg.depth.ratio <- 1;

	par(mar = c(1, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(2,1), ...);

	if (!add) {
		plot(x = c(min(max.end), max(max.end)), y = c(0,0.5), main = main, xlab = NA,
			ylab = "B allele frequency", type = "n", las = 1, xaxs = "i",
			yaxs = "i", xaxt = "n");
		}

	plotWindows(seqz.window = do.call(rbind, BAF.new), q.bg = error_col, m.col = dot_col, add = T);
	segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos,
		y0 = (segs.new$Bf), y1 = (segs.new$Bf), col = seg_col, lwd = 4, lend = 1);
	abline(v = max.end, lty = 1);

	if (!add) {
		plot(x = c(min(max.end), max(max.end)), y = ylim, main = "", xlab = NA,
			ylab = "Depth ratio", type = "n", las = 1, xaxs = "i", yaxs = "i", xaxt = "n");
		}

	plotWindows(seqz.window = do.call(rbind, ratio.new), q.bg = error_col, m.col = dot_col, add = T);
	segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos,
		y0 = (segs.new$depth.ratio), y1 = (segs.new$depth.ratio), col = seg_col, lwd = 4, lend = 1);

	if (!missing(ploidy) & !missing(cellularity)) {
		types <- types.matrix(CNt.min = 0, CNt.max = CNt.max, CNn = 2);
		depth.ratios <- model.points(
			cellularity = cellularity,
			ploidy = ploidy,
			avg.depth.ratio = avg.depth.ratio,
			types = types)[, "depth.ratio"];

		depth.ratios <- unique(data.frame(CNt = types$CNt, ratio = depth.ratios));
		abline(h = depth.ratios$ratio, lty = 2);
		axis(labels = as.character(depth.ratios$CNt), side = 4, line = 0, las = 1,
			at = depth.ratios$ratio);
		mtext(text = "Copy number", side = 4, line = 2, cex = par("cex.lab")*par("cex"));
		}

	abline(v = max.end, lty = 1);
	axis(labels = chrs, at = coords.names, side = 1, cex.axis = 1);
	}

#################
# PARSE OPTIONS #
#################

## For now, test if commands are in original, trailing format, or new opt-parse format
option_list <- list(
	make_option(c("-q", "--seqz_file"), type="character", default=NULL, 
		help="seqz.data file produced from Varscan input", 
		metavar="character"),
	make_option(c("-o", "--out_dir"), type="character", default=NULL,
		help="output directory [default= %default]", metavar="character"),
	make_option(c("-a", "--assembly"), type="character", default="hg38",
		help="Assembly the data were aligned to [default= %default]", metavar="character"),
	make_option(c("-g", "--gamma"), type="double", default=40,
		help="penalty for each discontinuity in the curve [default= %default]", metavar="double"),
	make_option(c("-k", "--kmin"), type="double", default=5,
		help="minimum number of probes in each segment [default= %default]", metavar="double"),
	make_option(c("-n", "--kmin.pcf"), type="double", default=5,
		help="minimum number of probes in each segment. kmin.pcf effective only when breaks.method is set to full [default= %default]", metavar="double"),
	make_option(c("-m", "--min_reads_normal"), type="double", default=10,
		help="min reads to make genotype call [default= %default]", metavar="double"),
	make_option(c("-w", "--window"), type="integer", default=1000000,
		help="window size for plotting [default= %default]", metavar="integer"),
	make_option(c("-r", "--ratio_priority"), type="logical", default=FALSE,
		help="Use only depth ratio for segmentation? [default= %default]", metavar="logical"),
	make_option(c("-f", "--min_reads_baf"), type="integer", default=1,
		help="Threshold on the depth of the positions included to calculate the average BAF for segment. Set to extreme value (ex. 1x10^6) to exclude BAF from calculations [default= %default]", metavar="integer")
	); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

seqz.file <- opt$seqz_file;
outdir <- opt$out_dir;

gamma <- opt$gamma;

min.reads.normal <- opt$min_reads_normal;
min.reads.baf <- opt$min_reads_baf;
window <- opt$window;
ratio_priority <- opt$ratio_priority;
kmin <- opt$kmin;
assembly <- opt$assembly;

########
# MAIN #
########

filestem <- gsub('.seqz','',basename(seqz.file));
if (grepl('.gz$', basename(seqz.file))) {
	filestem <- gsub('.seqz.gz','',basename(seqz.file));
	}

output_tune <- paste0(outdir,"/output_tuning/");
output_html <- paste0(outdir,"/output_tuning/html_selector/");

dir.create(output_html, showWarnings = FALSE, recursive=TRUE);

# test this gamma
dir.create(
	paste0(output_html, "/", "gamma_", gamma, "_chr_plots", "/"),
	showWarnings = FALSE, recursive = TRUE
	);

print(paste0("Reading sequenza.extract with gamma: ", gamma));

data <- sequenza.extract(
	seqz.file,
	breaks.method = 'het',
	window = window,
	gamma = gamma,
	min.reads.normal = min.reads.normal,
	min.reads.baf = min.reads.baf,
	kmin = kmin,
	assembly = assembly
	);

n_seg <- length(unlist(data$segments))/9;

out.data <- data.frame(
	n_seg = n_seg,
	gamma = gamma
	);

save(
	data,
	out.data,
	file = paste0(output_tune, "gamma_", gamma, ".RDS")
	);

# output plots for manual review
png(
	filename = paste0(output_html, "genome_view_gamma_", gamma, ".png"),
	height = 5,
	width = 15,
	units = "in",
	res = 300
	);

plotRawGenome(data);

dev.off();

# output all chromosomes as well
chroms <- data$chromosomes;

for (chr in chroms) {

	png(
		filename = paste0(output_html, "/gamma_", gamma, '_chr_plots/', chr, ".png"),
		height = 5,
		width = 8,
		units = "in",
		res = 300
		);

	chromosome.view(
		baf.windows = data$BAF[[chr]],
		ratio.windows = data$ratio[[chr]],
		segments = data$segments[[chr]],
		min.N.ratio = 1
		);

	dev.off();
	}
