### SequenzaSingleSample_optimized.R ###############################################################
# Script to run Sequenza (v2.1) on snps and cnv files from Varscan using optimized gamma parameter
#
# Modified from /cluster/projects/pughlab/src/sequenza_wrapper/SequenzaSingleSample_v2.1.R
# to incorporate reference type

### PREAMBLE #######################################################################################
library(sequenza);
library(optparse);
library(GenomicRanges);

### FUNCTIONS ######################################################################################
# function to extract alternative solutions
alternative.cp.solutions <- function(cp.table) {
	ci <- get.ci(cp.table);
	p.alt <- which(diff(sign(diff(ci$values.ploidy$y))) == -2) + 1;
	get.alt <- function(idx.p, cp.table) {
		idx.c <- which.max(cp.table$lpp[idx.p,]);
		c(
			cellularity = cp.table$cellularity[idx.c],
			ploidy = cp.table$ploidy[idx.p],
			SLPP = cp.table$lpp[idx.p, idx.c]
			);
		}
	res <- lapply(p.alt, FUN = function (x) get.alt(x, cp.table));
	res <- as.data.frame(do.call(rbind, res));
	if (nrow(res) > 0 ) {
		res[order(res$SLPP, decreasing = TRUE), ];
	} else {
		data.frame(
			cellularity = ci$max.cellularity, 
			 ploidy = ci$max.ploidy,
			SLPP = cp.table$lpp[which(cp.table$ploidy == ci$max.ploidy),
			which(cp.table$cellularity == ci$max.cellularity)]
			);
		}
	}

### PARSE OPTIONS ##################################################################################
# For now, test if commands are in original, trailing format, or new opt-parse format
option_list <- list(
	make_option(c("-q", "--seqz_file"), type="character", default=NULL,
		help="seqz.data file produced from Varscan input",
		metavar="character"),
	make_option(c("-o", "--out_dir"), type="character", default=NULL,
		help="output directory [default= %default]", metavar="character"),
	make_option(c("-m", "--min.reads.normal"), type="double", default=10,
		help="min reads to make genotype call [default= %default]", metavar="double"),
	make_option(c("-g", "--gamma"), type="character", default=40,
		help="penalty for each discontinuity in the curve [default= %default]", metavar="character"),
	make_option(c("-k", "--kmin"), type="double", default=5,
		help="minimum number of probes in each segment [default= %default]", metavar="double"),
	make_option(c("-n", "--remove_centromeres"), type="logical", default=TRUE,
		help="remove all snp and copy-number data from centromeric regions", metavar="logical"),
	make_option(c("-p", "--ploidy_priors"), type="character",
		default="/cluster/projects/pughlab/src/sequenza_wrapper/PANCAN_ASCAT_ploidy_prob.Rdata",
		help="Ploidy priors table [default= %default]", metavar="character"),
	make_option(c("-t", "--cancer_type"), type="character", default="all",
		help="Cancer type to define priors [default= %default]", metavar="character"),
	make_option(c("-w", "--window"), type="integer", default=1000000,
		help="window size for plotting [default= %default]", metavar="integer"),
	make_option(c("-r", "--ratio_priority"), type="logical", default=FALSE,
		help="Use only depth ratio for segmentation? [default= %default]", metavar="logical"),
	make_option(c("-a", "--assembly"), type="character", default='hg38', metavar="character",
		help="Reference type (one of hg19, GRCh37, hg38, GRCh38) [default= %default]"),
	make_option(c("-f", "--min_reads_baf"), type="integer", default=1,
		help="Threshold on the depth of the positions included to calculate the average BAF for segment. Set to extreme value (ex. 1x10^6) to exclude BAF from calculations [default= %default]", metavar="integer")
	); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# extract options
seqz.file <- opt$seqz_file;
outdir <- opt$out_dir;

usr_cancer_type <- opt$cancer_type;
purity <- 1;
ud_ploidy <- 2;

min.reads.normal <- opt$min.reads.normal;
min.reads.baf <- opt$min_reads_baf;

assembly <- opt$assembly;

print("Running Sequenza with the following options:");
print(opt);

### MAIN ###########################################################################################
# make a new directory so we don't overwrite old output
new.outdir <- paste0(outdir, "/optimized");

# add another output directory for user defined purity 
output_udp <- paste0(new.outdir,"/output_udp_",purity,"/");
dir.create(output_udp, showWarnings = FALSE, recursive = TRUE);

# find optimal gamma
print("Finding optimal gamma");

tuning_dir <- paste0(outdir,"/output_tuning/");
tuning.files <- list.files(path = tuning_dir, pattern = 'RDS', full.names = TRUE);

for (file in tuning.files) {

	load(file);

	if (file == tuning.files[1]) {
		gamma.frame <- out.data;
		} else {
		gamma.frame <- rbind(gamma.frame, out.data);
		}
	}

# Determine the optimal gamma value using previously run tuning step
gamma.frame <- gamma.frame[order(gamma.frame$gamma),];

n_seg <- gamma.frame$n_seg;
gamma_vect <- gamma.frame$gamma;

# calculate slope for each pair of points
diff_x <- diff(gamma_vect);
diff_y <- diff(n_seg);

f_slopes <- diff_x/diff_y;

# slope closest to -1
opt_slope <- which.min(abs(-1-f_slopes));

opt_gamma <- mean(gamma_vect[opt_slope+1],gamma_vect[opt_slope]);
opt_seg <- n_seg[gamma_vect==opt_gamma];

log_n_seg <- log10(n_seg);
log_opt_seg <- log10(opt_seg);

print(paste0(">>Optimal GAMMA: ", opt_gamma));

## plot gamma values vs. segment values
pdf(
        file = paste0(new.outdir, "/optimal_gamma_selection.pdf"),
        height = 5,
        width = 8
        );

plot(gamma_vect, n_seg/100, pch = 16, xlab = "Gamma Value", ylab = "Number of Segments (hundred)", log = "y", las = 1);

lines(gamma_vect,n_seg/100,col="red",lty=1,lwd=2);

arrows(
        x0 = opt_gamma,
        y0 = (opt_seg+(max(n_seg)*0.01))/100,
        x1 = opt_gamma,
        y1 = opt_seg/100,
        col = "grey",
        lwd = 2,
        length = 0.1
        );

text(
        x = opt_gamma,
        y = (opt_seg+(max(n_seg)*0.01))/100,
        labels = paste0("Optimal Gamma = ", opt_gamma, "\n", "Segments N = ", opt_seg), pos = 3
        );

dev.off();

### SEQUENZA #######################################################################################
filestem <- gsub('.seqz','',basename(seqz.file));

print("Running sequenza.extract with optimal gamma");

data <- sequenza.extract(
	seqz.file,
	breaks.method = 'het', 
	window = opt$window,
	gamma = opt_gamma,
	min.reads.normal = min.reads.normal,
	min.reads.baf = min.reads.baf,
	kmin = opt$kmin,
	assembly = assembly
	);

# feed in prior probabilities based on cancer type
print("Feeding in prior probabilities for ploidy");

load(opt$ploidy_priors);

priors <- subset(ploidy_table, cancer_type==toupper(usr_cancer_type))[c("CN","value")];

CP <- sequenza.fit(
	data,
	priors.table = priors,
	ratio.priority = opt$ratio_priority
	);

# output series of results files (images and text files)
sequenza.results(
	sequenza.extract = data,
	cp.table = CP,
	sample.id = filestem,
	out.dir = new.outdir
	);

### Extract all alternative solutions###
print("Getting alternative purity solutions");

alt <- alternative.cp.solutions(CP);

purities <- alt$cellularity;
ploidies <- alt$ploidy;

## output series of results files (images and text files) for the alternative purity solutions
## index starts at 2 because the first solution is already output as the top hit.
if (length(purities) > 1) { ## bug fix to not run alt solutions if there are no alt solutions.
	for (i in 2:length(purities)) {
		print(paste("Creating new directories and printing solution:", i));
		output_alt <- paste0(new.outdir,"/sol",i,"_",purities[i],"/");
		dir.create(output_udp, showWarnings = FALSE, recursive = TRUE);
		sequenza.results(
			sequenza.extract = data, 
			cp.table = CP,
			sample.id = filestem,
			out.dir = output_alt,
			cellularity = purities[i],
			ploidy = ploidies[i]
			);
		}
	}

## output results using forced purity parameter
## If you input purity alone, the ploidy result will still be selected from the previous
## fit, unless you re-fit with the user-defined purity as one of the (using the command below)
## options.	I have opted to instead [force the user to instead provide both the purity and ploidy
## values.
sequenza.results(
	sequenza.extract = data,
	cp.table = CP, 
	sample.id = filestem,
	out.dir = output_udp,
	cellularity = purity,
	ploidy = ud_ploidy
	);
