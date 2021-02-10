### SequenzaSingleSample_v2.1.R ####################################################################
# Script to run Sequenza (v2.1) on snps and cnv files from Varscan
#
# Modified from /cluster/projects/pughlab/src/sequenza_wrapper/SequenzaSingleSample_v2.1.R
# to incorporate reference type

### PREAMBLE #######################################################################################
library(sequenza);
library(optparse);
library(GenomicRanges);

### FUNCTIONS ######################################################################################
# function to find overlap among bed files
BedBedOverlap <- function(bed.frame1, bed.frame2, pad = 0) {

	colnames(bed.frame1)[1:3] <- c("Chromsome","Start_Position","End_Position");
	bed.frame1$Start_Position <- as.numeric(as.character(bed.frame1$Start_Position));
	bed.frame1$End_Position <- as.numeric(as.character(bed.frame1$End_Position));
	bed.frame1$Chromsome <- gsub("chr","",bed.frame1$Chromsome);

	colnames(bed.frame2)[1:3] <- c("Chromsome","Start_Position","End_Position");
	bed.frame2$Start_Position <- as.numeric(as.character(bed.frame2$Start_Position));
	bed.frame2$End_Position <- as.numeric(as.character(bed.frame2$End_Position));
	bed.frame2$Chromsome <- gsub("chr","",bed.frame2$Chromsome);

	## add padding to bed if desired
	bed.frame1$Start_Position <- bed.frame1$Start_Position - pad;
	bed.frame1$End_Position <- bed.frame1$End_Position + pad;

	bed.frame2$Start_Position <- bed.frame2$Start_Position - pad;
	bed.frame2$End_Position <- bed.frame2$End_Position + pad;

	bed_ranges1 <- with(
		bed.frame1,
		GRanges(Chromsome, IRanges(start=Start_Position, end=End_Position))
		);

	bed_ranges2 <- with(
		bed.frame2,
		GRanges(Chromsome, IRanges(start=Start_Position, end=End_Position))
		);

	hits <- overlapsAny(bed_ranges1,bed_ranges2);
	
	return(hits);
	}

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
	make_option(c("-s", "--snp_file"), type="character", default=NULL, 
		help="varscan snp file name", metavar="character"),
	make_option(c("-c", "--cnv_file"), type="character", default=NULL,
		help="varscan copy number file name [default= %default]", metavar="character"),
	make_option(c("-o", "--out_dir"), type="character", default=NULL,
		help="output directory [default= %default]", metavar="character"),
	make_option(c("-p", "--purity"), type="double", default=1,
		help="alternate purity [default= %default]", metavar="double"),
	make_option(c("-pl", "--ploidy"), type="double", default=2,
		help="alternate purity [default= %default]", metavar="double"),
	make_option(c("-m", "--min.reads.normal"), type="double", default=10,
		help="min reads to make genotype call [default= %default]", metavar="double"),
	make_option(c("-n", "--remove_centromeres"), type="logical", default=TRUE,
		help="remove all snp and copy-number data from centromeric regions", metavar="logical"),
	make_option(c("-t", "--cancer_type"), type="character", default="all",
		help="Cancer type to define priors [default= %default]", metavar="character"),
	make_option(c("-b", "--breaks_method"), type="character", default="het",
		help="breaks method for segmentation [default= %default]", metavar="character"),
	make_option(c("-w", "--window"), type="integer", default=1000000,
		help="window size for plotting [default= %default]", metavar="integer"),
	make_option(c("-r", "--ratio_priority"), type="logical", default=FALSE,
		help="Use only depth ratio for segmentation? [default= %default]", metavar="logical"),
	make_option(c("-y", "--remove_Y_vars"), type="logical", default=TRUE,
		help="Remove Y chromosome variants? [default= %default]", metavar="logical"),
	make_option(c("-g", "--ref_type"), type="character", default='hg38',
		help="Reference type (one of hg19, GRCh37, hg38, GRCh38) [default= %default]", metavar="character"),
	make_option(c("-f", "--min_reads_baf"), type="integer", default=1,
		help="Threshold on the depth of the positions included to calculate the average BAF for segment. Set to extreme value (ex. 1x10^6) to exclude BAF from calculations [default= %default]", metavar="integer")
	); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# extract options
snp.file <- opt$snp_file;
cnv.file <- opt$cnv_file;
outdir <- opt$out_dir;
purity <- opt$purity;

usr_cancer_type <- opt$cancer_type;
ud_ploidy <- opt$ploidy;

min.reads.normal <- opt$min.reads.normal;
min.reads.baf <- opt$min_reads_baf;

remove_Y <- opt$remove_Y_vars;

ref.type <- opt$ref_type;

print("Running Sequenza with the following options:");
print(opt);

### MAIN ###########################################################################################
# add another output directory for user defined purity 
output_udp <- paste0(outdir,"/output_udp_",purity,"/");
dir.create(output_udp, showWarnings = FALSE, recursive=TRUE);

# read in the snp file
snp <- read.delim(snp.file);

# read in the cnv file
cnv <- read.delim(cnv.file);

# remove any non-canonical chromosomes
chr_list <- c(paste0("chr",1:22),"chrX","chrY");
if (grepl('^GRCh', ref.type)) { chr_list <- c(1:22,"X","Y"); }

# remove chromosome Y because it seems to cause problems
if (remove_Y) {
	chr_list <- chr_list[!grepl('Y', chr_list)];
	}

snp <- subset(snp,chrom %in% chr_list);
cnv <- subset(cnv,chrom %in% chr_list);

# remove centromeric regions
if (opt$remove_centromeres) {

	if (ref.type %in% c('hg19','GRCh37')) {
		centromere.file <- '/cluster/projects/pughlab/references/ucsc/hg19/centromeres.txt';
		idx <- 1:3;
		} else if (ref.type %in% c('hg38','GRCh38')) {
		centromere.file <- '/cluster/projects/pughlab/references/ucsc/hg38/centromeres.txt';
		idx <- 2:4;
		}

	cent <- read.delim(centromere.file, header = FALSE);

	del_snp <- BedBedOverlap(snp[,c(1,2,2)], cent[,idx]);
	snp <- snp[!del_snp,];

	del_cnv <- BedBedOverlap(cnv[,c(1:3)], cent[,idx]);
	cnv <- cnv[!del_cnv,];
	 
	}

# bug where "chr" causes problems with "fast" breaks method
if (opt$breaks_method=="fast") {
	snp$chrom <- gsub("chr","",snp$chrom);
	cnv$chrom <- gsub("chr","",cnv$chrom);
	}

# prepare sequenza data file
seqz.data <- VarScan2seqz(
	varscan.somatic = snp, 
	varscan.copynumber = cnv
	);

file <- gsub(".snp", "",snp.file);
write.table(seqz.data, file , col.names = TRUE, row.names = FALSE, sep = "\t");

# full vs. het methods, applicable in different scenarios
data <- sequenza.extract(
	file, gz = FALSE,
	breaks.method = opt$breaks_method, 
	window = opt$window, 
	min.reads.normal = min.reads.normal,
	min.reads.baf = min.reads.baf
	);

### SEQUENZA #######################################################################################
# feed in prior probabilities based on cancer type
load(file="/cluster/projects/pughlab/src/sequenza_wrapper/PANCAN_ASCAT_ploidy_prob.Rdata")

priors <- subset(ploidy_table, cancer_type==toupper(usr_cancer_type))[c("CN","value")];

CP.example <- sequenza.fit(
	data,
	priors.table = priors,
	ratio.priority = opt$ratio_priority
	);

# output series of results files (images and text files)
sequenza.results(
	sequenza.extract = data,
	cp.table = CP.example,
	sample.id = basename(file),
	out.dir = outdir
	);

# output results for all original alternative solutions

### Extract all alternative solutions###
print ("Getting alternative purity solutions");

alt <- alternative.cp.solutions(CP.example);

purities <- alt$cellularity;
ploidies <- alt$ploidy;

## output series of results files (images and text files) for the alternative purity solutions
## index starts at 2 because the first solution is already output as the top hit.

if (length(purities) > 1) { ## bug fix to not run alt solutions if there are no alt solutions.
	for (i in 2:length(purities)) {
		print(paste("Creating new directories and printing solution:", i));
		output_alt <- paste0(outdir,"/sol",i,"_",purities[i],"/");
		dir.create(output_udp, showWarnings = FALSE, recursive = TRUE);
		sequenza.results(
			sequenza.extract = data, 
			cp.table = CP.example,
			sample.id = basename(file) , 
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
	cp.table = CP.example, 
	sample.id = basename(file), 
	out.dir = output_udp,
	cellularity = purity,
	ploidy = ud_ploidy
	);

## Ouput a total Copy-number seg file for viewing in IGV etc.
data.seg <- data$segments[[1]];
for (i in 2:length(data$segments)) {
	data.seg <- rbind(data.seg,data$segments[[i]]);
	}

colnames(data.seg) <- c("chrom","loc.start","loc.end","Bf","N.BAF","sd.BAF","seg.mean","num.mark","sd.ratio");
data.seg$seg.mean <- log2(data.seg$seg.mean);
data.seg$ID <- rep(basename(file),nrow(data.seg));
data.seg <- data.seg[,c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")];

setwd(outdir);

write.table(
	data.seg,
	paste0(basename(file),"_Total_CN.seg"),
	row.names = FALSE,
	quote = FALSE,
	sep = "\t"
	);
