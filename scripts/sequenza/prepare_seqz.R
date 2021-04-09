### SequenzaSingleSample_TuneA.R ###################################################################
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

### PARSE OPTIONS ##################################################################################
# For now, test if commands are in original, trailing format, or new opt-parse format
option_list <- list(
	make_option(c("-s", "--snp_file"), type="character", default=NULL, 
		help="varscan snp file name", metavar="character"),
	make_option(c("-c", "--cnv_file"), type="character", default=NULL,
		help="varscan copy number file name [default= %default]", metavar="character"),
	make_option(c("-o", "--out_dir"), type="character", default=NULL,
		help="output directory [default= %default]", metavar="character"),
	make_option(c("-m", "--min.reads.normal"), type="double", default=10,
		help="min reads to make genotype call [default= %default]", metavar="double"),
	make_option(c("-n", "--remove_centromeres"), type="logical", default=TRUE,
		help="remove all snp and copy-number data from centromeric regions", metavar="logical"),
	make_option(c("-w", "--window"), type="integer", default=1000000,
		help="window size for plotting [default= %default]", metavar="integer"),
	make_option(c("-r", "--ratio_priority"), type="logical", default=FALSE,
		help="Use only depth ratio for segmentation? [default= %default]", metavar="logical"),
	make_option(c("-a", "--assembly"), type="character", default='hg38',
		help="Reference type (one of hg19, GRCh37, hg38, GRCh38) [default= %default]",
		metavar="character"),
	make_option(c("-f", "--min_reads_baf"), type="integer", default=1,
		help="Threshold on the depth of the positions included to calculate the average BAF for segment. Set to extreme value (ex. 1x10^6) to exclude BAF from calculations [default= %default]", metavar="integer")
	); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# extract options
snp.file <- opt$snp_file;
cnv.file <- opt$cnv_file;
outdir <- opt$out_dir;

min.reads.normal <- opt$min.reads.normal;
min.reads.baf <- opt$min_reads_baf;

assembly <- opt$assembly;

print("Running Sequenza with the following options:");
print(opt);

### MAIN ###########################################################################################
# add output directory 
dir.create(outdir, showWarnings = FALSE, recursive = TRUE);

print("Reading in snp and cnv calls");

# read in the snp file
snp <- read.delim(snp.file);

# read in the cnv file
cnv <- read.delim(cnv.file);

# remove any non-canonical chromosomes
chr_list <- c(paste0("chr",1:22),"chrX");

snp <- subset(snp,chrom %in% chr_list);
cnv <- subset(cnv,chrom %in% chr_list);

# remove centromeric regions
if (opt$remove_centromeres) {

	print("Removing centromere regions");

	if (assembly %in% c('hg19','GRCh37')) {
		centromere.file <- '/cluster/projects/pughlab/references/ucsc/hg19/centromeres.txt';
		idx <- 1:3;
		} else if (assembly %in% c('hg38','GRCh38')) {
		centromere.file <- '/cluster/projects/pughlab/references/ucsc/hg38/centromeres.txt';
		idx <- 2:4;
		}

	cent <- read.delim(centromere.file, header = FALSE);

	del_snp <- BedBedOverlap(snp[,c(1,2,2)], cent[,idx]);
	snp <- snp[!del_snp,];

	del_cnv <- BedBedOverlap(cnv[,c(1:3)], cent[,idx]);
	cnv <- cnv[!del_cnv,];
	 
	}

# if using the 'called' input cnv file (this has filtered out low coverage regions)
# specifically, normal_depth < 20 & tumor_depth < 10
if (! 'log2_ratio' %in% colnames(cnv)) {
	colnames(cnv)[which(colnames(cnv) == 'raw_ratio')] <- 'log2_ratio';
	}

# prepare sequenza data file
print("Preparing SEQZ data");
seqz.data <- VarScan2seqz(
	varscan.somatic = snp, 
	varscan.copynumber = cnv
	);

filestem <- gsub(".snp", "", basename(snp.file));

print(paste0("Writing SEQZ to file: ", paste0(outdir, '/', filestem, '.seqz')));

write.table(
	seqz.data,
	file = paste0(outdir, '/', filestem, '.seqz'),
	col.names = TRUE,
	row.names = FALSE,
	sep = "\t"
	);

