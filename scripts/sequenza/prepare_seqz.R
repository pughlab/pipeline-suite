### SequenzaSingleSample_TuneA.R ###################################################################
# Script to run Sequenza (v2.1) on snps and cnv files from Varscan
#
# Modified from /cluster/projects/pughlab/src/sequenza_wrapper/SequenzaSingleSample_v2.1.R
# to incorporate reference type

### PREAMBLE #######################################################################################
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

# varscan2seqz function (source: sequenza v2.1.2)
VarScan2seqz <- function (varscan.somatic, varscan.copynumber = NULL, normal_var_freq = 0.25) {
	iupac.nucs <- setNames(c("A", "C", "G", "GT", "AC", "AG", "CG", "T", "AT", "CT"),
		c("A", "C", "G", "K", "M", "R", "S", "T", "W", "Y"));

	zygosity.vect <- setNames(c("hom", "hom", "hom", "het", "het", "het", "het", "hom", "het", "het"),
		c("A", "C", "G", "K", "M", "R", "S", "T", "W", "Y"));

	varscan.somatic <- varscan.somatic[varscan.somatic$somatic_status != "Unknown", ];
	varscan.somatic$normal_var_freq <- as.numeric(sub("%", "", varscan.somatic$normal_var_freq))/100;
	varscan.somatic$tumor_var_freq <- as.numeric(sub("%", "", varscan.somatic$tumor_var_freq))/100;

	zygosity.normal <- zygosity.vect[varscan.somatic$normal_gt];
	if (normal_var_freq > 0.5) { normal_var_freq <- 1 - normal_var_freq;  }
	idx <- which(zygosity.normal == "het" & (
		normal_var_freq >= varscan.somatic$normal_var_freq | 
		varscan.somatic$normal_var_freq >= (1 - normal_var_freq)));

	varscan.somatic <- varscan.somatic[-idx, ];
	zygosity.normal <- zygosity.vect[varscan.somatic$normal_gt];

	AB.normal <- iupac.nucs[varscan.somatic$normal_gt];
	AB.tumor <- rep(".", length(AB.normal));
	strand <- AB.tumor;

	depth.normal <- varscan.somatic$normal_reads1 + varscan.somatic$normal_reads2;
	depth.tumor <- varscan.somatic$tumor_reads1 + varscan.somatic$tumor_reads1;
	depth.ratio <- depth.tumor/depth.normal;

	Af <- 1 - varscan.somatic$tumor_var_freq;
	Bf <- rep(0, length(Af));
	idx <- zygosity.normal == "het" & Af < 0.5;
	Af[idx] <- 1 - Af[idx];
	idx <- zygosity.normal == "het";
	Bf[idx] <- 1 - Af[idx];

	genotype_1 <- varscan.somatic[idx, c("ref", "var")];
	genotype_1 <- apply(genotype_1, 1, paste, collapse = "");
	AB.normal[idx] <- genotype_1;
	idx <- idx & varscan.somatic$tumor_var_freq > 0.5;
	genotype_2 <- varscan.somatic[idx, c("var", "ref")];
	genotype_2 <- apply(genotype_2, 1, paste, collapse = "");
	AB.normal[idx] <- genotype_2;

	idx <- zygosity.normal == "hom" & varscan.somatic$somatic_status == "Somatic";
	if (sum(idx) > 0) {
		mut.b <- cbind(as.character(iupac.nucs[varscan.somatic$tumor_gt[idx]]), 
			as.character(varscan.somatic$normal_gt[idx]));

		mut.b <- sapply(X = 1:sum(idx), FUN = function(x) {
			gsub(x = mut.b[x,1], pattern = mut.b[x, 2], replacement = "");
			});

		mut <- paste0(mut.b, varscan.somatic$tumor_var_freq[idx]);
		strand[idx] <- paste0(mut.b,
			varscan.somatic$tumor_reads2_plus[idx]/(varscan.somatic$tumor_reads2_plus[idx] + 
			varscan.somatic$tumor_reads2_minus[idx]));
		AB.tumor[idx] <- mut;
		}

	res <- data.frame(
		chromosome = as.character(varscan.somatic$chrom),
		position = varscan.somatic$position,
		base.ref = as.character(varscan.somatic$ref),
		depth.normal = depth.normal,
		depth.tumor = depth.tumor,
		depth.ratio = depth.ratio,
		Af = round(Af, 3),
		Bf = round(Bf, 3),
		zygosity.normal = zygosity.normal,
		GC.percent = 50,
		good.reads = round(depth.tumor, 2),
		AB.normal = AB.normal,
		AB.tumor = AB.tumor,
		tumor.strand = strand,
		stringsAsFactors = FALSE
		);

	normal.pos <- res$zygosity.normal == "hom" & res$AB.tumor == ".";
	res <- res[res$depth.ratio > 0 & !is.infinite(res$depth.ratio) & !normal.pos, ];

	if (!is.null(varscan.copynumber)) {
		smart.id <- order(c(1:nrow(varscan.copynumber), 1:nrow(varscan.copynumber) + 0.5));
		varscan.copynumber$log2_ratio <- 2^(varscan.copynumber$log2_ratio);
		varscan.copynumber$normal_depth <- round(varscan.copynumber$normal_depth, 0);
		varscan.copynumber$tumor_depth <- round(varscan.copynumber$tumor_depth, 0);
		varscan.copynumber$chrom <- as.character(varscan.copynumber$chrom);

		mat.t <- data.frame(
			chromosome = c(varscan.copynumber$chrom, varscan.copynumber$chrom)[smart.id],
			position = c(varscan.copynumber$chr_start, varscan.copynumber$chr_stop)[smart.id],
			base.ref = "N",
                        depth.normal = c(varscan.copynumber$normal_depth, 
				varscan.copynumber$normal_depth)[smart.id],
			depth.tumor = c(varscan.copynumber$tumor_depth, 
				varscan.copynumber$tumor_depth)[smart.id],
			depth.ratio = c(varscan.copynumber$log2_ratio, 
				varscan.copynumber$log2_ratio)[smart.id],
			Af = 1,
			Bf = 0,
			zygosity.normal = "hom",
			GC.percent = c(varscan.copynumber$gc_content, 
				varscan.copynumber$gc_content)[smart.id],
			stringsAsFactors = FALSE
			);

		mat.t <- cbind(mat.t, good.reads = mat.t$depth.tumor, 
			AB.normal = "N", AB.tumor = ".", tumor.strand = ".");

		chrom.order <- unique(mat.t$chromosome);
		l.cnv <- split(mat.t, mat.t$chromosome);
		l.snp <- split(res, res$chromosome);

		for (i in names(l.snp)) {
			tab.i <- rbind(l.cnv[[i]], l.snp[[i]]);
			tab.i <- tab.i[order(tab.i$position), ];
			l.cnv[[i]] <- tab.i;
			}

		res <- do.call(rbind, l.cnv[chrom.order]);
		}
	return(res);
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
	make_option(c("-n", "--remove_centromeres"), type="logical", default=TRUE,
		help="remove all snp and copy-number data from centromeric regions", metavar="logical"),
	make_option(c("-a", "--assembly"), type="character", default='hg38',
		help="Reference type (one of hg19, GRCh37, hg38, GRCh38) [default= %default]",
		metavar="character"),
	make_option(c("-w", "--is_wgs"), type = "logical", default = FALSE, metavar = "logical",
		help="is the input data generated from WGS?")
	); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# extract options
snp.file <- opt$snp_file;
cnv.file <- opt$cnv_file;
outdir <- opt$out_dir;

assembly <- opt$assembly;

is.wgs <- opt$is_wgs;

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

#snp <- subset(snp,chrom %in% chr_list);
#cnv <- subset(cnv,chrom %in% chr_list);

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
	
filestem <- gsub(".snp", "", basename(snp.file));

if (! is.wgs) {

	snp_subset <- subset(snp,chrom %in% chr_list);
	cnv_subset <- subset(cnv,chrom %in% chr_list);

	# prepare sequenza data file
	print("Preparing SEQZ data");
	seqz.data <- VarScan2seqz(
		varscan.somatic = snp_subset, 
		varscan.copynumber = cnv_subset
		);

	write.table(
		seqz.data,
		file = paste0(outdir, '/', filestem, '.seqz'),
		col.names = TRUE,
		row.names = FALSE,
		sep = "\t"
		);

	} else {

	seqz.list <- list();

	for (chr in chr_list) {

		snp_subset <- subset(snp,chrom == chr);
		cnv_subset <- subset(cnv,chrom == chr);

		# prepare sequenza data file
		print(paste0("Preparing SEQZ data for ", chr, "..."));
		seqz.list[[chr]] <- VarScan2seqz(
			varscan.somatic = snp_subset, 
			varscan.copynumber = cnv_subset
			);
		gc();

		if (chr == chr_list[1]) {
			write.table(
				seqz.list[[chr]],
				file = paste0(outdir, '/', filestem, '.seqz'),
				col.names = TRUE,
				row.names = FALSE,
				sep = "\t"
				);
			} else {
			write.table(
				seqz.list[[chr]],
				file = paste0(outdir, '/', filestem, '.seqz'),
				col.names = FALSE,
				row.names = FALSE,
				sep = "\t",
				append = TRUE
				);
			}
		}
	}
	
print(paste0("Writing SEQZ to file: ", paste0(outdir, '/', filestem, '.seqz')));
