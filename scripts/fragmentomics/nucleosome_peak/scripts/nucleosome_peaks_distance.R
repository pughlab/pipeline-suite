# file: runFrag.R
# author: Derek Wong, Ph.D
# date: Aug 3rd, 2022
# modified: July 4th, 2023 by sprokopec

library(optparse);

## Set script variables
option_list <- list(
	make_option(c("--id"), type = "character", help = "sample id. Required"),
	make_option(c("--path"), type = "character", help = "Path to bedpe file. Required."),
	make_option(c("--outdir"), type = "character", help = "Path to output directory. Required."),
	make_option(c("--peaks"), type = "character", help = "Path to nucleosome peaks",
		default = '/cluster/projects/pughlab/projects/CHARM/nucleosome_peaks/hg38')
	);

opt <- parse_args(OptionParser(option_list = option_list));
print(opt);
options(scipen = 999, stringsAsFactors = FALSE);

## Load required packages
library(tidyverse);
library(data.table);
options(stringsAsFactors = FALSE);
options(bitmapType = 'cairo');

## Get variables from input script
id	<- opt$id;
outdir	<- opt$outdir;
bedpe_path	<- opt$path;
peaks_path	<- opt$peaks;

## Find files
files <- list.files(bedpe_path, id, full.names = TRUE);
peaks <- list.files(peaks_path, "chr", full.names = TRUE);
chrs <- files[grep("chr", files)];
chrs <- sapply(strsplit(chrs, "_\\s*"), tail, 1);
chrs <- gsub(".bedpe", "", chrs);

### Read in files for each chromosome
table <- data.frame(distance = c(-1000:1000));

for (i in c(1:length(chrs))) {

	### Set variables
	chr <- chrs[[i]];
	file <- files[grep(paste0(chr, ".bed"), files)];
	peak <- peaks[grep(paste0(chr, ".txt"), peaks)];

	### Read in files
	data <- read.delim(file, header = FALSE);
	nucleosomes <- read.delim(peak);

	### Format bedpe file
	data <- data[data$V1 == data$V4,];
	data <- data[data$V1 == chr,];
	data$V2 <- ifelse(data$V2 < data$V3, data$V2, data$V3);
	data$V5 <- ifelse(data$V5 < data$V6, data$V5, data$V6);
	data_f <- data[,1:2];
	data_r <- data[,c(4,6)];

	nucleosomes <- nucleosomes[nucleosomes$Chrom == chr,]
 
	### Find nearest peak
	setDT(data_f);
	setDT(data_r);
	setDT(nucleosomes);

	data_f[,merge:=V2];
	data_r[,merge:=V6];
	nucleosomes[,merge:=peak];

	setkeyv(data_f, c("merge"));
	setkeyv(data_r, c("merge"));
	setkeyv(nucleosomes, c("merge"));

	merged_f <- nucleosomes[data_f, roll = "nearest"];
	merged_f$distance <- merged_f$V2 - merged_f$peak;
	merged_f <- merged_f[abs(merged_f$distance) <= 1000,];

	merged_r <- nucleosomes[data_r, roll = "nearest"];
	merged_r$distance <- merged_r$V6 - merged_r$peak;
	merged_r <- merged_r[abs(merged_r$distance) <= 1000,];

	### Count frequencies
	distances <- data.frame(distance = c(-1000:1000));
	freq_f <- as.data.frame(table(merged_f$distance));
	freq_r <- as.data.frame(table(merged_r$distance));

	distances <- merge(distances, freq_f, by.x = "distance", by.y = "Var1", all = TRUE);
	distances <- merge(distances, freq_r, by.x = "distance", by.y = "Var1", all = TRUE);
	distances[is.na(distances)] <- 0;
	distances$frequency <- rowSums(distances[, 2:3]);

	distances <- distances[,c(1,4)];

	### Merge with table
	table <- merge(table, distances, by = "distance")
	}

colnames(table) <- c("distance", chrs);

write.table(
	table,
	file = file.path(outdir, paste0(id, "_peak_distance.txt")),
	row.names = FALSE,
	col.names = TRUE,
	sep = "\t"
	);

