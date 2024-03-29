# file: git_01-read_bam.R
# author: Derek Wong, Ph.D
# date: October 5th, 2021
# modified: June 28th, 2023 by sprokopec

## Read GAlignmentPairs
bamfile <- bam;
indexed.bam <- gsub("$", ".bai", bamfile);
if (!file.exists(indexed.bam)) {
	indexBam(bamfile);
	}

param <- ScanBamParam(
	flag = scanBamFlag(
		isPaired = TRUE,
		isProperPair = TRUE,
		isDuplicate = FALSE,
		isSecondaryAlignment = FALSE,
		isUnmappedQuery = FALSE
		),
	mapqFilter = 30
	);

galp <- readGAlignmentPairs(bamfile, param = param);
rm(param, indexed.bam);

## Filter reads: 90-220bp on Autosomes and mitochondrial reads
frags <- granges(
	keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode = "coarse"),
	on.discordant.seqnames = "drop"
	);

mito <- granges(
	keepSeqlevels(galp, paste0("chrM"), pruning.mode = "coarse"),
	on.discordant.seqnames = "drop"
	);

length <- length(galp);
rm(galp);

w.all <- width(frags);
frags <- frags[which(w.all >= 90 & w.all <= 220)];
rm(w.all);

gcs <- GCcontent(Hsapiens, unstrand(frags));
frags$gc <- gcs;
rm(gcs);

gc();
