#!/usr/bin/env Rscript

# will input command Rscript convert_pindel.R /path/to/input -o output.tab
args <- commandArgs(trailingOnly = TRUE);

n.args <- length(args);

input.dirs <- args[1:(n.args-2)];
output <- args[n.args];

# read in input file
input.files <- c();
for (input.path in input.dirs) {
	input.files <- c(
		input.files,
		list.files(path = input.path, full.name = TRUE)
		);
	}

input.files <- input.files[grepl('_D$|_INV$|_INT_final$|_LI$|_TD$', input.files)];

# initiate object to hold data
sv.data <- list();
total.variants <- 0;

# write function to read in data
read.pindel <- function(file) {
	tmp <- tryCatch(
		read.delim(file, comment.char = '#', header = F),
		error = function(e) { NA }
		);
	if (is.logical(tmp)) { return(NULL); } else { return(tmp); }
	}

# for deletions: only keep length > 50
deletions <- do.call(
	rbind,
	lapply(input.files[grep('_D$', input.files)], read.pindel)
	);
deletions <- deletions[!is.na(deletions$V8),];

total.variants <- total.variants + nrow(deletions);

# keep quality >= 60
deletions <- deletions[which(as.numeric(gsub('SUM_MS ','',as.character(deletions$V16))) >= 60),];
# keep read support >= 5
deletions <- deletions[which(deletions$V10 >= 5),];

if (nrow(deletions) > 0) {

	tmp <- data.frame(
		tracking_id = paste0('pindel-DEL','_',deletions$V1),
		event_type = rep('deletion', nrow(deletions)),
		break1_chromosome = as.character(sapply(deletions$V4, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break1_position_start = as.numeric(
			apply(deletions[,c('V5','V7')], 1, function(i) {
				min(
					unlist(strsplit(as.character(i[1]),' '))[2],
					unlist(strsplit(as.character(i[2]),' '))[2]
					)
				})
			),
		break1_position_end = as.numeric( 
			apply(deletions[,c('V5','V7')], 1, function(i) {
				max(
					unlist(strsplit(as.character(i[1]),' '))[2],
					unlist(strsplit(as.character(i[2]),' '))[2]
					)
				})
			),
		break2_chromosome = as.character(sapply(deletions$V4, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break2_position_start = apply(deletions[,c('V6','V8')],1,function(i) {
			min(as.numeric(as.character(i[1])), as.numeric(as.character(i[2]))) } ), 
		break2_position_end = apply(deletions[,c('V6','V8')],1,function(i) {
			max(as.numeric(as.character(i[1])), as.numeric(as.character(i[2]))) } ) 
		);

	tmp <- tmp[!duplicated(tmp[,-1]),];
	tmp <- tmp[which((tmp$break2_position_end - tmp$break1_position_start) >= 50),];

	if (nrow(tmp) > 0) { sv.data[['DEL']] <- tmp; }
	}
	
# for insertions (long):
insertions <- do.call(
	rbind,
	lapply(input.files[grep('_LI$', input.files)], read.pindel)
	);
insertions <- insertions[which(insertions$V8 != ''),];

total.variants <- total.variants + nrow(insertions);

# keep read support >= 5
insertions <- insertions[which(
	abs(as.numeric(gsub(' ','',insertions$V5))) >= 5 | abs(as.numeric(gsub(' ','',insertions$V7))) >= 5
	),];

if (nrow(insertions) > 0) {

	tmp <- data.frame(
		tracking_id = paste0('pindel-INS','_',insertions$V1),
		event_type = rep('insertion', nrow(insertions)),
		break1_chromosome = as.character(sapply(insertions$V3, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break1_position_start = as.numeric(as.character(insertions$V4)),
		break1_position_end = as.numeric(as.character(insertions$V4)),
		break2_chromosome = as.character(sapply(insertions$V3, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break2_position_start = as.numeric(as.character(insertions$V6)),
		break2_position_end = as.numeric(as.character(insertions$V6))
		);

	tmp <- tmp[!duplicated(tmp[,-1]),];
	tmp <- tmp[which(abs(tmp$break2_position_end - tmp$break1_position_start) >= 20),];

	if (nrow(tmp) > 0) { sv.data[['INS']] <- tmp; }
	}

# for duplications (tandom):
duplications <- do.call(
	rbind,
	lapply(input.files[grep('_TD$', input.files)], read.pindel)
	);
duplications <- duplications[!is.na(duplications$V8),];

total.variants <- total.variants + nrow(duplications);

# keep quality >= 60
duplications <- duplications[which(as.numeric(gsub('SUM_MS ','',as.character(duplications$V16))) >= 60),];
# keep read support >= 5
duplications <- duplications[which(duplications$V10 >= 5),];

if (nrow(duplications) > 0) {

	tmp <- data.frame(
		tracking_id = paste0('pindel-DUP','_',duplications$V1),
		event_type = rep('duplication', nrow(duplications)),
		break1_chromosome = as.character(sapply(duplications$V4, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break1_position_start = as.numeric(sapply(duplications$V5, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break1_position_end = as.numeric(sapply(duplications$V5, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break2_chromosome = as.character(sapply(duplications$V4, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break2_position_start = as.numeric(as.character(duplications$V6)),
		break2_position_end = as.numeric(as.character(duplications$V6))
		);

	tmp <- tmp[!duplicated(tmp[,-1]),];
	tmp <- tmp[which((tmp$break2_position_end - tmp$break1_position_start) >= 20),];

	if (nrow(tmp) > 0) { sv.data[['DUP']] <- tmp; }
	}

# for inversions:
inversions <- do.call(
	rbind,
	lapply(input.files[grep('_INV$', input.files)], read.pindel)
	);
inversions <- inversions[!is.na(inversions$V8),];

total.variants <- total.variants + nrow(inversions);

# keep quality >= 60
inversions <- inversions[which(as.numeric(gsub('SUM_MS ','',as.character(inversions$V16))) >= 60),];
# keep read support >= 5
inversions <- inversions[which(inversions$V10 >= 5),];

if (nrow(inversions) > 0) {

	tmp <- data.frame(
		tracking_id = paste0('pindel-INV','_',inversions$V1),
		event_type = rep('inversion', nrow(inversions)),
		break1_chromosome = as.character(sapply(inversions$V4, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break1_position_start = as.numeric(sapply(inversions$V5, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break1_position_end = as.numeric(sapply(inversions$V5, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break2_chromosome = as.character(sapply(inversions$V4, function(i) { unlist(strsplit(as.character(i),' '))[2] } )),
		break2_position_start = as.numeric(as.character(inversions$V6)),
		break2_position_end = as.numeric(as.character(inversions$V6))
		);

	tmp <- tmp[!duplicated(tmp[,-1]),];

	if (nrow(tmp) > 0) { sv.data[['INV']] <- tmp; }
	}

# for fusions (interchromosomal events; INT):
fusions <- do.call(
	rbind,
	lapply(input.files[grep('_INT_final$', input.files)], read.pindel)
	);

total.variants <- total.variants + nrow(fusions);

# keep read support >= 5
fusions <- fusions[which(fusions$V12 >= 5),];

if (nrow(fusions) > 0) {

	tmp <- data.frame(
		tracking_id = paste0('pindel-BND','_',1:nrow(fusions)),
		event_type = rep('translocation', nrow(fusions)),
		break1_chromosome = as.character(fusions$V2),
		break1_position_start = apply(fusions[,c('V16','V25')],1,function(i) { min(i, na.rm = T) } ),
		break1_position_end = apply(fusions[,c('V16','V25')],1,function(i) { max(i, na.rm = T) } ),
		break2_chromosome = as.character(fusions$V6),
		break2_position_start = apply(fusions[,c('V19','V28')],1,function(i) { min(i, na.rm = T) } ), 
		break2_position_end = apply(fusions[,c('V19','V28')],1,function(i) { max(i, na.rm = T) } )
		);

	tmp <- tmp[!duplicated(tmp[,-1]),];

	if (nrow(tmp) > 0) { sv.data[['BND']] <- tmp; }
	}

print(paste0('found ', length(input.files), ' input files consisting of ', total.variants, ' variants'));

# merge events
input.data <- do.call(rbind, sv.data);

# indicate chromosomes to keep
keep.chrs <- c(1:22,'X','Y');
if (any(grepl('chr',input.data$break1_chromosome))) { keep.chrs <- paste0('chr',keep.chrs); }

# filter data
output.data <- input.data[which(input.data$break1_chromosome %in% keep.chrs & 
	input.data$break2_chromosome %in% keep.chrs),];

print(paste0('collapsed to ', nrow(output.data), ' rows'));

# save to file
if (nrow(output.data) > 0) {

	colnames(output.data)[1] <- '#tracking_id';

	# save to file	
	write.table(
		output.data,
		file = output,
		row.names = FALSE,
		col.names = TRUE,
		quote = FALSE,
		sep = '\t'
		);
	} else {
	write(
		c('#break1_chromosome','break1_position_start','break1_position_end','break2_chromosome','break2_position_start','break2_position_end'),
		ncolumns = 6,
		file = output,
		sep = '\t'
		);	
	}
