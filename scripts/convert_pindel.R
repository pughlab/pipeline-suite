#!/usr/bin/env Rscript

# will input command Rscript convert_pindel.R /path/to/input -o output.tab
args <- commandArgs(trailingOnly = TRUE);

n.args <- length(args);

input.files <- args[1:(n.args-2)];
output <- args[n.args];

# write function to read in data
read.pindel <- function(file) {
	tmp <- tryCatch(
		read.delim(file, comment.char = '#', header = F),
		error = function(e) { NA }
		);
	if (is.logical(tmp)) { return(NULL); } else { return(tmp); }
	}

# read in input file
input.list <- lapply(
	input.files,
	read.pindel
	);

sv.data <- do.call(rbind, input.list);

if (is.null(sv.data)) {
	print(paste0('found ', length(input.files), ' input files consisting of 0 variants'));
	write(
		c('#break1_chromosome','break1_position_start','break1_position_end','break2_chromosome','break2_position_start','break2_position_end'),
		ncolumns = 6,
		file = output,
		sep = '\t'
		);	
	} else {

	# supply header
	colnames(sv.data) <- c('IDX','event_type','break1_chromosome','break1_position_start','break1_position_end','break2_chromosome','break2_position_start','break2_position_end','support_1','support_2','quality');

	# drop incomplete cases
	sv.data <- sv.data[!is.na(sv.data$break1_position_end),];
	sv.data <- sv.data[!is.na(sv.data$break2_position_end),];

	# initiate object for output
	total.variants <- nrow(sv.data);

	# print some metrics	
	print(paste0('found ', length(input.files), ' input files consisting of ', total.variants, ' variants'));

	gc();

	# do some filtering
	remove.idx <- c(
		which(sv.data$event_type == 'D' & sv.data$quality < 60),
		which(sv.data$event_type == 'TD' & sv.data$quality < 60),
		which(sv.data$event_type == 'INV' & sv.data$quality < 60),
		which(sv.data$event_type == 'D' & sv.data$support_1 < 5),
		which(sv.data$event_type == 'TD' & sv.data$support_1 < 5),
		which(sv.data$event_type == 'INV' & sv.data$support_1 < 5),
		which(sv.data$event_type == 'INT' & sv.data$support_1 < 5),
		which(sv.data$event_type == 'LI' & (sv.data$support_1 < 5 & sv.data$support_2 < 5))
		);
	sv.data <- sv.data[-remove.idx,];

	gc();

	# refactor event type
	sv.data$event_type <- factor(
		sv.data$event_type,
		levels = c('D','LI','INV','INT','TD'),
		labels = c('DEL','INS','INV','BND','DUP')
		);

	# calculate SV length
	sv.data$Length <- apply(
		sv.data[,grepl('position', colnames(sv.data))],
		1,
		function(i) {
			max(i) - min(i)
			}
		);

	# do some more filtering
	# for deletions: only keep length > 100 bases
	del.idx <- which(sv.data$event_type == 'DEL' & sv.data$Length >= 100);

	# for insertions: only keep length > 100 bases
	ins.idx <- which(sv.data$event_type == 'INS' & sv.data$Length >= 100);

	# for duplications: only keep length > 50 bases
	dup.idx <- which(sv.data$event_type == 'DUP' & sv.data$Length >= 50);

	# for inversions: keep all remaining
	inv.idx <- which(sv.data$event_type == 'INV');

	# for translocations: keep all remaining
	bnd.idx <- which(sv.data$event_type == 'BND');

	keep.idx <- c(del.idx, ins.idx, dup.idx, inv.idx, bnd.idx);

	if (length(keep.idx) > 0) {

		# fill in tracking ids
		sv.data$tracking_id <- NA;

		if (length(del.idx) > 0) {
			sv.data[del.idx,]$tracking_id <- paste0('pindel-DEL_',1:length(del.idx));
			}
		if (length(ins.idx) > 0) {
			sv.data[ins.idx,]$tracking_id <- paste0('pindel-INS_',1:length(ins.idx));
			}
		if (length(dup.idx) > 0) {
			sv.data[dup.idx,]$tracking_id <- paste0('pindel-DUP_',1:length(dup.idx));
			}
		if (length(inv.idx) > 0) {
			sv.data[inv.idx,]$tracking_id <- paste0('pindel-INV_',1:length(inv.idx));
			}
		if (length(bnd.idx) > 0) {
			sv.data[bnd.idx,]$tracking_id <- paste0('pindel-BND_',1:length(bnd.idx));
			}

		output.data <- sv.data[!is.na(sv.data$tracking_id),c(13, 2:8)];

		output.data <- output.data[!duplicated(output.data[,-1]),];

		gc();

		# indicate chromosomes to keep
		keep.chrs <- c(1:22,'X','Y');
		if (any(grepl('chr',output.data$break1_chromosome))) { keep.chrs <- paste0('chr',keep.chrs); }

		# filter data
		output.data <- output.data[which(output.data$break1_chromosome %in% keep.chrs & 
			output.data$break2_chromosome %in% keep.chrs),];

		if (nrow(output.data) > 0) {
			print(paste0('collapsed to ', nrow(output.data), ' rows'));

			colnames(output.data)[1] <- '#tracking_id';

			# refactor event type
			output.data$event_type <- factor(
				output.data$event_type,
				levels = c('DEL','INS','INV','BND','DUP'),
				labels = c('deletion','insertion','inversion','translocation','duplication')
				);

			# perform 1 final check to make sure all start < end
			if (any(output.data$break1_position_start > output.data$break1_position_end)) {
				idx <- which(output.data$break1_position_start > output.data$break1_position_end);
				tmp <- output.data;
				tmp[idx,'break1_position_start'] <- apply(
					output.data[idx,c('break1_position_start','break1_position_end')],
					1, function(i) { sort(i)[1]; } );
				tmp[idx,'break1_position_end'] <- apply(
					output.data[idx,c('break1_position_start','break1_position_end')],
					1, function(i) { sort(i)[2]; } );

				output.data <- tmp;
				}

			if (any(output.data$break2_position_start > output.data$break2_position_end)) {
				idx <- which(output.data$break2_position_start > output.data$break2_position_end);
				tmp <- output.data;
				tmp[idx,'break2_position_start'] <- apply(
					output.data[idx,c('break2_position_start','break2_position_end')],
					1, function(i) { sort(i)[1]; } );
				tmp[idx,'break2_position_end'] <- apply(
					output.data[idx,c('break2_position_start','break2_position_end')],
					1, function(i) { sort(i)[2]; } );

				output.data <- tmp;
				}

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

			print('No variants remaining after filtering.');

			write(
				c('#break1_chromosome','break1_position_start','break1_position_end','break2_chromosome','break2_position_start','break2_position_end'),
				ncolumns = 6,
				file = output,
				sep = '\t'
				);	
			}
		}
	}
