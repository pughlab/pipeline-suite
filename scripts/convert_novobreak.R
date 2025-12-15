#!/usr/bin/env Rscript

# will input command Rscript convert_novobreak.R input.txt -o output.tab
args <- commandArgs(trailingOnly = TRUE);

n.args <- length(args);

input.files <- args[1:(n.args-2)];
output <- args[n.args];

# read in input file
input.list <- list();
for (file in input.files) {
	input.list[[file]] <- read.delim(file, stringsAsFactors = FALSE);
	}
input.data <- do.call(rbind, input.list);

# quick input check
if (nrow(input.data) == 0) {
	print('No variants detected.');

	write(
		c('#break1_chromosome','break1_position_start','break1_position_end','break2_chromosome','break2_position_start','break2_position_end'),
		ncolumns = 6,
		file = output,
		sep = '\t'
                );

	} else {

	total.variants <- nrow(input.data);

	# print some metrics
	print(paste0('found ', length(input.files), ' input files consisting of ', 
		total.variants, ' variants'));

	# give unique ids
	input.data$ID <- paste0('novobreak-',input.data$sv_type,'_',1:nrow(input.data));

	# indicate chromosomes to keep
	keep.chrs <- c(1:22,'X','Y');
	if (any(grepl('chr',input.data$left_chrom))) { keep.chrs <- paste0('chr',keep.chrs); }

	# filter data
	input.data <- input.data[which(input.data$left_chrom %in% keep.chrs & 
		input.data$right_chrom %in% keep.chrs),];

	if (nrow(input.data) == 0) {
		print('No variants remaining after filtering.');
		write(
			c('#break1_chromosome','break1_position_start','break1_position_end','break2_chromosome','break2_position_start','break2_position_end'),
			ncolumns = 6,
			file = output,
			sep = '\t'
			);

		} else {

		# format output data
		output.data <- data.frame(
			tracking_id = input.data$ID,
			event_type = c('deletion','duplication','inversion','translocation')[match(input.data$sv_type, c('DEL','DUP','INV','TRA'))],
			break1_chromosome = paste0('chr', input.data$left_chrom),
			break1_position_start = as.integer(input.data$left_position - 10),
			break1_position_end = as.integer(input.data$left_position + 10),
			break2_chromosome = paste0('chr', input.data$right_chrom),
			break2_position_start = as.integer(input.data$right_position - 10),
			break2_position_end = as.integer(input.data$right_position + 10)
			);

		if (any(grepl('chrchr', output.data$break1_chromosome))) {
			output.data$break1_chromosome <- gsub('chrchr','chr', output.data$break1_chromosome);
			output.data$break2_chromosome <- gsub('chrchr','chr', output.data$break2_chromosome);
			}

		# if no variants remain
		if (nrow(output.data) == 0) {

			print('No variants remaining after filtering.');

			write(
				c('#break1_chromosome','break1_position_start','break1_position_end','break2_chromosome','break2_position_start','break2_position_end'),
				ncolumns = 6,
				file = output,
				sep = '\t'
				);

			} else {
			print(paste0('collapsed to ', nrow(output.data), ' rows'));

			# save to file
			write.table(
				output.data,
				file = output,
				row.names = FALSE,
				col.names = TRUE,
				quote = FALSE,
				sep = '\t'
				);
			}
		}
	}
