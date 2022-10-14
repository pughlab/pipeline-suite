#!/usr/bin/env Rscript

# will input command Rscript conver_arriba.R input.txt -o output.tab
args <- commandArgs(trailingOnly = TRUE);

n.args <- length(args);

input.files <- args[1:(length(args)-2)];
output <- args[length(args)];

# read in input file
input.list <- lapply(
	input.files,
	read.delim,
	stringsAsFactors = FALSE
	);

input.data <- do.call(rbind, input.list);

if ( (is.null(input.data)) | (nrow(input.data) == 0) ) {
	print(paste0('found ', length(input.files), ' input files consisting of 0 variants'));
	write(
		c('#break1_chromosome','break1_position_start','break1_position_end','break2_chromosome','break2_position_start','break2_position_end'),
		ncolumns = 6,
		file = output,
		sep = '\t'
		);	
	} else {

	input.data$ID <- paste0('arriba-',1:nrow(input.data));

	# format output data
	output.data <- NULL;

	for (i in 1:nrow(input.data)) {

		tracking.id <- input.data[i,]$ID;
		break1 <- unlist(strsplit(input.data[i,]$breakpoint1,':'));
		break2 <- unlist(strsplit(input.data[i,]$breakpoint2,':'));

		tmp <- data.frame(
			'#tracking_id'		= rep(tracking.id, 4),
			event_type		= 'None',
			break1_chromosome	= paste0('chr',break1[1]),
			break1_position_start	= break1[2],
			break1_position_end	= break1[2],
			break1_orientation	= c('L','L','R','R'),
			break1_strand		= '?',
			break1_seq		= 'None',
			break2_chromosome	= paste0('chr',break2[1]),
			break2_position_start	= break2[2],
			break2_position_end	= break2[2],
			break2_orientation	= c('L','R','L','R'),
			break2_strand		= '?',
			break2_seq		= 'None',
			opposing_strands	= c('True','False','False','True'),
			stranded		= 'False',
			tools			= 'arriba',
			untemplated_seq		= NA
			);

		output.data <- rbind(output.data,tmp);
		}

	if (any(grepl('chrchr', output.data$break1_chromosome))) {
		output.data$break1_chromosome <- gsub('chrchr','chr', output.data$break1_chromosome);
		output.data$break2_chromosome <- gsub('chrchr','chr', output.data$break2_chromosome);
		}

	write.table(
		output.data,
		file = output,
		row.names = FALSE,
		col.names = TRUE,
		quote = FALSE,
		sep = '\t'
		);
	}
