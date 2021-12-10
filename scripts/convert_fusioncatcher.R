#!/usr/bin/env Rscript

# will input command Rscript conver_fusioncatcher.R input.txt -o output.tab
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

	input.data$ID <- paste0('fusioncatcher-',1:nrow(input.data));

	unique.data <- aggregate(
		input.data$ID,
		by = list(
			'break1' = input.data$Fusion_point_for_gene_1.5end_fusion_partner.,
			'break2' = input.data$Fusion_point_for_gene_2.3end_fusion_partner.
			),
		FUN = function(i) { paste(i, collapse = ';') }
		);

	# format output data
	output.data <- NULL;

	for (i in 1:nrow(unique.data)) {

		tracking.id <- unique.data[i,]$x;
		break1 <- unlist(strsplit(unique.data[i,]$break1,':'));
		break2 <- unlist(strsplit(unique.data[i,]$break2,':'));

		tmp <- data.frame(
			'#tracking_id'		= rep(tracking.id, 4),
			event_type		= 'None',
			break1_chromosome	= paste0('chr',break1[1]),
			break1_position_start	= break1[2],
			break1_position_end	= break1[2],
			break1_orientation	= c('L','L','R','R'),
			break1_strand		= '?',#break1[3],
			break1_seq		= 'None',
			break2_chromosome	= paste0('chr',break2[1]),
			break2_position_start	= break2[2],
			break2_position_end	= break2[2],
			break2_orientation	= c('L','R','L','R'),
			break2_strand		= '?', #break2[3],
			break2_seq		= 'None',
			opposing_strands	= c('True','False','False','True'),
			stranded		= 'False',
			tools			= 'fusioncatcher',
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
