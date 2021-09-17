#!/usr/bin/env Rscript

# will input command Rscript convert_svict.R input.txt -o output.tab
args <- commandArgs(trailingOnly = TRUE);

n.args <- length(args);

input.files <- args[1:(n.args-2)];
output <- args[n.args];

# read in input file
input.list <- list();
total.variants <- 0;

for (file in input.files) {

	print(paste0("reading: ", file));

	input.data <- tryCatch(
		read.delim(file, stringsAsFactors = FALSE, comment.char = '#', header = F),
		error = function(e) { NA }
		);

	if (is.logical(input.data)) {
		print(paste0("found 0 rows"));
		next;
		} else {
		print(paste0("found ", nrow(input.data), " rows"));
		total.variants <- total.variants + nrow(input.data);
		}

	input.data$ID <- paste0('svict-',1:nrow(input.data));

	# indicate chromosomes to keep
	keep.chrs <- c(1:22,'X','Y');
	if (any(grepl('chr',input.data$V1))) { keep.chrs <- paste0('chr',keep.chrs); }

	# filter data
	input.data <- input.data[which(input.data$V1 %in% keep.chrs),];
	#input.data <- input.data[which(input.data$V7 == 'PASS'),];

	# extract key features from INFO field
	input.data$sv_type <- sapply(
		input.data$V8,
		function(i) { 
			parts <- unlist(strsplit(as.character(i),';|:'));
			type <- parts[grepl('SVTYPE=', parts)];
			return(sub('SVTYPE=','',type));
			}
		);

	# extract mate chr/pos
	input.data[,c('chr2','pos2','orientation')] <- NA;

	if (any(input.data$sv_type == 'INS')) {
		ins.idx <- which(input.data$sv_type == 'INS');
		input.data[ins.idx,]$chr2 <- input.data[ins.idx,]$V1;
		input.data[ins.idx,]$pos2 <- sapply(input.data[ins.idx,]$V8,
			function(i) {
				parts <- unlist(strsplit(as.character(i),';'));
				mate.pos <- sub('END=','',parts[grepl('^END=', parts)]);
				return(mate.pos);
				}
			);

		input.data[ins.idx,]$orientation <- sapply(input.data[ins.idx,]$V8,
			function(i) {
				parts <- unlist(strsplit(as.character(i),';'));
				direction <- unlist(strsplit(parts[grepl('SVTYPE=', parts)],':'))[2];
				return(direction);
				}
			);
		}

	if (any(input.data$sv_type == 'DEL')) {
		del.idx <- which(input.data$sv_type == 'DEL');
		input.data[del.idx,]$chr2 <- input.data[del.idx,]$V1;
		input.data[del.idx,]$pos2 <- sapply(input.data[del.idx,]$V8,
			function(i) {
				parts <- unlist(strsplit(as.character(i),';'));
				mate.pos <- sub('END=','',parts[grepl('^END=', parts)]);
				return(mate.pos);
				}
			);
		}

	if (any(input.data$sv_type == 'DUP')) {
		dup.idx <- which(input.data$sv_type == 'DUP');
		input.data[dup.idx,]$chr2 <- input.data[dup.idx,]$V1;
		input.data[dup.idx,]$pos2 <- sapply(input.data[dup.idx,]$V8,
			function(i) {
				parts <- unlist(strsplit(as.character(i),';|:'));
				mate.pos <- sub('END=','',parts[grepl('^END=', parts)]);
				return(mate.pos);
				}
			);
		}

	if (any(input.data$sv_type == 'BND')) {
		bnd.idx <- which(input.data$sv_type == 'BND');
		input.data[bnd.idx,]$chr2 <- sapply(
			input.data[bnd.idx,]$V5,
			function(i) { 
				parts <- unlist(strsplit(as.character(i),'\\[|\\]|:'));
				return(parts[2]);
				}
			);

		input.data[bnd.idx,]$pos2 <- sapply(
			input.data[bnd.idx,]$V5,
			function(i) { 
				parts <- unlist(strsplit(as.character(i),'\\[|:|\\]'));
				return(parts[3]);
				}
			);

		input.data[bnd.idx,]$orientation <- sapply(
			input.data[bnd.idx,]$V5,
			function(i) { 
				parts <- unlist(strsplit(as.character(i),''));
				bnd <- parts[which(parts %in% c('[',']','C','G','A','T'))];
				return(paste(gsub('C|G|A|T','r',bnd),collapse = ''));
				}
			);
		}

	# filter data
	input.data <- input.data[which(input.data$chr2 %in% keep.chrs),];

	# format output data
	output.data <- data.frame(
		tracking_id = input.data$ID,
		event_type = c('deletion','duplication','inversion','insertion','translocation')[match(input.data$sv_type, c('DEL','DUP','INV','INS','BND'))],
		break1_chromosome = paste0('chr', input.data$V1),
		break1_position_start = input.data$V2,
		break1_position_end = input.data$V2,
		break1_orientation = c('L','R','L','R','L','R')[match(input.data$orientation, c('L','R','r[[','[[r','r]]',']]r'))],
		break1_strand = '?',
		break1_seq = 'None',
		break2_chromosome = paste0('chr', input.data$chr2),
		break2_position_start = input.data$pos2,
		break2_position_end = input.data$pos2,
		break2_orientation = c('L','R','R','R','L','L')[match(input.data$orientation, c('L','R','r[[','[[r','r]]',']]r'))],
		break2_strand = '?',
		break2_seq = 'None',
		opposing_strands = 'False',
		stranded = 'False',
		tools = 'svict'
		);

	if (any(grepl('chrchr', output.data$break1_chromosome))) {
		output.data$break1_chromosome <- gsub('chrchr','chr', output.data$break1_chromosome);
		output.data$break2_chromosome <- gsub('chrchr','chr', output.data$break2_chromosome);
		}

	if (any(output.data$event_type == 'translocation') & any(output.data$break1_chromosome == output.data$break2_chromosome)) {
		output.data <- output.data[-which(output.data$event_type == 'translocation' & 
			(output.data$break1_chromosome == output.data$break2_chromosome)),];
		}

	if (nrow(output.data) > 0) {
		print(paste0("collapsed to ", nrow(output.data), " rows"));
		}

	input.list[[file]] <- output.data;
	}


output.combined <- do.call(rbind, input.list);

# if input file is empty
if ((0 == total.variants) | (is.null(output.combined))) {

	write(
		c('#break1_chromosome','break1_position_start','break1_position_end','break2_chromosome','break2_position_start','break2_position_end'),
		ncolumns = 6,
		file = output,
		sep = '\t'
		);	
	
	} else {

	rownames(output.combined) <- 1:nrow(output.combined);
	output.combined <- output.combined[!is.na(output.combined[,1]),];

	colnames(output.combined)[1] <- '#tracking_id';

	if (any(is.na(output.combined$break1_orientation) | is.na(output.combined$break2_orientation))) {
		output.combined <- output.combined[,c(1,2,3,4,5,9,10,11,17)];
		}

	write.table(
		output.combined,
		file = output,
		row.names = FALSE,
		col.names = TRUE,
		quote = FALSE,
		sep = '\t'
		);
	}
