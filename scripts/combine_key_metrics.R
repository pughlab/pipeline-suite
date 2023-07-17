# find files
fastqc_files <- list.files(pattern = 'metrics.txt');
md5_files <- list.files(pattern = 'md5$');

# read in and collect FASTQC metrics
metric.data <- list();

for (i in fastqc_files) {

	tmp <- read.delim(i, skip = 3, nrow = 7, header = FALSE, col.names = c('Measure','Value'));
	tmp$Filename <- tmp[1,2];
	metric.data[[length(metric.data) + 1]] <- tmp[4:7,];

	}

combined <- do.call(rbind, metric.data);

reshaped <- reshape(combined, direction = 'wide', idvar = 'Filename', timevar = 'Measure');
colnames(reshaped) <- gsub('Value.', '', colnames(reshaped));

# add in md5sum values
reshaped$md5sum <- NA;

for (i in md5_files) {

	tmp <- read.table(i, header = F, col.names = c('MD5','Filename'));
	if (nrow(tmp) > 0) {
		filepath <- tmp[1,2];
		parts <- unlist(strsplit(as.character(filepath),'/'));
		filename <- parts[length(parts)];

		reshaped[which(reshaped$Filename == filename),]$md5sum <- as.character(tmp[1,1]);
		}
	}

# write results to file
write.table(
	reshaped,
	file = 'fastqc_key_metrics.tsv',
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);
