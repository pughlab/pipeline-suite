### collect_snv_output.R ###########################################################################
# find, collate and format output from SNV caller

### FUNCTIONS ######################################################################################
# function to generate a standardized filename
generate.filename <- function(project.stem, file.core, extension, include.date = TRUE) {

	# build up the filename
	file.name <- paste(project.stem, file.core, sep = '_');
	file.name <- paste(file.name, extension, sep = '.');

	if (include.date) {
		file.name <- paste(Sys.Date(), file.name, sep = '_');
		}

	return(file.name);
	}

# function to write session profile to file
save.session.profile <- function(file.name) {

	# open the file
	sink(file = file.name, split = FALSE);

	# write memory usage to file
	cat('### MEMORY USAGE ###############################################################');
	print(proc.time());

	# write sessionInfo to file
	cat("\n### SESSION INFO ###############################################################");
	print(sessionInfo());

	# close the file
	sink();

	}

### PREPARE SESSION ################################################################################
# import libraries
library(argparse);
library(plyr);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-d', '--directory', type = 'character', help = 'path to data directory',
	default = getwd());
parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-t', '--t_depth', type = 'integer', 
	help = 'minimum tumour depth required (rna-seq only)', default = 20);

arguments <- parser$parse_args();

# what's the date?
date <- Sys.Date();

setwd(arguments$directory);


### MAF CODING
# ensure proper coding of maf columns; this should be consistent, however may change with versions
# required primarily for germline calls (very small; protects agains all 'T' allele being interpreted
# as logical TRUE)
maf.classes <- rep('character',132);
maf.classes[c(2,6,7,40:45,58,60,77:85,100:107,112:122,124:132)] <- 'numeric';

# determine type of experiment that was run
is.rnaseq <- grepl('RNASeq', arguments$directory);
is.exome <- grepl('Exome|WEX|WXS|ctDNA', arguments$directory);
is.wgs <- grepl('WGS', arguments$directory);
is.germline <- grepl('CPSR', arguments$directory);

# if this is RNA-Seq, only keep variants in coding regions
classes.to.keep <- c('RNA','Missense_Mutation','Splice_Region','Splice_Site','In_Frame_Del','In_Frame_Ins',
	'Frame_Shift_Del','Frame_Shift_Ins','Nonsense_Mutation','Nonstop_Mutation','Translation_Start_Site');

### MAIN ###########################################################################################
# find results files
maf.files <- list.files(pattern = '.maf$|.maf.gz', recursive = TRUE);

# read them in and store them
samples <- c();

maf.data <- list();

for (i in 1:length(maf.files)) {

	file <- maf.files[i];

	# extract sample ID
	smp <- unlist(strsplit(file,'\\/'))[2];
	samples <- unique(c(samples, smp));

	# read in data
	if (is.germline) {
		tmp <- read.delim(file, comment.char = "#", colClasses = maf.classes);
		} else {
		tmp <- read.delim(file, comment.char = "#", as.is = TRUE);
		}

	## do some filtering
	# for some reason WGS:SomaticSniper output has FILTER = '.'
	if (!(is.wgs & grepl('SomaticSniper', getwd()))) {
		tmp <- tmp[which(tmp$FILTER == 'PASS'),];
		}

	# add a depth filter (some 'PASS' calls (MuTect2 only) have depth = 0 
	# 	despite depth > 0 reported by other callers at the same positions)
	# we will only remove these low/no depth calls if the normal is also low/absent
	idx <- which(tmp$t_depth < 5 & (tmp$n_depth < 5 | is.na(tmp$n_depth)));
	if (length(idx) > 0) {
		tmp <- tmp[-idx,];
		}

	# if RNA-Seq, add a further depth filter
	if (is.rnaseq) {
		tmp <- tmp[which(tmp$t_depth > arguments$t_depth),];
		}

	# if this is Pindel, remove extra fields (if present)
	if (grepl('Pindel', file)) {
		tmp <- tmp[!grepl('ONP', tmp$Variant_Type),];
		tmp <- tmp[,setdiff(colnames(tmp),c('Fusion','Method','Frame','CONSENSUS'))];

		# and remove germline hits if present (should only ever provide somatic calls)
		tmp <- tmp[which(tmp$t_alt_count > 1 & tmp$t_depth >= 5 & (tmp$n_depth >= 5 | is.na(tmp$n_depth))),];

		if (!all(is.na(tmp$n_depth))) {
			n_vaf <- tmp$n_alt_count / tmp$n_depth;
			tmp <- tmp[which(n_vaf < 0.1 | is.na(n_vaf)),];
			}
		}

	maf.data[[i]] <- tmp;

	rm(tmp);
	}

# save full (combined) maf data to file
full.maf.data <- do.call(rbind, maf.data);

norm.idx <- which(full.maf.data$Matched_Norm_Sample_Barcode == 'NORMAL');
if (length(norm.idx) > 0) {
	full.maf.data[norm.idx,c('Matched_Norm_Sample_Barcode','Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2')] <- NA;
	}

if (!is.rnaseq & !is.germline) {
	norm.idx <- which(!is.na(full.maf.data$Matched_Norm_Sample_Barcode));
	if (length(norm.idx) > 0) { full.maf.data[norm.idx,]$Mutation_Status <- 'somatic'; }
	}

if (!is.rnaseq & is.germline) {
	norm.idx <- which(!is.na(full.maf.data$Matched_Norm_Sample_Barcode));
	if (length(norm.idx) > 0) { full.maf.data[norm.idx,]$Mutation_Status <- 'germline'; }
	}

# if this is germline data, make sure n_vaf > 0
# [there are cases with wonky alleles (genotype 3/3) that vcf2maf converts to reference?]
if (is.germline) {
	t_vaf <- full.maf.data$t_alt_count / full.maf.data$t_depth;
	n_vaf <- full.maf.data$n_alt_count / full.maf.data$n_depth;
	keep.idx <- which(
		(!is.na(full.maf.data$Matched_Norm_Sample_Barcode) & n_vaf > 0) |
		(is.na(full.maf.data$Matched_Norm_Sample_Barcode) & t_vaf > 0)
		);
	full.maf.data <- full.maf.data[keep.idx,];
	}

# replace NAs with blanks (required for cBioportal upload)
for (field in c('t_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count')) {
	if (any(is.na(full.maf.data[,field]))) {
		full.maf.data[is.na(full.maf.data[,field]),field] <- '';
		}
	}

# fill in missing Entrez_Gene_Ids where possible
if (any(full.maf.data$Entrez_Gene_Id == 0 & full.maf.data$SYMBOL_SOURCE == 'EntrezGene')) {
	full.maf.data$Entrez_Gene_Id <- as.character(full.maf.data$Entrez_Gene_Id);
	full.maf.data$Gene <- as.character(full.maf.data$Gene);
	idx <- which(full.maf.data$Entrez_Gene_Id == 0 & full.maf.data$SYMBOL_SOURCE == 'EntrezGene');
	full.maf.data[idx,]$Entrez_Gene_Id <- full.maf.data[idx,]$Gene;
	if (any(grepl('ENSG', full.maf.data$Entrez_Gene_Id))) {
		full.maf.data[grepl('ENSG', full.maf.data$Entrez_Gene_Id),]$Entrez_Gene_Id <- 0;
		}
	full.maf.data$Entrez_Gene_Id <- as.numeric(full.maf.data$Entrez_Gene_Id);
	}

# remove unnecessary columns
colnames(full.maf.data) <- gsub('vcf','variant',colnames(full.maf.data));
exclude.field <- grepl('variant', colnames(full.maf.data));


### fix some weird cases
# some cases where vcf2maf doesn't fill in Entrez_Gene_Id
if (any(full.maf.data$Hugo_Symbol == 'ATXN7' & full.maf.data$Chromosome == 'chr3')) {
        full.maf.data[which(full.maf.data$Hugo_Symbol == 'ATXN7' & full.maf.data$Chromosome == 'chr3'),]$Entrez_Gene_Id <- 6314;
        }

if (any(full.maf.data$Hugo_Symbol == 'EZHIP' & full.maf.data$Chromosome == 'chrX')) {
        full.maf.data[which(full.maf.data$Hugo_Symbol == 'EZHIP' & full.maf.data$Chromosome == 'chrX'),]$Entrez_Gene_Id <- 340602;
        }

if (any(full.maf.data$Hugo_Symbol == 'MAP3K14' & full.maf.data$Chromosome == 'chr17')) {
        full.maf.data[which(full.maf.data$Hugo_Symbol == 'MAP3K14' & full.maf.data$Chromosome == 'chr17'),]$Entrez_Gene_Id <- 9020;
        }

# some cases where vcf2maf fills in the wrong GeneID
# ie, EIF4A2 is given GeneID for SNORA81 (which is within EIF4A2) but not necessarily on the variant position
if (any(full.maf.data$Hugo_Symbol == 'EIF4A2' & full.maf.data$Chromosome == 'chr3')) {
        full.maf.data[which(full.maf.data$Hugo_Symbol == 'EIF4A2' & full.maf.data$Chromosome == 'chr3'),]$Entrez_Gene_Id <- 1974;
        }

if (any(full.maf.data$Hugo_Symbol == 'MEF2B' & full.maf.data$Chromosome == 'chr19')) {
        full.maf.data[which(full.maf.data$Hugo_Symbol == 'MEF2B' & full.maf.data$Chromosome == 'chr19'),]$Entrez_Gene_Id <- 100271849;
        }

if (any(full.maf.data$Hugo_Symbol == 'PRSS1' & full.maf.data$Chromosome == 'chr7')) {
        full.maf.data[which(full.maf.data$Hugo_Symbol == 'PRSS1' & full.maf.data$Chromosome == 'chr7'),]$Entrez_Gene_Id <- 5644;
        }

if (any(full.maf.data$Hugo_Symbol == 'SMG1' & full.maf.data$Chromosome == 'chr16')) {
        full.maf.data[which(full.maf.data$Hugo_Symbol == 'SMG1' & full.maf.data$Chromosome == 'chr16'),]$Entrez_Gene_Id <- 23049;
        }

if (any(full.maf.data$Hugo_Symbol == 'SRSF2' & full.maf.data$Chromosome == 'chr17')) {
        full.maf.data[which(full.maf.data$Hugo_Symbol == 'SRSF2' & full.maf.data$Chromosome == 'chr17'),]$Entrez_Gene_Id <- 6427;
        }

# vcf2maf mis-annotates a common downstream/regulartory_region TERC mutation (rs2293607) 
# as a downstream variant in ACTRT3; therefore we will manually update this for plotting
if (any(full.maf.data$dbSNP_RS == 'rs2293607' & full.maf.data$SYMBOL == 'ACTRT3')) {
        full.maf.data[which(full.maf.data$dbSNP_RS == 'rs2293607' & full.maf.data$SYMBOL == 'ACTRT3'),]$Hugo_Symbol <- 'TERC';
        full.maf.data[which(full.maf.data$dbSNP_RS == 'rs2293607' & full.maf.data$SYMBOL == 'ACTRT3'),]$Entrez_Gene_Id <- 7012;
        }

# there are most likely more wierd cases and I will add them as I find them
###

# save to file
write.table(
	full.maf.data[,!exclude.field],
	file = generate.filename(arguments$project, 'mutations_for_cbioportal', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CollectVariantData','SessionProfile','txt'));
