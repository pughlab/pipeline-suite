### create_fastq_yaml.pl ##########################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(any);
use YAML qw(DumpFile);
use Data::Dumper;

####################################################################################################
# version       author          comment
# 1.0		sprokopec	create config file in yaml format
# 2.0		sprokopec	default to TGL input format

### USAGE ##########################################################################################
# create_fastq_yaml.pl -i INPUT.samples -d DATA_DIR -o /path/to/OUTPUT_CONFIG.yaml -t dna|rna
#
# where INPUT.samples is a tab separated file containing sample info for use as YAML entries
# 	header = <Patient.ID\tSample.ID\tType> Sample.ID is also used to grab sample-specific files
#
# FASTQ files are expected to be paired end, with filenames ending with R1.fastq.gz and R2.fastq.gz

### GETOPTS PLUS ERROR CHECKING AND DEFAULT VALUES #################################################
# declare variables
my $sample_info;
my $data_directory;
my $output_config;

# indicate the data type (DNA or RNA; minor differences to output format)
my $data_type = 'dna';

# read in command line arguments
GetOptions(
	'i|input=s'		=> \$sample_info,
	'd|data_dir=s'		=> \$data_directory,
	't|data_type:s'		=> \$data_type,
	'o|output_config=s'	=> \$output_config
	);

# check for and remove trailing / from data dir
$data_directory =~ s/\/$//;

### HANDLING FILES #################################################################################
# find all fastq files
opendir(RAWFILES, $data_directory) or die "Cannot open '$data_directory' !";
my @fastqfiles = grep {/fastq.gz/} readdir(RAWFILES);
closedir(RAWFILES);

# open sample file
open(my $INPUT, '<', $sample_info) or die "Cannot open '$sample_info' file\n";
my $header = <$INPUT>;

my $smp_data = {};

# process each sample in SAMPLES
while (my $line = <$INPUT>) {

	chomp($line);
	
	# split line info
	my @file_info = split(/\t/, $line);
	my $subject = $file_info[0];
	my $sample = $file_info[1];
	my $type = $file_info[2];

	$smp_data->{$subject}->{$sample}->{type} = $type;

	# find sample fastq files
	# this assumes that the sample ID is present in the fastq filename
	my @subset = grep { /$sample/ } @fastqfiles;
	@subset = sort(@subset);
	my $size = scalar @subset;

	# Extract library and lane names
	# in general, the filename will be in the format libraryID_laneID_R1.fastq.gz
	my @parts;
	my ($lib,$lane,$read_dir);

	foreach my $fastq ( @subset ) {

		# a 'standard' library ID will end with the sequencing type
		# where WG = wgs, EX = whole-exome, WT = rna-seq
		if (grep /WG/, $fastq) {
			@parts = split(/WG_/, $fastq);
			$lib = $parts[0] . 'WG';
			$lane = $parts[1];
			$lane =~ s/_R1.fastq.gz|_R2.fastq.gz//;
			} elsif (grep /EX/, $fastq) {
			@parts = split(/EX_/, $fastq);
			$lib = $parts[0] . 'EX';
			$lane = $parts[1];
			$lane =~ s/_R1.fastq.gz|_R2.fastq.gz//;
			} elsif (grep /WT/, $fastq) {
			@parts = split(/WT_/, $fastq);
			$lib = $parts[0] . 'WT';
			$lane = $parts[1];
			$lane =~ s/_R1.fastq.gz|_R2.fastq.gz//;
			} else {
			# for a 'non-standard' case, let's assume the libraryID is the sampleID
			@parts = split(/_|\./, $fastq);
			$lib = $sample;

			# if non-unique lane identifier is used, it will generally be 
			# in the format L# (ie, L001)
			if (any { $_ =~ m/L00/ } @parts) {
				my $lane_idx = first_index { $_ eq 'L00' } @parts;
				$lane = $parts[$lane_idx];
				} else {
				$lane = 'lane_name';
				}
			}

		$read_dir = 'R1';
		if (grep /R2.fastq/, $fastq) { $read_dir = 'R2'; }

		if ('dna' eq $data_type) {
			$smp_data->{$subject}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{$read_dir} = $data_directory . '/' . $fastq;
			} else {
			$smp_data->{$subject}->{$sample}->{runlanes}->{$lane}->{$read_dir} = $data_directory . '/' . $fastq;
			}
		}
	}

close $INPUT;

local $YAML::Indent = 4;
DumpFile($output_config, $smp_data);


exit;
