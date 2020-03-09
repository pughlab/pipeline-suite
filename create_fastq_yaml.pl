### create_fastq_yaml.pl ##########################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);

####################################################################################################
# version       author          comment
# 1.0		sprokopec	create config file in yaml format

### USAGE ##########################################################################################
# create_fastq_yaml.pl -i INPUT.samples -d DATA_DIR -o /path/to/OUTPUT_CONFIG.yaml -t dna|rna
#
# where INPUT.samples is a tab separated file containing sample info for use as YAML entries
# 	header = <Patient.ID\tSample.ID>   Sample.ID is also used to grab sample-specific files
#
# Tumour/Normal are identified based on pattern matching (where any SK/BC/A = Normal)
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
# open output file
open (my $fh, '>', $output_config) or die "Cannot open '$output_config' !";
print $fh "---\n";

# identify pattern to use
my $pattern = 'R1.fastq.gz$';
	
# find all fastq files
opendir(RAWFILES, $data_directory) or die "Cannot open '$data_directory' !";
my @fastqfiles = grep {/$pattern/} readdir(RAWFILES);
closedir(RAWFILES);

# open sample file
open(my $INPUT, '<', $sample_info) or die "Cannot open '$sample_info' file\n";
my $header = <$INPUT>;

my @subject_names = ();

# process each sample in SAMPLES
while (my $line = <$INPUT>) {

	chomp($line);
	
	# split line info
	my @file_info = split(/\t/, $line);
	my $subject = $file_info[0];
	my $sample = $file_info[1];

	# find sample fastq files
	my @subset = grep {/$sample/} @fastqfiles;
	@subset = sort(@subset);
	my $size = scalar @subset;	

	my $type;
	if ($sample =~ m/BC|SK|A/) { $type = 'normal'; } else { $type = 'tumour'; }		# indicates if these are tumour or normal

	# if the subject has already been written (if multiple samples per patient)
	if (grep /$subject/, @subject_names) {
		print "$subject already exists; continuing with next sample for this patient\n";
		}
	else {
		push @subject_names, $subject;
		print $fh "$subject:\n";
		}
	
	print $fh "    $sample:\n";
	print $fh "        type: $type\n";

	if ('dna' eq $data_type) {
		print $fh "        libraries:\n";
		print $fh "            $sample:\n";  # for now we are assuming that only 1 library was used
		print $fh "                runlanes:\n";
	} else {
		print $fh "        runlanes:\n";
		}

	for (my $i = 1; $i <= $size; $i++) {

		my $r1 = join('/', $data_directory, $subset[$i-1]);
		my $r2 = join('/', $data_directory, $subset[$i-1]);
		$r2 =~ s/R1.fastq.gz/R2.fastq.gz/;

		my $lane = $subset[$i-1];
		$lane =~ s/$sample//;
		$lane =~ s/^_//;
		$lane =~ s/.R1.fastq.gz//;

		if ('dna' eq $data_type) {
			print $fh "                    $lane:\n";
			print $fh "                        fastq:\n";
			print $fh "                            R1: $r1\n";
			print $fh "                            R2: $r2\n";
		} else {
			print $fh "           $lane:\n";
			print $fh "               R1: $r1\n";
			print $fh "               R2: $r2\n";
			}
		}
	}

close $INPUT;
close $fh;
exit;
