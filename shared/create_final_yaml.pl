### create_final_yaml.pl ###########################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);

###################################################################################################
# version       author          comment
# 1.0           sprokopec       create config file in yaml format

### USAGE ##########################################################################################
# create_final_yaml.pl -d DATA_DIR -o /path/to/OUTPUT_CONFIG.yaml -p PATTERN
#
# where:
# 	DATA_DIR indicates top directory (ie, DATE_bwa_version)
# 	PATTERN indicates grep pattern (ie, markdup.bam$) for use with bwa output
#
# Tumour/Normal are identified based on pattern matching (where any SK/BC/A = Normal)

### GETOPTS PLUS ERROR CHECKING AND DEFAULT VALUES #################################################
# declare variables
my $data_directory;
my $output_config;

# identify pattern to use
my $pattern = '.bam$';

# read in command line arguments
GetOptions(
	'd|data_dir=s'		=> \$data_directory,
	'p|pattern=s'		=> \$pattern,
	'o|output_config=s'	=> \$output_config
	);

# check for and remove trailing / from data dir
$data_directory =~ s/\/$//;

### HANDLING FILES #################################################################################
# open output file
open (my $fh, '>', $output_config) or die "Cannot open '$output_config' !";
print $fh "---\n";

# find all patient/subject directories
opendir(SUBJECTS, $data_directory) or die "Cannot open '$data_directory' !";
my @subject_dirs = grep { !/logs|^\./ && -d "$data_directory/$_" } readdir(SUBJECTS);
closedir(SUBJECTS);

my @subject_names = ();

# process each patient in SUBJECTS
foreach my $subject (sort(@subject_dirs)) {

	my $subject_directory = join('/', $data_directory, $subject);

	# if the subject has already been written (if multiple samples per patient)
	if (grep /$subject/, @subject_names) {
		print "$subject already exists; continuing with next sample for this patient\n";
		}
	else {
		push @subject_names, $subject;
		print $fh "$subject:\n";
		}

	my @normals = ();
	my @tumours = ();

	# find sample directories
	opendir(OUTFILES, $subject_directory) or die "Cannot open '$subject_directory' !";
	my @out_files = grep {/$pattern/} readdir(OUTFILES);
	closedir(OUTFILES);

	my $sample;
	foreach my $file (sort(@out_files)) {

		my @sample_parts = split(/\_/, $file);
		$sample = $sample_parts[0];

		# if output comes from RSEM, trim off the file suffix
		$sample =~ s/.genes.results//;
		$sample =~ s/.isoforms.results//;

		if ($sample =~ m/BC|SK|A/) {
			push @normals, $sample . ": " . $subject_directory . "/" . $file;
			} else {
			push @tumours, $sample . ": " . $subject_directory . "/" . $file;
			}
		}

	if (scalar(@normals) > 0) {
		print $fh "    normal:\n";
		foreach my $normal (sort(@normals)) {
			print $fh "        $normal\n";
			}
		}

	if (scalar(@tumours) > 0) {
		print $fh "    tumour:\n";
		foreach my $tumour (sort(@tumours)) {
			print $fh "        $tumour\n";
			}
		}
	}

close $fh;

exit;
