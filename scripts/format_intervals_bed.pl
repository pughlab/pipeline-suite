#!/usr/bin/env perl
#### format_intervals_bed.pl #######################################################################
use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($intervals_file, $reference, $dictionary, $help);

GetOptions(
	'h|help'	=> \$help,
	'b|bed=s'	=> \$intervals_file,
	'r|reference=s'	=> \$reference
	);

if ($help) {
        my $help_msg = join("\n",
                "Options:",
                "\t--help|-h\tPrint this help message",
                "\t--bed|-b\tpath to file intervals bed file",
		"\t--reference|-r\tpath to reference fasta"
		);

	print "$help_msg\n";
        exit;
        }

$dictionary = $reference;
$dictionary =~ s/\.fa/\.dict/;

### FUNCTIONS ######################################################################################
# format command to convert intervals.bed to picard-style intervals.list
sub get_format_intervals_command {
	my %args = (
		input_bed	=> undef,
		picard_out	=> undef,
		@_
		);

	my $format_command .= "\n\n" . join(' ',
		'java -jar $picard_dir/picard.jar BedToIntervalList',
		'I=' . $args{input_bed},
		'SD=' . $dictionary,
		'O=' . $args{picard_out}
		);

	return($format_command);
	}

### RUN ############################################################################################
# read it in
open(my $INPUT, '<', $intervals_file) or die "Cannot open '$intervals_file' file\n";

# extract regions
my @regions;
while (my $line = <$INPUT>) {
	chomp($line);

	next if ($line =~ m/#/);

	my @parts = split /\t/, $line;

	# add padding
	my $chr = $parts[0];
	my $start = $parts[1] - 50;
	my $end = $parts[2] + 50;

	push @regions, join("\t", $chr, $start, $end) . "\n";
	}

close $INPUT;

# now write to new file
my $output_file = $intervals_file;
$output_file =~ s/\.bed/_padding100bp.bed/;

open (my $fh_output, '>', $output_file) or die "Could not open $output_file\n";

foreach my $interval (@regions) {
	print $fh_output $interval;
	}

close $fh_output;

# make sure there is a bgzipped/indexed version as well
my $compress_command = join("\n",
	'module load tabix',
	"bgzip -c -f $intervals_file > $intervals_file.gz",
	"tabix -p bed $intervals_file.gz",
	"bgzip -c -f $output_file > $output_file.gz",
	"tabix -p bed $output_file.gz"
	);

`$compress_command`;

# produce an accompanying picard-style interval list
my $picard_intervals = $intervals_file;
$picard_intervals =~ s/\.bed/\.interval_list/;

my $picard_command = 'module load picard';
$picard_command .= "\n" . get_format_intervals_command(
	input_bed	=> $output_file,
	picard_out	=> $picard_intervals
	);

`$picard_command`;

exit;
