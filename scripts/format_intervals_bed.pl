#!/usr/bin/env perl
#### format_intervals_bed.pl #######################################################################
use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($intervals_file, $help);

GetOptions(
	'h|help'	=> \$help,
	'b|bed=s'	=> \$intervals_file
	);

if ($help) {
        my $help_msg = join("\n",
                "Options:",
                "\t--help|-h\tPrint this help message",
                "\t--bed|-b\tpath to file intervals bed file"
		);

	print "$help_msg\n";
        exit;
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

exit;
