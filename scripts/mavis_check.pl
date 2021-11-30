#!/usr/bin/env perl
#### mavis_check.pl ################################################################################
use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(any all);

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($jobid_file, $help);

GetOptions(
	'h|help'	=> \$help,
	'j|jobs=s'	=> \$jobid_file
	);

if ($help) {
        my $help_msg = join("\n",
                "Options:",
                "\t--help|-h\tPrint this help message",
                "\t--jobs|-j\tpath to file containing job ids"
		);

	print "$help_msg\n";
        exit;
        }

### RUN ############################################################################################
# read it in
open(my $INPUT, '<', $jobid_file) or die "Cannot open '$jobid_file' file\n";

# extract ids
my @jobs;
while (my $line = <$INPUT>) {
	chomp($line);
	unless ('None' eq $line) { push @jobs, $line; }
	}
close $INPUT;

@jobs = sort(@jobs);
my $job_list = join(',', @jobs);
my $final_job = $jobs[-1];

# wait until it finishes
my $complete = 0;
my $timeouts = 0;

sleep(60);

while (!$complete && $timeouts < 20 ) {

	my $status = `sacct --format='State' -j $job_list`;

	# if we run into a server connection error (happens rarely with sacct)
	# increment timeouts (if we continue to repeatedly timeout, we will exit)
	if ($status =~ m/Connection timed out/) {
		$timeouts++;
		next;
		}

	my @job_statuses = split("\n", $status);
	@job_statuses = @job_statuses[2..$#job_statuses];

	print Dumper \@job_statuses;

	# if final job has finished successfully:
	if ( all { $_ =~ m/COMPLETED/ } @job_statuses ) {
		$complete = 1;
		print "All jobs completed successfully.";
		}

	# if none of the above, we will exit with an error
	elsif ( any { $_ =~ m/FAILED|TIMEOUT|CANCELLED/ } @job_statuses ) {
		die("Final MAVIS job: $final_job finished with errors.");
		}

	# if the job is still pending or running, try again in a bit
	# but also reset timeouts, because we only care about consecutive timeouts
	else {
		$timeouts = 0;
		}

	sleep(60);
	}

exit;
