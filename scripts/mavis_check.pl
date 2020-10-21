#!/usr/bin/env perl
#### mavis_check.pl ################################################################################
use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;

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
	push @jobs, $line;
	}
close $INPUT;

@jobs = sort(@jobs);
my $final_job = $jobs[-1];

# wait until it finishes
my $complete = 0;
my $timeouts = 0;

while (!$complete && $timeouts < 20 ) {
	sleep(30);
	my $status = `sacct --format='State' -j $final_job`;

	# if final job has finished successfully:
	if ($status =~ m/COMPLETED/s) { $complete = 1; }

	# if we run into a server connection error (happens rarely with sacct)
	# increment timeouts (if we continue to repeatedly timeout, we will exit)
	elsif ($status =~ m/Connection timed out/) {
		$timeouts++;
		}

	# if the job is still pending or running, try again in a bit
	# but also reset timeouts, because we only care about consecutive timeouts
	elsif ($status =~ m/PENDING|RUNNING/) {
		$timeouts = 0;
		}

	# if none of the above, we will exit with an error
	else {
		die("Final MAVIS job: $final_job finished with errors.");
		}
	}

exit;
