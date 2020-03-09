#!/usr/bin/env perl
### utilities.pl ###########################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use File::Basename;
use File::Path qw(make_path);

# function to check for missing/empty files
sub missing_file {
	my @bad = grep { ! -e $_ || -s $_ <= 1}
		map { ref($_) eq 'HASH' ? values %{$_} : $_ }
		@_;

	if (@bad) {
		my $s = scalar(@bad) > 1 ? 's' : '';
		print "Missing or empty file$s: ", join(',', @bad) . "\n";
		return('Y');
		}
	return('N');
	}

# save commands to shell script in log directory
sub write_script {
	my %args = (
		log_dir	=> undef,
		name	=> undef,
		cmd	=> undef,
		modules => [],
		@_
		);

	my $cmd_log_dir = join('/', $args{log_dir}, $args{name});
	unless(-e $cmd_log_dir) { make_path($cmd_log_dir); }

	my @modules_list;
	for (my $i=0; $i < scalar (@{$args{modules}}); $i++) {
		$modules_list[$i] = "module load $args{modules}->[$i]";
		}

	my $modules_to_load = join("\n", @modules_list);

	my $script = join('/', $cmd_log_dir, 'script.sh');
	open (my $fh_script, '>', $script) or Carp::croak("Cannot open file $script: $!");
	print $fh_script "#!/bin/bash\n"; 
	print $fh_script $modules_to_load . "\n";
	print $fh_script $args{cmd};
	close($fh_script);

	return($script);
	}

# execute commands/submit jobs
sub submit_job {
	my %args = (
		jobname		=> undef,
		shell_command	=> undef,
		dependencies	=> undef,
		max_time	=> undef,
		mem		=> undef,
		cpus_per_task	=> 1,
		dry_run		=> 'N',
		@_
		);

	my $jobdir = $args{shell_command};
	$jobdir =~ s/\/script.sh//;

	# make a 
	my $job_log_dir = join('/', $jobdir, 'slurm');
	unless(-e $job_log_dir) { make_path($job_log_dir); }

	my $job_command = join(' ',
		'sbatch',
		'--job-name', $args{jobname},
		#'--mail-type=BEGIN,END,FAIL', # doesn't seem to work?
		'-p all',
		'-t', $args{max_time},
		'--mem', $args{mem},
		'-D', $job_log_dir,
		'-c', $args{cpus_per_task}
		);

	if ('' ne $args{dependencies}) {
		$job_command = $job_command . ' --dependency=afterok:' . $args{dependencies};
		$job_command .= ' --kill-on-invalid-dep=yes';
		}

	$job_command = join(' ', $job_command, $args{shell_command});
	print "\nCOMMAND IS: " . $job_command . "\n";

	my $job_id = $args{jobname};
	if ('N' eq $args{dry_run}) {
		$job_id = `$job_command`;
		chop $job_id;
		if ($job_id =~ m/Submitted batch job ([0-9]+)/) {
			$job_id = $1;
		} else {
			die("Failed to submit job with command $job_command");
			}

		print "Job number $job_id submitted.\n\n";
		}

	return($job_id);
	}

# command to extract job status / metrics
sub collect_job_stats {
	my %args = (
		job_ids	=> undef,
		outfile	=> undef,
		@_
		);

	my $sacct_command = join(' ',
		'sacct -P --delimiter=","',
		'--format="User,JobID,Start,End,AllocCPUS,Elapsed,CPUTime,MaxRSS,ExitCode"',
		'-j', $args{job_ids},
		'>', $args{outfile} . ';',
		'sed -i "s/,/\t/g"', $args{outfile}
		);

	return($sacct_command);
	}
1;
