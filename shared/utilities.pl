#!/usr/bin/env perl
### utilities.pl ###################################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use File::Basename;
use File::Path qw(make_path);
use List::Util 'any';

### SHARED SUBROUTINES #############################################################################
# function to do standard error checking on input config files
sub error_checking {
	my %args = (
		tool_data	=> undef,
		pipeline	=> undef,
		@_
		);

	my $tool_data = $args{tool_data};
	my $pipeline  = $args{pipeline};

	# check to see if common arguments are supplied and/or are correct
	# for del_intermediates, if not Y then default to N
	if ( (!defined($tool_data->{del_intermediate})) || ('Y' ne $tool_data->{del_intermediate}) && ('N' ne $tool_data->{del_intermediate}) ) {
		print "Option del_intermediate is neither Y or N, defaulting to N\n";
		$tool_data->{del_intermediate} = 'N';
		}

	# for create_output_yaml, if not Y then default to N
	if ( ('Y' ne $tool_data->{create_output_yaml}) || (!defined($tool_data->{create_output_yaml})) ) {
		$tool_data->{create_output_yaml} = 'N';
		}

	# is ref_type either hg38 or hg19?
	if ( ('hg38' ne $tool_data->{ref_type}) && ('hg19' ne $tool_data->{ref_type}) ) {
		die("Unrecognized ref_type; must be one of hg19 or hg38.");
		}
	
	# check if this is a dry run or not
	if ((!defined($tool_data->{dry_run})) || ('Y' ne $tool_data->{dry_run})) {
		$tool_data->{dry_run} = 'N';
		}

	# check for compatible HPC driver; if not found, change dry_run to Y
	$tool_data->{HPC_driver} = lc $tool_data->{HPC_driver};
	my @compatible_drivers = qw(slurm);
	if ((!any { /$tool_data->{HPC_driver}/ } @compatible_drivers ) && ('N' eq $tool_data->{dry_run})) {
		print "Unrecognized HPC driver requested: setting dry_run to Y, jobs will not be submitted but commands will be written to file.\n";
		$tool_data->{dry_run} = 'Y';
		}

	my $is_ref_valid;

	# check for pipeline specific options
	if ('bwa' eq $pipeline) {
		if ('bwamem' ne $tool_data->{aligner}) { die("This pipeline is currently only compatible with BWA-MEM!"); }
		if (('Y' ne $tool_data->{mark_dup}) & ('N' ne $tool_data->{mark_dup})) {
			print "Option mark_dup must be either Y or N, defaulting to N\n";
			$tool_data->{mark_dup} = 'N';
			}

		if (!defined($tool_data->{reference}))  { die("Must supply path to reference genome!"); }
		$is_ref_valid = validate_ref(
			reference	=> $tool_data->{reference},
			pipeline	=> $pipeline,
			exts		=> [qw(.fa .fa.amb .fa.ann .fa.bwt .fa.fai .fa.pac .fa.sa)]
			);
		}

	if ('gatk' eq $pipeline) {

		if (!defined($tool_data->{intervals_bed})) {
			print "WARNING: no target intervals provided; if this is exome data, please provide the target regions!\n";
			}

		if (!defined($tool_data->{reference})) { die("Must supply path to reference genome!"); }
		$is_ref_valid = validate_ref(
			reference	=> $tool_data->{reference},
			pipeline	=> $pipeline,
			exts		=> [qw(.fa .dict .fa.fai)]
			);
		}

	if ('star' eq $pipeline) {
		if (!defined($tool_data->{reference_dir}))  { die("Must supply path to reference genome directory!"); }

		if (('Y' ne $tool_data->{mark_dup}) && ('N' ne $tool_data->{mark_dup})) {
			print "Option mark_dup is neither Y or N, defaulting to N\n";
			$tool_data->{mark_dup} = 'N';
			}
		}

	if ('star-fusion' eq $pipeline) {
		if (!defined($tool_data->{reference_dir}))  { die("Must supply path to reference genome directory!"); }
		}

	if ('rsem' eq $pipeline) {
		if (!defined($tool_data->{reference_dir}))  { die("Must supply path to reference genome directory!"); }

		my @strand_options = qw(none forward reverse);
		if (!defined($tool_data->{strandedness})) {
			print "No option provided for 'strandedness'; setting to default: none.\n";
			$tool_data->{strandedness} = 'none';
			}

		if ( !any /$tool_data->{strandedness}/ } @strand_options ) {
			print "Unrecognized 'strandedness' option: must be one of none, forward or reverse! Setting to default: none.\n";
			$tool_data->{strandedness} = 'none';
			}
		}

	return($tool_data);
	}

# function to ensure all necessary reference files are available
sub validate_ref {
	my %args = (
		reference	=> undef,
		pipeline	=> undef,
		exts		=> undef,
		@_
		);

	my $ref_file_base;

	if ( ('bwa' eq $args{pipeline}) || ('gatk' eq $args{pipeline}) ) {
		$ref_file_base = ($args{reference} =~ s/\.fa$//r);
		}

	foreach my $ext (@{$args{exts}}) {
		unless (-e $ref_file_base . $ext) { die("Missing reference file: " . $ref_file_base . $ext); }
		}

	# if nothing fails, return Y
	return('Y');
	}

# function to set output/resume directory
sub set_output_path {
	my %args = (
		tool_data => undef,
		@_
		);

	my $tool_data = $args{tool_data};

	my $path = $tool_data->{resume_dir};
	my $resume = 'N';

	if (defined($path)) {
		$resume = 'Y';
		$path =~ s/\/$//;
		} else {
		$tool_data->{output_dir} =~ s/\/$//;
		$path = join('/', $tool_data->{output_dir}, join('_', $tool_data->{date}, $tool_data->{tool}, $tool_data->{tool_version}));
		unless(-e $path) { make_path($path); }
		}

	my $log = join('/', $path, 'logs');
	unless(-e $log) { make_path($log); }

	return($resume, $path, $log);
	}

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
		hpc_driver	=> 'slurm',
		@_
		);

	my $jobdir = $args{shell_command};
	$jobdir =~ s/\/script.sh//;

	# make a directory for error/log output 
	my $job_log_dir = join('/', $jobdir, $args{hpc_driver});
	unless(-e $job_log_dir) { make_path($job_log_dir); }

	my $job_command;
	if ('slurm' eq $args{hpc_driver}) {
		$job_command = join(' ',
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
		} else {
		die("Unrecognized HPC driver: currently only compatible with slurm");
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

		print "Job number $job_id submitted.\n";
		}

	# print 1 extra blank line in the log to separate steps
	print "\n";

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
