#!/usr/bin/env perl
### methyldackel.pl ################################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use Getopt::Std;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(any first);
use YAML qw(LoadFile);
use IO::Handle;

my $cwd = dirname(__FILE__);
require "$cwd/utilities.pl";

# define some global variables
our ($reference) = undef;

####################################################################################################
# version       author		comment
# 1.0		sprokopec       script to run MethylDackel on EMSeq

### USAGE ##########################################################################################
# methyldackel.pl -t tool_config.yaml -d data_config.yaml -o /path/to/output/dir -c slurm --remove --dry_run
#
# where:
#	-t (tool.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data.yaml) contains sample information (YAML file containing paths to bwameth-aligned BAMs)
#	-o (/path/to/output/dir) indicates tool-specific output directory
#	-c indicates hpc driver (ie, slurm)
#	--remove indicates that intermediates will be removed
#	--dry_run indicates that this is a dry run

### DEFINE SUBROUTINES #############################################################################
# format command to calculate mbias
sub get_mbias_command {
	my %args = (
		bam		=> undef,
		output_stem	=> undef,
		@_
		);

	my $methylation_cmd = join(' ',
		'MethylDackel mbias',
		'--txt', # to output results to txt format
		$reference,
		$args{bam},
		$args{output_stem},
		'>', $args{output_stem} . '.txt'
		);
	
	return($methylation_cmd);
	}

# format command to run MethylDackel extract
sub get_extract_command {
	my %args = (
		bam		=> undef,
		output_stem	=> undef,
		n_cpus		=> undef,
		non_cpg		=> undef,
		@_
		);

	my $non_cpg = ('Y' eq $args{non_cpg}) ? '--CHH --CHG' : '';

	my $methylation_cmd = "mbias=''";
	$methylation_cmd .= "\n\n" . join(' ',
		'MethylDackel extract',
		'-@', $args{n_cpus},
		$non_cpg,
		'$mbias',
		'-o', $args{output_stem},
		$reference,
		$args{bam}
		);

	return($methylation_cmd);
	}

### MAIN ###########################################################################################
sub main {
	my %args = (
		tool_config		=> undef,
		data_config		=> undef,
		output_directory	=> undef,
		hpc_driver		=> undef,
		dry_run			=> undef,
		no_wait			=> undef,
		@_
		);

	my $tool_config = $args{tool_config};
	my $data_config = $args{data_config};

	### PREAMBLE ######################################################################################
	unless($args{dry_run}) {
		print "Initiating MethylDackel pipeline...\n";
		}

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'methylation');
	my $date = strftime "%F", localtime;

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_METHYLDACKEL_pipeline.log');

	# create a file to hold job metrics
	my (@files, $run_count, $outfile, $touch_exit_status);
	unless ($args{dry_run}) {
		# initiate a file to hold job metrics
		opendir(LOGFILES, $log_directory) or die "Cannot open $log_directory";
		@files = grep { /slurm_job_metrics/ } readdir(LOGFILES);
		$run_count = scalar(@files) + 1;
		closedir(LOGFILES);

		$outfile = $log_directory . '/slurm_job_metrics_' . $run_count . '.out';
		$touch_exit_status = system("touch $outfile");
		if (0 != $touch_exit_status) { Carp::croak("Cannot touch file $outfile"); }

		$log_file = join('/', $log_directory, 'run_METHYLDACKEL_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running MethylDackel pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n    Reference used: $tool_data->{reference}";

	$reference = $tool_data->{reference};

	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	# set tools and versions
	my $methyldackel = 'MethylDackel/' . $tool_data->{methyldackel_version};
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};
	my $r_version	= 'R/' . $tool_data->{r_version};

	# get user-specified tool parameters
	my $parameters = $tool_data->{methyldackel}->{parameters};

	# get optional HPC group
	my $hpc_group = defined($tool_data->{hpc_group}) ? "-A $tool_data->{hpc_group}" : undef;

	### RUN ###########################################################################################
	my ($run_script, $run_id, $link, $should_run_final, $prep_run_id);
	my (@all_jobs);

	# get sample data
	my $smp_data = LoadFile($data_config);

	unless($args{dry_run}) {
		print "Processing " . scalar(keys %{$smp_data}) . " patients.\n";
		}

	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient";

		# find bams
		my @normal_ids = keys %{$smp_data->{$patient}->{'normal'}};
		my @tumour_ids = keys %{$smp_data->{$patient}->{'tumour'}};

		my @sample_ids = @tumour_ids;
		push @sample_ids, @normal_ids;
		@sample_ids = sort(@sample_ids);

		# create some directories
		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my $link_directory = join('/', $patient_directory, 'bam_links');
		unless(-e $link_directory) { make_path($link_directory); }

		# create some symlinks
		foreach my $normal (@normal_ids) {
			my @tmp = split /\//, $smp_data->{$patient}->{normal}->{$normal};
			$link = join('/', $link_directory, $tmp[-1]);
			symlink($smp_data->{$patient}->{normal}->{$normal}, $link);
			}
		foreach my $tumour (@tumour_ids) {
			my @tmp = split /\//, $smp_data->{$patient}->{tumour}->{$tumour};
			$link = join('/', $link_directory, $tmp[-1]);
			symlink($smp_data->{$patient}->{tumour}->{$tumour}, $link);
			}

		my (@final_outputs, @patient_jobs, $cleanup_cmd);

		# for each tumour sample
		foreach my $sample (@sample_ids) {

			# if there are any samples to run, we will run the final combine job
			$should_run_final = 1;

			# what type of sample is this?
			my $type = 'tumour';
			if ( (any { $_ =~ m/$sample/ } @normal_ids) ) {
				$type = 'normal';
				print $log "\n  NORMAL: $sample\n";
				} else {
				print $log "\n  TUMOUR: $sample\n";
				}

			my $mbias_directory = join('/', $patient_directory, 'mbias');
			unless(-e $mbias_directory) { make_path($mbias_directory); }

			$run_id = '';

			my $output_stem = join('/', $patient_directory, $sample);

			# start with mbias command
			# note: that this information is not capturable and must be input manually 
			# if desired!
			my $mbias_command = get_mbias_command(
				bam		=> $smp_data->{$patient}->{$type}->{$sample},
				output_stem	=> join('/', $mbias_directory, $sample . '_mbias')
				);

			my $mbias_output = join('/', $mbias_directory, $sample . '_mbias.COMPLETE');
	
			$mbias_command .= "\n\n" . join(' ',
				'echo methyldackel mbias complete',
				'>', $mbias_output
				);

			# check if this should be run
			if ('Y' eq missing_file($mbias_output)) {

				# record command (in log directory) and then run job
				print $log "  >> Submitting job for MethylDackel MBIAS...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_methyldackel_mbias_' . $sample,
					cmd	=> $mbias_command,
					modules	=> [$methyldackel],
					max_time	=> $parameters->{mbias}->{time},
					mem		=> $parameters->{mbias}->{mem},
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'run_methyldackel_mbias_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @all_jobs, $run_id;
				push @patient_jobs, $run_id;
				} else {
				print $log "  >> Skipping MethylDackel MBIAS because this has already been completed!\n";
				}

			# extract methylation levels
			if (!defined($parameters->{extract}->{n_cpus})) {
				$parameters->{extract}->{n_cpus} = 1;
				}
			
			my $extract_command = get_extract_command(
				bam		=> $smp_data->{$patient}->{$type}->{$sample},
				output_stem	=> $output_stem,
				non_cpg		=> $parameters->{extract}->{non_cpg},
				n_cpus		=> $parameters->{extract}->{n_cpus}
				);

			my $methyldackel_output = join('/',
				$patient_directory, $sample . '_methyldackel.COMPLETE');

			$extract_command .= "\n\n" . join(' ',
				'echo methyldackel extraction complete',
				'>', $methyldackel_output
				);

			# check if this should be run
			if ('Y' eq missing_file($methyldackel_output)) {

				# record command (in log directory) and then run job
				print $log "  >> Submitting job for MethylDackel EXTRACT...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_methyldackel_extract_' . $sample,
					cmd	=> $extract_command,
					modules	=> [$methyldackel],
					max_time	=> $parameters->{extract}->{time},
					mem		=> $parameters->{extract}->{mem},
					cpus_per_task	=> $parameters->{extract}->{n_cpus},
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'run_methyldackel_extract_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @all_jobs, $run_id;
				push @patient_jobs, $run_id;
				} else {
				print $log "  >> Skipping MethylDackel EXTRACT because this has already been completed!\n";
				}

			push @final_outputs, $output_stem . '_CpG.bedGraph';
			if ('Y' eq $parameters->{extract}->{non_cpg}) {
				push @final_outputs, $output_stem . '_CHG.bedGraph';
				push @final_outputs, $output_stem . '_CHH.bedGraph';
				}
			}

		print $log "\nFINAL OUTPUT:\n  " . join("\n  ", @final_outputs) . "\n";
		print $log "---\n";
		}

	# collate results
	my $collect_output = join(' ',
		"Rscript $cwd/collect_methyldackel_output.R",
		'-d', $output_directory,
		'-p', $tool_data->{project_name},
		'-t', $tool_data->{targets_bed},
		'-r', $tool_data->{ref_type},
		'-s', $data_config
		);

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'combine_methyldackel_output',
		cmd	=> $collect_output,
		modules	=> [$r_version],
		dependencies	=> join(':', @all_jobs),
		mem		=> '16G',
		max_time	=> '24:00:00',
		hpc_driver	=> $args{hpc_driver},
		extra_args	=> [$hpc_group]
		);

	$run_id = submit_job(
		jobname		=> 'combine_methyldackel_output',
		shell_command	=> $run_script,
		hpc_driver	=> $args{hpc_driver},
		dry_run		=> $args{dry_run},
		log_file	=> $log
		);

	push @all_jobs, $run_id;

	# if this is not a dry run OR there are jobs to assess (run or resumed with jobs submitted) then
	# collect job metrics (exit status, mem, run time)
	unless ( ($args{dry_run}) || (scalar(@all_jobs) == 0) ) {

		# collect job stats
		my $collect_metrics = collect_job_stats(
			job_ids		=> join(',', @all_jobs),
			outfile		=> $outfile,
			hpc_driver	=> $args{hpc_driver}
			);

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'output_job_metrics_' . $run_count,
			cmd	=> $collect_metrics,
			dependencies	=> join(':', @all_jobs),
			mem		=> '256M',
			hpc_driver	=> $args{hpc_driver},
			kill_on_error	=> 0,
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'output_job_metrics',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @all_jobs, $run_id;

		# do some logging
		print "Number of jobs submitted: " . scalar(@all_jobs) . "\n";

		my $n_queued = `squeue -r | wc -l`;
		print "Total number of jobs in queue: " . $n_queued . "\n";

		# wait until it finishes
		unless ($args{no_wait}) {
			check_final_status(job_id => $run_id);
			}
		}

	# finish up
	print $log "\nProgramming terminated successfully.\n\n";
	close $log;
	}

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($tool_config, $data_config, $output_directory);
my $hpc_driver = 'slurm';
my ($dry_run, $help, $no_wait);

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	'd|data=s'	=> \$data_config,
	't|tool=s'	=> \$tool_config,
	'o|out_dir=s'	=> \$output_directory,
	'c|cluster=s'	=> \$hpc_driver,
	'dry-run'	=> \$dry_run,
	'no-wait'	=> \$no_wait
	);

if ($help) {
	my $help_msg = join("\n",
		"Options:",
		"\t--help|-h\tPrint this help message",
		"\t--data|-d\t<string> data config (yaml format)",
		"\t--tool|-t\t<string> tool config (yaml format)",
		"\t--out_dir|-o\t<string> path to output directory",
		"\t--cluster|-c\t<string> cluster scheduler (default: slurm)",
		"\t--dry-run\t<boolean> should jobs be submitted? (default: false)",
		"\t--no-wait\t<boolean> should we exit after job submission (true) or wait until all jobs have completed (false)? (default: false)"
		);

	print "$help_msg\n";
	exit;
	}

# do some quick error checks to confirm valid arguments	
if (!defined($tool_config)) { die("No tool config file defined; please provide -t | --tool (ie, tool_config.yaml)"); }
if (!defined($data_config)) { die("No data config file defined; please provide -d | --data (ie, sample_config.yaml)"); }
if (!defined($output_directory)) { die("No output directory defined; please provide -o | --out_dir"); }

main(
	tool_config		=> $tool_config,
	data_config		=> $data_config,
	output_directory	=> $output_directory,
	hpc_driver		=> $hpc_driver,
	dry_run			=> $dry_run,
	no_wait			=> $no_wait
	);
