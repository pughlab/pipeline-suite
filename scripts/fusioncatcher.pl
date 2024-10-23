#!/usr/bin/env perl
### fusioncatcher.pl ###############################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use Getopt::Std;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use File::Path qw(make_path);
use YAML qw(LoadFile);
use IO::Handle;

my $cwd = dirname(__FILE__);
require "$cwd/utilities.pl";

####################################################################################################
# version	author	  	comment
# 1.0		sprokopec	script to run FusionCatcher on RNA-Seq data
# 1.1		sprokopec	minor updates for compatibility with larger pipeline	
# 1.2		sprokopec	added help message and cleaned up code
# 1.3		sprokopec	minor updates for tool config

### USAGE ##########################################################################################
# fusioncatcher.pl -t tool_config.yaml -d data_config.yaml -o /path/to/output/dir -c slurm --remove --dry_run
#
# where:
#	-t (tool_config.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data_config.yaml) contains sample information (YAML file containing paths to FASTQ files)
#	-o (/path/to/output/dir) indicates tool-specific output directory
#	-c indicates hpc driver (ie, slurm)
#	--remove remove intermediate files?
#	--dry_run is this a dry run?

### SUBROUTINES ####################################################################################
# format command to run FusionCatcher
sub get_fusion_command {
	my %args = (
		ref_dir		=> undef,
		input		=> undef,
		output_dir	=> undef,
		tmp_dir		=> undef,
		ref_type	=> 'hg38',
		java_mem	=> undef,
		@_
		);

	my $fusion_command = join(' ',
		'fusioncatcher.py',
		'-i', $args{input},
		'-o', $args{output_dir},
		'-d', $args{ref_dir},
		'--keep-viruses-alignments',
		'--limitOutSJcollapsed 2500000',
		'--limitSjdbInsertNsj 2000000',
		'--tmp=' . $args{tmp_dir}
		);

	if ('hg38' eq $args{ref_type}) {
		$fusion_command .= ' --skip-conversion-grch37';
		}

	if (defined($args{java_mem})) {
		$fusion_command .= ' --Xmx=' . $args{java_mem};
		}

	return($fusion_command);
	}

### MAIN ##########################################################################################
sub main {
	my %args = (
		tool_config		=> undef,
		data_config		=> undef,
		output_directory	=> undef,
		hpc_driver		=> undef,
		del_intermediates	=> undef,
		dry_run			=> undef,
		no_wait			=> undef,
		@_
		);

	my $tool_config = $args{tool_config};
	my $data_config = $args{data_config};

	### PREAMBLE ######################################################################################
	unless($args{dry_run}) {
		print "Initiating FusionCatcher pipeline...\n";
		}

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'fusioncatcher');

	# clean up reference_dir (aesthetic reasons only)
	$tool_data->{fusioncatcher_reference_dir} =~ s/\/$//;

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_FusionCatcher_pipeline.log');

	# create a file to hold job metrics
	my (@files, $run_count, $outfile, $touch_exit_status);
	unless ($args{dry_run}) {
		# initiate a file to hold job metrics (ensures that an existing file isn't overwritten by concurrent jobs)
		opendir(LOGFILES, $log_directory) or die "Cannot open $log_directory";
		@files = grep { /slurm_job_metrics/ } readdir(LOGFILES);
		$run_count = scalar(@files) + 1;
		closedir(LOGFILES);

		$outfile = $log_directory . '/slurm_job_metrics_' . $run_count . '.out';
		$touch_exit_status = system("touch $outfile");
		if (0 != $touch_exit_status) { Carp::croak("Cannot touch file $outfile"); }

		$log_file = join('/', $log_directory, 'run_FusionCatcher_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	# start logging
	print $log "---\n";
	print $log "Running FusionCatcher (+ViralAlignment) pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n    FusionCatcher reference directory: $tool_data->{fusioncatcher_reference_dir}";
	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	# set tools and versions
	my $fusioncatcher	= 'fusioncatcher/' . $tool_data->{fusioncatcher_version};
	my $perl		= 'perl/5.30.0';
	my $r_version		= 'R/' . $tool_data->{r_version};

	# get user-specified tool parameters
	my $parameters = $tool_data->{fusioncatcher}->{parameters};

	# get optional HPC group
	my $hpc_group = defined($tool_data->{hpc_group}) ? "-A $tool_data->{hpc_group}" : undef;

	### RUN ############################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	unless($args{dry_run}) {
		print "Processing " . scalar(keys %{$smp_data}) . " patients.\n";
		}

	my ($run_script, $run_id, $raw_link, $should_run_final);
	my @all_jobs;

	# process each patient in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";

		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my (@final_outputs, @patient_jobs);
		my $cleanup_cmd = '';

		# process each separate sample for this patient
		foreach my $sample (sort keys %{$smp_data->{$patient}}) {

			if ('normal' eq $smp_data->{$patient}->{$sample}->{type}) {
				print $log "\n>> SAMPLE: $sample is labelled as normal - will be skipped!\n";
				next;
				}

			# if there are any samples to run, we will run the final combine job
			$should_run_final = 1;

			print $log "  SAMPLE: $sample\n";

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			my $link_directory = join('/', $sample_directory, 'input_links');
			unless(-e $link_directory) { make_path($link_directory); }

			my $tmp_directory = join('/', $sample_directory, 'TEMP');
			unless(-e $tmp_directory) { make_path($tmp_directory); }
			$cleanup_cmd .= "\nrm -rf $tmp_directory";

			my @fastqs;
			my @libraries = keys %{$smp_data->{$patient}->{$sample}->{libraries}};
			foreach my $lib ( @libraries ) {
				print $log "\n    LIBRARY: $lib\n";
				my $lib_data = $smp_data->{$patient}->{$sample}->{libraries}->{$lib};
				my @lanes = keys %{$lib_data->{runlanes}};
				foreach my $lane ( sort(@lanes) ) {
					print $log "      LANE: $lane\n";
					my $r1 = $lib_data->{runlanes}->{$lane}->{fastq}->{R1};
					my $r2 = $lib_data->{runlanes}->{$lane}->{fastq}->{R2};

					print $log "        R1: $r1\n";
					print $log "        R2: $r2\n\n";

					# create a symlink for this file
					my @tmp = split /\//, $r1;
					$raw_link = join('/', $link_directory, $tmp[-1]);
					symlink($r1, $raw_link);

					@tmp = split /\//, $r2;
					$raw_link = join('/', $link_directory, $tmp[-1]);
					symlink($r2, $raw_link);

					push @fastqs, $r1;
					push @fastqs, $r2;
					}
				}

			# reset run_id for this sample
			$run_id = '';

			## run FusionCatcher on these fastqs
			my $fusion_output = join('/', $sample_directory, 'final-list_candidate-fusion-genes.txt');

			my $fusion_cmd = get_fusion_command(
				ref_dir		=> $tool_data->{fusioncatcher_reference_dir},
				input		=> join(',', @fastqs),
				output_dir	=> $sample_directory,
				tmp_dir		=> $tmp_directory,
				ref_type	=> $tool_data->{ref_type},
				java_mem	=> $parameters->{fusioncatcher}->{java_mem}
				);

			# check if this should be run
			if ('Y' eq missing_file($fusion_output)) {

				# record command (in log directory) and then run job
				print $log "  >> Submitting job to run FusionCatcher...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_fusionCatcher_' . $sample,
					cmd	=> $fusion_cmd,
					modules => [$fusioncatcher],
					max_time	=> $parameters->{fusioncatcher}->{time},
					mem		=> $parameters->{fusioncatcher}->{mem},
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'run_fusionCatcher_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;

				} else {
				print $log "  >> Skipping FusionCatcher because output already exists...\n";
				}

			# add output from STAR-Fusion to final_outputs
			push @final_outputs, $fusion_output;
			}

		# remove temporary directories (once per patient)
		if ($args{del_intermediates}) {

			if (scalar(@patient_jobs) == 0) {
				`$cleanup_cmd`;
				} else {

				print $log "Submitting job to clean up temporary/intermediate files...\n";

				# make sure final output exists before removing intermediate files!
				$cleanup_cmd = join("\n",
					"if [ -s " . join(" ] && [ -s ", @final_outputs) . " ]; then",
					"  $cleanup_cmd",
					"else",
					'  echo "One or more FINAL OUTPUT FILES is missing; not removing intermediates"',
					"fi"
					);

				$run_script = write_script(
					log_dir => $log_directory,
					name    => 'run_cleanup_' . $patient,
					cmd     => $cleanup_cmd,
					dependencies	=> join(':', @patient_jobs),
					mem		=> '256M',
					hpc_driver	=> $args{hpc_driver},
					kill_on_error	=> 0,
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'run_cleanup_' . $patient,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);
				}
			}

		print $log "\nFINAL OUTPUT:\n  " . join("\n  ", @final_outputs) . "\n";
		print $log "---\n";
		}

	# collect final results
	if ($should_run_final) {

		my $collect_results = join(' ',
			"Rscript $cwd/collect_fusioncatcher_output.R",
			'-d', $output_directory,
			'-p', $tool_data->{project_name}
			);

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'combine_and_format_results',
			cmd	=> $collect_results,
			modules		=> [$r_version],
			dependencies	=> join(':', @all_jobs),
			max_time	=> $parameters->{combine_results}->{time},
			mem		=> $parameters->{combine_results}->{mem},
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'combine_and_format_results',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @all_jobs, $run_id;
		}

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
my ($data_config, $tool_config, $output_directory);
my $hpc_driver = 'slurm';
my ($remove_junk, $dry_run, $help, $no_wait);

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	'd|data=s'	=> \$data_config,
	't|tool=s'	=> \$tool_config,
	'o|out_dir=s'	=> \$output_directory,
	'c|cluster=s'	=> \$hpc_driver,
	'remove'	=> \$remove_junk,
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
		"\t--remove\t<boolean> should intermediates be removed? (default: false)",
		"\t--dry-run\t<boolean> should jobs be submitted? (default: false)",
		"\t--no-wait\t<boolean> should we exit after job submission (true) or wait until all jobs have completed (false)? (default: false)"
		);

	print "$help_msg\n";
	exit;
	}

# do some quick error checks to ensure valid input
if (!defined($tool_config)) { die("No tool config file defined; please provide -t | --tool (ie, tool_config.yaml)"); }
if (!defined($data_config)) { die("No data config file defined; please provide -d | --data (ie, sample_config.yaml)"); }
if (!defined($output_directory)) { die("No output directory defined; please provide -o | --out_dir"); }

main(
	tool_config		=> $tool_config,
	data_config		=> $data_config,
	output_directory 	=> $output_directory,
	hpc_driver		=> $hpc_driver,
	del_intermediates	=> $remove_junk,
	dry_run			=> $dry_run,
	no_wait			=> $no_wait
	);
