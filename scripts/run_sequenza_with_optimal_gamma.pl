#!/usr/bin/env perl
### run_sequenza_with_optimal_gamma.pl #############################################################
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
use List::Util qw(any);
use File::Find;

use IO::Handle;

my $cwd = dirname($0);
require "$cwd/utilities.pl";

# define some global variables
our ($reference) = undef;

####################################################################################################
# version       author		comment
# 1.0		sprokopec	run sequenza with gamma tuning

### DEFINE SUBROUTINES #############################################################################
# find files recursively
sub _get_files {
	my ($dirs, $exten) = @_;

	my @files;
	my $want = sub {
		-e && /\Q$exten\E$/ && push @files, $File::Find::name
		};

	find($want, $dirs);
	return(@files);
	}

# run sequenza: create seqz file
sub create_seqz_command {
	my %args = (
		snp		=> undef,
		cnv		=> undef,
		out_dir		=> undef,
		ref_type	=> undef,
		@_
		);

	my $sequenza_command = join(' ',
		"Rscript $cwd/sequenza/prepare_seqz.R",
		'--snp_file', $args{snp},
		'--cnv_file', $args{cnv},
		'--out_dir', $args{out_dir},
		'--assembly', $args{ref_type}
		);

	return($sequenza_command);
	}

# run sequenza: gamma tuning
sub get_gamma_tuning_command {
	my %args = (
		seqz		=> undef,
		out_dir		=> undef,
		ref_type	=> undef,
		@_
		);

	my $sequenza_command = 'for i in {0..20..5} {30..150..10} {200..500..100} {1000..3000..500}; do';

	$sequenza_command .= "\n  " . join(' ',
		"Rscript $cwd/sequenza/test_gamma.R",
		'--seqz_file', $args{seqz},
		'--out_dir', $args{out_dir},
		'--assembly', $args{ref_type},
		'--gamma $i'
		);

	$sequenza_command .= "\ndone";

	return($sequenza_command);
	}

# run sequenza with optimal gamma
sub get_final_sequenza_command {
	my %args = (
		seqz		=> undef,
		out_dir		=> undef,
		cancer_type	=> undef,
		ref_type	=> undef,
		ploidy_priors	=> undef,
		@_
		);

	if (!defined($args{cancer_type})) { $args{cancer_type} = 'all'; }

	my $sequenza_command = join(' ',
		"Rscript $cwd/sequenza/run_optimized_sequenza.R",
		'--seqz_file', $args{seqz},
		'--out_dir', $args{out_dir},
		'--assembly', $args{ref_type},
		'--cancer_type', $args{cancer_type},
		'--ploidy_priors', $args{ploidy_priors}
		);

	return($sequenza_command);
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

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'varscan');

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs/SEQUENZA');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_Sequenza_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_Sequenza_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running Sequenza tuning pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n  Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	# get user-specified tool parameters
	my $parameters = $tool_data->{varscan}->{parameters};

	### RUN ###########################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	# find .snp and .copynumber (or .copynumber.called)
	my @extensions = qw(VarScan__exome.snp VarScan.snp VarScan.copynumber VarScan.copynumber.called);

	# initialize objects
	my ($run_script, $run_id);
	my @all_jobs;

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";

		# find bams
		my @tumour_ids = sort keys %{$smp_data->{$patient}->{'tumour'}};
		my @normal_ids = sort keys %{$smp_data->{$patient}->{'normal'}};

		if (scalar(@tumour_ids) == 0) {
			print $log "\n>> No tumour BAM provided, skipping patient.\n";
			next;
			} elsif (scalar(@normal_ids) == 0) {
			print $log "\n>> No normal BAM provided, skipping patient.\n";
			next;
			}

		# create some directories
		my $patient_directory = join('/', $output_directory, $patient);

		foreach my $tumour (@tumour_ids) {

			my $sample_directory = join('/', $patient_directory, $tumour);
			my $sequenza_directory = join('/', $sample_directory, 'Sequenza');
			my $tuning_directory = join('/', $sequenza_directory, 'output_tuning');
			my $optimized_directory = join('/', $sequenza_directory, 'optimized');

			my @snp_files = _get_files($sequenza_directory, '.snp');
			my @cnv_files = _get_files($sequenza_directory, '.copynumber.called');

			my $seqz_cmd = create_seqz_command(
				snp		=> $snp_files[0],
				cnv		=> $cnv_files[0],
				out_dir		=> $sequenza_directory,
				ref_type	=> $tool_data->{ref_type}
				);

			my $seqz_file = $snp_files[0];
			$seqz_file =~ s/.snp/.seqz/;

			$run_id = '';

			# check if this should be run
			if ('Y' eq missing_file($seqz_file)) {

				# record command (in log directory) and then run job
				print $log "Submitting job for Sequenza (SEQZ)...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_sequenza_prep_' . $tumour,
					cmd	=> $seqz_cmd,
					modules	=> ['R/3.3.0'],
					max_time	=> $parameters->{sequenza}->{time},
					mem		=> $parameters->{sequenza}->{mem},
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id = submit_job(
					jobname		=> 'run_sequenza_prep_' . $tumour,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @all_jobs, $run_id;
				} else {
				print $log "Skipping Sequenza (SEQZ) because this has already been completed!\n";
				}

			my $tuning_cmd = get_gamma_tuning_command(
				seqz		=> $seqz_file,
				out_dir		=> $sequenza_directory,
				ref_type	=> $tool_data->{ref_type}
				);

			my $tuning_status = join('/', $tuning_directory, 'gamma_tuning_status');

			$tuning_cmd .= "\n\n" . 'echo "gamma tuning complete" > ' . $tuning_status;

			# set time/memory requirements
			my ($time_req, $mem_req);
			if ('exome' eq $tool_data->{seq_type}) {
				$time_req = '24:00:00';
				$mem_req = '4G';
				} elsif ('wgs' eq $tool_data->{seq_type}) {
				$time_req = '72:00:00';
				$mem_req = '8G';
				}

			# check if this should be run
			if ('Y' eq missing_file($tuning_status)) {

				# record command (in log directory) and then run job
				print $log "Submitting job for Sequenza gamma tuning...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_sequenza_gamma_tuning_' . $tumour,
					cmd	=> $tuning_cmd,
					modules	=> ['R/3.6.1'],
					dependencies	=> $run_id,
					max_time	=> $time_req,
					mem		=> $mem_req,
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id = submit_job(
					jobname		=> 'run_sequenza_gamma_tuning_' . $tumour,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @all_jobs, $run_id;
				} else {
				print $log "Skipping Sequenza gamma tuning because this has already been completed!\n";
				}

			my $run_sequenza_cmd = get_final_sequenza_command(
				seqz		=> $seqz_file,
				out_dir		=> $sequenza_directory,
				ref_type	=> $tool_data->{ref_type},
				cancer_type	=> $parameters->{sequenza}->{cancer_type_prior},
				ploidy_priors	=> $parameters->{sequenza}->{ploidy_priors}
				);

			my $final_calls = join('/', $optimized_directory, basename($seqz_file));
			$final_calls =~ s/.seqz/_Total_CN.seg/;

			# check if this should be run
			if ('Y' eq missing_file($final_calls)) {

				# record command (in log directory) and then run job
				print $log "Submitting job for Sequenza (OPTIMIZED)...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_sequenza_optimized_' . $tumour,
					cmd	=> $run_sequenza_cmd,
					modules	=> ['R/3.6.1'],
					dependencies	=> $run_id,
					max_time	=> $parameters->{sequenza}->{time},
					mem		=> $parameters->{sequenza}->{mem},
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id = submit_job(
					jobname		=> 'run_sequenza_optimized_' . $tumour,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @all_jobs, $run_id;
				} else {
				print $log "Skipping Sequenza (OPTIMIZED) because this has already been completed!\n";
				}
			}
		}

	# collate results
	my $collect_output = join(' ',
		"Rscript $cwd/collect_sequenza_output.R",
		'-d', $output_directory,
		'-p', $tool_data->{project_name},
		'-g', $tool_data->{gtf}
		);

	if ( ('exome' eq $tool_data->{seq_type}) && (defined($tool_data->{intervals_bed})) ) {
		$collect_output .= " -t $tool_data->{intervals_bed}";
		}

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'combine_sequenza_segment_calls',
		cmd	=> $collect_output,
		modules	=> ['R/3.6.1'],
		dependencies	=> join(':', @all_jobs),
		mem		=> '6G',
		max_time	=> '24:00:00',
		hpc_driver	=> $args{hpc_driver}
		);

	$run_id = submit_job(
		jobname		=> 'combine_sequenza_segment_calls',
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
			job_ids	=> join(',', @all_jobs),
			outfile	=> $outfile
			);

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'output_job_metrics_' . $run_count,
			cmd	=> $collect_metrics,
			dependencies	=> join(':', @all_jobs),
			mem		=> '256M',
			hpc_driver	=> $args{hpc_driver}
			);

		$run_id = submit_job(
			jobname		=> 'output_job_metrics',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		# wait until it finishes
		unless ($args{no_wait}) {

			my $complete = 0;
			my $timeouts = 0;

			while (!$complete && $timeouts < 20 ) {
				sleep(30);
				my $status = `sacct --format='State' -j $run_id`;

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
					die("Final SEQUENZA accounting job: $run_id finished with errors.");
					}
				}
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