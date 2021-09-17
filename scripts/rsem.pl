#!/usr/bin/env perl
### rsem.pl ########################################################################################
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
use List::Util 'any';
use IO::Handle;

my $cwd = dirname(__FILE__);
require "$cwd/utilities.pl";

####################################################################################################
# version       author		comment
# 1.0		sprokopec       script to run RSEM on STAR-aligned RNASeq data
# 1.1		sprokopec	minor updates for compatibility with larger pipeline
# 1.2		sprokopec	added help message and cleaned up code
# 1.3           sprokopec       minor updates for tool config

### USAGE ##########################################################################################
# rsem.pl -t tool_config.yaml -d data_config.yaml -o /path/to/output/dir -c slurm --remove --dry_run
#
# where:
#	-t (tool_config.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data_config.yaml) contains sample information (YAML file containing paths to FASTQ files)
#	-o (/path/to/output/dir) indicates tool-specific output directory
#	-c indicates hpc driver (ie, slurm)
#	--remove remove intermediate files?
#	--dry_run is this a dry run?

### DEFINE SUBROUTINES #############################################################################
# format command to run RSEM
sub get_rsem_command {
	my %args = (
		input		=> undef,
		output_stem	=> undef,
		ref_dir		=> undef,
		tmp_dir		=> undef,
		strand		=> undef,
		@_
		);

	my $rsem_command = join(' ',
		'rsem-calculate-expression',
		'--paired-end --bam --estimate-rspd --output-genome-bam',
		'--temporary-folder', $args{tmp_dir},
		'-p 8'
		);

	if (defined($args{strand})) {
		$rsem_command = join(' ',
			$rsem_command,
			'--strandedness', $args{strand}
			);
		}

	$rsem_command = join(' ',
		$rsem_command,
		$args{input},
		$args{ref_dir},
		$args{output_stem}
		);

	return($rsem_command);
	}

### MAIN ###########################################################################################
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
		print "Initiating RSEM pipeline...\n";
		}

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'rsem');

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_RSEM_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_RSEM_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running RNASeq RSEM pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n    Reference used: $tool_data->{rsem_reference}";
	print $log "\n    Output directory: $output_directory";
	print $log "\n    Strandedness: $tool_data->{rsem}->{strandedness}";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	# set tools and versions
	my $rsem = 'rsem/' . $tool_data->{rsem_version};
	my $r_version = 'R/' . $tool_data->{r_version};

	# get user-specified tool parameters
	my $parameters = $tool_data->{rsem}->{parameters};

	### RUN ###########################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	unless($args{dry_run}) {
		print "Processing " . scalar(keys %{$smp_data}) . " patients.\n";
		}

	my ($run_script, $run_id, $raw_link);
	my @all_jobs;

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";

		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my @normal_ids = keys %{$smp_data->{$patient}->{'normal'}};
		my @tumour_ids = keys %{$smp_data->{$patient}->{'tumour'}};

		my @samples = @tumour_ids;
		if (scalar(@normal_ids) > 0) { push @samples, @normal_ids; }

		my (@patient_jobs, @final_outputs);
		my $cleanup_cmd = '';

		foreach my $sample (@samples) {

			print $log "  SAMPLE: $sample\n";

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			my $link_directory = join('/', $sample_directory, 'bam_links');
			unless(-e $link_directory) { make_path($link_directory); }

			my $tmp_directory = join('/', $sample_directory, 'TEMP');
			unless(-e $tmp_directory) { make_path($tmp_directory); }
			$cleanup_cmd .= "\nrm -rf $tmp_directory";

			my $type;
			if ( (any { $_ =~ m/$sample/ } @normal_ids) ) {
				$type = 'normal';
				} else {
				$type = 'tumour';
				}

			# because smp_data points to STAR-aligned, picard merged/markdup BAMs,
			# we need to first, get the parent directory
			my $data_dir = $smp_data->{$patient}->{$type}->{$sample};
			my @parts = split /\//, $data_dir;
			$data_dir =~ s/$parts[-1]//;
			$data_dir =~ s/\/$//;

			# now direct it to the specific sample of interest
			my $aligned_bam = join('/', $data_dir, $sample, 'Aligned.toTranscriptome.out.bam');

			print $log "    STAR aligned BAM: $aligned_bam\n\n";
			
			# create a symlink for the bam
			my @tmp = split /\//, $aligned_bam;
			$raw_link = join('/', $link_directory, $tmp[-1]);
			symlink($aligned_bam, $raw_link);

			# reset run_id for this sample
			$run_id = '';
			my @smp_jobs;

			## RUN RSEM
			my $rsem_cmd = "cd $sample_directory;\n";
			$rsem_cmd .= get_rsem_command( 
				input		=> $aligned_bam,
				output_stem	=> $sample,
				ref_dir		=> $tool_data->{rsem_reference},
				tmp_dir		=> $tmp_directory,
				strand		=> $tool_data->{rsem}->{strandedness}
				);

			my $genes_file = join('/', $sample_directory, $sample . '.genes.results');
			my $isoforms_file = join('/', $sample_directory, $sample . '.isoforms.results');

			# add files to cleanup
			$cleanup_cmd .= "\nrm " . join('/', $sample_directory, '*.bam');

			# check if this should be run
			if ('Y' eq missing_file($genes_file)) {

				# record command (in log directory) and then run job
				print $log "Submitting job for RSEM...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_RSEM_' . $sample,
					cmd	=> $rsem_cmd,
					modules	=> ['perl', $rsem],
					max_time	=> $parameters->{rsem}->{time},
					mem		=> $parameters->{rsem}->{mem},
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id = submit_job(
					jobname		=> 'run_RSEM_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;
				}
			else {
				print $log "Skipping rsem-calculate-expression because this has already been completed!\n";
				}

			push @final_outputs, $genes_file;
			push @final_outputs, $isoforms_file;
			}

		# run cleanup, once per patient
		if ($args{del_intermediates}) {

			if (scalar(@patient_jobs) == 0) {
				`$cleanup_cmd`;
				} else {

				print $log "Submitting job to clean up temporary/intermediate files...\n";

				# make sure final output exists before removing intermediate files!
				$cleanup_cmd = join("\n",
					"if [ -s " . join(" ] && [ -s ", @final_outputs) . " ]; then",
					$cleanup_cmd,
					"else",
					'echo "One or more FINAL OUTPUT FILES is missing; not removing intermediates"',
					"fi"
					);

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_cleanup_' . $patient,
					cmd	=> $cleanup_cmd,
					dependencies	=> join(':', @patient_jobs),
					mem		=> '256M',
					hpc_driver	=> $args{hpc_driver},
					kill_on_error	=> 0
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

		print $log "\nFINAL OUTPUT:\n" . join("\n  ", @final_outputs) . "\n";
		print $log "---\n";
		}

	# collect and combine results
	my $collect_results = join(' ',
		"Rscript $cwd/collect_rsem_output.R",
		'-d', $output_directory,
		'-p', $tool_data->{project_name}
		);

	if (defined($tool_data->{reference_gtf})) {
		$collect_results .= " -g $tool_data->{reference_gtf}";
		}

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'combine_and_format_results',
		cmd	=> $collect_results,
		modules		=> [$r_version],
		dependencies	=> join(':', @all_jobs),
		max_time	=> $parameters->{combine_results}->{time},
		mem		=> $parameters->{combine_results}->{mem},
		hpc_driver	=> $args{hpc_driver}
		);

	$run_id = submit_job(
		jobname		=> 'combine_and_format_results',
		shell_command	=> $run_script,
		hpc_driver	=> $args{hpc_driver},
		dry_run		=> $args{dry_run},
		log_file	=> $log
		);

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
			hpc_driver	=> $args{hpc_driver},
			kill_on_error	=> 0
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
					die("Final RSEM accounting job: $run_id finished with errors.");
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

# quick error checks to confirm valid arguments
if (!defined($tool_config)) { die("No tool config file defined; please provide -t | --tool (ie, tool_config.yaml)"); }
if (!defined($data_config)) { die("No data config file defined; please provide -d | --data (ie, sample_config.yaml)"); }
if (!defined($output_directory)) { die("No output directory defined; please provide -o | --out_dir"); }

# run it!
main(
	tool_config		=> $tool_config,
	data_config		=> $data_config,
	output_directory	=> $output_directory,
	hpc_driver		=> $hpc_driver,
	del_intermediates	=> $remove_junk,
	dry_run			=> $dry_run,
	no_wait			=> $no_wait
	);
