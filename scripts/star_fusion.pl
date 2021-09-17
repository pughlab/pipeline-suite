#/!/usr/bin/env perl
### star_fusion.pl #################################################################################
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
# version       author	  	comment
# 1.0		sprokopec       script to run STAR-fusion on RNASeq data aligned with STAR
# 1.1		sprokopec	minor updates for compatibility with larger pipeline
# 1.2		sprokopec	added help message and cleaned up code
# 1.3           sprokopec       minor updates for tool config

### USAGE ##########################################################################################
# star_fusion.pl -t tool_config.yaml -c data_config.yaml -o /path/to/output/dir -c slurm --remove --dry_run
#
# where:
#	-t (tool_config.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data_config.yaml) contains sample information (YAML file containing paths to FASTQ files)
#	-o (/path/to/output/dir) indicates tool-specific output directory
#	-c indicates hpc driver (ie, slurm)
#	--remove remove intermediate files?
#	--dry_run is this a dry run?

### SUBROUTINES ####################################################################################
# format command to run STAR-Fusion
sub get_star_fusion_command {
	my %args = (
		fusion_call	=> undef,
		reference	=> undef,
		input		=> undef,
		output_dir	=> undef,
		tmp_dir		=> undef,
		fusion_inspect	=> undef,
		r1		=> undef,
		r2		=> undef,
		@_
		);

	if (!defined($args{fusion_call})) {
		$args{fusion_call} = 'STAR-Fusion';
		}

	my $star_command = join(' ',
		$args{fusion_call},
		'--genome_lib_dir', $args{reference},
		'--chimeric_junction', $args{input}, # Chimeric.out.junction
		'--output_dir', $args{output_dir},
		'--CPU 4',
		'--tmpdir', $args{tmp_dir}
		);

	if (defined($args{fusion_inspect})) {
		$star_command .= join(' ',
			' --FusionInspector', $args{fusion_inspect},
			'--left_fq', $args{r1},
			'--right_fq', $args{r2}
			);

		if ('inspect' eq $args{fusion_inspect}) {
			$star_command .= ' --extract_fusion_reads';
			}
		}

	return($star_command);
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
		print "Initiating STAR-Fusion pipeline...\n";
		}

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'star-fusion');

	# clean up reference_dir (aesthetic reasons only)
	$tool_data->{star_fusion_reference_dir} =~ s/\/$//;

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_STAR-Fusion_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_STAR-Fusion_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running STAR-fusion pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	if (defined($tool_data->{star_fusion_path})) {
		print $log "\n    STAR-fusion path: $tool_data->{star_fusion_path}";
		}
	print $log "\n    STAR-fusion reference directory: $tool_data->{star_fusion_reference_dir}";
	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	# set tools and versions
	my $star_fusion = 'STAR-Fusion/' . $tool_data->{star_fusion_version};
	if (defined($tool_data->{star_fusion_path})) {
		$star_fusion = '';
		}
	my $star	= 'STAR/' . $tool_data->{star_version};
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};
	my $perl	= 'perl/5.30.0';
	my $r_version	= 'R/' . $tool_data->{r_version};

	# get user-specified tool parameters
	my $parameters = $tool_data->{star_fusion}->{parameters};

	### RUN ############################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	unless($args{dry_run}) {
		print "Processing " . scalar(keys %{$smp_data}) . " patients.\n";
		}

	my ($run_script, $run_id, $raw_link);
	my @all_jobs;

	# process each patient in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";

		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my @normal_ids = keys %{$smp_data->{$patient}->{'normal'}};
		my @tumour_ids = keys %{$smp_data->{$patient}->{'tumour'}};

		my @samples = @tumour_ids;
		if (scalar(@normal_ids) > 0) { push @samples, @normal_ids; }

		my (@final_outputs, @patient_jobs);
		my $cleanup_cmd = '';

		# process each separate sample for this patient
		foreach my $sample (@samples) {

			print $log "  SAMPLE: $sample\n";

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			my $link_directory = join('/', $sample_directory, 'input_links');
			unless(-e $link_directory) { make_path($link_directory); }

			my $tmp_directory = join('/', $sample_directory, 'TEMP');
			unless(-e $tmp_directory) { make_path($tmp_directory); }
			$cleanup_cmd .= "\n  rm -rf $tmp_directory";

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
			my $junctions_file = join('/', $data_dir, $sample, 'Chimeric.out.junction');
			print $log "    STAR chimeric junctions: $junctions_file\n\n";

			# create a symlink for this file
			my @tmp = split /\//, $junctions_file;
			$raw_link = join('/', $link_directory, 'Chimeric.out.junction');
			symlink($junctions_file, $raw_link);

			# reset run_id for this sample
			$run_id = '';
			my @smp_jobs;

			# if we will be running FusionInspect validate, we will need the fastq files as well
			my ($r1, $r2);
			my (@r1_fastq_files, @r2_fastq_files);
			if (defined($parameters->{FusionInspect})) {

				$data_dir = join('/', $data_dir, $sample, 'fastq_links');
				opendir(RAWFILES, $data_dir) or die "Cannot open $data_dir !";
				my @dir_files = readdir(RAWFILES);
				@r1_fastq_files = grep {/R1.fastq.gz|R1_001.fastq.gz/} @dir_files;
				@r2_fastq_files = grep {/R2.fastq.gz|R2_001.fastq.gz/} @dir_files;
				closedir(RAWFILES);

				# if there were multiple lanes, we will need to combine them
				$r1 = join('/', $data_dir, $r1_fastq_files[0]);
				$r2 = join('/', $data_dir, $r2_fastq_files[0]);

				if (scalar(@r1_fastq_files) > 1) {

					$r1 = join('/', $link_directory, $sample . '_combined_lanes.R1.fastq');
					$r2 = join('/', $link_directory, $sample . '_combined_lanes.R2.fastq');

					my $cat_fq_cmd = "cd $data_dir";
					$cat_fq_cmd .= "\necho Concatenating fastq files...";

					$cat_fq_cmd .= join(' ',
						"\nzcat",
						join(' ', sort(@r1_fastq_files)),
						'>', $r1
						);

					$cat_fq_cmd .= join(' ',
						"\nzcat",
						join(' ', sort(@r2_fastq_files)),
						'>', $r2
						);

					$cat_fq_cmd .= "\necho Finished concatenating files...";

					my $gzip_fq_cmd = "echo Compressing concatenated fastq files...";

					$gzip_fq_cmd .= "\ngzip " . $r1;
					$gzip_fq_cmd .= "\ngzip " . $r2;
					$gzip_fq_cmd .= "\necho Finished compressing fastq files. Now ready for FusionInspector";

					$r1 .= '.gz';
					$r2 .= '.gz';

					if ('Y' eq missing_file($r2)) {

						# record command (in log directory) and then run job
						print $log "Submitting job to prepare input for FusionInspector...\n";

						$run_script = write_script(
							log_dir	=> $log_directory,
							name	=> 'run_prepare_fastq_for_FusionInspector_' . $sample,
							cmd	=> $cat_fq_cmd . "\n" . $gzip_fq_cmd,
							max_time	=> '08:00:00',
							mem		=> '1G',
							hpc_driver	=> $args{hpc_driver}
							);

						$run_id = submit_job(
							jobname		=> 'run_prepare_fastq_for_FusionInspector_' . $sample,
							shell_command	=> $run_script,
							hpc_driver	=> $args{hpc_driver},
							dry_run		=> $args{dry_run},
							log_file	=> $log
							);

						push @smp_jobs, $run_id;
						push @all_jobs, $run_id;
						}
					}
				}

			## run STAR-Fusion on these junctions
			my $fusion_output = join('/',
				$sample_directory,
				'star-fusion.fusion_predictions.abridged.tsv'
				);

			my $fusion_cmd = get_star_fusion_command(
				fusion_call	=> $tool_data->{star_fusion_path},
				reference	=> $tool_data->{star_fusion_reference_dir},
				input		=> $junctions_file,
				output_dir	=> $sample_directory,
				tmp_dir		=> $tmp_directory,
				fusion_inspect	=> $parameters->{FusionInspect},
				r1		=> $r1,
				r2		=> $r2
				);

			# check if this should be run
			if ('Y' eq missing_file($fusion_output)) {

				# record command (in log directory) and then run job
				print $log "Submitting job to run STAR-Fusion...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_STAR-fusion_' . $sample,
					cmd	=> $fusion_cmd,
					modules => [$star, $perl, $samtools, 'tabix', 'python', $star_fusion],
					# requires python igv-reports
					dependencies	=> $run_id,
					max_time	=> $parameters->{star_fusion}->{time},
					mem		=> $parameters->{star_fusion}->{mem},
					hpc_driver	=> $args{hpc_driver},
					cpus_per_task	=> 4
					);

				$run_id = submit_job(
					jobname		=> 'run_STAR-Fusion_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @smp_jobs, $run_id;
				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;

				# add some stuff to the cleanup
				$cleanup_cmd .= "\n  rm $sample_directory/pipeliner*";
				$cleanup_cmd .= "\n  rm -rf $sample_directory/_starF_checkpoints";
				$cleanup_cmd .= "\n  rm $sample_directory/*.ok";

				push @patient_jobs, $run_id;
				}
			else {
				print $log "Skipping STAR-Fusion because output already exists...\n";
				}

			# add output from STAR-Fusion to final_outputs
			push @final_outputs, $fusion_output;
			$cleanup_cmd .= "\n  rm -rf " . join('/', 
				$sample_directory,
				'star-fusion.preliminary',
				'star-fusion.filter.intermediates_dir'
				);

			$cleanup_cmd .= "\n  tar -czvf $sample_directory/star-fusion.preliminary.tar.gz " .
				"$sample_directory/star-fusion.preliminary/ --remove-files";

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
					$cleanup_cmd,
					"else",
					'echo "One or more FINAL OUTPUT FILES is missing; not removing intermediates"',
					"fi"
					);

				$run_script = write_script(
					log_dir => $log_directory,
					name    => 'run_cleanup_' . $patient,
					cmd     => $cleanup_cmd,
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
		"Rscript $cwd/collect_star-fusion_output.R",
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
			job_ids => join(',', @all_jobs),
			outfile => $outfile
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
					die("Final STAR-FUSION accounting job: $run_id finished with errors.");
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
	't|tool=s'	=> \$tool_config,
	'd|data=s'	=> \$data_config,
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
	output_directory	=> $output_directory,
	hpc_driver		=> $hpc_driver,
	del_intermediates	=> $remove_junk,
	dry_run			=> $dry_run,
	no_wait			=> $no_wait
	);
