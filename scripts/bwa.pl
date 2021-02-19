#!/usr/bin/env perl
### bwa.pl #########################################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use YAML qw(LoadFile);
use IO::Handle;

my $cwd = dirname($0);
require "$cwd/utilities.pl";

####################################################################################################
# version       author		comment
# 1.0           sprokopec	run BWA-MEM on raw FASTQ files (paired end only!)
# 1.1		sprokopec	minor updates for compatibility with larger pipeline
# 1.2		sprokopec	added help msg and cleaned up code
# 1.3		sprokopec	minor updates for tool config

#### USAGE ##########################################################################################
# bwa.pl -t tool_config.yaml -d data_config.yaml -o /path/to/output/dir -b /path/to/output/yaml -c slurm --remove --dry_run
#
# where:
# 	-t (tool_config.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data_config.yaml) contains sample information (YAML file containing paths to FASTQ files)
#	-o (/path/to/output/dir) indicates tool-specific output directory
#	-b (/path/to/output/yaml) indicates output yaml for aligned BAMs
#	-c indicates hpc driver (ie, slurm)
#	--remove indicates whether or not to remove intermediate files
#	--dry_run indicates whether or not this is a dry run

### SUBROUTINES ####################################################################################
# format command to generate readgroup
sub set_readgroup {
	my %args = (
		subject		=> undef,
		sample		=> undef,
		lane		=> undef,
		library		=> undef,
		platform	=> undef,
		center		=> undef,
		@_
		);

	my $readgroup = '"@RG\tID:id\tSM:smp\tCN:center\tPL:platform\tLB:library\tPU:unit"';
	$readgroup =~ s/id/"$args{library}\_$args{lane}"/;
	$readgroup =~ s/smp/$args{sample}/;
	$readgroup =~ s/unit/$args{lane}/;
	$readgroup =~ s/library/$args{library}/;
	$readgroup =~ s/platform/$args{platform}/;
	$readgroup =~ s/center/$args{center}/;

	return($readgroup);
	}

# format command to run bwa then sort and index the output and remove tmp files
sub get_bwa_command {
	my %args = (
		r1		=> undef,
		r2		=> undef,
		stem		=> undef,
		readgroup	=> undef,
		reference	=> undef,
		@_
		);

	my $bwa_command = join(' ',
		'bwa mem -M -t4',  # 4 threads
		'-R', $args{readgroup},
		$args{reference},
		$args{r1}, $args{r2},
		'>', $args{stem} . '.sam'
		);

	return($bwa_command);
	}	

sub get_bwa_sort_command {
	my %args = (
		stem => undef,
		@_
		);

	my $sort_command = join(' ',
		'samtools sort -@4 -O bam',
		'-o', $args{stem} . '.bam',
		'-T', $args{stem}, $args{stem} . '.sam'
		);

	return($sort_command);
	}

sub get_index_bam_command {
	my %args = (
		output => undef,
		@_
		);

	my $index_command = join(' ',
		'samtools index',
		$args{output}
		);

	my $index_command2 = join(' ',
		'md5sum',
		$args{output},
		'>', $args{output} . '.md5'
		);
	
	my $command = join(";\n", $index_command, $index_command2);
	return($command);
	}

# format command to merge bams
sub get_merge_command {
	my %args = (
		input		=> undef,
		output		=> undef,
		tmp_dir		=> undef,
		java_mem	=> undef,
		@_
		);

	my $merge_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-jar $picard_dir/picard.jar MergeSamFiles',
		'INPUT=' . $args{input},
		'OUTPUT=' . $args{output},
		'ASSUME_SORTED=true USE_THREADING=true',
		'CREATE_INDEX=true CREATE_MD5_FILE=true',
		'TMP_DIR=' . $args{tmp_dir}
		);

	return($merge_command);
	}

sub get_merge_markdup_command {
	my %args = (
		input		=> undef,
		output		=> undef,
		tmp_dir		=> undef,
		java_mem	=> undef,
		@_
		);

	my $merge_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-jar $picard_dir/picard.jar MarkDuplicates',
		'INPUT=' . $args{input},
		'OUTPUT=' . $args{output},
		'METRICS_FILE=' . $args{output} . '.metrics',
		'TMP_DIR=' . $args{tmp_dir},
		'ASSUME_SORTED=true CREATE_INDEX=true CREATE_MD5_FILE=true',
		'MAX_RECORDS_IN_RAM=100000 VALIDATION_STRINGENCY=SILENT'
		);

	return($merge_command);
	}

### MAIN ##########################################################################################
sub main {
	my %args = (
		tool_config		=> undef,
		data_config		=> undef,
		output_directory	=> undef,
		output_config		=> undef,
		hpc_driver		=> undef,
		del_intermediates	=> undef,
		dry_run			=> undef,
		no_wait			=> undef,
		@_
		);

	my $tool_config = $args{tool_config};
	my $data_config = $args{data_config};

	### PREAMBLE ######################################################################################

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'bwa');

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_BWA_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_BWA_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running BWA-MEM alignment pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n    Reference: $tool_data->{reference}";
	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config\n";
	print $log "\n---";

	# set tools and versions
	my $bwa_version	= 'bwa/' . $tool_data->{bwa_version};
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};
	my $picard	= 'picard/' . $tool_data->{picard_version};

	# get user-specified tool parameters
	my $parameters = $tool_data->{bwa}->{parameters};

	### HANDLING FILES #################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	my ($run_script, $run_id, $run_id_extra, $raw_link);
	my @all_jobs;

	# initiate final output yaml file
	my $output_yaml = join('/', $output_directory, 'bwa_bam_config.yaml');
	if (defined($args{output_config})) {
		$output_yaml = $args{output_config};
		}
	open (my $yaml, '>', $output_yaml) or die "Cannot open '$output_yaml' !";
	print $yaml "---\n";

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";
		print $yaml "$patient:\n";

		# make a sample-specific directory
		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my @final_outputs;
		my $type = '';

		my (%tumours, %normals);

		foreach my $sample (sort keys %{$smp_data->{$patient}}) {

			print $log "  SAMPLE: $sample\n";

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			my $tmp_directory = join('/', $sample_directory, 'TEMP');
			unless(-e $tmp_directory) { make_path($tmp_directory); }
			my $cleanup_cmd = "rm -rf $tmp_directory";

			my @libraries = sort keys %{$smp_data->{$patient}->{$sample}->{libraries}};
			my (@lane_intermediates, @lane_holds, @smp_jobs);

			# determine sample type
			$type = $smp_data->{$patient}->{$sample}->{type};

			foreach my $lib ( @libraries ) {

				print $log "    LIBRARY: $lib\n";

				my @lanes = sort keys %{$smp_data->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}};

				foreach my $lane ( @lanes ) {

					print $log "      LANE: $lane\n";

					my $lane_directory = join('/', $sample_directory, $lane);
					unless(-e $lane_directory) { make_path($lane_directory); }

					my $raw_directory = join('/', $lane_directory, 'fastq_links');
					unless(-e $raw_directory) { make_path($raw_directory); }

					# collect input files
					my $r1 = $smp_data->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R1};
					my $r2 = $smp_data->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R2};

					print $log "	R1: $r1\n";
					print $log "	R2: $r2\n\n";

					my @tmp = split /\//, $r1;
					$raw_link = join('/', $raw_directory, $tmp[-1]);
					symlink($r1, $raw_link);
					@tmp = split /\//, $r2;
					$raw_link = join('/', $raw_directory, $tmp[-1]);
					symlink($r2, $raw_link);

					my $filestem = join('_', $sample, $lane);
					($run_id, $run_id_extra) = '';

					# run BWA-MEM on these fastqs
					my $readgroup = set_readgroup(
						subject		=> $patient,
						sample		=> $sample,
						lane		=> $lane,
						library 	=> $lib,
						platform 	=> $tool_data->{platform},
						center		=> $tool_data->{seq_center}
						);

					my $bwa = get_bwa_command(
						r1		=> $r1,
						r2		=> $r2,
						stem		=> $filestem,
						readgroup	=> $readgroup,
						reference	=> $tool_data->{bwa}->{reference}
						);

					$bwa = 'cd ' . $lane_directory . ";\n" . $bwa ;

					# check if this should be run
					if ( 
						('Y' eq missing_file(join('/', $lane_directory, $filestem . '.sam'))) &&
						# if the index is missing, then this has been run before and the 
						# resulting sam deleted, (so don't do it again)
						('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam.bai')))
						) {
						# record command (in log directory) and then run job
						print $log "Submitting job for bwa...\n";
						$run_script = write_script(
							log_dir => $log_directory,
							name	=> 'run_bwa_mem_' . $filestem,
							cmd	=> $bwa,
							modules => [$bwa_version, $samtools],
							max_time	=> $parameters->{bwa}->{time}->{$type},
							mem		=> $parameters->{bwa}->{mem}->{$type},
							cpus_per_task	=> 4,
							hpc_driver	=> $args{hpc_driver}
							);

						$run_id = submit_job(
							jobname		=> 'run_bwa_mem_' . $filestem,
							shell_command	=> $run_script,
							hpc_driver	=> $args{hpc_driver},
							dry_run		=> $args{dry_run},
							log_file	=> $log
							);

						push @smp_jobs, $run_id;
						push @all_jobs, $run_id;
						}
					else {
						print $log "Skipping alignment because output already exists...\n";
						}

					# sort the resulting SAM and output BAM
					my $sort_cmd = "cd $lane_directory;\n";
					$sort_cmd .= get_bwa_sort_command(
						stem => $filestem
						);

					# check if this should be run
					if (
						('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam'))) &&
						# if the index is missing, then this has been run before and
						# the resulting sam deleted (so don't do it again)
						('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam.bai')))
						) {
						# record command (in log directory) and then run job
						print $log "Submitting job to sort bam...\n";
						$run_script = write_script(
							log_dir	=> $log_directory,
							name	=> 'run_sort_sam_' . $filestem,
							cmd	=> $sort_cmd,
							modules	=> [$samtools],
							dependencies	=> $run_id,
							max_time	=> $parameters->{sort}->{time}->{$type},
							mem		=> $parameters->{sort}->{mem}->{$type},
							hpc_driver	=> $args{hpc_driver}
							);

						$run_id = submit_job(
							jobname		=> 'run_sort_sam_' . $filestem,
							shell_command	=> $run_script,
							hpc_driver	=> $args{hpc_driver},
							dry_run		=> $args{dry_run},
							log_file	=> $log
							);

						push @smp_jobs, $run_id;
						push @all_jobs, $run_id;
						}
					else {
						print $log "Skipping sort step because this was performed previously...\n";
						}

					# index the resulting BAM and remove intermediate SAM
					my $index_lane_cmd = "cd $lane_directory;\n";
					$index_lane_cmd .= get_index_bam_command(
						output => $filestem . '.bam'
						);

					$cleanup_cmd .= "\nrm $lane_directory/$filestem.sam";

					# check if this should be run
					if (
						('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam.bai'))) ||
						('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam.md5')))
						) {
						# record command (in log directory) and then run job
						print $log "Submitting job to index bam...\n";
						$run_script = write_script(
							log_dir	=> $log_directory,
							name	=> 'run_bam_index_' . $filestem,
							cmd	=> $index_lane_cmd,
							modules	=> [$samtools],
							dependencies	=> $run_id,
							max_time	=> $parameters->{index}->{time}->{$type},
							mem		=> $parameters->{index}->{mem}->{$type},
							hpc_driver	=> $args{hpc_driver}
							);

						$run_id = submit_job(
							jobname		=> 'run_bam_index_' . $filestem,
							shell_command	=> $run_script,
							hpc_driver	=> $args{hpc_driver},
							dry_run		=> $args{dry_run},
							log_file	=> $log
							);

						push @smp_jobs, $run_id;
						push @all_jobs, $run_id;
						}
					else {
						print $log "Skipping index step because this was performed previously...\n";
						}

					push @lane_intermediates, $lane_directory . '/' . $filestem . '.bam';
					push @lane_holds, $run_id;
					}
				}

			# if there are multiple libraries for a single sample, merge them
			my ($smp_output, $jobname, $merge_cmd);
			my $input_files = join(' INPUT=', @lane_intermediates);

			if ('N' eq $parameters->{merge}->{mark_dup}) {

				$smp_output = $patient_directory . '/' . $sample . '_bwamem_aligned_sorted_merged.bam';
				$jobname = 'run_merge_lanes_' . $sample;

				$merge_cmd = get_merge_command(
					input		=> $input_files,
					output		=> $smp_output,
					tmp_dir		=> $tmp_directory,
					java_mem	=> $parameters->{merge}->{java_mem}->{$type}
					);
				}

			elsif ('Y' eq $parameters->{merge}->{mark_dup}) {

				$smp_output = join('/', $patient_directory, $sample . '_bwamem_aligned_sorted_merged_markdup.bam');
				$jobname = 'run_merge_markdup_' . $sample;

				$merge_cmd = get_merge_markdup_command(
					input		=> $input_files,
					output		=> $smp_output,
					tmp_dir		=> $tmp_directory,
					java_mem	=> $parameters->{merge}->{java_mem}->{$type}
					);
				}

			$cleanup_cmd .= "\nrm " . join(";\nrm ", @lane_intermediates) . "\n";

			# add output file to list
			if ('normal' eq $type) { $normals{$sample} = $smp_output; }
			if ('tumour' eq $type) { $tumours{$sample} = $smp_output; }

			# check if this should be run
			if ('Y' eq missing_file("$smp_output.md5")) {

				# record command (in log directory) and then run job
				print $log "Submitting job to merge lanes and mark duplicates...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> $jobname,
					cmd	=> $merge_cmd,
					modules	=> [$picard],
					dependencies	=> join(':', @lane_holds),
					mem 		=> $parameters->{merge}->{mem}->{$type},
					max_time 	=> $parameters->{merge}->{time}->{$type},
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id = submit_job(
					jobname		=> $jobname,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @smp_jobs, $run_id;
				push @all_jobs, $run_id;
				}
			else {
				print $log "Skipping mark duplicate step because this was performed previously...\n";	
				}

			push @final_outputs, $smp_output;

			# clean up/remove intermediate files (once per sample)
			if ($args{del_intermediates}) {

				if (scalar(@smp_jobs) == 0) {
					`rm -rf $tmp_directory`;
					} else {

					print $log "Submitting job to clean up temporary/intermediate files...\n";

					# make sure final output exists before removing intermediate files!
					$cleanup_cmd = join("\n",
						"if [ -s $smp_output ]; then",
						$cleanup_cmd,
						"else",
						"  echo " . '"FINAL OUTPUT: ' . $smp_output . ' is missing; not removing intermediates"',
						"fi"
						);

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_cleanup_' . $sample,
						cmd	=> $cleanup_cmd,
						dependencies	=> join(':', @smp_jobs),
						mem		=> '256M',
						hpc_driver	=> $args{hpc_driver},
						kill_on_error	=> 0
						);

					$run_id_extra = submit_job(
						jobname		=> 'run_cleanup_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);
					}
				}
			}

		# and finally, add the final files to the output yaml
		my $key;
		if (scalar(keys %tumours) > 0) {
			print $yaml "    tumour:\n";
			foreach $key (keys %tumours) { print $yaml "        $key: $tumours{$key}\n"; }
			}
		if (scalar(keys %normals) > 0) {
			print $yaml "    normal:\n";
			foreach $key (keys %normals) { print $yaml "        $key: $normals{$key}\n"; }
			}

		print $log "FINAL OUTPUT:\n" . join("\n  ", @final_outputs) . "\n";
		print $log "---\n";
		}

	close $yaml;

	# if this is not a dry run OR there are jobs to assess (run or resumed with jobs submitted) then
	# collect job metrics (exit status, mem, run time)
	unless ( ($args{dry_run}) || (scalar(@all_jobs) == 0) ) {

		# collect job stats
		my $collect_metrics = collect_job_stats(
			job_ids	=> join(',', @all_jobs),
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
					die("Final BWA accounting job: $run_id finished with errors.");
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
my ($data_config, $tool_config, $output_directory, $output_config) = undef;
my $hpc_driver = 'slurm';
my ($remove_junk, $dry_run, $no_wait);
my $help;

# read in command line arguments
GetOptions(
	'h|help'	=> \$help,
	'd|data=s'	=> \$data_config,
	't|tool=s'	=> \$tool_config,
	'o|out_dir=s'	=> \$output_directory,
	'b|out_yaml=s'	=> \$output_config,
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
		"\t--out_yaml|-b\t<string> path to output yaml (listing BWA-aligned BAMs)",
		"\t--cluster|-c\t<string> cluster scheduler (default: slurm)",
		"\t--remove\t<boolean> should intermediates be removed? (default: false)",
		"\t--dry-run\t<boolean> should jobs be submitted? (default: false)",
		"\t--no-wait\t<boolean> should we exit after job submission (true) or wait until all jobs have completed (false)? (default: false)"
		);

	print "$help_msg\n";
	exit;
	}

if (!defined($tool_config)) { die("No tool config file defined; please provide -t | --tool (ie, tool_config.yaml)"); }
if (!defined($data_config)) { die("No data config file defined; please provide -d | --data (ie, sample_config.yaml)"); }
if (!defined($output_directory)) { die("No output directory defined; please provide -o | --out_dir"); }

main(
	tool_config		=> $tool_config,
	data_config		=> $data_config,
	output_directory	=> $output_directory,
	output_config		=> $output_config,
	hpc_driver		=> $hpc_driver,
	del_intermediates	=> $remove_junk,
	dry_run			=> $dry_run,
	no_wait			=> $no_wait
	);
