#!/usr/bin/env perl
### bwa.pl ###########################################################################
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

my $cwd = dirname($0);
require "$cwd/utilities.pl";

### SUBROUTINES ####################################################################################
# format command to generate readgroup
sub set_readgroup {
	my %args = (
		subject		=> undef,
		sample		=> undef,
		lane		=> undef,
		library		=> undef,
		platform	=> undef,
		@_
		);

	my $readgroup = '"@RG\tID:patient\tSM:smp\tPL:platform\tPU:unit\tLB:library"';
	$readgroup =~ s/patient/$args{subject}/;
	$readgroup =~ s/smp/$args{sample}/;
	$readgroup =~ s/unit/$args{lane}/;
	$readgroup =~ s/library/$args{library}/;
	$readgroup =~ s/platform/$args{platform}/;

	return($readgroup);
	}

# format command to run bwa then sort and index the output and remove tmp files
sub get_bwa_command {
	my %args = (
		reference	=> undef,
		r1		=> undef,
		r2		=> undef,
		stem		=> undef,
		readgroup	=> undef,
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

sub main {
	my %args = (
		tool_config => undef,
		data_config => undef,
		@_
		);

	my $tool_config = $args{tool_config};
	my $data_config = $args{data_config};

	my $date = strftime "%F", localtime;

	### PREAMBLE ######################################################################################

	if (!defined($tool_config)) { die("No tool config file defined; please provide -t | --tool (ie, tool_config.yaml)"); }
	if (!defined($data_config)) { die("No data config file defined; please provide -c | --config (ie, sample_config.yaml)"); }

	# load tool config
	my $tool_data = LoadFile($tool_config);

	# start logging
	print "---\n";
	print "Running BWA-MEM alignment pipeline for NGS data.\n";
	print "\n  Tool config used: $tool_config";
	print "\n    Reference: $tool_data->{reference}";
	print "\n    Output directory: $tool_data->{output_dir}";
	print "\n  Sample config used: $data_config\n";
	print "\n---";

	# set tools and versions
	my $bwa_version = $tool_data->{tool} . '/' . $tool_data->{tool_version};
	my $samtools = 'samtools/' . $tool_data->{samtools_version};
	my $picard = 'picard/' . $tool_data->{picard_version};

	if ('bwamem' ne $tool_data->{aligner}) { die("This pipeline is currently only compatible with BWA-MEM!"); }
	if (('Y' ne $tool_data->{mark_dup}) & ('N' ne $tool_data->{mark_dup})) { die("Option mark_dup must be either Y or N !"); }
	if (('Y' ne $tool_data->{del_intermediate}) & ('N' ne $tool_data->{del_intermediate})) { die("Option del_intermediate must be either Y or N !"); }

	if ((!defined($tool_data->{dry_run})) || ('Y' ne $tool_data->{dry_run})) {
		$tool_data->{dry_run} = 'N';
		}

	### CREATE DIRECTORY STRUCTURE ####################################################################
	my $output_directory;
	my $resume = 'N';

	# check for RESUME directory
	if (defined ($tool_data->{resume_dir})) {
		$resume = 'Y';
		$tool_data->{resume_dir} =~ s/\/$//;
		$output_directory = $tool_data->{resume_dir};
		}

	# otherwise, make a directory for this batch
	else {
		# make sure output directory exists
		unless(-e $tool_data->{output_dir}) { make_path($tool_data->{output_dir}); }

		$tool_data->{output_dir} =~ s/\/$//;
		$output_directory = join('/', $tool_data->{output_dir},  join('_', $date, $tool_data->{tool}, $tool_data->{tool_version}));
		unless(-e $output_directory) { make_path($output_directory); }
		}

	# and directory for log files
	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	### HANDLING FILES #################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	my ($run_script, $run_id, $link);

	my (@all_jobs);

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print "\nInitiating process for PATIENT: $patient\n";

		# make a sample-specific directory
		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my @final_outputs;
		my $type;

		foreach my $sample (sort keys %{$smp_data->{$patient}}) {

			print "  SAMPLE: $sample\n";

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			my $tmp_directory = join('/', $sample_directory, 'TEMP');
			unless(-e $tmp_directory) { make_path($tmp_directory); }

			my @libraries = keys %{$smp_data->{$patient}->{$sample}->{libraries}};
			my (@lane_intermediates, @lane_holds);

			# determine sample type
			$type = $smp_data->{$patient}->{$sample}->{type};

			foreach my $lib ( @libraries ) {

				print "    LIBRARY: $lib\n";

				my @lanes = keys %{$smp_data->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}};

				foreach my $lane ( @lanes ) {

					print "      LANE: $lane\n";

					my $lane_directory = join('/', $sample_directory, $lane);
					unless(-e $lane_directory) { make_path($lane_directory); }

					my $raw_directory = join('/', $lane_directory, 'fastq_links');
					unless(-e $raw_directory) { make_path($raw_directory); }

					# collect input files
					my $r1 = $smp_data->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R1};
					my $r2 = $smp_data->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R2};

					print "	R1: $r1\n";
					print "	R2: $r2\n";

					my @tmp = split /\//, $r1;
					$link = join('/', $raw_directory, $tmp[-1]);
					symlink($r1, $link);
					$link =~ s/R1/R2/;
					symlink($r2, $link);

					my $filestem = join('_', $sample, $lane);
					$run_id = '';

					# run BWA-MEM on these fastqs
					my $readgroup = set_readgroup(
						subject		=> $patient,
						sample		=> $sample,
						lane		=> $lane,
						library 	=> $lib,
						platform 	=> $tool_data->{platform}
						);

					my $bwa = get_bwa_command(
						reference	=> $tool_data->{reference},
						r1		=> $r1,
						r2		=> $r2,
						stem		=> $filestem,
						readgroup	=> $readgroup
						);

					$bwa = 'cd ' . $lane_directory . ";\n" . $bwa ;

					# check if this should be run
					if ( ('N' eq $resume) ||
						( 
							('Y' eq missing_file(join('/', $lane_directory, $filestem . '.sam'))) &&
							# if the index is missing, then this has been run before and the resulting sam deleted
							# (so don't do it again)
							('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam.bai')))
							)
						) {
						# record command (in log directory) and then run job
						print "Submitting job for bwa...\n";
						$run_script = write_script(
							log_dir => $log_directory,
							name	=> 'run_bwa_mem_' . $filestem,
							cmd	=> $bwa,
							modules => [$bwa_version, $samtools]
							);

						$run_id = submit_job(
							jobname		=> 'run_bwa_mem_' . $filestem,
							shell_command	=> $run_script,
							dependencies	=> $run_id,
							max_time	=> $tool_data->{parameters}->{bwa}->{time}->{$type},
							mem		=> $tool_data->{parameters}->{bwa}->{mem}->{$type},
							cpus_per_task	=> 4,
							dry_run		=> $tool_data->{dry_run}
							);

						push @all_jobs, $run_id;
						}
					else {
						print "Skipping alignment because this step was performed previously...\n";
						}

					# sort the resulting SAM and output BAM
					my $sort_cmd = get_bwa_sort_command(
						stem => $filestem
						);
					$sort_cmd = 'cd ' . $lane_directory . ";\n" . $sort_cmd;

					# check if this should be run
					if ( ('N' eq $resume) ||
						(
							('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam'))) &&
							# if the index is missing, then this has been run before and the resulting sam deleted (so don't do it again)
							('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam.bai')))
							)
						) {
						# record command (in log directory) and then run job
						print "Submitting job to sort bam...\n";
						$run_script = write_script(
							log_dir	=> $log_directory,
							name	=> 'run_sort_sam_' . $filestem,
							cmd	=> $sort_cmd,
							modules	=> [$samtools]
							);

						$run_id = submit_job(
							jobname		=> 'run_sort_sam_' . $filestem,
							shell_command	=> $run_script,
							dependencies	=> $run_id,
							max_time	=> $tool_data->{parameters}->{sort}->{time}->{$type},
							mem		=> $tool_data->{parameters}->{sort}->{mem}->{$type},
							cpus_per_task 	=> 1,
							dry_run		=> $tool_data->{dry_run}
							);

						push @all_jobs, $run_id;
						}
					else {
						print "Skipping sort step because this was performed previously...\n";
						}

					# index the resulting BAM and remove intermediate SAM
					my $index_lane_cmd = get_index_bam_command(
						output => $filestem . '.bam'
						);
					$index_lane_cmd = 'cd ' . $lane_directory . ";\n" . $index_lane_cmd;

					if ('Y' eq $tool_data->{del_intermediate}) {
						$index_lane_cmd = $index_lane_cmd . ";\nrm $filestem.sam";
						}

					# check if this should be run
					if ( ('N' eq $resume) ||
						(
							('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam.bai'))) ||
							('Y' eq missing_file(join('/', $lane_directory, $filestem . '.bam.md5')))
							)
						) {
						# record command (in log directory) and then run job
						print "Submitting job to index bam...\n";
						$run_script = write_script(
							log_dir	=> $log_directory,
							name	=> 'run_bam_index_' . $filestem,
							cmd	=> $index_lane_cmd,
							modules	=> [$samtools]
							);

						$run_id = submit_job(
							jobname		=> 'run_bam_index_' . $filestem,
							shell_command	=> $run_script,
							dependencies	=> $run_id,
							max_time	=> $tool_data->{parameters}->{index}->{time}->{$type},
							mem		=> $tool_data->{parameters}->{index}->{mem}->{$type},
							cpus_per_task	=> 1,
							dry_run		=> $tool_data->{dry_run}
							);

						push @all_jobs, $run_id;
						}
					else {
						print "Skipping index step because this was performed previously...\n";
						}

					push @lane_intermediates, $lane_directory . '/' . $filestem . '.bam';
					push @lane_holds, $run_id;
					}
				}

			# if there are multiple libraries for a single sample, merge them
			my ($smp_output, $jobname, $merge_cmd);
			my $input_files = join(' INPUT=', @lane_intermediates);

			if ('N' eq $tool_data->{mark_dup}) {

				$smp_output = $patient_directory . '/' . $sample . '_bwamem_aligned_sorted_merged.bam';
				$jobname = 'run_merge_lanes_' . $sample;

				$merge_cmd = get_merge_command(
					input		=> $input_files,
					output		=> $smp_output,
					tmp_dir		=> $tmp_directory,
					java_mem	=> $tool_data->{parameters}->{merge}->{java_mem}->{$type}
					);
				}

			elsif ('Y' eq $tool_data->{mark_dup}) {

				$smp_output = $patient_directory . '/' . $sample . '_bwamem_aligned_sorted_merged_markdup.bam';
				$jobname = 'run_merge_markdup_' . $sample;

				$merge_cmd = get_merge_markdup_command(
					input		=> $input_files,
					output		=> $smp_output,
					tmp_dir		=> $tmp_directory,
					java_mem	=> $tool_data->{parameters}->{merge}->{java_mem}->{$type}
					);
				}

			if ('Y' eq $tool_data->{del_intermediate}) {
				$merge_cmd .= ";\nrm " . join(";\nrm ", @lane_intermediates) . "\n";
				}

			# check if this should be run
			if ( ('N' eq $resume) || ('Y' eq missing_file($smp_output)) ) {
				# record command (in log directory) and then run job
				print "Submitting job to merge lanes and mark dupilcates...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> $jobname,
					cmd	=> $merge_cmd,
					modules	=> [$picard]
					);

				$run_id = submit_job(
					jobname		=> $jobname,
					shell_command	=> $run_script,
					dependencies	=> join(',', @lane_holds),
					mem 		=> $tool_data->{parameters}->{merge}->{mem}->{$type},
					max_time 	=> $tool_data->{parameters}->{merge}->{time}->{$type},
					cpus_per_task	=> 1,
					dry_run		=> $tool_data->{dry_run}
					);

				push @all_jobs, $run_id;
				}
			else {
				print "Skipping mark dupicate step because this was performed previously...\n";	
				}

			push @final_outputs, $smp_output;
			}

		print "\nFINAL OUTPUT:\n" . join("\n  ", @final_outputs) . "\n";
		print "---\n";
		}

	# collect job stats
	my $collect_metrics = collect_job_stats(
		job_ids	=> join(',', @all_jobs),
		outfile	=> $log_directory . '/slurm_job_metrics.out'
		);

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'output_job_metrics',
		cmd	=> $collect_metrics
		);

	$run_id = submit_job(
		jobname		=> 'output_job_metrics',
		shell_command	=> $run_script,
		dependencies	=> join(',', @all_jobs),
		max_time	=> '0:10:00',
		mem		=> '1G',
		cpus_per_task	=> 1,
		dry_run		=> $tool_data->{dry_run}
		);

	# final job to output a BAM config for downstream stuff
	my $output_yaml_cmd = join(' ',
		'perl ./create_bam_yaml.pl',
		'-d', $output_directory,
		'-o', $output_directory . '/bam_config.yaml'
		);

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'output_final_yaml',
		cmd	=> $output_yaml_cmd,
		modules	=> ['perl']
		);

	$run_id = submit_job(
		jobname		=> 'output_final_yaml',
		shell_command	=> $run_script,
		dependencies	=> join(',', @all_jobs),
		max_time	=> '0:10:00',
		mem		=> '1G',
		cpus_per_tasks	=> 1,
		dry_run		=> $tool_data->{dry_run}
		);

	# finish up
	print "\nProgramming terminated successfully.\n\n";

	}

### GETOPTS PLUS ERROR CHECKING AND DEFAULT VALUES #################################################
# declare variables
my $tool_config;
my $data_config;

# read in command line arguments
GetOptions(
	't|tool=s'	=> \$tool_config,
	'c|config=s'    => \$data_config
	);

main(tool_config => $tool_config, data_config => $data_config);
