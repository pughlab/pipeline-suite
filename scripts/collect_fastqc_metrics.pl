#!/usr/bin/env perl
### collect_fastqc_metrics.pl ######################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use Getopt::Std;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use File::Path qw(make_path);
use List::Util 'any';
use YAML qw(LoadFile);

my $cwd = dirname(__FILE__);
require "$cwd/utilities.pl";

####################################################################################################
# version       author		comment
# 1.0		sprokopec       script to collect quality metrics from fastq or BAM files

### USAGE ##########################################################################################
# collect_fastqc_metrics.pl -t tool.yaml -d data.yaml { --rna } -c slurm --dry_run
#
# where:
#	-t (tool.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data.yaml) contains sample information (YAML file containing paths to fastq files)
#	--rna binary indicating if input is RNA (because the data.yaml will be slightly different)
#	-c indicates hpc driver (ie, slurm)
#	--dry_run indicates that this is a dry run

### DEFINE SUBROUTINES #############################################################################
# format command to run FASTQC
sub get_fastqc_cmd {
	my %args = (
		input		=> undef,
		output_dir	=> undef,
		@_
		);

	my $fastqc_command = join(' ',
		'fastqc -o', $args{output_dir},
		$args{input}
		);

	return($fastqc_command);
	}

# format command to extract metrics
sub get_fastqc_metrics {
	my %args = (
		output_dir	=> undef,
		@_
		);

	my $extract_cmd = join("\n",
		"cd $args{output_dir}",
		'for i in *.zip',
		'  do',
		'    TMP=${i//.zip/};',
		'    OUT=${i//.zip/_metrics.txt};',
		'    unzip $i $TMP/fastqc_data.txt;',
		'    cp $TMP/fastqc_data.txt $OUT;',
		'    rm -rf $TMP/;',
		'  done'
		);

	return($extract_cmd);
	}

### MAIN ###########################################################################################
sub main{
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
		print "Initiating FASTQC pipeline...\n";
		}

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'fastqc');

	# organize output directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	# start logging
	my $log_file = join('/', $log_directory, 'run_FASTQC_pipeline.log');

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
		if (0 != $touch_exit_status) { Carp::croak("Cannot touch file $outfile"); 

		$log_file = join('/', $log_directory, 'run_FASTQC_pipeline_' . $run_count . '.log');}
		}

	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running FASTQC pipeline to collect quality metrics from fastq files...\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---\n";

	# set tools and versions
	my $fastqc	= 'fastqc';
	my $r_version   = 'R/' . $tool_data->{r_version};

	# get user-specified tool parameters
	my $parameters = $tool_data->{fastqc}->{parameters};

	# get optional HPC group
	my $hpc_group = defined($tool_data->{hpc_group}) ? "-A $tool_data->{hpc_group}" : undef;

	### RUN ###########################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	unless($args{dry_run}) {
		print "Processing fastq files for " . scalar(keys %{$smp_data}) . " patients.\n";
		}

	my ($run_script, $run_id);
	my @all_jobs;

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient";

		foreach my $sample (sort keys %{$smp_data->{$patient}}) {

			my @fastqs;

			# is this DNA or RNA?
#			unless ( 'rna' eq $tool_data->{seq_type} ) {

				my @libraries = keys %{$smp_data->{$patient}->{$sample}->{libraries}};

				foreach my $lib ( @libraries ) {

					my @lanes = keys %{$smp_data->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}};
					foreach my $lane ( @lanes ) {

						my $r1 = $smp_data->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R1};
						my $r2 = $smp_data->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R2};

						push @fastqs, $r1;
						push @fastqs, $r2;
						}
					}

#				} else {

#				my @lanes = keys %{$smp_data->{$patient}->{$sample}->{runlanes}};

#				foreach my $lane ( @lanes ) {

#					my $r1 = $smp_data->{$patient}->{$sample}->{runlanes}->{$lane}->{R1};
#					my $r2 = $smp_data->{$patient}->{$sample}->{runlanes}->{$lane}->{R2};

#					push @fastqs, $r1;
#					push @fastqs, $r2;
#					}
#				}

			my $n_fastqs = scalar(@fastqs);
			print $log "  SAMPLE: $sample has $n_fastqs fastqs to QC\n";

			# Run FASTQC
			my $out_file = join('/', $output_directory, $sample . '_fastqc_status.COMPLETE');
			my $fastqc_cmd = get_fastqc_cmd( 
				input		=> join(' ', @fastqs),
				output_dir	=> $output_directory
				);

			$fastqc_cmd .= "\n\n" . "echo 'fastqc finished successfully' > $out_file";

			# check if this should be run
			if ('Y' eq missing_file($out_file)) {

				# record command (in log directory) and then run job
				print $log "  >> Submitting job for FASTQC...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_fastqc_' . $sample,
					cmd	=> $fastqc_cmd,
					modules	=> [$fastqc],
					max_time	=> $parameters->{fastqc}->{time},
					mem		=> $parameters->{fastqc}->{mem},
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'run_fastqc_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @all_jobs, $run_id;
				} else {
				print $log "  >> Skipping FASTQC because output already exists...\n";
				}

			# get MD5SUMs (if necessary)
			$out_file = join('/', $output_directory, $sample . '_md5sums.COMPLETE');
			my $md5_cmd = "cd $output_directory\n";
			foreach my $file ( @fastqs ) {
				my $stem = basename($file);
				$md5_cmd .= "md5sum $file > $stem.md5\n";
				}

			$md5_cmd .= "\n\n" . "echo 'md5sum collection finished successfully' > $out_file";

			# check if this should be run
			if ('Y' eq missing_file($out_file)) {

				# record command (in log directory) and then run job
				print $log "  >> Submitting job for MD5SUM...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_md5sums_' . $sample,
					cmd	=> $md5_cmd,
					max_time	=> '24:00:00',
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'run_md5sums_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @all_jobs, $run_id;
				} else {
				print $log "  >> Skipping MD5SUM collection because output already exists...\n";
				}
			}
		}

	# format command to extract info from output
	my $extract_cmd = get_fastqc_metrics(output_dir => $output_directory);

	# format command to collate results
	my $combine_cmd = "Rscript $cwd/combine_key_metrics.R";

	my $collate_cmd = $extract_cmd . "\n" . $combine_cmd;
	$collate_cmd .= "\n\n" . join("\n",
		"if [ -s fastqc_key_metrics.tsv ]; then",
		"  rm *metrics.txt",
		"  rm *md5",
		"else",
		"  echo Error in creating final output: fastqc_key_metrics.tsv",
		"fi"
		);

	print $log "Submitting job to extract and collate results...\n";

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'run_collate_results',
		cmd	=> $collate_cmd,
		modules	=> [$r_version],
		dependencies	=> join(',', @all_jobs),
		max_time	=> '01:00:00',
		mem		=> '1G',
		hpc_driver	=> $args{hpc_driver}
		);

	$run_id = submit_job(
		jobname		=> 'run_collate_results',
		shell_command	=> $run_script,
		hpc_driver	=> $args{hpc_driver},
		dry_run		=> $args{dry_run},
		log_file	=> $log
		);

	push @all_jobs, $run_id;

	print $log "\nFINAL OUTPUT: $output_directory/fastqc_key_metrics.tsv\n";
	print $log "---\n";

	# collect job metrics if not dry_run
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
			kill_on_error	=> 0,
			hpc_driver	=> $args{hpc_driver},
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
my ($help, $no_wait, $dry_run);

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	't|tool=s'	=> \$tool_config,
	'd|data=s'	=> \$data_config,
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

main(
	tool_config		=> $tool_config,
	data_config		=> $data_config,
	output_directory	=> $output_directory,
	hpc_driver		=> $hpc_driver,
	dry_run			=> $dry_run,
	no_wait			=> $no_wait
	);
