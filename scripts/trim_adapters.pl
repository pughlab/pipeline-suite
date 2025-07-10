#!/usr/bin/env perl
### trim_adapters.pl ###############################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use version;
use Carp;
use POSIX qw(strftime);
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use IO::Handle;
use YAML qw(LoadFile DumpFile);
use Data::Dumper;

my $cwd = dirname(__FILE__);
require "$cwd/utilities.pl";

# define some global variables
our $tool_version = undef;

####################################################################################################
# version       author	  	comment
# 1.0		sprokopec       script to run trim_galore or fastp on EM-Seq fastq files

### USAGE ##########################################################################################
# trim_adapters.pl -t tool_config.yaml -d data_config.yaml -o /path/to/output/dir -c slurm --remove --dry_run
#
# where:
#	-t (tool_config.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data_config.yaml) contains sample information (YAML file containing paths to FASTQ files)
#	-o (/path/to/output/dir) indicates tool-specific output directory
#	-b (/path/to/output/yaml) indicates YAML containing processed BAMs
#	-c indicates hpc driver (ie, slurm)
#	--remove remove intermediate files?
#	--dry_run is this a dry run?

### SUBROUTINES ####################################################################################
# format command to run trim_galore
sub get_trim_galore_command {
	my %args = (
		r1	=> undef,
		r2	=> undef,
		out_dir	=> undef,
		@_
		);

	my $trim_command = join(' ',
		'trim_galore',
		'--output_dir', $args{out_dir},
		'--2colour 20',
		'--paired', $args{r1}, $args{r2}
		);

	return($trim_command);
	}

# format command to run fastp
sub get_fastp_command {
	my %args = (
		r1	=> undef,
		r2	=> undef,
		stem	=> undef,
		@_
		);

	my $trim_command = join(' ',
		'fastp',
		'--in1', $args{r1},
		'--in2', $args{r2},
		'--out1', $args{stem} . '_R1.fastq.gz',
		'--out2', $args{stem} . '_R2.fastq.gz',
		'-l 2 -Q -p --trim_poly_g',
		'--html', $args{stem} . '.html',
		'--json', $args{stem} . '.json'
		);

	return($trim_command);
	}

# format command to run custom UMI trimming (ConcensusCruncher)
sub get_extractUMI_command {
	my %args = (
		r1	=> undef,
		r2	=> undef,
		stem	=> undef,
		@_
		);

	my $trim_command = join(' ',
		'python3 /cluster/projects/pughlab/src/ConsensusCruncher/scripts/extract_barcodes.py',
		'--blist /cluster/projects/pughlab/src/UMI-tools/idt_dual_index_barcodes.txt',
		'--read1', $args{r1},
		'--read2', $args{r2},
		'--outfile', $args{stem}
		);

	return($trim_command);
	}

### MAIN ###########################################################################################
sub main {
	my %args = (
		tool_config		=> undef,
		data_config		=> undef,
		output_directory	=> undef,
		output_config		=> undef,
		hpc_driver		=> undef,
		dry_run			=> undef,
		no_wait			=> undef,
		@_
		);

	my $tool_config = $args{tool_config};
	my $data_config = $args{data_config};

	### PREAMBLE ######################################################################################
	unless($args{dry_run}) {
		print "Initiating AdapterTrim pipeline...\n";
		}

	# load tool config
	my $tool_data = LoadFile($tool_config);

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_TRIM_ADAPTERS_pipeline.log');

	# create a file to hold job metrics
	my (@files, $outfile, $touch_exit_status);
	my $run_count = 0;
	unless ($args{dry_run}) {
		# initiate a file to hold job metrics (ensures that an existing file isn't overwritten by concurrent jobs)
		opendir(LOGFILES, $log_directory) or die "Cannot open $log_directory";
		@files = grep { /slurm_job_metrics/ } readdir(LOGFILES);
		$run_count = scalar(@files) + 1;
		closedir(LOGFILES);

		$outfile = $log_directory . '/slurm_job_metrics_' . $run_count . '.out';
		$touch_exit_status = system("touch $outfile");
		if (0 != $touch_exit_status) { Carp::croak("Cannot touch file $outfile"); }

		$log_file = join('/', $log_directory, 'run_TRIM_ADAPTERS_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running AdapterTrim pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n  Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";

	# set tools and versions
	my $trim_tool = '';
	if ('caphla' eq $tool_data->{seq_type}) {
		print $log "\n\n  Trimming using tool: ConsensusCruncher";
		$trim_tool = 'python3/3.7.2';
		} elsif (defined($tool_data->{baits_bed})) {
		print $log "\n\n  Trimming using tool: trim_galore";
		$trim_tool = 'trim_galore/' . $tool_data->{trim_galore_version};
		} else {
		print $log "\n\n  Trimming using tool: fastq";
		$trim_tool = 'fastp/' . $tool_data->{fastp_version};
		}

	print $log "\n---";

	# get user-specified tool parameters
	my $parameters = $tool_data->{trim_adapters}->{parameters};

	# get optional HPC group
	my $hpc_group = defined($tool_data->{hpc_group}) ? "-A $tool_data->{hpc_group}" : undef;

	### MAIN ###########################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	unless($args{dry_run}) {
		print "Processing " . scalar(keys %{$smp_data}) . " patients.\n";
		}

	my ($run_script, $run_id);
	my @all_jobs;

	# initiate final output yaml file
	my $output_yaml = join('/', $output_directory, 'emseq_fastq_trimmed_config.yaml');
	if (defined($args{output_config})) {
		$output_yaml = $args{output_config};
		}

	# use input config as template
	my $smp_data_out = $smp_data;

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";

		foreach my $sample (sort keys %{$smp_data->{$patient}}) {

			print $log "  SAMPLE: $sample\n";

			# determine sample type
			my $type = $smp_data->{$patient}->{$sample}->{type};

			# get list of libraries
			my @libraries = keys %{$smp_data->{$patient}->{$sample}->{libraries}};

			foreach my $lib ( @libraries ) {
				print $log "\n    LIBRARY: $lib\n";
				my $lib_data = $smp_data->{$patient}->{$sample}->{libraries}->{$lib};

				# get list of lanes
				my @lanes = keys %{$lib_data->{runlanes}};

				foreach my $lane ( sort(@lanes) ) {
					print $log "      LANE: $lane\n";

					my $r1 = $lib_data->{runlanes}->{$lane}->{fastq}->{R1};
					my $r2 = $lib_data->{runlanes}->{$lane}->{fastq}->{R2};

					print $log "      R1: $r1\n";
					print $log "      R2: $r2\n\n";

					my $stem = join('_', $lib, $lane);
					my $output_stem = join('/', $output_directory, $stem);
					my ($trim_command, $new_r1, $new_r2);

					if ('caphla' eq $tool_data->{seq_type}) {

						$trim_command = get_extractUMI_command(
							r1	=> $r1,
							r2	=> $r2,
							stem	=> $output_stem
							);

						$new_r1 = $output_stem . '_barcode_R1.fastq';
						$new_r2 = $output_stem . '_barcode_R2.fastq';

						} elsif (defined($tool_data->{baits_bed})) {

						$trim_command = get_trim_galore_command(
							r1	=> $r1,
							r2	=> $r2,
							out_dir	=> $output_directory
							);

						$new_r1 = join('/', $output_directory, basename($r1));
						$new_r1 =~ s/.fastq.gz/_val_1.fq.gz/;

						$new_r2 = join('/', $output_directory, basename($r2));
						$new_r2 =~ s/.fastq.gz/_val_2.fq.gz/;


						} else {

						$trim_command = get_fastp_command(
							r1	=> $r1,
							r2	=> $r2,
							stem	=> $output_stem . '_trimmed'
							);

						$new_r1 = $output_stem . '_trimmed_R1.fastq.gz';
						$new_r2 = $output_stem . '_trimmed_R2.fastq.gz';
						}

					$smp_data_out->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R1} = $new_r1;
					$smp_data_out->{$patient}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R2} = $new_r2;

					# check if this should be run
					if ('Y' eq missing_file($new_r2)) {

						# record command (in log directory) and then run job
						print $log "      >> Submitting job to run AdapterTrim...\n";

						$run_script = write_script(
							log_dir	=> $log_directory,
							name	=> 'run_trim_adapters_' . $stem,
							cmd	=> $trim_command,
							modules	=> [$trim_tool],
							max_time	=> $parameters->{trim}->{time},
							mem		=> $parameters->{trim}->{mem},
							hpc_driver	=> $args{hpc_driver},
							extra_args	=> [$hpc_group]
							);

						$run_id = submit_job(
							jobname		=>'run_trim_adapters_' . $stem,
							shell_command	=> $run_script,
							hpc_driver	=> $args{hpc_driver},
							dry_run		=> $args{dry_run},
							log_file	=> $log
							);

						push @all_jobs, $run_id;
						} else {
						print $log "      >> Skipping trim step because output already exists.\n";
						}
					}
				}
			}
		}

	local $YAML::Indent = 4;
	DumpFile($output_yaml, $smp_data_out);

	# if this is not a dry run OR there are jobs to assess (run or resumed with jobs submitted) then
	# collect job metrics (exit status, mem, run time)
	unless ( ($args{dry_run}) || (scalar(@all_jobs) == 0) ) {

		# collect job metrics
		my $collect_metrics = collect_job_stats(
			job_ids		=> join(',', @all_jobs),
			outfile		=> $outfile,
			hpc_driver	=> $args{hpc_driver}
			);

		$run_script = write_script(
			log_dir => $log_directory,
			name    => 'output_job_metrics_' . $run_count,
			cmd     => $collect_metrics,
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
my ($data_config, $tool_config, $output_directory, $output_config);
my $hpc_driver = 'slurm';
my ($dry_run, $help, $no_wait);

# read in command line arguments
GetOptions(
	'h|help'	=> \$help,
	'd|data=s'	=> \$data_config,
	'o|out_dir=s'	=> \$output_directory,
	'b|out_yaml=s'	=> \$output_config,
	't|tool=s'	=> \$tool_config,
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
		"\t--out_yaml|-b\t<string> path to output yaml (listing trimmed fastq files)",
		"\t--cluster|-c\t<string> cluster scheduler (default: slurm)",
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
	dry_run			=> $dry_run,
	no_wait			=> $no_wait
	);
