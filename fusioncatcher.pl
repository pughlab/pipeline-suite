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

my $cwd = dirname($0);
require "$cwd/shared/utilities.pl";

####################################################################################################
# version       author	  	comment
# 1.0		sprokopec       script to run FusionCatcher on RNA-Seq data

### USAGE ##########################################################################################
# fusioncatcher.pl -t tool_config.yaml -c data_config.yaml
#
# where:
#	- tool_config.yaml contains tool versions and parameters, output directory,
#	reference information, etc.
#	- data_config.yaml contains sample information (YAML file containing paths to FASTQ files)

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
		tool_config => undef,
		data_config => undef,
		@_
		);

	my $tool_config = $args{tool_config};
	my $data_config = $args{data_config};

	### PREAMBLE ######################################################################################

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'fusioncatcher');
	$tool_data->{date} = strftime "%F", localtime;

	# clean up reference_dir (aesthetic reasons only)
	$tool_data->{reference_dir} =~ s/\/$//;

	# check for resume and confirm output directories
	my ($resume, $output_directory, $log_directory) = set_output_path(tool_data => $tool_data);

	# start logging
	print "---\n";
	print "Running FusionCatcher (+ViralAlignment) pipeline.\n";
	print "\n  Tool config used: $tool_config";
	print "\n    FusionCatcher reference directory: $tool_data->{reference_dir}";
	print "\n    Output directory: $output_directory";
	print "\n  Sample config used: $data_config";
	print "\n---";

	# set tools and versions
	my $fusioncatcher = $tool_data->{tool} . '/' . $tool_data->{tool_version};
	my $perl	= 'perl/5.30.0';

	# create a file to hold job metrics
	my (@files, $run_count, $outfile, $touch_exit_status);
	if ('N' eq $tool_data->{dry_run}) {
		# initiate a file to hold job metrics (ensures that an existing file isn't overwritten by concurrent jobs)
		opendir(LOGFILES, $log_directory) or die "Cannot open $log_directory";
		@files = grep { /slurm_job_metrics/ } readdir(LOGFILES);
		$run_count = scalar(@files) + 1;
		closedir(LOGFILES);

		$outfile = $log_directory . '/slurm_job_metrics_' . $run_count . '.out';
		$touch_exit_status = system("touch $outfile");
		if (0 != $touch_exit_status) { Carp::croak("Cannot touch file $outfile"); }
		}

	### RUN ############################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	my ($run_script, $run_id, $raw_link, $final_link);
	my @all_jobs;

	# process each patient in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print "\nInitiating process for PATIENT: $patient\n";

		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my (@final_outputs, @patient_jobs);
		my $cleanup_cmd = '';

		# process each separate sample for this patient
		foreach my $sample (sort keys %{$smp_data->{$patient}}) {

			if ($sample =~ m/BC$|SK$|A$/) {
				print "\n>> SAMPLE: $sample is labelled as normal - will be skipped!\n";
				next;
				}

			print "  SAMPLE: $sample\n";

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			my $link_directory = join('/', $sample_directory, 'input_links');
			unless(-e $link_directory) { make_path($link_directory); }

			my $tmp_directory = join('/', $sample_directory, 'TEMP');
			unless(-e $tmp_directory) { make_path($tmp_directory); }
			$cleanup_cmd .= "\nrm -rf $tmp_directory";

			my @lanes = keys %{$smp_data->{$patient}->{$sample}->{runlanes}};
			my @fastqs;

			foreach my $lane ( @lanes ) {

				print "    LANE: $lane\n";

				# collect input files
				my $r1 = $smp_data->{$patient}->{$sample}->{runlanes}->{$lane}->{R1};
				my $r2 = $smp_data->{$patient}->{$sample}->{runlanes}->{$lane}->{R2};

				print "      R1: $r1\n";
				print "      R2: $r2\n\n";

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

			# reset run_id for this sample
			$run_id = '';

			## run FusionCatcher on these fastqs
			my $fusion_output = join('/', $sample_directory, 'final-list_candidate-fusion-genes.txt');

			my $fusion_cmd = get_fusion_command(
				ref_dir		=> $tool_data->{reference_dir},
				input		=> join(',', @fastqs),
				output_dir	=> $sample_directory,
				tmp_dir		=> $tmp_directory,
				ref_type	=> $tool_data->{ref_type},
				java_mem	=> $tool_data->{parameters}->{fusioncatcher}->{java_mem}
				);

			# IF THIS STEP IS SUCCESSFULLY RUN,
			# create a symlink for the final output in the TOP directory
			my $smp_output = $sample . "_candidate-fusion-genes.txt";
			my $links_cmd = join("\n",
				"cd $patient_directory",
				"ln -s $fusion_output $smp_output",
				"cd $tool_data->{output_dir}"
				);

			if (-l $smp_output) { unlink $smp_output or die "Failed to remove previous symlink: $smp_output;\n"; }
			$links_cmd .= "\nln -s $fusion_output $smp_output";

			$fusion_cmd .= "\n$links_cmd";

			# check if this should be run
			if ( ('N' eq $resume) || ('Y' eq missing_file($fusion_output))) {

				# record command (in log directory) and then run job
				print "Submitting job to run FusionCatcher...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_FusionCatcher_' . $sample,
					cmd	=> $fusion_cmd,
					modules => [$fusioncatcher],
					max_time	=> $tool_data->{parameters}->{fusioncatcher}->{time},
					mem		=> $tool_data->{parameters}->{fusioncatcher}->{mem},
					hpc_driver	=> $tool_data->{HPC_driver},
					extra_args	=> '-p himem'
					);

				$run_id = submit_job(
					jobname		=> 'run_FusionCatcher_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $tool_data->{HPC_driver},
					dry_run		=> $tool_data->{dry_run}
					);

				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;

				}
			else {
				print "Skipping FusionCatcher because output already exists...\n";
				}

			# add output from STAR-Fusion to final_outputs
			push @final_outputs, $fusion_output;
			}

		# remove temporary directories (once per patient)
		if ('Y' eq $tool_data->{del_intermediate}) {

			print "Submitting job to clean up temporary/intermediate files...\n";

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
				dependencies	=> join(',', @patient_jobs),
				max_time	=> '00:05:00',
				mem		=> '256M',
				hpc_driver	=> $tool_data->{HPC_driver}
				);

			$run_id = submit_job(
				jobname		=> 'run_cleanup_' . $patient,
				shell_command	=> $run_script,
				hpc_driver	=> $tool_data->{HPC_driver},
				dry_run		=> $tool_data->{dry_run}
				);
			}

		print "\nFINAL OUTPUT:\n" . join("\n  ", @final_outputs) . "\n";
		print "---\n";
		}

	# collect job metrics if any were run
	if ('N' eq $tool_data->{dry_run}) {

		# collect job stats
		my $collect_metrics = collect_job_stats(
			job_ids => join(',', @all_jobs),
			outfile => $outfile
			);

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'output_job_metrics_' . $run_count,
			cmd	=> $collect_metrics,
			dependencies	=> join(',', @all_jobs),
			max_time	=> '0:05:00',
			mem		=> '256M',
			hpc_driver	=> $tool_data->{HPC_driver}
			);

		$run_id = submit_job(
			jobname		=> 'output_job_metrics',
			shell_command	=> $run_script,
			hpc_driver	=> $tool_data->{HPC_driver},
			dry_run		=> 'N'
			);
		}

	# final job to output a BAM config for downstream stuff
	if ('Y' eq $tool_data->{create_output_yaml}) {

		print "Creating config yaml for output fusion calls...\n";

		my $output_yaml_cmd = join(' ',
			"perl $cwd/shared/create_final_yaml.pl",
			'-d', $output_directory,
			'-o', $output_directory . '/fusions_config.yaml',
			'-p', 'final-list_candidate-fusion-genes.txt'
			);

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'output_final_yaml',
			cmd	=> $output_yaml_cmd,,
			modules	=> ['perl'],
			dependencies	=> join(',', @all_jobs),
			max_time	=> '00:10:00',
			mem		=> '1G',
			hpc_driver	=> $tool_data->{HPC_driver}
			);

		$run_id = submit_job(
			jobname		=> 'output_final_yaml',
			shell_command	=> $run_script,
			hpc_driver	=> $tool_data->{HPC_driver},
			dry_run		=> $tool_data->{dry_run}
			);

		} else {
			print "Not creating output config yaml as requested...\n";
		}

	# finish up
	print "\nProgramming terminated successfully.\n\n";

	}

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my $tool_config;
my $data_config;

# get command line arguments
GetOptions(
	't|tool=s'	=> \$tool_config,
	'c|config=s'	=> \$data_config
	);

# do some quick error checks to ensure valid input	
if (!defined($tool_config)) { die("No tool config file defined; please provide -t | --tool (ie, tool_config.yaml)"); }
if (!defined($data_config)) { die("No data config file defined; please provide -c | --config (ie, sample_config.yaml)"); }

main(tool_config => $tool_config, data_config => $data_config);
