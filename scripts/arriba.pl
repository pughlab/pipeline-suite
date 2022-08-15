#/!/usr/bin/env perl
### arriba.pl ######################################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use version;
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

# define some global variables
our ($reference, $gtf, $blacklist, $known_fusions, $protein_domains) = undef;

####################################################################################################
# version       author	  	comment
# 1.0		sprokopec       script to run Arriba on RNASeq data

### USAGE ##########################################################################################
# arriba.pl -t tool_config.yaml -c data_config.yaml -o /path/to/output/dir -c slurm --remove --dry_run
#
# where:
#	-t (tool_config.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data_config.yaml) contains sample information (YAML file containing paths to FASTQ files)
#	-o (/path/to/output/dir) indicates tool-specific output directory
#	-c indicates hpc driver (ie, slurm)
#	--remove remove intermediate files?
#	--dry_run is this a dry run?

### SUBROUTINES ####################################################################################
# format command for STAR
sub get_star_command {
	my %args = (
		r1		=> undef,
		r2		=> undef,
		reference_dir	=> undef,
		@_
		);

	my $star_command = join(' ',
		'STAR --runMode alignReads',
		# basic options
		'--genomeDir', $args{reference_dir},
		'--genomeLoad NoSharedMemory',
		'--runThreadN 1',
		'--readFilesCommand zcat',
		'--readFilesIn', $args{r1}, $args{r2},
		# output options
		'--outSAMtype BAM Unsorted',
		'--outSAMunmapped Within',
		'--outFilterMultimapNmax 50',
		'--peOverlapNbasesMin 10',
		'--alignSplicedMateMapLminOverLmate 0.5',
		'--alignSJstitchMismatchNmax 5 -1 5 5',
		'--chimSegmentMin 10',
		'--chimOutType WithinBAM HardClip',
		'--chimJunctionOverhangMin 10',
		'--chimScoreDropMax 30',
		'--chimScoreJunctionNonGTAG 0',
		'--chimScoreSeparation 1',
		'--chimSegmentReadGapMax 3',
		'--chimMultimapNmax 50'
		);

	return($star_command);
	}

# format command to run Arriba:
sub get_arriba_virus_command {
	my %args = (
		output		=> undef,
		input_bam	=> undef,
		@_
		);

	my $arriba_command = 'DIRNAME=$(which arriba | xargs dirname)';

	$arriba_command .= "\n\n" . join(' ',
		'$DIRNAME/scripts/quantify_virus_expression.sh',
		$args{input_bam},	 # Aligned.out.bam
		$args{output}
		);

	return($arriba_command);
	}

# format command to run Arriba
sub get_arriba_command {
	my %args = (
		output_stem	=> undef,
		input_bam	=> undef,
		chimeric	=> undef,
		@_
		);

	my $arriba_command = 'DIRNAME=$(which arriba | xargs dirname)';

	$arriba_command .= "\n\n" . join(' ',
		'arriba',
		'-x', $args{input_bam},
		'-a', $reference,
		'-g', $gtf,
		'-o', $args{output_stem} . '.tsv',
		'-O', $args{output_stem} . '_discarded.tsv',
		'-b', '$DIRNAME/database/' . $blacklist,
		'-p', '$DIRNAME/database/' . $protein_domains
		);

	if ('' ne $known_fusions) {
		$arriba_command .= ' ' . join(' ',
			'-k', '$DIRNAME/database/' . $known_fusions,
			'-t', '$DIRNAME/database/' . $known_fusions
			);
		}

	if (defined($args{chimeric})) {
		$arriba_command .= " -c $args{chimeric}"; # Chimeric.out.sam
		}

	return($arriba_command);
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
		print "Initiating Arriba pipeline...\n";
		}

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'arriba');

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_Arriba_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_Arriba_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running Arriba pipeline.\n";
	print $log "\n  Tool config used: $tool_config";

	my $arriba_ref;
	if (defined($tool_data->{arriba}->{star_reference})) {
		$arriba_ref = $tool_data->{arriba}->{star_reference};
		} else {
		$arriba_ref = $tool_data->{star_reference_dir};
		}

	print $log "\n    STAR reference directory: $arriba_ref";
	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	if (defined($tool_data->{arriba}->{reference})) {
		$reference	= $tool_data->{arriba}->{reference};
		} else {
		$reference	= $tool_data->{reference};
		}

	if (defined($tool_data->{arriba}->{gtf})) {
		$gtf 		= $tool_data->{arriba}->{gtf};
		} else {
		$gtf 		= $tool_data->{reference_gtf};
		}

	if ( ('hg38' eq $tool_data->{ref_type}) || ('GRCh38' eq $tool_data->{ref_type}) ) {
		$blacklist = 'blacklist_hg38_GRCh38_*.tsv.gz';
		$known_fusions = 'known_fusions_hg38_GRCh38_*.tsv.gz';
		$protein_domains = 'protein_domains_hg38_GRCh38_*.gff3';
		} elsif ( ('hg19' eq $tool_data->{ref_type}) || ('GRCh37' eq $tool_data->{ref_type}) ) {
		$blacklist = 'blacklist_hg19_hs37d5_GRCh37_*.tsv.gz';
		$known_fusions = 'known_fusions_hg19_hs37d5_GRCh37_*.tsv.gz';
		$protein_domains = 'protein_domains_hg19_hs37d5_GRCh37_*.gff3';
		}

	# set tools and versions
	my $star;
	if (defined($tool_data->{arriba}->{star_version})) {
		$star	= 'STAR/' . $tool_data->{arriba}->{star_version};
		} else {
		$star	= 'STAR/' . $tool_data->{star_version};
		}
	
	my $arriba	= 'arriba/' . $tool_data->{arriba_version};
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};
	my $perl	= 'perl/5.30.0';
	my $r_version	= 'R/' . $tool_data->{r_version};

	# if this is an older version, not all annotation files are provided
	my $needed = version->declare('2.3.0')->numify;
	my $given = version->declare($tool_data->{arriba_version})->numify;

	if ($given < $needed) { $known_fusions = ''; }

	# get user-specified tool parameters
	my $parameters = $tool_data->{arriba}->{parameters};

	# get optional HPC group
	my $hpc_group = defined($tool_data->{hpc_group}) ? "-A $tool_data->{hpc_group}" : undef;

	### RUN ############################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	unless($args{dry_run}) {
		print "Processing " . scalar(keys %{$smp_data}) . " patients.\n";
		}

	my ($run_script, $run_id, $star_run_id, $raw_link, $should_run_final);
	my @all_jobs;

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";

		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my $cleanup_cmd = '';
		my (@final_outputs, @patient_jobs);

		my (%tumours, %normals);

		foreach my $sample (sort keys %{$smp_data->{$patient}}) {

			# if there are any samples at all, run the final step
			$should_run_final = 1;

			# clear out run_id for this sample
			my @sample_jobs;
			$star_run_id = '';
			$run_id = '';

			# process this sample
			print $log "  SAMPLE: $sample\n";

			# determine sample type
			my $type = $smp_data->{$patient}->{$sample}->{type};

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			my $raw_directory = join('/', $sample_directory, 'fastq_links');
			unless(-e $raw_directory) { make_path($raw_directory); }

			my $tmp_directory = join('/', $sample_directory, 'intermediate_files');
			unless(-e $tmp_directory) { make_path($tmp_directory); }
			$cleanup_cmd .= "\n  rm -rf $tmp_directory";

			my @lanes = keys %{$smp_data->{$patient}->{$sample}->{runlanes}};
			my (@r1_fastqs, @r2_fastqs);

			foreach my $lane ( @lanes ) {

				print $log "    LANE: $lane\n";

				# collect input files
				my $r1 = $smp_data->{$patient}->{$sample}->{runlanes}->{$lane}->{R1};
				my $r2 = $smp_data->{$patient}->{$sample}->{runlanes}->{$lane}->{R2};

				print $log "      R1: $r1\n";
				print $log "      R2: $r2\n\n";

				my @tmp = split /\//, $r1;
				my $raw_link = join('/', $raw_directory, $tmp[-1]);
				symlink($r1, $raw_link);
				$raw_link =~ s/R1/R2/;
				symlink($r2, $raw_link);

				# add to respective lists
				push @r1_fastqs, $r1;
				push @r2_fastqs, $r2;
				}

			# set up commands and output
			my $bam = join('/', $tmp_directory, 'Aligned.out.bam');
			my $output_stem = join('/', $sample_directory, 'fusions');
			my $virus_output = join('/', $sample_directory, 'virus_expression.tsv');

			my $star_cmd = "cd $tmp_directory\n\n";
			$star_cmd .= get_star_command(
				r1		=> join(',', @r1_fastqs),
				r2		=> join(',', @r2_fastqs),
				reference_dir	=> $arriba_ref
				);

			# check if this should be run
			if ('Y' eq missing_file($bam)) {

				# record command (in log directory) and then run job
				print $log "Submitting job to run STAR...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_star_' . $sample,
					cmd	=> $star_cmd,
					modules => [$star],
					max_time	=> $parameters->{star}->{time},
					mem		=> $parameters->{star}->{mem},
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				$star_run_id = submit_job(
					jobname		=> 'run_star_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @sample_jobs, $run_id;
				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;
				} else {
				print $log "Skipping STAR because output already exists...\n";
				}

			# run quantify viral expression
			my $virus_cmd = get_arriba_virus_command(
				input_bam	=> $bam,
				output		=> $virus_output
				);

			# check if this should be run
			if ('Y' eq missing_file($virus_output)) {

				# record command (in log directory) and then run job
				print $log "Submitting job to run Arriba: VirusQuantification...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_arriba_quantify_viral_expression_' . $sample,
					cmd	=> $virus_cmd,
					modules => [$star, $arriba],
					dependencies	=> $star_run_id,
					max_time	=> $parameters->{quantify_virus}->{time},
					mem		=> $parameters->{quantify_virus}->{mem},
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'run_arriba_quantify_viral_expression_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @sample_jobs, $run_id;
				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;
				} else {
				print $log "Skipping Arriba: VirusQuantification because output already exists...\n";
				}

			push @final_outputs, $virus_output;

			# run fusion detection
			my $arriba_cmd = get_arriba_command(
				input_bam	=> $bam,
				output_stem	=> $output_stem
				);

			# check if this should be run
			if ('Y' eq missing_file($output_stem . '.tsv')) {

				# record command (in log directory) and then run job
				print $log "Submitting job to run Arriba...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_arriba_' . $sample,
					cmd	=> $arriba_cmd,
					modules => [$star, $arriba],
					dependencies	=> $star_run_id,
					max_time	=> $parameters->{arriba}->{time},
					mem		=> $parameters->{arriba}->{mem},
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'run_arriba_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @sample_jobs, $run_id;
				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;
				} else {
				print $log "Skipping Arriba because output already exists...\n";
				}

			push @final_outputs, "$output_stem.tsv";

			# remove temporary directories (once per patient)
			if ($args{del_intermediates}) {

				if (scalar(@sample_jobs) == 0) {
					`$cleanup_cmd`;
					} else {

					print $log "Submitting job to clean up temporary/intermediate files...\n";

					# make sure final output exists before removing intermediate files!
					$cleanup_cmd = join("\n",
						"if [ -s " . join(" ] && [ -s ", @final_outputs) . " ]; then",
						$cleanup_cmd,
						"else",
						'echo "FINAL OUTPUT FILE is missing; not removing intermediates"',
						"fi"
						);

					$run_script = write_script(
						log_dir => $log_directory,
						name    => 'run_cleanup_' . $sample,
						cmd     => $cleanup_cmd,
						dependencies	=> join(':', @sample_jobs),
						mem		=> '256M',
						hpc_driver	=> $args{hpc_driver},
						kill_on_error	=> 0,
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_cleanup_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);
					}
				}
			}

		print $log "\nFINAL OUTPUT:\n" . join("\n  ", @final_outputs) . "\n";
		print $log "---\n";
		}

	# collect and combine results
	if ($should_run_final) {

		my $collect_results = join(' ',
			"Rscript $cwd/collect_arriba_output.R",
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
		"\t--data|-d\t<string> data (fastq) config (yaml format)",
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
