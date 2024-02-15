#!/usr/bin/env perl
### star.pl ########################################################################################
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
use YAML qw(LoadFile DumpFile);
use Data::Dumper;
use IO::Handle;

my $cwd = dirname(__FILE__);
require "$cwd/utilities.pl";

# define some global variables
our $tool_version = undef;

####################################################################################################
# version       author	  	comment
# 1.0		sprokopec       script to run STAR alignment on RNASeq data
# 1.1		sprokopec	minor updates for compatibility with larger pipeline
# 1.2		sprokopec	added help message and cleaned up code
# 1.3           sprokopec       minor updates for tool config

### USAGE ##########################################################################################
# star.pl -t tool_config.yaml -d data_config.yaml -o /path/to/output/dir -c slurm --remove --dry_run
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
# format command to run alignment using STAR
sub get_star_command {
	my %args = (
		r1		=> undef,
		r2		=> undef,
		reference_dir	=> undef,
		readgroup	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $star_command = join(' ',
		'STAR --runMode alignReads',
		# basic options
		'--genomeDir', $args{reference_dir},
		'--runThreadN 1',
		'--readFilesCommand zcat',
		'--readFilesIn', $args{r1}, $args{r2},
		'--twopassMode Basic',
		# output options
	#	'--outFileNamePrefix', $stem,
		'--outTmpDir', $args{tmp_dir},
		'--outSAMtype BAM SortedByCoordinate',
		'--outSAMunmapped Within',
		'--outSAMprimaryFlag AllBestScore',
		'--outBAMsortingThreadN 1',
		'--outFilterIntronMotifs RemoveNoncanonical',
		# Note: --chimOutType SeparateSAMold has been deprecated; as of the new version of STAR (2.7.2), these are the recommended args:
		# https://github.com/STAR-Fusion/STAR-Fusion/wiki
		'--chimSegmentMin 10',
	#	'--chimOutType WithinBAM',
		'--chimJunctionOverhangMin 10',
		'--chimOutJunctionFormat 1',
		'--alignSJDBoverhangMin 10',
		'--alignMatesGapMax 100000',
		'--alignIntronMax 100000',
		'--alignSJstitchMismatchNmax 5 -1 5 5',
		'--outSAMattrRGline', $args{readgroup},
		'--chimMultimapScoreRange 3',
		'--chimScoreJunctionNonGTAG -4',
		'--chimMultimapNmax 20',
		'--chimNonchimScoreDropMin 10',
		'--peOverlapNbasesMin 12',
		'--peOverlapMMp 0.1',
		# transcript quant mode
		'--quantMode GeneCounts TranscriptomeSAM',
		# limit options
#		'--limitIObufferSize 250000000',
		'--limitBAMsortRAM 29000000000'
		);

	if ($tool_version < version->declare('2.7.9')->numify) {
		$star_command .= " --limitIObufferSize 250000000";
		} else {
		$star_command .= " --limitIObufferSize 30000000 50000000";
		}

	return($star_command);
	}

# format command to add readgroup and sort BAM
sub format_readgroup {
	my %args = (
		subject		=> undef,
		sample		=> undef,
		lane		=> undef,
		lib		=> undef,
		platform	=> 'Illumina',
		center		=> undef,
		@_
		);

	my $readgroup = "ID:patient\tSM:smp\tCN:center\tPL:platform\tLB:library\tPU:unit";
	$readgroup =~ s/patient/$args{subject}/;
	$readgroup =~ s/smp/$args{sample}/;
	$readgroup =~ s/unit/$args{lane}/;
	$readgroup =~ s/library/$args{lib}/;
	$readgroup =~ s/platform/$args{platform}/;
	$readgroup =~ s/center/$args{center}/;

	return($readgroup);
	}

# format command to index genome-sorted bam
sub get_index_bam_command {
	my %args = (
		output => undef,
		@_
		);

	my $index_command = join(' ',
		'samtools index',
		$args{output}
		);

	return($index_command);
	}

# format command to mark duplicates
sub get_markdup_command {
	my %args = (
		input		=> undef,
		output		=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		flowcell	=> 'random',
		@_
		);

	my $markdup_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $picard_dir/picard.jar MarkDuplicates',
		'INPUT=' . $args{input},
		'OUTPUT=' . $args{output},
		'METRICS_FILE=' . $args{output} . '.metrics',
		'ASSUME_SORTED=true CREATE_INDEX=true CREATE_MD5_FILE=true',
		'MAX_RECORDS_IN_RAM=100000 VALIDATION_STRINGENCY=SILENT'
		);

	if ('patterned' eq $args{flowcell}) {
		$markdup_command .= " OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500";
		}

	return($markdup_command);
	}

# format command to run RNA-SeQC
sub get_rnaseqc_cmd {
	my %args = (
		input		=> undef,
		output_dir	=> undef,
		bwa		=> '/cluster/tools/software/bwa/0.7.15/',
		reference	=> undef,
		gtf		=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $format_gtf_command = "cat $args{gtf}";
	if ($args{gtf} =~ m/.gz$/) { $format_gtf_command = "zcat $args{gtf}"; }

	$format_gtf_command .= " | grep -e '#' -e 'transcript_id'";
	$format_gtf_command .= " > $args{tmp_dir}/transcript_ids.gtf";

	my $qc_command = $format_gtf_command . "\n\n" . join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $rnaseqc_dir/RNA-SeQC.jar',
		'-bwa', $args{bwa},
		'-o', $args{output_dir},
		'-t', "$args{tmp_dir}/transcript_ids.gtf",
		'-r', $args{reference},
		'-singleEnd no',
		'-s', $args{input}
		);

	return($qc_command);
	}

### MAIN ###########################################################################################
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
	unless($args{dry_run}) {
		print "Initiating STAR pipeline...\n";
		}

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(
		tool_data => $tool_data_orig, pipeline => 'star', data_type => 'rna');

	# check version
	$tool_version = $tool_data->{star_version};
	$tool_version =~ s/\D$//;
	$tool_version = version->declare($tool_version)->numify;

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_STAR_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_STAR_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running STAR alignment pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n    STAR reference directory: $tool_data->{star_reference_dir}";
	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	# set tools and versions
	my $star_version	= 'STAR/' . $tool_data->{star_version};
	my $samtools		= 'samtools/' . $tool_data->{samtools_version};
	my $picard		= 'picard/' . $tool_data->{picard_version};
	my $rnaseqc		= 'rna_seqc/' . $tool_data->{rna_seqc_version};
	my $r_version		= 'R/' . $tool_data->{r_version};

	# get user-specified tool parameters
	my $parameters = $tool_data->{star}->{parameters};

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
	my $output_yaml = join('/', $output_directory, 'star_bam_config.yaml');
	if (defined($args{output_config})) {
		$output_yaml = $args{output_config};
		}

	# initiate output yaml objects
	my ($smp_data_out, $smp_data_markdup);

	# create sample sheet (tab-delim file with id/path/group)
	my $qc_directory = join('/', $output_directory, 'RNASeQC_' . $run_count);
	unless ( -e $qc_directory ) { make_path($qc_directory); }

	my $sample_sheet = join('/', $qc_directory, 'sample_sheet.tsv');
	open(my $fh, '>', $sample_sheet) or die "Cannot open '$sample_sheet' !";

	# add header
	print $fh "SampleID\tLocation\tGroup\n";

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";

		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my $cleanup_cmd = '';
		my (@final_outputs, @patient_jobs);

		foreach my $sample (sort keys %{$smp_data->{$patient}}) {

			print $log "\n  SAMPLE: $sample\n";

			# determine sample type
			my $type = $smp_data->{$patient}->{$sample}->{type};

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			my $raw_directory = join('/', $sample_directory, 'fastq_links');
			unless(-e $raw_directory) { make_path($raw_directory); }

			my $temp_star = join('/', $sample_directory, 'intermediate_files');
			$cleanup_cmd .= "\n  rm -rf $temp_star";

			my $tmp_directory = join('/', $sample_directory, 'TEMP');
			unless(-e $tmp_directory) { make_path($tmp_directory); }
			$cleanup_cmd .= "\n  rm -rf $tmp_directory";

			my (@libraries, @all_lanes, @r1_fastqs, @r2_fastqs);
			@libraries = keys %{$smp_data->{$patient}->{$sample}->{libraries}};

			foreach my $lib ( @libraries ) {
				print $log "\n    LIBRARY: $lib\n";
				my $lib_data = $smp_data->{$patient}->{$sample}->{libraries}->{$lib};
				my @lanes = keys %{$lib_data->{runlanes}};
				push @all_lanes, @lanes;
				foreach my $lane ( sort(@lanes) ) {
					print $log "      LANE: $lane\n";
					my $r1 = $lib_data->{runlanes}->{$lane}->{fastq}->{R1};
					my $r2 = $lib_data->{runlanes}->{$lane}->{fastq}->{R2};

					print $log "        R1: $r1\n";
					print $log "        R2: $r2\n\n";

					# create a symlink for this file
					my @tmp = split /\//, $r1;
					my $raw_link = join('/', $raw_directory, $tmp[-1]);
					symlink($r1, $raw_link);

					@tmp = split /\//, $r2;
					$raw_link = join('/', $raw_directory, $tmp[-1]);
					symlink($r2, $raw_link);

					push @r1_fastqs, $r1;
					push @r2_fastqs, $r2;
					}
				}

			# clear out run_id for this sample
			$run_id = '';

			## run STAR on these fastqs
			my $readgroup = format_readgroup(
				subject		=> $patient,
				sample		=> $sample,
				lane		=> join(',', @all_lanes),
				lib		=> join(',', @libraries),
				platform	=> $tool_data->{platform},
				center		=> $tool_data->{seq_center}
				);

			my $star = "cd $sample_directory\n";

			if ( -e $temp_star ) {
				`rm -rf $temp_star`;
				}

			$star .= get_star_command(
				r1		=> join(',', @r1_fastqs),
				r2		=> join(',', @r2_fastqs),
				reference_dir	=> $tool_data->{star_reference_dir},
				readgroup	=> $readgroup,
				tmp_dir		=> $temp_star
				);

			my $genome_bam = join('/', $sample_directory, 'Aligned.sortedByCoord.out.bam');

			# add output file to list
			if ('normal' eq $type) {
				$smp_data_out->{$patient}->{normal}->{$sample} = $genome_bam;
				} elsif ('tumour' eq $type) {
				$smp_data_out->{$patient}->{tumour}->{$sample} = $genome_bam;
				}

			# check if this should be run
			if ('Y' eq missing_file($genome_bam)) {

				# record command (in log directory) and then run job
				print $log "  >> Submitting job to run STAR...\n";
				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_star_alignment_' . $sample,
					cmd	=> $star,
					modules	=> [$star_version],
					dependencies	=> $run_id,
					max_time	=> $parameters->{star}->{time},
					mem		=> $parameters->{star}->{mem},
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				# initial test took 15 hours with 22G and 1 node
				$run_id = submit_job(
					jobname		=>'run_star_alignment_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;
				} else {
				print $log "  >> Skipping alignment step because output already exists...\n";
				}

			# index the resulting BAM and remove intermediate SAM
			my $index_cmd = get_index_bam_command(
				output => $genome_bam
				);

			# check if this should be run
			if ( ('N' eq missing_file($genome_bam)) &&  
				('Y' eq missing_file("$genome_bam.bai")) ) {

				# record command (in log directory) and then run job
				print $log "  >> Submitting job to run BAM INDEX...\n";
				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_bam_index_' . $sample,
					cmd	=> $index_cmd,
					modules	=> [$samtools],
					dependencies	=> $run_id,
					max_time	=> '12:00:00',
					mem		=> '2G',
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				my $idx_run_id = submit_job(
					jobname		=>'run_bam_index_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $idx_run_id;
				push @all_jobs, $idx_run_id;
				} else {
				print $log "  >> Skipping index step because output already exists...\n";
				}

			## mark duplicates
			my $dedup_bam = join('/', $patient_directory, $sample . '_sorted_markdup.bam');

			if ('N' eq $parameters->{markdup}->{run}) {

				print $fh "$sample\t$genome_bam\tRNASeq\n";
				push @final_outputs, $genome_bam;

				} else {

				print $fh "$sample\t$dedup_bam\tRNASeq\n";

				my $markdup_cmd = get_markdup_command(
					input		=> $genome_bam,
					output		=> $dedup_bam,
					java_mem	=> $parameters->{markdup}->{java_mem},
					tmp_dir		=> $tmp_directory,
					flowcell	=> $tool_data->{flowcell_type}
					);
		
				# add output file to list
				if ('normal' eq $type) {
					$smp_data_out->{$patient}->{normal}->{$sample} = $dedup_bam;
					} elsif ('tumour' eq $type) {
					$smp_data_out->{$patient}->{tumour}->{$sample} = $dedup_bam;
					}

				# check if this should be run
				if ('Y' eq missing_file($dedup_bam . '.md5')) {

					# record command (in log directory) and then run job
					print $log "  >> Submitting job to merge lanes and mark duplicates...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_picard_markdup_' . $sample,
						cmd	=> $markdup_cmd,
						modules	=> [$picard, $samtools],
						dependencies	=> $run_id,
						max_time	=> $parameters->{markdup}->{time},
						mem		=> $parameters->{markdup}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_picard_markdup_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log "  >> Skipping mark duplicate step because output already exists...\n";
					}

				push @final_outputs, $dedup_bam;
				}
			}

		# once per patient, run cleanup
		if ( ($args{del_intermediates}) && (scalar(@patient_jobs) > 0) ) {

			print $log ">> Submitting job to clean up temporary/intermediate files...\n";

			# make sure final output exists before removing intermediate files!
			$cleanup_cmd = join("\n",
				"if [ -s " . join(" ] && [ -s ", @final_outputs) . " ]; then",
				$cleanup_cmd,
				"else",
				'echo "One or more FINAL OUTPUT FILES is missing; not removing intermediates"',
				"fi"
				);

			# if all lane alignments + mark dup are successful, clean up tmp directories
			$run_script = write_script(
				log_dir	=> $log_directory,
				name	=> 'run_cleanup_' . $patient,
				cmd	=> $cleanup_cmd,
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

		print $log "\nFINAL OUTPUT:\n  " . join("\n  ", @final_outputs) . "\n";
		print $log "---\n";
		}

	close $fh;

	# get command for RNASeQC
	my $qc_cmd = get_rnaseqc_cmd(
		input		=> $sample_sheet,
		output_dir	=> $qc_directory,
		bwa		=> $tool_data->{bwa_path},
		reference	=> $tool_data->{reference},
		gtf		=> $tool_data->{reference_gtf},
		java_mem	=> $parameters->{rna_seqc}->{java_mem},
		tmp_dir		=> $qc_directory
		);

	$qc_cmd .= ";\n\ncd $qc_directory;";
	$qc_cmd .= "\nif [ -s metrics.tsv ]; then";
	$qc_cmd .= "\n  find . -name '*tmp.txt*' -exec rm {} ". '\;';
	$qc_cmd .= "\n  rm transcript_ids.gtf;";
	$qc_cmd .= "\nfi";

	# record command (in log directory) and then run job
	print $log "\n>> Submitting job for RNA-SeQC...\n";

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'run_rna_seqc_cohort',
		cmd	=> $qc_cmd,
		modules	=> [$rnaseqc],
		dependencies	=> join(':', @all_jobs),
		max_time	=> $parameters->{rna_seqc}->{time},
		mem		=> $parameters->{rna_seqc}->{mem},
		hpc_driver	=> $args{hpc_driver},
		extra_args	=> [$hpc_group]
		);

	$run_id = submit_job(
		jobname		=> 'run_rna_seqc_cohort',
		shell_command	=> $run_script,
		hpc_driver	=> $args{hpc_driver},
		dry_run		=> $args{dry_run},
		log_file	=> $log
		);

	push @all_jobs, $run_id;

	# collect and combine results
	my $collect_results = join(' ',
		"Rscript $cwd/collect_rnaseqc_output.R",
		'-d', $output_directory,
		'-p', $tool_data->{project_name}
		);

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'combine_and_format_qc_results',
		cmd	=> $collect_results,
		modules		=> [$r_version],
		dependencies	=> $run_id,
		max_time	=> $parameters->{combine_results}->{time},
		mem		=> $parameters->{combine_results}->{mem},
		hpc_driver	=> $args{hpc_driver},
		extra_args	=> [$hpc_group]
		);

	$run_id = submit_job(
		jobname		=> 'combine_and_format_qc_results',
		shell_command	=> $run_script,
		hpc_driver	=> $args{hpc_driver},
		dry_run		=> $args{dry_run},
		log_file	=> $log
		);

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

	# output final data config
	local $YAML::Indent = 4;
	DumpFile($output_yaml, $smp_data_out);

	# finish up
	print $log "\nProgramming terminated successfully.\n\n";
	close $log;
	}

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($data_config, $tool_config, $output_directory, $output_config);
my $hpc_driver = 'slurm';
my ($remove_junk, $dry_run, $help, $no_wait);

# read in command line arguments
GetOptions(
	'h|help'	=> \$help,
	'd|data=s'	=> \$data_config,
	'o|out_dir=s'	=> \$output_directory,
	'b|out_yaml=s'	=> \$output_config,
	't|tool=s'	=> \$tool_config,
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
