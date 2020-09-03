#!/usr/bin/env perl
### gatk.pl ########################################################################################
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
use List::Util qw(any);
use IO::Handle;

my $cwd = dirname($0);
require "$cwd/utilities.pl";

# define some global variables
our ($data_type, $reference, $known_1000g, $known_indels, $known_mills, $dbsnp);

####################################################################################################
# version	author	  	comment
# 1.1		sprokopec	run GATKs indel realignment and recalibration on BWA aligned bams
# 1.2		sprokopec	added functionality for processing of STAR-aligned RNA-Seq bams
# 1.3		sprokopec	minor updates for compatibility with larger pipeline
# 1.4		sprokopec	added help msg and cleaned up code
# 1.5		sprokopec	minor updates for tool config
#
### USAGE ##########################################################################################
# gatk.pl -t tool.yaml -d data.yaml -o /path/to/output/dir -b /path/to/output/yaml -c slurm --remove --dry_ryn --rna
#
# where:
#	-t (tool.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data.yaml) contains sample information (YAML file containing paths to BWA-aligned or STAR-aligned BAMs)
#	-o (/path/to/output/dir) indicates tool-specific output directory
#	-b (/path/to/output/yaml) indicates output yaml for processed BAMs
#	-c indicates hpc driver (ie, slurm)
#	--remove indicates that intermediates will be removed
#	--dry_run indicates that this is a dry run
#	--rna to indicate whether input is RNA (STAR-aligned BAMs)

### SUBROUTINES ####################################################################################
# format command to split Cigar Reads
sub get_split_command {
	my %args = (
		input		=> undef,
		output		=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $split_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $gatk_dir/GenomeAnalysisTK.jar -T SplitNCigarReads',
		'-R', $reference,
		'-I', $args{input},
		'-o', $args{output},
		'-rf ReassignOneMappingQuality -rf UnmappedRead',
		'-RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS',
		'--generate_md5'
		);

	return($split_command);
	}

# format command to run GATK RealignerTargetCreator
sub get_target_intervals_command {
	my %args = (
		input		=> undef,
		n_Ramples	=> 1,
		output		=> undef,
		intervals	=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $target_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $gatk_dir/GenomeAnalysisTK.jar -T RealignerTargetCreator',
		'-R', $reference,
		'-I', $args{input},
		'-o', $args{output},
		'-known', $known_indels,
		'-known', $known_mills
		);

	# if this is DNA data, add additional options
	if ('dna' eq $data_type) {
		$target_command = join(' ',
			$target_command,
			'--disable_auto_index_creation_and_locking_when_reading_rods -nt', $args{n_samples},
			'-dt None'
			);

		if (defined($args{intervals})) {
			$target_command = join(' ',
				$target_command,
				'--intervals', $args{intervals},
				'--interval_padding 100'
				);
			}
		}

	return($target_command);
	}

# format command to run GATK IndelRealigner
sub get_indelrealign_command {
	my %args = (
		input		=> undef,
		intervals	=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		output		=> undef,
		@_
		);

	my $realign_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $gatk_dir/GenomeAnalysisTK.jar -T IndelRealigner',
		'-I', $args{input},
		'-R', $reference,
		'-targetIntervals', $args{intervals}
		);

	# if this is DNA data, add additional options
	if ('dna' eq $data_type) {

		$realign_command = join(' ',
			$realign_command,
			'--disable_auto_index_creation_and_locking_when_reading_rods',
			'-nWayOut _realigned.bam',
			'-known', $known_indels,
			'-known', $known_mills,
			'-compress 0 -dt None'
			);

		# otherwise, it is RNA, so add these options
		} elsif ('rna' eq $data_type) {

		$realign_command = join(' ',
			$realign_command,
			'-o', $args{output},
			'--generate_md5'
			);
		}

	return($realign_command);
	}

# format command to run GATK BaseRecalibrator
sub create_recalibration_table {
	my %args = (
		input		=> undef,
		output		=> undef,
		intervals	=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $bqsr_command = join(' ', 
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $gatk_dir/GenomeAnalysisTK.jar -T BaseRecalibrator',
		'-I', $args{input},
		'-R', $reference,
		'-knownSites', $known_1000g,
		'-knownSites', $dbsnp,
		'-o', $args{output}
		);

	# if this is DNA data, add additional options
	if ('dna' eq $data_type) {

		$bqsr_command = join(' ',
			$bqsr_command,
			'--disable_auto_index_creation_and_locking_when_reading_rods -nct 8',
			'-rf BadCigar',
			'--covariate ReadGroupCovariate',
			'--covariate QualityScoreCovariate',
			'--covariate CycleCovariate',
			'--covariate ContextCovariate',
			'-dt None'
			);

		if (defined($args{intervals})) {
			$bqsr_command = join(' ',
				$bqsr_command,
				'--intervals', $args{intervals},
				'--interval_padding 100'
				);
			}
		}

	return($bqsr_command);
	}

# format command to run GATK PrintReads
sub create_recalibrated_bam {
	my %args = (
		input		=> undef,
		bqsr		=> undef,
		output		=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $recal_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $gatk_dir/GenomeAnalysisTK.jar -T PrintReads',
		'-I', $args{input},
		'-R', $reference,
		'-BQSR', $args{bqsr},
		'-o', $args{output},
		'--generate_md5'
		);

	# if this is DNA data, add additional options
	if ('dna' eq $data_type) {

		$recal_command = join(' ',
			$recal_command,
			'--disable_auto_index_creation_and_locking_when_reading_rods -nct 8',
			'-rf BadCigar -dt None'
			);
		}

	return($recal_command);
	}

### RUN ############################################################################################
sub main {
	my %args = (
		tool_config		=> undef,
		data_config		=> undef,
		data_type		=> undef,
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
	$data_type = $args{data_type};

	### PREAMBLE ######################################################################################

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(
		tool_data	=> $tool_data_orig,
		pipeline	=> 'gatk',
		data_type	=> $data_type
		);

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_GATK_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_GATK_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	if ('dna' eq $data_type) {
		print $log "Running GATK pipeline for BWA-aligned DNA data.\n";
		} elsif ('rna' eq $data_type) {
		print $log "Running GATK pipeline for STAR-aligned RNA-Seq data.\n";
		}
	print $log "\n  Tool config used: $tool_config";
	print $log "\n    Reference: $tool_data->{reference}";

	$reference = $tool_data->{reference};
	if ('hg38' eq $tool_data->{ref_type}) {

		print $log "\n      Using GATK's hg38bundle files: /cluster/tools/data/genomes/human/hg38/hg38bundle/";
		$known_1000g	= '/cluster/tools/data/genomes/human/hg38/hg38bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz';
		$known_indels	= '/cluster/tools/data/genomes/human/hg38/hg38bundle/Homo_sapiens_assembly38.known_indels.vcf.gz';
		$known_mills	= '/cluster/tools/data/genomes/human/hg38/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz';
		$dbsnp		= '/cluster/tools/data/genomes/human/hg38/hg38bundle/dbsnp_144.hg38.vcf.gz';

		} elsif ('hg19' eq $tool_data->{ref_type}) {

		print $log "\n      Using hg19 variant calling files: /cluster/tools/data/genomes/human/hg19/variantcallingdata/";
		$known_1000g	= '/cluster/tools/data/genomes/human/hg19/variantcallingdata/1000G_phase1.snps.high_confidence.hg19.vcf';
		$known_indels	= '/cluster/tools/data/genomes/human/hg19/variantcallingdata/1000G_phase1.indels.hg19.vcf';
		$known_mills	= '/cluster/tools/data/genomes/human/hg19/variantcallingdata/Mills_and_1000G_gold_standard.indels.hg19.vcf';
		$dbsnp		= '/cluster/tools/data/genomes/human/hg19/variantcallingdata/dbsnp_138.hg19.vcf';

		}

	if (defined($tool_data->{dbsnp})) {
		print $log "\n      dbSNP: $tool_data->{dbsnp}";
		$dbsnp = $tool_data->{dbsnp};
		} else {
		print $log "\n      dbSNP: using default provided in above path (v138 for hg19; v144 for hg38)";
		}

	if (defined($tool_data->{intervals_bed})) {
		print $log "\n    Target intervals (exome): $tool_data->{intervals_bed}";
		}

	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	# set tools and versions
	my $gatk	= 'gatk/' . $tool_data->{gatk_version};
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};
	my $picard	= 'picard/' . $tool_data->{picard_version};

	# get user-specified tool parameters
	my $parameters = $tool_data->{gatk}->{parameters};

	### HANDLING FILES #################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	my ($run_script, $run_id_patient, $run_id_sample, $run_id_extra, $raw_link);
	my @all_jobs;

	# initiate final output yaml file
	my $output_yaml = join('/', $output_directory, 'gatk_bam_config.yaml');
	if (defined($args{output_config})) {
		$output_yaml = $args{output_config};
		}
	open (my $yaml, '>', $output_yaml) or die "Cannot open '$output_yaml' !";
	print $yaml "---\n";

	# process each patient in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";
		print $yaml "$patient:\n";

		# make a sample-specific directory
		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		# find sample IDs and paths to BAM input files
		my @normal_ids = keys %{$smp_data->{$patient}->{'normal'}};
		my @tumour_ids = keys %{$smp_data->{$patient}->{'tumour'}};
		my @normal_paths = values %{$smp_data->{$patient}->{'normal'}};
		my @tumour_paths = values %{$smp_data->{$patient}->{'tumour'}};

		my @samples = @tumour_ids;
		if (scalar(@normal_ids) > 0) { push @samples, @normal_ids; }

		# initiate some variables
		my (@final_outputs, @patient_jobs);
		$run_id_patient = '';

		# make a directory for intermediate files
		my $intermediate_directory = join('/', $patient_directory, 'intermediate_files');
		unless(-e $intermediate_directory) { make_path($intermediate_directory); }

		# make a TEMP directory
		my $tmp_directory = join('/', $patient_directory, 'TEMP');
		unless(-e $tmp_directory) { make_path($tmp_directory); }
		my $cleanup_cmd = "rm -rf $tmp_directory";

		# make a directory to link input files
		my $raw_directory = join('/', $patient_directory, 'bam_links');
		unless(-e $raw_directory) { make_path($raw_directory); }

		# create symlinks for the input files
		my @input_bams;
		foreach my $bam (@normal_paths) {
			my @tmp = split /\//, $bam;
			$raw_link = join('/', $raw_directory, $tmp[-1]);
			symlink($bam, $raw_link);
			push @input_bams, $tmp[-1];
			}

		foreach my $bam (@tumour_paths) {
			my @tmp = split /\//, $bam;
			$raw_link = join('/', $raw_directory, $tmp[-1]);
			symlink($bam, $raw_link);
			push @input_bams, $tmp[-1];
			}

		## for DNA, indel realigner target creation and indel realigner use all patient input files
		my ($input_string, $target_intervals, $stage1_cmd, $stage2_cmd, $java_check);
		my @realign_bams_dna;

		if ('dna' eq $data_type) {

			# combine input paths to a single string
			if (scalar(@normal_paths) > 0) {
				$input_string .=  join(' -I ', @normal_paths);
				}
			if ( (scalar(@normal_paths) > 0) & (scalar(@tumour_paths) > 0) ) {
				$input_string .= ' -I ';
				}
			if (scalar(@tumour_paths) > 0) {
				$input_string .= join(' -I ', @tumour_paths);
				}

			## RealignerTargetCreator
			$target_intervals = join('/', $intermediate_directory, $patient . '_target.intervals');
			$stage1_cmd = get_target_intervals_command(
				input		=> $input_string,
				n_samples	=> scalar(@input_bams),
				output		=> $target_intervals,
				intervals	=> $tool_data->{intervals_bed},
				java_mem	=> $parameters->{target_creator}->{java_mem},
				tmp_dir		=> $tmp_directory
				);

			$stage1_cmd .= "\n" . check_java_output(
				extra_cmd => "md5sum $target_intervals > $target_intervals.md5"
				);

			# check if this should be run
			if ('Y' eq missing_file($target_intervals . '.md5')) {

				# record command (in log directory) and then run job
				print $log "Submitting job for RealignerTargetCreator...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_indel_realigner_target_creator_' . $patient,
					cmd	=> $stage1_cmd,
					modules	=> [$gatk],
					max_time	=> $parameters->{target_creator}->{time},
					mem		=> $parameters->{target_creator}->{mem},
					cpus_per_task	=> scalar(@input_bams),
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id_patient = submit_job(
					jobname		=> 'run_indel_realigner_target_creator_' . $patient,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $run_id_patient;
				push @all_jobs, $run_id_patient;
				}
			else {
				print $log "Skipping RealignerTargetCreator because this has already been completed!\n";
				}

			## IndelRealigner
			$stage2_cmd = get_indelrealign_command(
				input		=> $input_string,
				intervals	=> $target_intervals,
				java_mem	=> $parameters->{realign}->{java_mem},
				tmp_dir		=> $tmp_directory
				);

			$stage2_cmd = "cd $intermediate_directory;\n$stage2_cmd;";

			my ($outbam, $md5_cmd);
			foreach my $inbam (@input_bams) {
				$outbam = $inbam; 
				$outbam =~ s/.bam/_realigned.bam/; 
				$md5_cmd .= "\n  md5sum $outbam > $outbam.md5;";
				push @realign_bams_dna, join('/', $intermediate_directory, $outbam);

				$outbam = join('/', $intermediate_directory, $outbam);
				my $outbai = $outbam;
				$outbai =~ s/bam$/bai/;
				$cleanup_cmd .= ";\nrm " . $outbam;
				$cleanup_cmd .= ";\nrm " . $outbai;
				}

			# this is a java-based command, so run a final check
			$java_check = "samtools quickcheck $realign_bams_dna[0]";
			$java_check .= "\n" . check_java_output(
				extra_cmd => $md5_cmd
				);

			$stage2_cmd .= "\n$java_check";

			# check if this should be run
			if ('Y' eq missing_file($realign_bams_dna[-1] . '.md5')) {

				# record command (in log directory) and then run job
				print $log "Submitting job for IndelRealigner...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_indel_realigner_' . $patient,
					cmd	=> $stage2_cmd,
					modules	=> [$gatk, $samtools],
					dependencies	=> $run_id_patient,
					max_time	=> $parameters->{realign}->{time},
					mem		=> $parameters->{realign}->{mem},
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id_patient = submit_job(
					jobname		=> 'run_indel_realigner_' . $patient,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $run_id_patient;
				push @all_jobs, $run_id_patient;
				}
			else {
				print $log "Skipping IndelRealigner because this has already been completed!\n";
				}
			}

		# Run per-sample steps (BQSR for DNA, all for RNA)
		my (%tumours, %normals);

		foreach my $sample (@samples) {

			$run_id_sample = '';

			# determine sample type
			my $type;
			if ( (any { $_ =~ m/$sample/ } @normal_ids) ) {
				$type = 'normal';
				} else {
				$type = 'tumour';
				}

			# initiate some variables
			my ($realigned_bam, $realigned_bai);

			if ('rna' eq $data_type) {

				print $log "  SAMPLE: $sample\n\n";

				my $aligned_bam = $smp_data->{$patient}->{$type}->{$sample};

				## first, split cigar reads
				my $split_bam = join('/', $intermediate_directory, $sample . '_split.bam');
				my $split_bai = join('/', $intermediate_directory, $sample . '_split.bai');

				my $split_cmd = get_split_command(
					input		=> $aligned_bam,
					output		=> $split_bam,
					java_mem	=> $parameters->{split_cigar}->{java_mem},
					tmp_dir		=> $tmp_directory
					);

				# this is a java-based command, so run a final check
				$java_check = "samtools quickcheck $split_bam";
				$java_check .= "\n" . check_java_output();

				$split_cmd .= "\n$java_check";

				$cleanup_cmd .= ";\nrm " . $split_bam;
				$cleanup_cmd .= ";\nrm " . $split_bai;

				# check if this should be run
				if ('Y' eq missing_file($split_bam . '.md5')) {

					# record command (in log directory) and then run job
					print $log "Submitting job for SplitNCigarReads...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_split_cigar_' . $sample,
						cmd	=> $split_cmd,
						modules	=> [$gatk, $samtools],
						max_time	=> $parameters->{split_cigar}->{time},
						mem		=> $parameters->{split_cigar}->{mem},
						hpc_driver	=> $args{hpc_driver}
						);

					$run_id_sample = submit_job(
						jobname		=> 'run_split_cigar_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id_sample;
					push @all_jobs, $run_id_sample;
					}
				else {
					print $log "Skipping SplitNCigarReads because this has already been completed!\n";
					}

				# create target intervals
				$target_intervals = join('/', $intermediate_directory, $sample . '_target.intervals');

				$stage1_cmd = get_target_intervals_command(
					input		=> $split_bam,
					output		=> $target_intervals,
					java_mem	=> $parameters->{target_creator}->{java_mem},
					tmp_dir		=> $tmp_directory
					);

				$stage1_cmd .= ";\nmd5sum " . join(' ', $target_intervals, '>', $target_intervals . '.md5');

				# check if this should be run
				if ('Y' eq missing_file($target_intervals . '.md5')) {

					# record command (in log directory) and then run job
					print $log "Submitting job for RealignerTargetCreator...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_indel_realigner_target_creator_' . $sample,
						cmd	=> $stage1_cmd,
						modules	=> [$gatk],
						dependencies	=> $run_id_sample,
						max_time	=> $parameters->{target_creator}->{time},
						mem		=> $parameters->{target_creator}->{mem},
						hpc_driver	=> $args{hpc_driver}
						);

					$run_id_sample = submit_job(
						jobname		=> 'run_indel_realigner_target_creator_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id_sample;
					push @all_jobs, $run_id_sample;
					}
				else {
					print $log "Skipping RealignerTargetCreator because this has already been completed!\n";
					}

				# perform indel realignment
				$realigned_bam = join('/', $intermediate_directory, $sample . '_split_realigned.bam');
				$realigned_bai = join('/', $intermediate_directory, $sample . '_split_realigned.bai');

				$stage2_cmd = get_indelrealign_command(
					input		=> $split_bam,
					output		=> $realigned_bam,
					intervals	=> $target_intervals,
					java_mem	=> $parameters->{realign}->{java_mem},
					tmp_dir		=> $tmp_directory
					);

				# this is a java-based command, so run a final check
				$java_check = "samtools quickcheck $realigned_bam";
				$java_check .= "\n" . check_java_output();

				$stage2_cmd .= "\n$java_check";

				$cleanup_cmd .= ";\nrm " . $realigned_bam;
				$cleanup_cmd .= ";\nrm " . $realigned_bai;

				# check if this should be run
				if ('Y' eq missing_file($realigned_bam . '.md5')) {

					# record command (in log directory) and then run job
					print $log "Submitting job for IndelRealigner...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_indel_realigner_' . $sample,
						cmd	=> $stage2_cmd,
						modules	=> [$gatk, $samtools],
						dependencies	=> $run_id_sample,
						max_time	=> $parameters->{realign}->{time},
						mem		=> $parameters->{realign}->{mem},
						hpc_driver	=> $args{hpc_driver}
						);

					$run_id_sample = submit_job(
						jobname		=> 'run_indel_realigner_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id_sample;
					push @all_jobs, $run_id_sample;
					}
				else {
					print $log "Skipping IndelRealigner because this has already been completed!\n";
					}
				}

			## BaseRecalibrator
			print $log "Performing base recalibration steps for: $sample\n";

			my $bqsr_file = join('/', $intermediate_directory, $sample . '.recal_data.grp');

			if ('dna' eq $data_type) {
				my @tmp = grep { /$sample/ } @realign_bams_dna;
				$realigned_bam = $tmp[0];
				}
 
			my $stage3_cmd = create_recalibration_table(
				input		=> $realigned_bam,
				output		=> $bqsr_file,
				intervals	=> $tool_data->{intervals_bed},
				java_mem	=> $parameters->{bqsr}->{java_mem},
				tmp_dir		=> $tmp_directory
				);

			$stage3_cmd .= "\n" . check_java_output();

			# check if this should be run
			if ('Y' eq missing_file($bqsr_file)) {

				# record command (in log directory) and then run job
				print $log "\nSubmitting job for BaseRecalibrator...";

				if ('dna' eq $data_type) {

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_base_quality_score_recalibrator_' . $sample,
						cmd	=> $stage3_cmd,
						modules	=> [$gatk],
						dependencies	=> $run_id_patient,
						max_time	=> $parameters->{bqsr}->{time}->{$type},
						mem		=> $parameters->{bqsr}->{mem},
						cpus_per_task	=> 8,
						hpc_driver	=> $args{hpc_driver}
						);
					}

				elsif ('rna' eq $data_type) {

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_base_quality_score_recalibrator_' . $sample,
						cmd	=> $stage3_cmd,
						modules	=> [$gatk],
						dependencies	=> $run_id_sample,
						max_time	=> $parameters->{bqsr}->{time}->{$type},
						mem		=> $parameters->{bqsr}->{mem},
						cpus_per_task	=> 1,
						hpc_driver	=> $args{hpc_driver}
						);
					}

				$run_id_sample = submit_job(
					jobname		=> 'run_base_quality_score_recalibrator_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $run_id_sample;
				push @all_jobs, $run_id_sample;
				}
			else {
				print $log "Skipping BaseRecalibrator because this has already been completed!\n";
				}

			## PrintReads
			my $recal_bam = join('/', $patient_directory, $sample . '_realigned_recalibrated.bam');

			my $stage4_cmd = create_recalibrated_bam(
				input		=> $realigned_bam,
				bqsr		=> $bqsr_file,
				output		=> $recal_bam,
				java_mem	=> $parameters->{recalibrate}->{java_mem},
				tmp_dir		=> $tmp_directory
				);

			if ('normal' eq $type) { $normals{$sample} = $recal_bam; }
			if ('tumour' eq $type) { $tumours{$sample} = $recal_bam; } 

			# check if this should be run
			if ('Y' eq missing_file($recal_bam . '.md5')) {

				# IF THIS FINAL STEP IS SUCCESSFULLY RUN,
				$java_check = "samtools quickcheck $recal_bam";

				$stage4_cmd .= "\n" . $java_check;

				# record command (in log directory) and then run job
				print $log "Submitting job for PrintReads (applying base recalibration)...\n";

				# determine number of cpus to request
				my $n_cpus;
				if ('dna' eq $data_type) {
					$n_cpus = 8;
					} elsif ('rna' eq $data_type) {
					$n_cpus = 1;
					}

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_apply_base_recalibration_' . $sample,
					cmd	=> $stage4_cmd,
					modules	=> [$gatk, $samtools],
					dependencies	=> $run_id_sample,
					max_time	=> $parameters->{recalibrate}->{time}->{$type},
					mem		=> $parameters->{recalibrate}->{mem},
					cpus_per_task	=> $n_cpus,
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id_sample = submit_job(
					jobname		=> 'run_apply_base_recalibration_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $run_id_sample;
				push @all_jobs, $run_id_sample;
				}
			else {
				print $log "Skipping PrintReads (apply base recalibration) because this has already been completed!\n";
				}

			push @final_outputs, $recal_bam;
			}
	
		# and finally, add the final files to the output yaml
		my $key;
		if (scalar(@tumour_ids) > 0) {
			print $yaml "    tumour:\n";
			foreach $key (keys %tumours) { print $yaml "        $key: $tumours{$key}\n"; }
			}
		if (scalar(@normal_ids) > 0) {
			print $yaml "    normal:\n";
			foreach $key (keys %normals) { print $yaml "        $key: $normals{$key}\n"; }
			}

		# clean up/remove intermediate files
		if ($args{del_intermediates}) {

			if (scalar(@patient_jobs) == 0) {
				`rm -rf $tmp_directory`;
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
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id_patient = submit_job(
					jobname		=> 'run_cleanup_' . $patient,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);
				}
			}

		print $log "FINAL OUTPUT:\n" . join("\n  ", @final_outputs) . "\n";
		print $log "---\n";
		}

	# if this is not a dry run OR there are jobs to assess (run or resumed with jobs submitted) then
	# collect job metrics (exit status, mem, run time)
	unless ( ($args{dry_run}) || (scalar(@all_jobs) == 0) ) {

		# collect job metrics
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

		$run_id_extra = submit_job(
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
				my $status = `sacct --format='State' -j $run_id_extra`;

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
					die("Final GATK accounting job: $run_id_extra finished with errors.");
					}
				}
			}
		}

	# finish up
	close $yaml;

	print $log "\nProgramming terminated successfully.\n\n";
	close $log;
	}

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($data_config, $tool_config, $output_directory, $output_config) = undef;
my $hpc_driver = 'slurm';
my ($remove_junk, $dry_run, $rna, $help, $no_wait);

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
	'no-wait'	=> \$no_wait,
	'rna'		=> \$rna
	);

if ($help) {
	my $help_msg = join("\n",
		"Options:",
		"\t--help|-h\tPrint this help message",
		"\t--data|-d\t<string> data config (yaml format)",
		"\t--tool|-t\t<string> tool config (yaml format)",
		"\t--out_dir|-o\t<string> path to output directory",
		"\t--out_yaml|-b\t<string> path to output yaml (listing GATK-processed BAMs)",
		"\t--cluster|-c\t<string> cluster scheduler (default: slurm)",
		"\t--remove\t<boolean> should intermediates be removed? (default: false)",
		"\t--dry-run\t<boolean> should jobs be submitted? (default: false)",
		"\t--no-wait\t<boolean> should we exit after job submission (true) or wait until all jobs have completed (false)? (default: false)",
		"\t--rna\t<boolean> is the input RNA (STAR-aligned BAMs)? (default: false)"
		);

	print "$help_msg\n";
	exit;
	}

# quick error checks to confirm valid arguments
if (!defined($tool_config)) { die("No tool config file defined; please provide -t | --tool (ie, tool_config.yaml)"); }
if (!defined($data_config)) { die("No data config file defined; please provide -d | --data (ie, sample_config.yaml)"); }
if (!defined($output_directory)) { die("No output directory defined; please provide -o | --out_dir"); }

$data_type = 'dna';
if ($rna) { $data_type = 'rna'; }

# run it!
main(
	tool_config 		=> $tool_config,
	data_config 		=> $data_config,
	output_directory	=> $output_directory,
	output_config		=> $output_config,
	hpc_driver		=> $hpc_driver,
	del_intermediates	=> $remove_junk,
	dry_run			=> $dry_run,
	data_type		=> $data_type,
	no_wait			=> $no_wait
	);
