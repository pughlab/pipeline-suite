#!/usr/bin/env perl
### haplotype_caller.pl ############################################################################
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
use List::Util qw(any);
use IO::Handle;

my $cwd = dirname(__FILE__);
require "$cwd/utilities.pl";

# define some global variables
our ($reference, $dbsnp);

####################################################################################################
# version       author		comment
# 1.1		sprokopec       script to run HaplotypeCaller and Oncotator on RNASeq data
# 1.2		sprokopec	updates to move indel realign/recalibration to gatk.pl AND
# 				change oncotator/funcotator to VEP (vcf2maf)
# 1.3		sprokopec	added DNA-compatible options
# 1.4		sprokopec	minor updates for compatibility with larger pipeline
# 1.5		sprokopec	added help msg and cleaned up code
# 1.6           sprokopec       minor updates for tool config

### USAGE ##########################################################################################
# haplotype_caller.pl -t tool.yaml -d data.yaml -o /path/to/output/dir -h slurm -r Y -n Y --rna
#
# where:
# 	-t (tool.yaml) contains tool versions and parameters, reference information, etc.
# 	-d (data.yaml) contains sample information (YAML file containing paths to BWA-aligned,
# 	or STAR-aligned BAMs, post gatk processing)
# 	-o (/path/to/output/dir) indicates tool-specific output directory
# 	-c indicates hpc driver (ie, slurm)
# 	--remove indicates that intermediates will be removed
# 	--dry_run indicates that this is a dry run
# 	--rna to indicate whether input is RNA (STAR-aligned BAMs)

### DEFINE SUBROUTINES #############################################################################
# format command to run GATK HaplotypeCaller
sub get_haplotype_command {
	my %args = (
		bam		=> undef,
		output		=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		data_type	=> undef,
		intervals	=> undef,
		@_
		);

	my $gatk_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $gatk_dir/GenomeAnalysisTK.jar -T HaplotypeCaller',
		'-R', $reference,
		'-I', $args{bam},
		'-o', $args{output}
		);

	if ('rna' eq $args{data_type}) {
		$gatk_command .= ' ' . join(' ',
			'-dontUseSoftClippedBases',
			'-stand_call_conf 20.0'
			);
		} elsif ('dna' eq $args{data_type}) {
		$gatk_command .= ' ' . join(' ',
			'-stand_call_conf 30.0',
			'--emitRefConfidence GVCF',
			'-variant_index_type LINEAR -variant_index_parameter 128000',
			'--dbsnp', $dbsnp
			);

		if (defined($args{intervals})) {
			$gatk_command .= ' ' . join(' ',
				'--intervals', $args{intervals},
				'--interval_padding 100'
				);
			}
		}

	return($gatk_command);
	}

# format command to run variant filter
sub get_filter_command {
	my %args = (
		input		=> undef,
		reference	=> undef,
		output		=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $gatk_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $gatk_dir/GenomeAnalysisTK.jar -T VariantFiltration',
		'-R', $reference,
		'-V', $args{input},
		'-window 35',
		'-cluster 3',
		'-filterName FS -filter "FS > 30.0"',
		'-filterName QD -filter "QD < 2.0"',
		'-o', $args{output}
		);

	return($gatk_command);
	}

# format command to combine gVCFs (per batch)
sub get_combine_gvcf_command {
	my %args = (
		input		=> undef,
		output		=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $gatk_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $gatk_dir/GenomeAnalysisTK.jar -T CombineGVCFs',
		'-R', $reference,
		 $args{input},
		'-o', $args{output}
		);

	return($gatk_command);
	}

### MAIN ###########################################################################################
sub main{
	my %args = (
		tool_config		=> undef,
		data_config		=> undef,
		data_type		=> undef,
		output_directory	=> undef,
		hpc_driver		=> undef,
		del_intermediates	=> undef,
		dry_run			=> undef,
		no_wait			=> undef,
		@_
		);

	my $tool_config = $args{tool_config};
	my $data_config = $args{data_config};
	my $data_type	= $args{data_type};

	### PREAMBLE ######################################################################################
	unless($args{dry_run}) {
		print "Initiating HaplotypeCaller pipeline...\n";
		}

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

	my $log_directory = join('/', $output_directory, 'logs', 'RUN_HAPLOTYPE_CALLER');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_VariantCalling_pipeline.log');

	# create a file to hold job metrics
	my (@files, $run_count, $outfile, $touch_exit_status);
	unless ($args{dry_run}) {
		# initiate a file to hold job metrics
		opendir(LOGFILES, $log_directory) or die "Cannot open $log_directory";
		@files = grep { /slurm_job_metrics/ } readdir(LOGFILES);
		$run_count = scalar(@files) + 1;
		closedir(LOGFILES);

		$outfile = $log_directory . '/slurm_job_metrics_' . $run_count . '.out';
		$touch_exit_status = system("touch $outfile");
		if (0 != $touch_exit_status) { Carp::croak("Cannot touch file $outfile"); }

		$log_file = join('/', $log_directory, 'run_VariantCalling_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running Variant Call (HaplotypeCaller) pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n    Reference used: $tool_data->{reference}";

	$reference = $tool_data->{reference};
	
	if (defined($tool_data->{dbsnp})) {
		print $log "\n      dbSNP: $tool_data->{dbsnp}";
		$dbsnp = $tool_data->{dbsnp};
		} else {
		print $log "\n      No dbSNP provided.";
		}

	# is this WGS (will split by chromosome)
	my $is_wgs = 0;
	if ('wgs' eq $tool_data->{seq_type}) {
		$is_wgs = 1;
		$tool_data->{targets_bed} = undef;
		}	

	if (defined($tool_data->{targets_bed})) {
		print $log "\n    Target intervals: $tool_data->{targets_bed}";
		}

	# over which intervals should this be run?
	my $string;
	if ('targeted' eq $tool_data->{seq_type}) {
		$string = 'panel';
		} elsif (defined($tool_data->{haplotypecaller}->{chromosomes})) {
		$string = $tool_data->{haplotypecaller}->{chromosomes}; 
		} elsif ( ('hg38' eq $tool_data->{ref_type}) || ('hg19' eq $tool_data->{ref_type})) {
		$string = 'chr' . join(',chr', 1..22) . ',chrX,chrY';
		} elsif ( ('GRCh37' eq $tool_data->{ref_type}) || ('GRCh38' eq $tool_data->{ref_type})) {
		$string = join(',', 1..22) . ',X,Y';
		} else {
		# if no chromosomes can be determined, run as a whole (very very slow!)
		print $log "  >> Could not determine chromosomes to run\n";
		print $log "  >> Will run full genome, however this will be very very slow!\n";
		$string = 'genome';
		}

	my @chroms = split(',', $string);

	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	# get user-specified tool parameters
	my $parameters = $tool_data->{haplotype_caller}->{parameters};

	# set tools and versions
	my $perl	= 'perl/' . $tool_data->{perl_version};
	my $gatk	= 'gatk/' . $tool_data->{gatk_version};
	my $picard	= 'picard/' . $tool_data->{picard_version};
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};
	my $vcftools	= 'vcftools/' . $tool_data->{vcftools_version};
	my $r_version	= 'R/' . $tool_data->{r_version};

	# only used for RNA-seq
	my $vcf2maf = undef;
	if (defined($tool_data->{vcf2maf_version})) {
		$vcf2maf = 'vcf2maf/' . $tool_data->{vcf2maf_version};
		$parameters->{annotate}->{vcf2maf_path} = undef;
		}

	# get optional HPC group
	my $hpc_group = defined($tool_data->{hpc_group}) ? "-A $tool_data->{hpc_group}" : undef;

	### RUN ###########################################################################################
	# get sample data
	my $smp_data = LoadFile($data_config);

	unless($args{dry_run}) {
		print "Processing " . scalar(keys %{$smp_data}) . " patients.\n";
		}

	# if multiple chromosomes are to be run (separately):
	my $chr_file = join('/', $output_directory, 'chromosome_list.txt');
	if (scalar(@chroms) > 1) {
		open (my $chr_list, '>', $chr_file) or die "Could not open $chr_file for writing.";	
		foreach my $chrom ( @chroms ) {
			print $chr_list "$chrom\n";
			}
		}

	my ($run_script, $run_id, $link, $java_check, $cleanup_cmd_dna, $cleanup_cmd_rna);
	my (@all_jobs, @gvcfs);

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient";

		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my $tmp_directory = join('/', $patient_directory, 'TEMP');
		unless(-e $tmp_directory) { make_path($tmp_directory); }

		# indicate this should be removed at the end
		$cleanup_cmd_dna .= "\nrm -rf $tmp_directory";
		$cleanup_cmd_rna = "rm -rf $tmp_directory";

		my $link_directory = join('/', $patient_directory, 'bam_links');
		unless(-e $link_directory) { make_path($link_directory); }

		my @normal_ids = keys %{$smp_data->{$patient}->{'normal'}};
		my @tumour_ids = keys %{$smp_data->{$patient}->{'tumour'}};

		my @sample_ids = @tumour_ids;
		push @sample_ids, @normal_ids;

		# create some symlinks
		foreach my $normal (@normal_ids) {
			my @tmp = split /\//, $smp_data->{$patient}->{normal}->{$normal};
			$link = join('/', $link_directory, $tmp[-1]);
			symlink($smp_data->{$patient}->{normal}->{$normal}, $link);
			}
		foreach my $tumour (@tumour_ids) {
			my @tmp = split /\//, $smp_data->{$patient}->{tumour}->{$tumour};
			$link = join('/', $link_directory, $tmp[-1]);
			symlink($smp_data->{$patient}->{tumour}->{$tumour}, $link);
			}

		# create an array to hold final outputs and all patient job ids
		my (@final_outputs, @patient_jobs);
		my ($run_id, $java_check) = '';

		foreach my $sample (@sample_ids) {

			print $log "\n  SAMPLE: $sample\n";

			my $type;
			if ( (any { $_ =~ m/$sample/ } @normal_ids) ) {
				$type = 'normal';
				} else {
				$type = 'tumour';
				}

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			# run HaplotypeCaller
			my $hc_vcf = join('/', $sample_directory, $sample . '_HaplotypeCaller.vcf');
			if ('dna' eq $data_type) {
				$hc_vcf = join('/', $sample_directory, $sample . '_HaplotypeCaller.g.vcf');
				push @gvcfs, " -V:$sample $hc_vcf.gz";
				}

			# set up HC command
			my $call_variants_cmd = get_haplotype_command(
				bam		=> $smp_data->{$patient}->{$type}->{$sample},
				output		=> $hc_vcf,
				data_type	=> $data_type,
				intervals	=> $tool_data->{targets_bed},
				java_mem	=> $parameters->{haplotype_call}->{java_mem},
				tmp_dir		=> $tmp_directory
				);

			# if samples were aligned elsewhere, the @RG SM tag may differ from the sample ID
			# this causes problems if there are multiple files with the same SM tag (ie, when
			# running combined output from ConsensusCruncher)
			$call_variants_cmd .= "\n\n" . join("\n",
				"echo $sample | bcftools reheader -s - $hc_vcf -o $hc_vcf.reheadered",
				"mv $hc_vcf.reheadered $hc_vcf",
				"md5sum $hc_vcf > $hc_vcf.md5",
				"bgzip $hc_vcf",
				"tabix -p vcf $hc_vcf.gz"
				);

			# special case if multiple chromosomes and SLURM HPC driver
			if ( (scalar(@chroms) > 1) && ('slurm' eq $args{hpc_driver}) && ($is_wgs) ) {

				my $split_hc_vcf = join('/', $tmp_directory, $sample . '_HC_${CHROM}.g.vcf');

				$call_variants_cmd = 'CHROM=$(sed -n "$SLURM_ARRAY_TASK_ID"p ' . $chr_file . ')';
				$call_variants_cmd .= "\necho Running chromosome: " . '$CHROM';

				$call_variants_cmd .= "\n\n" . join("\n",
					"if [ -s $split_hc_vcf.md5 ]; then",
					'  echo Output file for $CHROM already exists',
					'  exit',
					'fi'
					);

				$call_variants_cmd .= "\n\n" . get_haplotype_command(
					bam		=> $smp_data->{$patient}->{$type}->{$sample},
					output		=> $split_hc_vcf,
					data_type	=> $data_type,
					intervals	=> '$CHROM',
					java_mem	=> $parameters->{haplotype_call}->{java_mem},
					tmp_dir		=> $tmp_directory
					);

				$call_variants_cmd .= "\n\n" . join("\n",
					"md5sum $split_hc_vcf > $split_hc_vcf.md5",
					"bgzip $split_hc_vcf",
					"tabix -p vcf $split_hc_vcf.gz"
					);
				}

			# check if this should be run
			if ('Y' eq missing_file($hc_vcf . '.md5')) {

				# record command (in log directory) and then run job
				print $log "  >> Submitting job for HaplotypeCaller...\n";

				if ( (scalar(@chroms) > 1) && ('slurm' eq $args{hpc_driver}) && ($is_wgs) ) {
					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_haplotype_caller_' . $sample,
						cmd	=> $call_variants_cmd,
						modules	=> [$gatk, $samtools],
						max_time	=> $parameters->{haplotype_call}->{time},
						mem		=> $parameters->{haplotype_call}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group, '--array=1-'. scalar(@chroms)]
						);
					} else {
					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_haplotype_caller_' . $sample,
						cmd	=> $call_variants_cmd,
						modules	=> [$gatk, $samtools],
						max_time	=> $parameters->{haplotype_call}->{time},
						mem		=> $parameters->{haplotype_call}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);
					}

				$run_id = submit_job(
					jobname		=> 'run_haplotype_caller_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;
				} else {
				print $log "  >> Skipping HaplotypeCaller because this has already been completed!\n";
				}

			# special case if multiple chromosomes and SLURM HPC driver
			if ( (scalar(@chroms) > 1) && ('slurm' eq $args{hpc_driver}) && ($is_wgs) ) {

				my $merge_command = "cd $tmp_directory";
				$merge_command .= "\n\n" . join(' ',
					'vcf-concat',
					$sample . '_HC*.g.vcf.gz',
					'>', $sample . '_HC_merged.g.vcf'
					);
				$merge_command .= "\n" . join(' ',
					'vcf-sort -c -t', $tmp_directory,
					$sample . '_HC_merged.g.vcf',
					'>', $hc_vcf					
					);

				# if samples were aligned elsewhere, the @RG SM tag may differ from the
				# sample ID this causes problems if there are multiple files with the same 
				# SM tag (ie, when running combined output from ConsensusCruncher)
				$merge_command .= "\n\n" . join("",
					"echo $sample | ",
					"bcftools reheader -s - $hc_vcf -o $hc_vcf.reheadered",
					"\nmv $hc_vcf.reheadered $hc_vcf",
					"\nmd5sum $hc_vcf > $hc_vcf.md5",
					"\nbgzip $hc_vcf",
					"\ntabix -p vcf $hc_vcf.gz"
					);

				# check if this should be run
				if ('Y' eq missing_file($hc_vcf . '.md5')) {

					# record command (in log directory) and then run job
					print $log "  >> Submitting job for MERGE (per-chrom)...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_merge_per_chromosome_vcfs_' . $sample,
						cmd	=> $merge_command,
						modules	=> [$vcftools, $samtools],
						dependencies	=> $run_id,
						max_time	=> '24:00:00',
						mem		=> '2G',
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_merge_per_chromosome_vcfs_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log "  >> Skipping MERGE because this has already been completed!\n";
					}
				}

			if ('dna' eq $data_type) { push @final_outputs, $hc_vcf; }

			# if this is RNA-Seq, run filter + vcf2maf
			if ('rna' eq $data_type) {

				$cleanup_cmd_rna .= "\nrm $hc_vcf.gz";

				# run filter variants
				my $filtered_vcf = join('/',
					$sample_directory,
					$sample . '_HaplotypeCaller_filtered.vcf'
					);

				$cleanup_cmd_rna .= "\nrm $filtered_vcf";

				my $filter_cmd = get_filter_command(
					input		=> "$hc_vcf.gz",
					output		=> $filtered_vcf,
					reference	=> $tool_data->{reference},
					java_mem	=> $parameters->{filter_raw}->{java_mem},
					tmp_dir		=> $tmp_directory
					);

				$filter_cmd .= "\nmd5sum $filtered_vcf > $filtered_vcf.md5";

				# check if this should be run
				if ('Y' eq missing_file($filtered_vcf . '.md5')) {

					# record command (in log directory) and then run job
					print $log "  >> Submitting job for Variant Filtration...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_filter_variants_' . $sample,
						cmd	=> $filter_cmd,
						modules	=> [$gatk],
						dependencies	=> $run_id,
						max_time	=> $parameters->{filter_raw}->{time},
						mem		=> $parameters->{filter_raw}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_filter_variants_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log "  >> Skipping Variant Filtration because this has already been completed!\n";
					}

				### Run variant annotation (VEP + vcf2maf)
				my $final_vcf = join('/',
					$sample_directory,
					$sample . '_HaplotypeCaller_filtered_annotated.vcf'
					);
				my $final_maf = join('/',
					$sample_directory,
					$sample . '_HaplotypeCaller_filtered_annotated.maf'
					);

				# currently assumes that no RNA-Seq was run on normals
				my $vcf2maf_cmd = get_vcf2maf_command(
					input		=> $filtered_vcf,
					tumour_id	=> $sample,
					reference	=> $tool_data->{reference},
					ref_type	=> $tool_data->{ref_type},
					output		=> $final_maf,
					tmp_dir		=> $sample_directory,
					parameters	=> $parameters->{annotate}
					);

				# check if this should be run
				if ('Y' eq missing_file($final_maf . '.md5')) {

					# IF THIS FINAL STEP IS SUCCESSFULLY RUN,
					$vcf2maf_cmd .= "\n\n" . join("\n",
						"if [ -s $final_maf ]; then",
						"  md5sum $final_maf > $final_maf.md5",
						"  mv $sample_directory/$sample" . "_HaplotypeCaller_filtered.vep.vcf $final_vcf",
						"  md5sum $final_vcf > $final_vcf.md5",
						"  bgzip $final_vcf",
						"  tabix -p vcf $final_vcf.gz",
						"else",
						'  echo "FINAL OUTPUT MAF is missing; not running md5sum/bgzip/tabix..."',
						"fi"
						);

					# record command (in log directory) and then run job
					print $log "  >> Submitting job for vcf2maf...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_vcf2maf_and_VEP_' . $sample,
						cmd	=> $vcf2maf_cmd,
						modules => [$perl, $samtools, $vcf2maf],
						dependencies	=> $run_id,
						cpus_per_task	=> $parameters->{annotate}->{n_cpus},
						max_time	=> $parameters->{annotate}->{time},
						mem		=> $parameters->{annotate}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_vcf2maf_and_VEP_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log "  >> Skipping vcf2maf because this has already been completed!\n";
					}

				push @final_outputs, $final_maf;
				}
			}

		# should intermediate files be removed
		# run per patient
		if ( ($args{del_intermediates}) && ('rna' eq $data_type) ) {

			if (scalar(@patient_jobs) == 0) {
				`rm -rf $tmp_directory`;
				} else {

				print $log ">> Submitting job to clean up temporary/intermediate files...\n";

				# make sure final output exists before removing intermediate files!
				my @files_to_check;
				foreach my $tmp ( @final_outputs ) {
					push @files_to_check, $tmp . '.md5';
					}

				$cleanup_cmd_rna = join("\n",
					"if [ -s " . join(" ] && [ -s ", @files_to_check) . " ]; then",
					$cleanup_cmd_rna,
					"else",
					'echo "One or more FINAL OUTPUT FILES is missing; not removing intermediates"',
					"fi"
					);

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_cleanup_' . $patient,
					cmd	=> $cleanup_cmd_rna,
					dependencies	=> join(':', @patient_jobs),
					max_time	=> '00:05:00',
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
			}

		print $log "\nFINAL OUTPUT:\n" . join("\n  ", @final_outputs) . "\n";
		print $log "---\n";
		}

	# if this is DNA-Seq (WXS or WGS), combine above gvcfs (per-batch)
	if ('dna' eq $data_type) {

		my @batches = ([]);
		my $batch_idx = 0;
		my $file_count = 0;

		foreach my $gvcf ( @gvcfs ) {

			if ($file_count > 15) {
				$batch_idx++;
				@batches[$batch_idx] = [];
				$file_count = 0;
				}
			push @{$batches[$batch_idx]}, $gvcf;

			$file_count++;

			}

		$batch_idx = 0;

		my (@combined_gvcfs, @batch_jobs);

		foreach my $batch ( @batches ) {

			$batch_idx++;

			# run CombineGVCFs
			my $combined_gvcf = join('/', $output_directory, 'haplotype_caller_' . $batch_idx . '.g.vcf');

			my $combine_cmd = get_combine_gvcf_command(
				input		=> join(' ', @{$batch}),
				output		=> $combined_gvcf,
				java_mem	=> $parameters->{combine_gvcfs}->{java_mem},
				tmp_dir		=> $output_directory
				);

			# this is a java-based command, so run a final check
			my $extra_cmds = "\n  " . join("\n  ",
				"md5sum $combined_gvcf > $combined_gvcf.md5",
				"bgzip $combined_gvcf",
				"tabix -p vcf $combined_gvcf.gz\n"
				);

			$java_check = "\n" . check_java_output(
				extra_cmd => $extra_cmds
				);

			$combine_cmd .= "\n$java_check";

			# check if this should be run
			if ('Y' eq missing_file($combined_gvcf . '.md5')) {

				# record command (in log directory) and then run job
				print $log ">> Submitting job for CombineGVCFs...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_combine_gvcfs_batch_' . $batch_idx,
					cmd	=> $combine_cmd,
					modules	=> [$gatk, 'tabix'],
					dependencies	=> join(':', @all_jobs),
					max_time	=> $parameters->{combine_gvcfs}->{time},
					mem		=> $parameters->{combine_gvcfs}->{mem},
					hpc_driver	=> $args{hpc_driver},
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'run_combine_gvcfs_batch_' . $batch_idx,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @batch_jobs, $run_id;
				}

			push @combined_gvcfs, "$combined_gvcf.gz";
			}

		push @all_jobs, @batch_jobs;

		# should intermediate files be removed
		if ($args{del_intermediates}) {

			print $log ">> Submitting job to clean up temporary/intermediate files...\n";

			# make sure final output exists before removing intermediate files!
			$cleanup_cmd_dna = join("\n",
				"if [ -s " . join(" ] && [ -s ", @combined_gvcfs) . " ]; then",
				"  $cleanup_cmd_dna",
				"else",
				'  echo "One or more FINAL OUTPUT FILES is missing; not removing intermediates"',
				"fi"
				);

			$run_script = write_script(
				log_dir	=> $log_directory,
				name	=> 'run_final_cleanup',
				cmd	=> $cleanup_cmd_dna,
				dependencies	=> join(':', @all_jobs),
				mem		=> '256M',
				hpc_driver	=> $args{hpc_driver},
				kill_on_error	=> 0,
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'run_final_cleanup',
				shell_command	=> $run_script,
				hpc_driver	=> $args{hpc_driver},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);
			}
		}

	if ( ('rna' eq $data_type) ) {

		# collect and combine results
		my $collect_results = join(' ',
			"Rscript $cwd/collect_snv_output.R",
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
			extra_args	=> [$hpc_group],
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
my ($remove_junk, $dry_run, $rna, $help, $no_wait);

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	'd|data=s'	=> \$data_config,
	't|tool=s'	=> \$tool_config,
	'o|out_dir=s'	=> \$output_directory,
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
		"\t--cluster|-c\t<string> cluster scheduler (default: slurm)",
		"\t--remove\t<boolean> should intermediates be removed? (default: false)",
		"\t--dry-run\t<boolean> should jobs be submitted? (default: false)",
		"\t--no-wait\t<boolean> should we exit after job submission (true) or wait until all jobs have completed (false)? (default: false)",
		"\t--rna\t<boolean> is the input RNA (STAR-aligned BAMs)? (default: false)"
		);

	print "$help_msg\n";
	exit;
	}

# do some quick error checks to confirm valid arguments	
if (!defined($tool_config)) { die("No tool config file defined; please provide -t | --tool (ie, tool_config.yaml)"); }
if (!defined($data_config)) { die("No data config file defined; please provide -d | --data (ie, sample_config.yaml)"); }
if (!defined($output_directory)) { die("No output directory defined; please provide -o | --out_dir"); }

my $data_type = 'dna';
if ($rna) { $data_type = 'rna'; }

main(
	tool_config		=> $tool_config,
	data_config		=> $data_config,
	output_directory	=> $output_directory,
	hpc_driver		=> $hpc_driver,
	del_intermediates	=> $remove_junk,
	dry_run			=> $dry_run,
	no_wait			=> $no_wait,
	data_type		=> $data_type
	);
