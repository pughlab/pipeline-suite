#!/usr/bin/env perl
### pughlab_pipeline_summarize_output.pl ###########################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use Getopt::Std;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(any all first);
use List::MoreUtils qw(first_index);
use YAML qw(LoadFile);

my $cwd = dirname(__FILE__);
require "$cwd/utilities.pl";

####################################################################################################
# version       author		comment
# 1.0		sprokopec       tool to automatically generate reports

### MAIN ###########################################################################################
sub main {
	my %args = (
		config		=> undef,
		cluster		=> undef,
		report		=> undef,
		dry_run		=> undef,
		run_date	=> undef,
		no_wait		=> undef,
		@_
		);

	### PREAMBLE ######################################################################################
	unless($args{dry_run}) {
		print "Initiating SUMMARIZE pipeline...\n";
		}

	my $tool_data = LoadFile($args{config});
	my $run_date = $args{run_date};
	my $current_date = strftime "%F", localtime;

	# check for and/or create output directories
	my $output_directory = $tool_data->{output_dir};
	$output_directory =~ s/\/$//;

	my $summary_directory = join('/', $output_directory, 'SUMMARY_RESULTS');
	unless(-e $summary_directory) { make_path($summary_directory); }

	my $log_directory = join('/', $summary_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	# Initiate logging
	my $log_file = join('/', $log_directory, 'run_SUMMARIZE_pipeline.log');

	# create a file to hold job metrics
	my (@files, $run_count, $outfile, $touch_exit_status);
	unless ($args{dry_run}) {
		# initiate a file to hold job metrics 
		# (ensures that an existing file isn't overwritten by concurrent jobs)
		opendir(LOGFILES, $log_directory) or die "Cannot open $log_directory";
		@files = grep { /slurm_job_metrics/ } readdir(LOGFILES);
		$run_count = scalar(@files) + 1;
		closedir(LOGFILES);

		$outfile = $log_directory . '/slurm_job_metrics_' . $run_count . '.out';
		$touch_exit_status = system("touch $outfile");
		if (0 != $touch_exit_status) { Carp::croak("Cannot touch file $outfile"); }

		$log_file = join('/', $log_directory, 'run_SUMMARIZE_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running Summarize Results Pipeline...\n";
	print $log "\n  Tool config used: $args{config}";
	print $log "\n  Output directory: $output_directory";
	print $log "\n---\n";

	# get tool versions
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};
	my $vcftools	= 'vcftools/' . $tool_data->{vcftools_version};
	my $r_version	= 'R/' . $tool_data->{r_version};
	my $vcf2maf	= undef;
	if (defined($tool_data->{vcf2maf_version})) {
		$vcf2maf = 'vcf2maf/' . $tool_data->{vcf2maf_version};
		$tool_data->{annotate}->{vcf2maf_path} = undef;
		}

	# get optional HPC group
	my $hpc_group = defined($tool_data->{hpc_group}) ? "-A $tool_data->{hpc_group}" : undef;

	### RUN ####################################################################################
	# set up directory structure
	my $data_directory = join('/', $summary_directory, 'data');
	unless(-e $data_directory) { make_path($data_directory); }

	my $report_directory = join('/', $summary_directory, 'Report');
	unless(-e $report_directory) { make_path($report_directory); }

	my $plot_directory = join('/', $report_directory, 'plots');
	unless(-e $plot_directory) { make_path($plot_directory); }

	# check which tools have been requested
	my %tool_set = (
		'bwa'	=> defined($tool_data->{bwa}->{run}) ? $tool_data->{bwa}->{run} : 'N',
		'gatk'	=> defined($tool_data->{gatk}->{run}) ? $tool_data->{gatk}->{run} : 'N',
		'bamqc'	=> defined($tool_data->{bamqc}->{run}) ? $tool_data->{bamqc}->{run} : 'N',
		'haplotype_caller' => defined($tool_data->{haplotype_caller}->{run}) ? $tool_data->{haplotype_caller}->{run} : 'N',
		'mutect'	=> defined($tool_data->{mutect}->{run}) ? $tool_data->{mutect}->{run} : 'N',
		'mutect2'	=> defined($tool_data->{mutect2}->{run}) ? $tool_data->{mutect2}->{run} : 'N',
		'somaticsniper'	=> defined($tool_data->{somaticsniper}->{run}) ? $tool_data->{somaticsniper}->{run} : 'N',
		'strelka'	=> defined($tool_data->{strelka}->{run}) ? $tool_data->{strelka}->{run} : 'N',
		'varscan'	=> defined($tool_data->{varscan}->{run}) ? $tool_data->{varscan}->{run} : 'N',
		'vardict'	=> defined($tool_data->{vardict}->{run}) ? $tool_data->{vardict}->{run} : 'N',
		'pindel'	=> defined($tool_data->{pindel}->{run}) ? $tool_data->{pindel}->{run} : 'N',
		'gatk_cnv'	=> defined($tool_data->{gatk_cnv}->{run}) ? $tool_data->{gatk_cnv}->{run} : 'N',
		'cnvkit'	=> defined($tool_data->{cnvkit}->{run}) ? $tool_data->{cnvkit}->{run} : 'N',
		'gatk_gcnv'	=> defined($tool_data->{gatk_gcnv}->{run}) ? $tool_data->{gatk_gcnv}->{run} : 'N',
		'erds_gcnv'	=> defined($tool_data->{erds_gcnv}->{run}) ? $tool_data->{erds_gcnv}->{run} : 'N',
		'ichor_cna'	=> defined($tool_data->{ichor_cna}->{run}) ? $tool_data->{ichor_cna}->{run} : 'N',
		'mops'		=> defined($tool_data->{panelcn_mops}->{run}) ? $tool_data->{panelcn_mops}->{run} : 'N',
		'ascat'		=> defined($tool_data->{ascat}->{run}) ? $tool_data->{ascat}->{run} : 'N',
		'novobreak'	=> defined($tool_data->{novobreak}->{run}) ? $tool_data->{novobreak}->{run} : 'N',
		'delly'		=> defined($tool_data->{delly}->{run}) ? $tool_data->{delly}->{run} : 'N',
		'svict'		=> defined($tool_data->{svict}->{run}) ? $tool_data->{svict}->{run} : 'N',
		'mavis'		=> defined($tool_data->{mavis}->{run}) ? $tool_data->{mavis}->{run} : 'N',
		'msi'		=> defined($tool_data->{msi_sensor}->{run}) ? $tool_data->{msi_sensor}->{run} : 'N',
		'star'		=> defined($tool_data->{star}->{run}) ? $tool_data->{star}->{run} : 'N',
		'rsem'		=> defined($tool_data->{rsem}->{run}) ? $tool_data->{rsem}->{run} : 'N',
		'arriba'	=> defined($tool_data->{arriba}->{run}) ? $tool_data->{arriba}->{run} : 'N',
		'star_fusion'	=> defined($tool_data->{star_fusion}->{run}) ? $tool_data->{star_fusion}->{run} : 'N',
		'fusioncatcher'	=> defined($tool_data->{fusioncatcher}->{run}) ? $tool_data->{fusioncatcher}->{run} : 'N',
		'cosmic_sbs'	=> defined($tool_data->{summarize_steps}->{run_cosmic_sbs}) ? $tool_data->{summarize_steps}->{run_cosmic_sbs} : 'N',
		'chord'		=> defined($tool_data->{summarize_steps}->{run_chord}) ? $tool_data->{summarize_steps}->{run_chord} : 'N',
		'hrdetect'	=> defined($tool_data->{summarize_steps}->{run_hrdetect}) ? $tool_data->{summarize_steps}->{run_hrdetect} : 'N',
		'mutsig'	=> defined($tool_data->{summarize_steps}->{run_mutsig}) ? $tool_data->{summarize_steps}->{run_mutsig} : 'N'
		);

	# initiate objects
	my @job_ids;
	my ($qc_dir, $germ_dir, $cpsr_dir);
	my ($correlations, $qc_data, $cb_data, $seqqc_data, $callability_data, $contest_data, $cpsr_calls) = undef;
	my ($run_script, $run_id);

	### find the required input files
	if ('rna' eq $tool_data->{seq_type}) {

		$qc_dir = join('/', $output_directory, 'STAR');

		# coverage summary metrics
		opendir(RNASEQC, $qc_dir) or die "Cannot open '$qc_dir' !";
		my @rnaseqc_files = grep { /.tsv/ } readdir(RNASEQC);
		my @coverage_files = grep { /rnaseqc_output.tsv/ } @rnaseqc_files;
		@coverage_files = sort @coverage_files;
		my @correlation_files = grep { /Pearson_correlations.tsv/ } @rnaseqc_files;
		@correlation_files = sort @correlation_files;
		closedir(RNASEQC);

		$qc_data = join('/', $qc_dir, $coverage_files[-1]);
		$correlations = join('/', $qc_dir, $correlation_files[-1]);

		if ( -l join('/', $data_directory, 'rnaseqc_output.tsv')) {
			unlink join('/', $data_directory, 'rnaseqc_output.tsv');
			}
		if ( -l join('/', $data_directory, 'genes.rpkm.gct_correlations.tsv')) {
			unlink join('/', $data_directory, 'genes.rpkm.gct_correlations.tsv');
			}

		symlink($qc_data, join('/', $data_directory, 'rnaseqc_output.tsv'));
		symlink($correlations, join('/', $data_directory, 'genes.rpkm_correlations.tsv'));

		# create some QC plots
		my $qc_command = "cd $output_directory\n";
		$qc_command .= "Rscript $cwd/report/plot_qc_metrics.R";
		$qc_command .= " " . join(' ',
			'-g', $correlations,
			'-q', $qc_data,
			'-o', $plot_directory,
			'-p', $tool_data->{project_name},
			'-t', $tool_data->{seq_type}
			);

		# run command
		print $log "Submitting job to create QC plots...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'create_qc_plots',
			cmd		=> $qc_command,
			modules		=> [$r_version],
			mem		=> '2G',
			hpc_driver	=> $args{cluster},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'create_qc_plots',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;

		# rna expression values
		if ('Y' eq $tool_set{'rsem'}) {
			my $rsem_dir = join('/', $output_directory, 'RSEM');
			opendir(RSEM, $rsem_dir) or die "Cannot open '$rsem_dir' !";
			my @rsem_calls = grep { /expression_TPM_for_cbioportal.tsv/ } readdir(RSEM);
			@rsem_calls = sort @rsem_calls;
			closedir(RSEM);

			my $rsem_data = join('/', $rsem_dir, $rsem_calls[-1]);

			if ( -l join('/', $data_directory, 'rsem_expression.tsv')) {
				unlink join('/', $data_directory, 'rsem_expression.tsv');
				}

			symlink($rsem_data, join('/', $data_directory, 'rsem_expression.tsv'));

			# plot RNA expression profile
			my $rna_plot_command = join(' ',
				"Rscript $cwd/report/plot_expression_summary.R",
				'-p', $tool_data->{project_name},
				'-o', $plot_directory,
				'-r', $rsem_data
				);

			# run command
			print $log "Submitting job to plot RNA expression values...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_rna_expression_summary',
				cmd		=> $rna_plot_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'plot_rna_expression_summary',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;
			}

		# if variant calling was performed
		if ('Y' eq $tool_set{'haplotype_caller'}) {
			my $snv_dir = join('/', $output_directory, 'HaplotypeCaller');
			opendir(VARIANTS, $snv_dir) or die "Cannot open '$snv_dir' !";
			my @snv_calls = grep { /_mutations_for_cbioportal.tsv/ } readdir(VARIANTS);
			@snv_calls = sort @snv_calls;
			closedir(VARIANTS);

			my $snv_data = join('/', $snv_dir, $snv_calls[-1]);

			if ( -l join('/', $data_directory, 'haplotype_caller_snv_calls.tsv')) {
				unlink join('/', $data_directory, 'haplotype_caller_snv_calls.tsv');
				}

			symlink($snv_data, join('/', $data_directory, 'haplotype_caller_snv_calls.tsv'));

			# plot SNV summary
			my $snv_plot_command = join(' ',
				"Rscript $cwd/report/plot_rna_snv_summary.R",
				'-p', $tool_data->{project_name},
				'-o', $plot_directory,
				'-m', $snv_data
				);

			# run command
			print $log "Submitting job to plot SNV Summary...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_rna_snv_summary',
				cmd		=> $snv_plot_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'plot_rna_snv_summary',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;

			}

		# rna fusion values
		my $rna_fusion_command = join(' ',
			"Rscript $cwd/report/write_rna_fusion_summary.R",
			'-p', $tool_data->{project_name},
			'-o', $plot_directory
			);

		my $tool_count = 0;
		if ('Y' eq $tool_set{'star_fusion'}) {
			my $starfus_dir = join('/', $output_directory, 'STAR-Fusion');
			opendir(FUSION, $starfus_dir) or die "Cannot open '$starfus_dir' !";
			my @starfus_calls = grep { /_for_cbioportal.tsv/ } readdir(FUSION);
			@starfus_calls = sort @starfus_calls;
			closedir(FUSION);

			my $starfus_data = join('/', $starfus_dir, $starfus_calls[-1]);

			if ( -l join('/', $data_directory, 'starfusion_calls.tsv')) {
				unlink join('/', $data_directory, 'starfusion_calls.tsv');
				}

			symlink($starfus_data, join('/', $data_directory, 'starfusion_calls.tsv'));

			$rna_fusion_command .= " -s $starfus_data";
			$tool_count++;
			}

		if ('Y' eq $tool_set{'fusioncatcher'}) {
			my $fuscatch_dir = join('/', $output_directory, 'FusionCatcher');
			opendir(FUSION, $fuscatch_dir) or die "Cannot open '$fuscatch_dir' !";
			my @fuscatch_files = readdir(FUSION);
			closedir(FUSION);

			my @fuscatch_calls = grep { /_for_cbioportal.tsv/ } @fuscatch_files;
			@fuscatch_calls = sort @fuscatch_calls;

			my @viral_counts = grep { /viral_counts.tsv/ } @fuscatch_files;
			@viral_counts = sort @viral_counts;

			my $fuscatch_data = join('/', $fuscatch_dir, $fuscatch_calls[-1]);
			my $viral_data = join('/', $fuscatch_dir, $viral_counts[-1]);

			if ( -l join('/', $data_directory, 'fusioncatcher_calls.tsv')) {
				unlink join('/', $data_directory, 'fusioncatcher_calls.tsv');
				}

			if ( -l join('/', $data_directory, 'fusioncatcher_viral_counts.tsv')) {
				unlink join('/', $data_directory, 'fusioncatcher_viral_counts.tsv');
				}

			symlink($fuscatch_data, join('/', $data_directory, 'fusioncatcher_calls.tsv'));
			symlink($viral_data, join('/', $data_directory, 'fusioncatcher_viral_counts.tsv'));

			# for fusion summary
			$rna_fusion_command .= " -f $fuscatch_data";
			$tool_count++;

			# for viral summary
			my $rna_virus_command = join(' ',
				"Rscript $cwd/report/plot_viral_counts.R",
				'-p', $tool_data->{project_name},
				'-o', $plot_directory,
				'-v', $viral_data
				);			

			# run command
			print $log "Submitting job to summarize viral content...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'summarize_viral_counts',
				cmd		=> $rna_virus_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'summary_viral_counts',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;

			}

		# get SVs calls from MAVIS
		if ('Y' eq $tool_set{'mavis'}) {
			my $mavis_dir = join('/', $output_directory, 'Mavis');

			opendir(MAVIS, $mavis_dir) or die "Cannot open '$mavis_dir' !";
			my @mavis_files = grep { /mavis_output.tsv/ } readdir(MAVIS);
			@mavis_files = sort @mavis_files;
			closedir(MAVIS);

			my $sv_data = join('/', $mavis_dir, $mavis_files[-1]);

			if ( -l join('/', $data_directory, 'mavis_sv_data.tsv')) {
				unlink join('/', $data_directory, 'mavis_sv_data.tsv');
				}

			symlink($sv_data, join('/', $data_directory, 'mavis_sv_data.tsv'));

			# for fusion summary
			$rna_fusion_command .= " -m $sv_data";
			}

		if ($tool_count > 0) {

			# run command
			print $log "Submitting job to summarize RNA fusions...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'summarize_rna_fusions',
				cmd		=> $rna_fusion_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'summarize_rna_fusions',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;
			}

		} else {

		$qc_dir		= join('/', $output_directory, 'BAMQC');
		$germ_dir	= join('/', $output_directory, 'HaplotypeCaller/cohort/germline_variants');
		$cpsr_dir	= join('/', $output_directory, 'HaplotypeCaller/CPSR');

		# contamination estimates (T/N, meta + readgroup)
		if (-e "$qc_dir/ContEst") {
			opendir(CONTEST, $qc_dir . "/ContEst") or die "Cannot open '$qc_dir/ContEst' !";
			my @contest_files = grep { /ContEst_output.tsv/ } readdir(CONTEST);
			@contest_files = sort @contest_files;
			closedir(CONTEST);

			$contest_data = join('/', $qc_dir, 'ContEst', $contest_files[-1]);

			if ( -l join('/', $data_directory, 'ContEst_output.tsv')) {
				unlink join('/', $data_directory, 'ContEst_output.tsv');
				}

			symlink($contest_data, join('/', $data_directory, 'ContEst_output.tsv'));
			}

		# contamination estimates (all samples)
		if (-e "$qc_dir/SequenceMetrics") {
			opendir(QC, $qc_dir . "/SequenceMetrics");
			my @seqqc_files = grep { /tsv/ } readdir(QC);
			@seqqc_files = sort @seqqc_files;

			my @qc_files = grep { /ContaminationEstimates.tsv/ } @seqqc_files;
			my @cov_files = grep { /WGSMetrics.tsv/ } @seqqc_files;

			closedir(QC);

			$seqqc_data = join('/', $qc_dir, 'SequenceMetrics', $qc_files[-1]);

			if ( -l join('/', $data_directory, 'ContaminationEstimates.tsv')) {
				unlink join('/', $data_directory, 'ContaminationEstimates.tsv');
				}

			symlink($seqqc_data, join('/', $data_directory, 'ContaminationEstimates.tsv'));

			if ('wgs' eq $tool_data->{seq_type}) {
				$callability_data = join('/', $qc_dir, 'SequenceMetrics', $cov_files[-1]);
				if ( -l join('/', $data_directory, 'wgs_callability.tsv')) {
					unlink join('/', $data_directory, 'wgs_callability.tsv');
					}

				symlink($callability_data, join('/', $data_directory, 'wgs_callability.tsv'));
				}
			}

		# coverage summary metrics
		if (-e "$qc_dir/Coverage") {
			opendir(COVERAGE, $qc_dir . "/Coverage") or die "Cannot open '$qc_dir/Coverage' !";
			my @all_coverage_files = grep { /tsv/ } readdir(COVERAGE);

			my @coverage_files = grep { /Coverage_summary.tsv/ } @all_coverage_files;
			@coverage_files = sort @coverage_files;

			my @cb_files = grep { /total_bases_covered.tsv/ } @all_coverage_files;
			@cb_files = sort @cb_files;

			closedir(COVERAGE);

			$qc_data = join('/', $qc_dir, 'Coverage', $coverage_files[-1]);
			$cb_data = join('/', $qc_dir, 'Coverage', $cb_files[-1]);

			if ( -l join('/', $data_directory, 'Coverage_summary.tsv')) {
				unlink join('/', $data_directory, 'Coverage_summary.tsv');
				}

			if ( -l join('/', $data_directory, 'total_bases_covered.tsv')) {
				unlink join('/', $data_directory, 'total_bases_covered.tsv');
				}		

			symlink($qc_data, join('/', $data_directory, 'Coverage_summary.tsv'));
			symlink($cb_data, join('/', $data_directory, 'total_bases_covered.tsv'));
			}

		# germline variant correlations
		if (-e $germ_dir) {
			opendir(GERMLINE, $germ_dir) or die "Cannot open '$germ_dir' !";
			my @correlation_files = grep { /germline_correlation.tsv/ } readdir(GERMLINE);
			@correlation_files = sort @correlation_files;
			closedir(GERMLINE);

			$correlations = join('/', $germ_dir, $correlation_files[-1]);

			if ( -l join('/', $data_directory, 'germline_correlation.tsv')) {
				unlink join('/', $data_directory, 'germline_correlation.tsv');
				}

			symlink($correlations, join('/', $data_directory, 'germline_correlation.tsv'));
			}

		# create some QC plots
		my $qc_command = "cd $output_directory\n";
		$qc_command .= "Rscript $cwd/report/plot_qc_metrics.R";
		$qc_command .= " " . join(' ',
			'-g', $correlations,
			'-q', $qc_data,
			'-o', $plot_directory,
			'-p', $tool_data->{project_name},
			'-t', $tool_data->{seq_type}
			);

		if (defined($contest_data)) { $qc_command .= " -c $contest_data"; }
		if (defined($seqqc_data)) { $qc_command .= " -m $seqqc_data"; }

		my $qc_run_id;
		if ('Y' eq $tool_set{'bamqc'}) {

			# run command
			print $log "Submitting job to create QC plots...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'create_qc_plots',
				cmd		=> $qc_command,
				modules		=> [$r_version],
				mem		=> '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$qc_run_id = submit_job(
				jobname		=> 'create_qc_plots',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $qc_run_id;
			}

		# make mutation-type directories
		my $snv_directory = join('/', $summary_directory, 'MutationSummary');
		unless(-e $snv_directory) { make_path($snv_directory); }

		my $cnv_directory = join('/', $summary_directory, 'CNA_Summary');
		unless(-e $cnv_directory) { make_path($cnv_directory); }

		my $sv_directory = join('/', $summary_directory, 'SV_Summary');
		unless(-e $sv_directory) { make_path($sv_directory); }

		my $sig_directory = join('/', $summary_directory, 'MutationSignatures');
		if ( ('Y' eq $tool_set{'cosmic_sbs'}) || ('Y' eq $tool_set{'chord'}) ||
			('Y' eq $tool_set{'hrdetect'}) ) {
			unless(-e $sig_directory) { make_path($sig_directory); }
			}

		# significant germline variants
		if (-e $cpsr_dir) {
			opendir(CPSR, $cpsr_dir) or die "Cannot open '$cpsr_dir' !";
			my @cpsr_files = grep { /mutations_for_cbioportal.tsv/ } readdir(CPSR);
			@cpsr_files = sort @cpsr_files;
			closedir(CPSR);

			$cpsr_calls = join('/', $cpsr_dir, $cpsr_files[-1]);

			if ( -l join('/', $data_directory, 'cpsr_germline_variants.tsv')) {
				unlink join('/', $data_directory, 'cpsr_germline_variants.tsv');
				}

			symlink($cpsr_calls, join('/', $data_directory, 'cpsr_germline_variants.tsv'));

			# summarize/plot CPSR output
			my $plot_command = join(' ',
				"Rscript $cwd/report/plot_germline_snv_summary.R",
				'-p', $tool_data->{project_name},
				'-m', join('/', $data_directory, 'cpsr_germline_variants.tsv'),
				'-o', $snv_directory,
				'--report', $plot_directory
				);

			# run command
			print $log "Submitting job to create germline SNV plots...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_germline_snv_summary',
				cmd		=> $plot_command,
				modules		=> [$r_version],
				mem		=> '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'plot_germline_snv_summary',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;
			}

		# somatic variants
		my ($mutect_data, $mutect2_data, $strelka_data, $varscan_data, $sniper_data, $vardict_data);
		my ($pindel_data, $sequenza_data, $segments, $ploidy_data, $gatk_cnv, $gatk_pga, $msi_data);
		my ($gatk_gcnv, $cnvkit, $erds_gcnv, $mops_cnv, $ascat, $ichor_cnv);

		# find ENSEMBLE mutations
		my $ensemble_command .= join(' ',
			"Rscript $cwd/report/format_ensemble_mutations.R",
			'-p', $tool_data->{project_name},
			'-o', join('/', $data_directory, 'ENSEMBLE'),
			'-t', $tool_data->{seq_type}
			);

		# get mutation calls from MuTect
		if ('Y' eq $tool_set{'mutect'}) {
			my $mutect_dir = join('/', $output_directory, 'MuTect');
			opendir(MUTECT, $mutect_dir) or die "Cannot open '$mutect_dir' !";
			my @mutect_calls = grep { /mutations_for_cbioportal.tsv/ } readdir(MUTECT);
			@mutect_calls = sort @mutect_calls;
			closedir(MUTECT);

			$mutect_data = join('/', $mutect_dir, $mutect_calls[-1]);

			if ( -l join('/', $data_directory, 'mutect_somatic_variants.tsv')) {
				unlink join('/', $data_directory, 'mutect_somatic_variants.tsv');
				}

			symlink($mutect_data, join('/', $data_directory, 'mutect_somatic_variants.tsv'));

			$ensemble_command .= ' --mutect ' . join('/', $data_directory, 'mutect_somatic_variants.tsv');
			}

		# get mutation calls from MuTect2
		if ('Y' eq $tool_set{'mutect2'}) {
			my $mutect2_dir = join('/', $output_directory, 'MuTect2');
			opendir(MUTECT2, $mutect2_dir) or die "Cannot open '$mutect2_dir' !";
			my @mutect2_calls = grep { /mutations_for_cbioportal.tsv/ } readdir(MUTECT2);
			@mutect2_calls = sort @mutect2_calls;
			closedir(MUTECT2);

			$mutect2_data = join('/', $mutect2_dir, $mutect2_calls[-1]);

			if ( -l join('/', $data_directory, 'mutect2_somatic_variants.tsv')) {
				unlink join('/', $data_directory, 'mutect2_somatic_variants.tsv');
				}

			symlink($mutect2_data, join('/', $data_directory, 'mutect2_somatic_variants.tsv'));

			$ensemble_command .= " --mutect2 " . join('/', $data_directory, 'mutect2_somatic_variants.tsv');
			}

		# get mutation calls from Pindel
		if ('Y' eq $tool_set{'pindel'}) {
			my $pindel_dir = join('/', $output_directory, 'Pindel');
			opendir(PINDEL, $pindel_dir) or die "Cannot open '$pindel_dir' !";
			my @pindel_calls = grep { /mutations_for_cbioportal.tsv/ } readdir(PINDEL);
			@pindel_calls = sort @pindel_calls;
			closedir(PINDEL);

			$pindel_data = join('/', $pindel_dir, $pindel_calls[-1]);

			if ( -l join('/', $data_directory, 'pindel_somatic_variants.tsv')) {
				unlink join('/', $data_directory, 'pindel_somatic_variants.tsv');
				}

			symlink($pindel_data, join('/', $data_directory, 'pindel_somatic_variants.tsv'));

			$ensemble_command .= " --pindel " . join('/', $data_directory, 'pindel_somatic_variants.tsv');
			}

		# get mutation calls from Strelka
		if ('Y' eq $tool_set{'strelka'}) {
			my $strelka_dir = join('/', $output_directory, 'Strelka');
			opendir(STRELKA, $strelka_dir) or die "Cannot open '$strelka_dir' !";
			my @strelka_calls = grep { /mutations_for_cbioportal.tsv/ } readdir(STRELKA);
			@strelka_calls = sort @strelka_calls;
			closedir(STRELKA);

			$strelka_data = join('/', $strelka_dir, $strelka_calls[-1]);

			if ( -l join('/', $data_directory, 'strelka_somatic_variants.tsv')) {
				unlink join('/', $data_directory, 'strelka_somatic_variants.tsv');
				}

			symlink($strelka_data, join('/', $data_directory, 'strelka_somatic_variants.tsv'));

			$ensemble_command .= " --strelka " . join('/', $data_directory, 'strelka_somatic_variants.tsv');
			}

		# get mutation calls from SomaticSniper
		if ('Y' eq $tool_set{'somaticsniper'}) {
			my $sniper_dir = join('/', $output_directory, 'SomaticSniper');
			opendir(SNIPER, $sniper_dir) or die "Cannot open '$sniper_dir' !";
			my @sniper_calls = grep { /mutations_for_cbioportal.tsv/ } readdir(SNIPER);
			@sniper_calls = sort @sniper_calls;
			closedir(SNIPER);

			$sniper_data = join('/', $sniper_dir, $sniper_calls[-1]);

			if ( -l join('/', $data_directory, 'somaticsniper_somatic_variants.tsv')) {
				unlink join('/', $data_directory, 'somaticsniper_somatic_variants.tsv');
				}

			symlink($sniper_data, join('/', $data_directory, 'somaticsniper_somatic_variants.tsv'));

			$ensemble_command .= " --somaticsniper " . join('/', $data_directory, 'somaticsniper_somatic_variants.tsv');
			}

		# get mutation calls from VarDict
		if ('Y' eq $tool_set{'vardict'}) {
			my $vardict_dir = join('/', $output_directory, 'VarDict');
			opendir(VARDICT, $vardict_dir) or die "Cannot open '$vardict_dir' !";
			my @vardict_calls = grep { /mutations_for_cbioportal.tsv/ } readdir(VARDICT);
			@vardict_calls = sort @vardict_calls;
			closedir(VARDICT);

			$vardict_data = join('/', $vardict_dir, $vardict_calls[-1]);

			if ( -l join('/', $data_directory, 'vardict_somatic_variants.tsv')) {
				unlink join('/', $data_directory, 'vardict_somatic_variants.tsv');
				}

			symlink($vardict_data, join('/', $data_directory, 'vardict_somatic_variants.tsv'));

			$ensemble_command .= " --vardict " . join('/', $data_directory, 'vardict_somatic_variants.tsv');
			}

		# get mutation calls from VarScan
		my (@cna_calls, @seg_calls, @ploidy);

		if ('Y' eq $tool_set{'varscan'}) {
			my $varscan_dir = join('/', $output_directory, 'VarScan');

			opendir(VARSCAN, $varscan_dir) or die "Cannot open '$varscan_dir' !";
			my @varscan_files = readdir(VARSCAN);

			my @varscan_calls = grep { /mutations_for_cbioportal.tsv/ } @varscan_files;
			@varscan_calls = sort @varscan_calls;

			@cna_calls = grep { /Sequenza_ratio_gene_matrix.tsv/ } @varscan_files;
			@cna_calls = sort @cna_calls;

			@seg_calls = grep { /Sequenza_segments.tsv/ } @varscan_files;
			@seg_calls = sort @seg_calls;

			@ploidy = grep { /Sequenza_ploidy_purity.tsv/ } @varscan_files;
			@ploidy = sort @ploidy;

			closedir(VARSCAN);

			$varscan_data = join('/', $varscan_dir, $varscan_calls[-1]);

			if ( -l join('/', $data_directory, 'varscan_somatic_variants.tsv')) {
				unlink join('/', $data_directory, 'varscan_somatic_variants.tsv');
				}

			symlink($varscan_data, join('/', $data_directory, 'varscan_somatic_variants.tsv'));

			$ensemble_command .= " --varscan " . join('/', $data_directory, 'varscan_somatic_variants.tsv');

			# get cna calls from sequenza
			unless (scalar(@cna_calls) == 0) {

				$sequenza_data = join('/', $varscan_dir, $cna_calls[-1]);
				$segments = join('/', $varscan_dir, $seg_calls[-1]);
				$ploidy_data = join('/', $varscan_dir, $ploidy[-1]);

				if ( -l join('/', $data_directory, 'sequenza_ratio_matrix.tsv')) {
					unlink join('/', $data_directory, 'sequenza_ratio_matrix.tsv');
					}
				if ( -l join('/', $data_directory, 'sequenza_cna_metrics.tsv')) {
					unlink join('/', $data_directory, 'sequenza_cna_metrics.tsv');
					}

				symlink($sequenza_data, join('/', $data_directory, 'sequenza_ratio_matrix.tsv'));
				symlink($ploidy_data, join('/', $data_directory, 'sequenza_cna_metrics.tsv'));

				# plot CNA summary
				my $cna_plot_command = join(' ',
					"Rscript $cwd/report/plot_seqz_cna_summary.R",
					'-p', $tool_data->{project_name},
					'-o', $cnv_directory,
					'--report', $plot_directory,
					'-c', join('/', $data_directory, 'sequenza_ratio_matrix.tsv'),
					'-m', join('/', $data_directory, 'sequenza_cna_metrics.tsv'),
					'-s', 'ratio'
					);

				# run command
				print $log "Submitting job to plot somatic copy-number variants...\n";
				$run_script = write_script(
					log_dir		=> $log_directory,
					name		=> 'plot_seqz_cna_summary',
					cmd		=> $cna_plot_command,
					modules		=> [$r_version],
					max_time	=> '04:00:00',
					mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
					hpc_driver	=> $args{cluster},
					extra_args	=> [$hpc_group]
					);

				$run_id = submit_job(
					jobname		=> 'plot_seqz_cna_summary',
					shell_command	=> $run_script,
					hpc_driver	=> $args{cluster},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @job_ids, $run_id;
				}
			}

		# get CNAs calls from GATK:CNV
		if ('Y' eq $tool_set{'gatk_cnv'}) {
			my $gatk_cnv_dir = join('/', $output_directory, 'GATK_CNV');

			opendir(GATKCNV, $gatk_cnv_dir) or die "Cannot open '$gatk_cnv_dir' !";
			my @gatk_cnv_files = readdir(GATKCNV);
			my @cnv_files = grep { /ratio_gene_matrix.tsv/ } @gatk_cnv_files;
			@cnv_files = sort @cnv_files;
			my @pga_files = grep { /pga_estimates.tsv/ } @gatk_cnv_files;
			@pga_files = sort @pga_files;
			closedir(GATKCNV);

			my $cnv_data = join('/', $gatk_cnv_dir, $cnv_files[-1]);
			my $pga_data = join('/', $gatk_cnv_dir, $pga_files[-1]);

			if ( -l join('/', $data_directory, 'gatk_cnv_ratio_matrix.tsv')) {
				unlink join('/', $data_directory, 'gatk_cnv_ratio_matrix.tsv');
				}

			if ( -l join('/', $data_directory, 'gatk_cnv_pga_estimates.tsv')) {
				unlink join('/', $data_directory, 'gatk_cnv_pga_estimates.tsv');
				}

			symlink($cnv_data, join('/', $data_directory, 'gatk_cnv_ratio_matrix.tsv'));
			symlink($pga_data, join('/', $data_directory, 'gatk_cnv_pga_estimates.tsv'));

			# plot CNV summary
			my $cnv_plot_command = join(' ',
				"Rscript $cwd/report/plot_gatk_cna_summary.R",
				'-p', $tool_data->{project_name},
				'-o', $cnv_directory,
				'--report', $plot_directory,
				'-c', join('/', $data_directory, 'gatk_cnv_ratio_matrix.tsv'),
				'-m', join('/', $data_directory, 'gatk_cnv_pga_estimates.tsv'),
				'-s', 'ratio'
				);

			# run command
			print $log "Submitting job to plot GATK CNVs...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_gatk_cnv_summary',
				cmd		=> $cnv_plot_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'plot_gatk_cnv_summary',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;
			}

		# get CNAs calls from ASCAT
		if ('Y' eq $tool_set{'ascat'}) {
			my $ascat_dir = join('/', $output_directory, 'ASCAT');

			opendir(ASCATCNV, $ascat_dir) or die "Cannot open '$ascat_dir' !";
			my @ascat_files = readdir(ASCATCNV);
			my @cnv_files = grep { /ratio_gene_matrix.tsv/ } @ascat_files;
			@cnv_files = sort @cnv_files;
			my @pga_files = grep { /purity_ploidy.tsv/ } @ascat_files;
			@pga_files = sort @pga_files;
			closedir(ASCATCNV);

			my $cnv_data = join('/', $ascat_dir, $cnv_files[-1]);
			my $pga_data = join('/', $ascat_dir, $pga_files[-1]);

			if ( -l join('/', $data_directory, 'ascat_cnv_ratio_matrix.tsv')) {
				unlink join('/', $data_directory, 'ascat_cnv_ratio_matrix.tsv');
				}

			if ( -l join('/', $data_directory, 'ascat_purity_estimates.tsv')) {
				unlink join('/', $data_directory, 'ascat_purity_estimates.tsv');
				}

			symlink($cnv_data, join('/', $data_directory, 'ascat_cnv_ratio_matrix.tsv'));
			symlink($pga_data, join('/', $data_directory, 'ascat_purity_estimates.tsv'));

			# plot CNV summary
			my $cnv_plot_command = join(' ',
				"Rscript $cwd/report/plot_ascat_summary.R",
				'-p', $tool_data->{project_name},
				'-o', $cnv_directory,
				'--report', $plot_directory,
				'-c', join('/', $data_directory, 'ascat_cnv_ratio_matrix.tsv'),
				'-m', join('/', $data_directory, 'ascat_purity_estimates.tsv'),
				'-s', 'ratio'
				);

			# run command
			print $log "Submitting job to plot ASCAT CNVs...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_ascat_summary',
				cmd		=> $cnv_plot_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'plot_ascat_summary',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;
			}

		# get CNAs calls from ichorCNA
		if ('Y' eq $tool_set{'ichor_cna'}) {
			my $ichor_cna_dir = join('/', $output_directory, 'IchorCNA');

			opendir(ICHOR, $ichor_cna_dir) or die "Cannot open '$ichor_cna_dir' !";
			my @ichor_cna_files = readdir(ICHOR);
			my @cnv_files = grep { /perbin_cna_status.tsv/ } @ichor_cna_files;
			@cnv_files = sort @cnv_files;
			my @metric_files = grep { /ichorCNA_estimates.tsv/ } @ichor_cna_files;
			@metric_files = sort @metric_files;
			closedir(ICHOR);

			my $cnv_data = join('/', $ichor_cna_dir, $cnv_files[-1]);
			my $metric_data = join('/', $ichor_cna_dir, $metric_files[-1]);

			if ( -l join('/', $data_directory, 'ichor_cna_output.tsv')) {
				unlink join('/', $data_directory, 'ichor_cna_output.tsv');
				}

			if ( -l join('/', $data_directory, 'ichor_cna_estimates.tsv')) {
				unlink join('/', $data_directory, 'ichor_cna_estimates.tsv');
				}

			symlink($cnv_data, join('/', $data_directory, 'ichor_cna_output.tsv'));
			symlink($metric_data, join('/', $data_directory, 'ichor_cna_estimates.tsv'));

			# plot CNV summary
			my $cnv_plot_command = join(' ',
				"Rscript $cwd/report/plot_ichor_cna_summary.R",
				'-p', $tool_data->{project_name},
				'-o', $cnv_directory,
				'--report', $plot_directory,
				'-c', join('/', $data_directory, 'ichor_cna_output.tsv'),
				'-m', join('/', $data_directory, 'ichor_cna_estimates.tsv'),
				'-s', 'ratio'
				);

			# run command
			print $log "Submitting job to plot Ichor CNAs...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_ichor_cna_summary',
				cmd		=> $cnv_plot_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'plot_ichor_cna_summary',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;
			}

		# get CNAs calls from PanelCN.mops
		if ('Y' eq $tool_set{'mops'}) {
			my $mops_cnv_dir = join('/', $output_directory, 'panelCNmops');

			opendir(MOPSCNV, $mops_cnv_dir) or die "Cannot open '$mops_cnv_dir' !";
			my @mops_cnv_files = readdir(MOPSCNV);
			my @cnv_files = grep { /panelCN.mops_output.tsv/ } @mops_cnv_files;
			@cnv_files = sort @cnv_files;
			closedir(MOPSCNV);

			my $cnv_data = join('/', $mops_cnv_dir, $cnv_files[-1]);

			if ( -l join('/', $data_directory, 'mops_combined_output.tsv')) {
				unlink join('/', $data_directory, 'mops_combined_output.tsv');
				}

			symlink($cnv_data, join('/', $data_directory, 'mops_combined_output.tsv'));

			# plot SV summary
			my $cnv_plot_command = join(' ',
				"Rscript $cwd/report/plot_cn_mops_summary.R",
				'-p', $tool_data->{project_name},
				'-o', $cnv_directory,
				'--report', $plot_directory,
				'-c', join('/', $data_directory, 'mops_combined_output.tsv'),
				'-s', 'ratio'
				);

			# run command
			print $log "Submitting job to plot PanelCN.mops CNVs...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_mops_cnv_summary',
				cmd		=> $cnv_plot_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'plot_mops_cnv_summary',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;
			}

		# get SVs calls from MAVIS
		if ('Y' eq $tool_set{'mavis'}) {
			my $mavis_dir = join('/', $output_directory, 'Mavis');

			opendir(MAVIS, $mavis_dir) or die "Cannot open '$mavis_dir' !";
			my @mavis_files;
			if ('targeted' eq $tool_data->{seq_type}) {
				@mavis_files = grep { /mavis_output_filtered.tsv/ } readdir(MAVIS);
				} else {
				@mavis_files = grep { /mavis_output.tsv/ } readdir(MAVIS);
				}
			@mavis_files = sort @mavis_files;
			closedir(MAVIS);

			my $sv_data = join('/', $mavis_dir, $mavis_files[-1]);

			if ( -l join('/', $data_directory, 'mavis_sv_data.tsv')) {
				unlink join('/', $data_directory, 'mavis_sv_data.tsv');
				}

			symlink($sv_data, join('/', $data_directory, 'mavis_sv_data.tsv'));

			# plot SV summary
			my $sv_plot_command = join(' ',
				"Rscript $cwd/report/plot_sv_summary.R",
				'-p', $tool_data->{project_name},
				'-o', $sv_directory,
				'-r', $tool_data->{ref_type},
				'--report', $plot_directory,
				'-m', join('/', $data_directory, 'mavis_sv_data.tsv')
				);

			# run command
			print $log "Submitting job to plot SVs...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_sv_summary',
				cmd		=> $sv_plot_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$run_id = submit_job(
				jobname		=> 'plot_sv_summary',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $run_id;
			}

		# get MSI estimates
		if ('Y' eq $tool_set{'msi'}) {
			my $msi_dir = join('/', $output_directory, 'MSI');
			opendir(MSI, $msi_dir) or die "Cannot open '$msi_dir' !";
			my @msi_estimates = grep { /msi_estimates.tsv/ } readdir(MSI);
			@msi_estimates = sort @msi_estimates;
			closedir(MSI);

			$msi_data = join('/', $msi_dir, $msi_estimates[-1]);

			if ( -l join('/', $data_directory, 'msi_estimates.tsv')) {
				unlink join('/', $data_directory, 'msi_estimates.tsv');
				}

			symlink($msi_data, join('/', $data_directory, 'msi_estimates.tsv'));
			}

		# find ENSEMBLE mutations
		# run command
		print $log "Submitting job to collect somatic variants...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'collect_somatic_variant_calls',
			cmd		=> $ensemble_command,
			modules		=> [$r_version],
			max_time	=> '48:00:00',
			mem		=> ('wgs' eq $tool_data->{seq_type}) ? '48G' : '8G',
			hpc_driver	=> $args{cluster},
			extra_args	=> [$hpc_group]
			);

		my $ensemble_run_id = submit_job(
			jobname		=> 'collect_somatic_variant_calls',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $ensemble_run_id;

		# split ensemble calls to produce ensemble VCFs
		my $maf2vcf_command = 'maf2vcf.pl';
		if (defined($tool_data->{annotate}->{vcf2maf_path})) {
			$maf2vcf_command = "perl $tool_data->{annotate}->{vcf2maf_path}";
			$maf2vcf_command =~ s/vcf2maf.pl/maf2vcf.pl/;
			}

		$maf2vcf_command .= ' ' . join(' ',
			'--input-maf', join('/', $data_directory, 'ENSEMBLE', 'ensemble_mutation_data.tsv'),
			'--output-dir', join('/', $data_directory, 'ENSEMBLE', 'vcfs'),
			'--ref-fasta', $tool_data->{reference},
			'--per-tn-vcfs'
			);

		$maf2vcf_command .= "\n\ncd " . join('/', $data_directory, 'ENSEMBLE', 'vcfs');
		$maf2vcf_command .= "\n" . 'for VCF in *vcf; do';
		$maf2vcf_command .= "\n" . join("\n",
			'  vcf-sort -c $VCF | bgzip -f > $VCF.gz',
			'  tabix $VCF.gz',
			'done'
			);

		# run command
		print $log "Submitting job to run maf2vcf on ENSEMBLE calls...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'produce_ensemble_vcfs',
			cmd		=> $maf2vcf_command,
			modules		=> ['perl',$vcf2maf,$vcftools,$samtools,'tabix'],
			dependencies	=> $ensemble_run_id,
			max_time	=> '08:00:00',
			mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
			hpc_driver	=> $args{cluster},
			extra_args	=> [$hpc_group]
			);

		my $maf2vcf_run_id = submit_job(
			jobname		=> 'produce_ensemble_vcfs',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $maf2vcf_run_id;

		my ($toolsummary_run_id, $tmb_run_id, $mutsig_run_id) = '';
		my ($cosmicsbs_run_id, $chord_run_id, $hrdetect_run_id) = '';

		# plot mutation overlap (by tool)
		my $snv_overlap_command = join(' ',
			"Rscript $cwd/report/plot_snv_tool_summary.R",
			'-p', $tool_data->{project_name},
			'-o', $snv_directory,
			'--report', $plot_directory,
			'-i', join('/', $data_directory, 'ENSEMBLE', 'CombinedMutationData.RData')
			);

		# run command
		print $log "Submitting job to plot tool overlap for somatic variants...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'plot_snv_tool_overlap',
			cmd		=> $snv_overlap_command,
			modules		=> [$r_version],
			dependencies	=> $ensemble_run_id,
			max_time	=> '24:00:00',
			mem		=> ('wgs' eq $tool_data->{seq_type}) ? '8G' : '2G',
			hpc_driver	=> $args{cluster},
			extra_args	=> [$hpc_group]
			);

		$toolsummary_run_id = submit_job(
			jobname		=> 'plot_snv_tool_summary',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $toolsummary_run_id;

		# calculate TMB
		my $tmb_command = join(' ',
			"Rscript $cwd/report/calculate_tmb.R",
			'-p', $tool_data->{project_name},
			'-o', $snv_directory,
			'-r', $tool_data->{ref_type},
			'--maf', join('/', $data_directory, 'ENSEMBLE', 'ensemble_mutation_data.tsv'),
			'--callable', join('/', $data_directory, 'total_bases_covered.tsv'),
			'--method both'
			);

		my $tmb_data = join('/', $snv_directory, 'tumour_mutation_burden.tsv');

		# run command
		print $log "Submitting job to calculate tumour mutation burden...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'calculate_tmb',
			cmd		=> $tmb_command,
			modules		=> [$r_version],
			dependencies	=> $ensemble_run_id,
			max_time	=> '24:00:00',
			mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
			hpc_driver	=> $args{cluster},
			extra_args	=> [$hpc_group]
			);

		$tmb_run_id = submit_job(
			jobname		=> 'calculate_tmb',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $tmb_run_id;

		# run MutSigCV
		my $mutsig_command = join(' ',
			'sh /cluster/tools/software/MutSigCV/1.4/run_MutSigCV.sh',
			'/cluster/tools/software/MCR/8.1/v81',
			join('/', $data_directory, 'ENSEMBLE', 'ensemble_mutation_data.tsv'),
			'/cluster/projects/pughlab/references/MutSigCV/exome_full192.coverage.txt',
			'/cluster/projects/pughlab/references/MutSigCV/gene.covariates.txt',
			join('/', $plot_directory, $tool_data->{project_name} . '_MutSigCV'),
			'/cluster/projects/pughlab/references/MutSigCV/mutation_type_dictionary_file.txt'
			);

		if ( ('hg19' eq $tool_data->{ref_type}) || ('GRCh37' eq $tool_data->{ref_type}) ) {
			$mutsig_command .= ' /cluster/projects/pughlab/references/MutSigCV/chr_files_hg19';
			} elsif (('hg38' eq $tool_data->{ref_type}) || ('GRCh38' eq $tool_data->{ref_type})) {
			$mutsig_command .= ' /cluster/projects/pughlab/references/MutSigCV/chr_files_hg38';
			} else {
			print $log "MutSigCV requested but unknown ref_type provided; not running.\n";
			$tool_set{'mutsig'} = 'N';
			}

		my $significance_data = join('/',
			$plot_directory,
			$tool_data->{project_name} . '_MutSigCV.sig_genes.txt'
			);

		if ('Y' eq $tool_set{'mutsig'}) {

			# run command
			print $log "Submitting job to check mutation significance...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'run_mutsigcv',
				cmd		=> $mutsig_command,
				modules		=> ['MutSigCV/1.4'],
				dependencies	=> $ensemble_run_id,
				max_time	=> '06:00:00',
				mem		=> '8G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$mutsig_run_id = submit_job(
				jobname		=> 'run_mutsigcv',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $mutsig_run_id;
			}

		# plot mutation signatures
		if ('Y' eq $tool_set{'cosmic_sbs'}) {

			my $sbs_directory = join('/', $sig_directory, 'COSMIC_SBS');
			unless(-e $sbs_directory) { make_path($sbs_directory); }

			my $sig_plot_command = join(' ',
				"Rscript $cwd/report/apply_cosmic_mutation_signatures.R",
				'-p', $tool_data->{project_name},
				'-o', $sbs_directory,
				'--report', $plot_directory,
				'--maf', join('/', $data_directory, 'ENSEMBLE', 'ensemble_mutation_data.tsv'),
				'-r', $tool_data->{ref_type},
				'-t', $tool_data->{seq_type}
				);

			if (defined($tool_data->{mutation_signatures})) {
				$sig_plot_command .= " -s $tool_data->{mutation_signatures}";
				}

			# run command
			print $log "Submitting job to check mutation signatures...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_snv_mutation_signatures',
				cmd		=> $sig_plot_command,
				modules		=> [$r_version],
				dependencies	=> $ensemble_run_id,
				max_time	=> '12:00:00',
				mem		=> ('wgs' eq $tool_data->{seq_type}) ? '8G' : '4G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$cosmicsbs_run_id = submit_job(
				jobname		=> 'plot_snv_mutation_signatures',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $cosmicsbs_run_id;
			}

		# plot HRD signatures
		my ($chord_data, $hrdetect_data) = '';
		if ('Y' eq $tool_set{'chord'}) {

			my $chord_directory = join('/', $sig_directory, 'CHORD');
			unless(-e $chord_directory) { make_path($chord_directory); }

			$chord_data = join('/', $chord_directory, 'CHORD_predictions.tsv');

			my $sig_plot_command = join(' ',
				"Rscript $cwd/report/get_chord_signatures.R",
				'-p', $tool_data->{project_name},
				'-o', $chord_directory,
				'--report', $plot_directory,
				'-m', join('/', $data_directory, 'ENSEMBLE', 'ensemble_mutation_data.tsv'),
				'-s', join('/', $data_directory, 'mavis_sv_data.tsv'), 
				'-r', $tool_data->{ref_type}
				);

			if (defined($tool_data->{summarize_steps}->{lib_path})) {
				$sig_plot_command .= ' -l ' . $tool_data->{summarize_steps}->{lib_path};
				}

			# run command
			print $log "Submitting job to check CHORD signatures...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_chord_signatures',
				cmd		=> $sig_plot_command,
				modules		=> [$r_version],
				dependencies	=> $ensemble_run_id,
				max_time	=> '12:00:00',
				mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$chord_run_id = submit_job(
				jobname		=> 'plot_chord_signatures',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $chord_run_id;
			}

		if ('Y' eq $tool_set{'hrdetect'}) {

			my $hrdetect_directory = join('/', $sig_directory, 'HRDetect');
			unless(-e $hrdetect_directory) { make_path($hrdetect_directory); }

			$hrdetect_data = join('/', $hrdetect_directory, 'HRDetect_scores.tsv');

			my $sig_plot_command = join(' ',
				"Rscript $cwd/report/get_hrdetect_signatures.R",
				'-p', $tool_data->{project_name},
				'-o', $hrdetect_directory,
				'--report', $plot_directory,
				'--vcf_dir', join('/', $data_directory, 'ENSEMBLE', 'vcfs'),
				'-s', join('/', $data_directory, 'mavis_sv_data.tsv'),
				'-c', $segments,
				'-r', $tool_data->{ref_type}
				);

			if (defined($tool_data->{summarize_steps}->{lib_path})) {
				$sig_plot_command .= ' -l ' . $tool_data->{summarize_steps}->{lib_path};
				}

			# run command
			print $log "Submitting job to check HRDetect signatures...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_hrdetect_signatures',
				cmd		=> $sig_plot_command,
				modules		=> [$r_version],
				dependencies	=> join(':', $ensemble_run_id, $maf2vcf_run_id),
				max_time	=> '12:00:00',
				mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
				hpc_driver	=> $args{cluster},
				extra_args	=> [$hpc_group]
				);

			$hrdetect_run_id = submit_job(
				jobname		=> 'plot_hrdetect_signatures',
				shell_command	=> $run_script,
				hpc_driver	=> $args{cluster},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @job_ids, $hrdetect_run_id;
			}

		# plot mutation summary
		my @snv_dependencies = ( $ensemble_run_id );

		my $snv_plot_command = join(' ',
			"Rscript $cwd/report/summarize_mutation_data.R",
			'-p', $tool_data->{project_name},
			'-o', $snv_directory,
			'--report', $plot_directory,
			'--maf', join('/', $data_directory, 'ENSEMBLE', 'ensemble_mutation_data.tsv')
			);

		if ('Y' eq $tool_set{'msi'}) {
			$snv_plot_command .= " --msi " . join('/', $data_directory, 'msi_estimates.tsv');
			}

		if ('Y' eq $tool_set{'chord'}) {
			$snv_plot_command .= " --chord $chord_data";
			push @snv_dependencies, $chord_run_id;
			}

		if ('Y' eq $tool_set{'hrdetect'}) {
			$snv_plot_command .= " --hrdetect $hrdetect_data";
			push @snv_dependencies, $hrdetect_run_id;
			}

		# run command
		print $log "Submitting job to summarize (short) somatic variants...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'summarize_mutations',
			cmd		=> $snv_plot_command,
			modules		=> [$r_version],
			dependencies	=> join(':', @snv_dependencies),
			max_time	=> '12:00:00',
			mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
			hpc_driver	=> $args{cluster},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'summarize_mutations',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;

		# plot recurrently mutated geneset
		@snv_dependencies = ( $ensemble_run_id, $tmb_run_id );

		$snv_plot_command = join(' ',
			"Rscript $cwd/report/plot_recurrent_mutations.R",
			'-p', $tool_data->{project_name},
			'-o', $snv_directory,
			'--report', $plot_directory,
			'--maf', join('/', $data_directory, 'ENSEMBLE', 'ensemble_mutation_data.tsv'),
			'--tmb', $tmb_data,
			'-t', $tool_data->{seq_type}
			);

		if ('Y' eq $tool_set{'mutsig'}) {
			$snv_plot_command .= " --mutsig $significance_data";
			push @snv_dependencies, $mutsig_run_id;
			}

		# run command
		print $log "Submitting job to plot (short) somatic variants...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'plot_recurrent_geneset',
			cmd		=> $snv_plot_command,
			modules		=> [$r_version],
			dependencies	=> join(':', @snv_dependencies),
			max_time	=> '12:00:00',
			mem		=> ('wgs' eq $tool_data->{seq_type}) ? '4G' : '2G',
			hpc_driver	=> $args{cluster},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'plot_recurrent_geneset',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;

		}

	# write some methods
	my $methods_command = "perl $cwd/report/write_wxs_methods.pl";
	if ('wgs' eq $tool_data->{seq_type}) {
		$methods_command = "perl $cwd/report/write_wgs_methods.pl";
		} elsif ('rna' eq $tool_data->{seq_type}) {
		$methods_command = "perl $cwd/report/write_rna_methods.pl";
		} elsif ('targeted' eq $tool_data->{seq_type}) {
		$methods_command = "perl $cwd/report/write_targetseq_methods.pl";
		}

	$methods_command .= " -t $args{config} -d $plot_directory";

	# run command
	print $log "Submitting job to write methods...\n";
	$run_script = write_script(
		log_dir		=> $log_directory,
		name		=> 'write_pipeline_methods',
		cmd		=> $methods_command,
		modules		=> ['perl'],
		mem		=> '256M',
		hpc_driver	=> $args{cluster},
		extra_args	=> [$hpc_group]
		);

	$run_id = submit_job(
		jobname		=> 'write_pipeline_methods',
		shell_command	=> $run_script,
		hpc_driver	=> $args{cluster},
		dry_run		=> $args{dry_run},
		log_file	=> $log
		);

	push @job_ids, $run_id;

	# generate report.tex
	my $report_command = "cd $report_directory\n\n";
	$report_command .= join(' ',
		"perl $cwd/report/generate_report.pl",
		'-o', $report_directory,
		'-i', $plot_directory,
		'-t', $tool_data->{project_name},
		'-d', $run_date
		);

	# run pdflatex (twice, to ensure proper page numbering)
	$report_command .= "\n\n" . join("\n",
		"pdflatex Report.tex",
		"pdflatex Report.tex"
		);

	$report_command .= "\n\n" . join("\n",
		'if [ 0 == $? ]; then',
		'  echo Generated and compiled LaTeX file with pdflatex',
		'else',
		'  echo Error running pdflatex: $?',
		'fi'
		);
	
	# run command
	print $log "Submitting job to generate report file...\n";
	$run_script = write_script(
		log_dir		=> $log_directory,
		name		=> 'create_report',
		cmd		=> $report_command,
		dependencies	=> join(':', @job_ids),
		modules		=> ['perl','texlive'],
		mem		=> '256M',
		hpc_driver	=> $args{cluster},
		extra_args	=> [$hpc_group]
		);

	if ($args{report}) {
		$run_id = submit_job(
			jobname		=> 'create_report',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);
		}

	# if this is not a dry run OR there are jobs to assess (run or resumed with jobs submitted) then
	# collect job metrics (exit status, mem, run time)
	unless ( ($args{dry_run}) || (scalar(@job_ids) == 0) ) {

		# collect job stats
		my $collect_metrics = collect_job_stats(
			job_ids		=> join(',', @job_ids),
			outfile		=> $outfile,
			hpc_driver	=> $args{cluster}
			);

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'output_job_metrics_' . $run_count,
			cmd	=> $collect_metrics,
			dependencies	=> join(':', @job_ids),
			mem		=> '256M',
			hpc_driver	=> $args{cluster},
			kill_on_error	=> 0,
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'output_job_metrics',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;

		# do some logging
		print "Number of jobs submitted: " . scalar(@job_ids) . "\n";

		my $n_queued = `squeue -r | wc -l`;
		print "Total number of jobs in queue: " . $n_queued . "\n";

		# wait until it finishes
		unless ($args{no_wait}) {
			check_final_status(job_id => $run_id);
			}
		}
	}
 
### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($help, $tool_config, $hpc_driver, $dry_run, $run_date, $no_wait, $create_report);

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	't|tools=s'	=> \$tool_config,
	'd|date=s'	=> \$run_date,
	'create_report'	=> \$create_report,
	'c|cluster=s'	=> \$hpc_driver,
	'dry-run'	=> \$dry_run,
	'no-wait'	=> \$no_wait
	);

if ($help) {
	my $help_msg = join("\n",
		"Options:",
		"\t--help|-h\tPrint this help message",
		"\t--tools|-t\t<string> master tool config (yaml format)",
		"\t--date|-d\t<string> Date the pipeline was initiated",
		"\t--create_report\t<boolean> should the final report (pdf) be created? (default: false)",
		"\t--cluster|-c\t<string> cluster scheduler (default: slurm)",
		"\t--dry-run\t<boolean> should jobs be submitted? (default: false)",
		"\t--no-wait\t<boolean> should we exit after job submission (true) or wait until all jobs have completed (false)? (default: false)"
		);

	print "$help_msg\n";
	exit;
	}

# do some quick error checks to confirm valid arguments	
main(
	config		=> $tool_config,
	run_date	=> $run_date,
	report		=> $create_report,
	cluster		=> $hpc_driver,
	dry_run		=> $dry_run,
	no_wait		=> $no_wait
	);
