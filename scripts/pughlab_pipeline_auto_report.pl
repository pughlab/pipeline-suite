#!/usr/bin/env perl
### pughlab_pipeline_auto_report.pl ################################################################
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

my $cwd = dirname($0);
require "$cwd/utilities.pl";

####################################################################################################
# version       author		comment
# 1.0		sprokopec       tool to automatically generate reports

### MAIN ###########################################################################################
sub main {
	my %args = (
		config		=> undef,
		cluster		=> undef,
		dry_run		=> undef,
		run_date	=> undef,
		no_wait		=> undef,
		@_
		);

	my $tool_data = LoadFile($args{config});
	my $run_date = $args{run_date};
	my $current_date = strftime "%F", localtime;

	# check for and/or create output directories
	my $output_directory = $tool_data->{output_dir};
	$output_directory =~ s/\/$//;

	my $report_directory = join('/', $output_directory, 'Report');
	unless(-e $report_directory) { make_path($report_directory); }

	my $log_directory = join('/', $report_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $data_directory = join('/', $report_directory, 'data');
	unless(-e $data_directory) { make_path($data_directory); }

	my $plot_directory = join('/', $report_directory, 'plots');
	unless(-e $plot_directory) { make_path($plot_directory); }

	# Initiate logging
	my $log_file = join('/', $log_directory, 'run_Report_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_Report_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running Report Pipeline...\n";

	# get tool versions
	my $r_version = 'R/' . $tool_data->{r_version};

	### RUN ####################################################################################

	my @job_ids;
	my ($qc_dir, $germ_dir, $cpsr_dir);
	my ($correlations, $qc_data, $cb_data, $seqqc_data, $contest_data, $cpsr_calls) = undef;
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
			max_time	=> '01:00:00',
			mem		=> '2G',
			hpc_driver	=> $args{cluster}
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
		if ('Y' eq $tool_data->{rsem}->{run}) {
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
				dependencies	=> $run_id,
				max_time	=> '04:00:00',
				mem		=> '2G',
				hpc_driver	=> $args{cluster}
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
		if ('Y' eq $tool_data->{haplotype_caller}->{run}) {
			my $snv_dir = join('/', $output_directory, 'HaplotypeCaller');
			opendir(VARIANTS, $snv_dir) or die "Cannot open '$snv_dir' !";
			my @snv_calls = grep { /_variant_by_patient.tsv/ } readdir(VARIANTS);
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
				hpc_driver	=> $args{cluster}
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
		if ('Y' eq $tool_data->{star_fusion}->{run}) {
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

		if ('Y' eq $tool_data->{fusioncatcher}->{run}) {
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
				hpc_driver	=> $args{cluster}
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
				hpc_driver	=> $args{cluster}
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
			my @seqqc_files = grep { /ContaminationEstimates.tsv/ } readdir(QC);
			@seqqc_files = sort @seqqc_files;
			closedir(QC);

			$seqqc_data = join('/', $qc_dir, 'SequenceMetrics', $seqqc_files[-1]);

			if ( -l join('/', $data_directory, 'ContaminationEstimates.tsv')) {
				unlink join('/', $data_directory, 'ContaminationEstimates.tsv');
				}

			symlink($seqqc_data, join('/', $data_directory, 'ContaminationEstimates.tsv'));
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
			'-t', $tool_data->{seq_type},
			'-c', $contest_data,
			'-m', $seqqc_data
			);

		# run command
		print $log "Submitting job to create QC plots...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'create_qc_plots',
			cmd		=> $qc_command,
			modules		=> [$r_version],
			max_time	=> '01:00:00',
			mem		=> '2G',
			hpc_driver	=> $args{cluster}
			);

		my $qc_run_id = submit_job(
			jobname		=> 'create_qc_plots',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $qc_run_id;

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
			my $plot_command = "cd $output_directory\n";
			$plot_command .= "Rscript $cwd/report/plot_germline_snv_summary.R";
			$plot_command .= " " . join(' ',
				'-m', $cpsr_calls,
				'-o', $plot_directory,
				'-p', $tool_data->{project_name}
				);

			# run command
			print $log "Submitting job to create germline SNV plots...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_germline_snv_summary',
				cmd		=> $plot_command,
				modules		=> [$r_version],
				max_time	=> '01:00:00',
				mem		=> '2G',
				hpc_driver	=> $args{cluster}
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
		my ($sequenza_data, $ploidy_data, $gatk_cnv, $gatk_pga, $msi_data);
		my $n_tools = 0;

		# find ENSEMBLE mutations
		my $ensemble_command .= join(' ',
			"Rscript $cwd/report/format_ensemble_mutations.R",
			'-p', $tool_data->{project_name},
			'-o', $plot_directory
			);

		# get mutation calls from MuTect
		if ('Y' eq $tool_data->{mutect}->{run}) {
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

			$ensemble_command .= " --mutect $mutect_data";
			$n_tools++;
			}

		# get mutation calls from MuTect2
		if ('Y' eq $tool_data->{mutect2}->{run}) {
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

			$ensemble_command .= " --mutect2 $mutect2_data";
			$n_tools++;
			}

		# get mutation calls from Strelka
		if ('Y' eq $tool_data->{strelka}->{run}) {
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

			$ensemble_command .= " --strelka $strelka_data";
			$n_tools++;
			}

		# get mutation calls from SomaticSniper
		if ('Y' eq $tool_data->{somaticsniper}->{run}) {
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

			$ensemble_command .= " --somaticsniper $sniper_data";
			$n_tools++;
			}

		# get mutation calls from VarDict
		if ('Y' eq $tool_data->{vardict}->{run}) {
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

			$ensemble_command .= " --vardict $vardict_data";
			$n_tools++;
			}


		# get mutation calls from VarScan
		my (@cna_calls, @ploidy);
		if ('Y' eq $tool_data->{varscan}->{run}) {
			my $varscan_dir = join('/', $output_directory, 'VarScan');

			opendir(VARSCAN, $varscan_dir) or die "Cannot open '$varscan_dir' !";
			my @varscan_files = readdir(VARSCAN);

			my @varscan_calls = grep { /mutations_for_cbioportal.tsv/ } @varscan_files;
			@varscan_calls = sort @varscan_calls;

			@cna_calls = grep { /Sequenza_ratio_gene_matrix.tsv/ } @varscan_files;
			@cna_calls = sort @cna_calls;

			@ploidy = grep { /Sequenza_ploidy_purity.tsv/ } @varscan_files;
			@ploidy = sort @ploidy;

			closedir(VARSCAN);

			$varscan_data = join('/', $varscan_dir, $varscan_calls[-1]);

			if ( -l join('/', $data_directory, 'varscan_somatic_variants.tsv')) {
				unlink join('/', $data_directory, 'varscan_somatic_variants.tsv');
				}

			symlink($varscan_data, join('/', $data_directory, 'varscan_somatic_variants.tsv'));

			$ensemble_command .= " --varscan $varscan_data";
			$n_tools++;

			# get cna calls from sequenza
			unless (scalar(@cna_calls) == 0) {

				$sequenza_data = join('/', $varscan_dir, $cna_calls[-1]);
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
					"Rscript $cwd/report/plot_cna_summary.R",
					'-p', $tool_data->{project_name},
					'-o', $plot_directory,
					'-c', $sequenza_data,
					'-m', $ploidy_data,
					'-s', 'ratio'
					);

				# run command
				print $log "Submitting job to plot somatic copy-number variants...\n";
				$run_script = write_script(
					log_dir		=> $log_directory,
					name		=> 'plot_cna_summary',
					cmd		=> $cna_plot_command,
					modules		=> [$r_version],
					max_time	=> '04:00:00',
					mem		=> '2G',
					hpc_driver	=> $args{cluster}
					);

				$run_id = submit_job(
					jobname		=> 'plot_cna_summary',
					shell_command	=> $run_script,
					hpc_driver	=> $args{cluster},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @job_ids, $run_id;
				}
			}

		# get CNAs calls from GATK:CNV
		if ('Y' eq $tool_data->{gatk_cnv}->{run}) {
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

			# plot SV summary
			my $cnv_plot_command = join(' ',
				"Rscript $cwd/report/plot_gatk_cna_summary.R",
				'-p', $tool_data->{project_name},
				'-o', $plot_directory,
				'-c', $cnv_data,
				'-m', $pga_data,
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
				mem		=> '2G',
				hpc_driver	=> $args{cluster}
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

		# get SVs calls from MAVIS
		if ('Y' eq $tool_data->{mavis}->{run}) {
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

			# plot SV summary
			my $sv_plot_command = join(' ',
				"Rscript $cwd/report/plot_sv_summary.R",
				'-p', $tool_data->{project_name},
				'-o', $plot_directory,
				'-m', $sv_data
				);

			# run command
			print $log "Submitting job to plot SVs...\n";
			$run_script = write_script(
				log_dir		=> $log_directory,
				name		=> 'plot_sv_summary',
				cmd		=> $sv_plot_command,
				modules		=> [$r_version],
				max_time	=> '04:00:00',
				mem		=> '2G',
				hpc_driver	=> $args{cluster}
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
		if ('Y' eq $tool_data->{other_tools}->{msi_run}) {
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
			dependencies	=> $qc_run_id,
			modules		=> [$r_version],
			max_time	=> '48:00:00',
			mem		=> '4G',
			hpc_driver	=> $args{cluster}
			);

		$run_id = submit_job(
			jobname		=> 'collect_somatic_variant_calls',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;

		# plot mutation summary
		my $snv_plot_command = join(' ',
			"Rscript $cwd/report/plot_snv_summary.R",
			'-p', $tool_data->{project_name},
			'-o', $plot_directory,
			'-i', join('/', $plot_directory, 'ensemble_mutation_data.tsv'),
			'-c', $cb_data,
			'-t', $tool_data->{seq_type}
			);

		if ('Y' eq $tool_data->{other_tools}->{msi_run}) {
			$snv_plot_command .= " -m $msi_data";
			}

		# run command
		print $log "Submitting job to plot (short) somatic variants...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'plot_snv_summary',
			cmd		=> $snv_plot_command,
			modules		=> [$r_version],
			dependencies	=> $run_id,
			max_time	=> '04:00:00',
			mem		=> '2G',
			hpc_driver	=> $args{cluster}
			);

		$run_id = submit_job(
			jobname		=> 'plot_snv_summary',
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
		}

	$methods_command .= " -t $args{config} -d $plot_directory";

	# run command
	print $log "Submitting job to write methods...\n";
	$run_script = write_script(
		log_dir		=> $log_directory,
		name		=> 'write_pipeline_methods',
		cmd		=> $methods_command,
		modules		=> ['perl'],
		max_time	=> '00:30:00',
		mem		=> '256M',
		hpc_driver	=> $args{cluster}
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
		'-i', './plots',
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
		max_time	=> '0:30:00',
		mem		=> '256M',
		hpc_driver	=> $args{cluster}
		);

	$run_id = submit_job(
		jobname		=> 'create_report',
		shell_command	=> $run_script,
		hpc_driver	=> $args{cluster},
		dry_run		=> $args{dry_run},
		log_file	=> $log
		);

	# if this is not a dry run OR there are jobs to assess (run or resumed with jobs submitted) then
	# collect job metrics (exit status, mem, run time)
	unless ( ($args{dry_run}) || (scalar(@job_ids) == 0) ) {

		# collect job stats
		my $collect_metrics = collect_job_stats(
			job_ids => join(',', @job_ids),
			outfile => $outfile
			);

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'output_job_metrics_' . $run_count,
			cmd	=> $collect_metrics,
			dependencies	=> join(':', @job_ids),
			mem		=> '256M',
			hpc_driver	=> $args{cluster},
			kill_on_error	=> 0
			);

		$run_id = submit_job(
			jobname		=> 'output_job_metrics',
			shell_command	=> $run_script,
			hpc_driver	=> $args{cluster},
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
					die("Final REPORT accounting job: $run_id finished with errors.");
					}
				}
			}
		}
	}
 
### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($help, $tool_config, $hpc_driver, $dry_run, $run_date, $no_wait);

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	't|tools=s'	=> \$tool_config,
	'd|date=s'	=> \$run_date,
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
	cluster		=> $hpc_driver,
	dry_run		=> $dry_run,
	no_wait		=> $no_wait
	);
