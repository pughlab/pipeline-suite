#!/usr/bin/env perl
### multimmr_reports.pl ############################################################################
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
use List::Util qw(any none);
use Data::Dumper;

my $cwd = dirname(__FILE__);
require "$cwd/../../utilities.pl";

####################################################################################################
# version       author	  	comment
# 1.0		sprokopec       script to run PughLab DNASeq + EMSeq pipeline

### USAGE ##########################################################################################
# pughlab_dnaseq_pipeline.pl -c tool_config.yaml -d data.yaml
#
# where:
#	- dnaseq_tool_config.yaml contains tool versions and parameters, output directory, reference
#	information, etc. for the DNA-seq pipeline
#	- dnaseq_data_config.yaml contains sample information (YAML file containing paths to FASTQ 
#	files) for the DNA-seq pipeline
#	- emseq_tool_config.yaml contains tool versions and parameters, output directory, reference
#	information, etc. for the EM-seq pipeline
#	- emseq_data_config.yaml contains sample information (YAML file containing paths to FASTQ 
#	files) for the EM-seq pipeline

### MAIN ###########################################################################################
sub main {
	my %args = (
		wgs_tool_config		=> undef,
		wgs_data_config		=> undef,
		dna_tool_config		=> undef,
		dna_data_config		=> undef,
		em_tool_config		=> undef,
		em_data_config		=> undef,
		output_directory	=> undef,
		hpc_driver		=> undef,
		dry_run			=> undef,
		no_wait			=> undef,
		@_
		);

	### PREAMBLE ######################################################################################
	my $date = strftime "%F", localtime;
	my $timestamp = strftime "%F_%H-%M-%S", localtime;

	# get analysis directories
	my ($tool_data, $wgs_directory, $em_directory, $dna_directory);

	if (defined($args{wgs_tool_config})) {
		$tool_data = LoadFile($args{wgs_tool_config});
		$wgs_directory = $tool_data->{output_dir};
		}

	if (defined($args{em_tool_config})) {
		$tool_data = LoadFile($args{em_tool_config});
		$em_directory = $tool_data->{output_dir};
		}

	if (defined($args{dna_tool_config})) {
		$tool_data = LoadFile($args{dna_tool_config});
		$dna_directory = $tool_data->{output_dir};
		}

	# organize output and log directories
	my $output_directory = $args{output_directory};
	my $log_directory = join('/', $output_directory, 'logs');

	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_MultiMMR_report_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_MultiMMR_report_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running MultiMMR pipeline.\n";

	my ($run_wgs, $run_dna, $run_em);

	if ( (defined($args{wgs_tool_config})) & (defined($args{wgs_data_config})) ) {
		$run_wgs = 1;
		print $log "\n  Config used for sWGS: $args{wgs_tool_config}";
		print $log "\n  Sample config used for sWGS: $args{wgs_data_config}\n";
		}

	if ( (defined($args{dna_tool_config})) & (defined($args{dna_data_config})) ) {
		$run_dna = 1;
		print $log "\n  Config used for DNA-Seq: $args{dna_tool_config}";
		print $log "\n  Sample config used for DNA-Seq: $args{dna_data_config}\n";
		}

	if ( (defined($args{em_tool_config})) & (defined($args{em_data_config})) ) {
		$run_em = 1;
		print $log "\n  Config used for EM-Seq: $args{em_tool_config}";
		print $log "\n  Sample config used for EM-Seq: $args{em_data_config}";
		}

	print $log "\n  Output directory: $output_directory";
	print $log "\n---\n\n";

	# set tools and versions
	my $perl	= 'perl/' . $tool_data->{perl_version};
	my $r_version	= 'R/' . $tool_data->{r_version};

        # indicate maximum time limit for parent jobs to wait
	my $max_time = '12:00:00';
	my $max_mem = '2G';

	# get optional HPC group
	my $hpc_group = defined($tool_data->{hpc_group}) ? "-A $tool_data->{hpc_group}" : undef;

	### RUN ####################################################################################
	# set up directory structure
	my $data_directory = join('/', $output_directory, 'data');
	unless(-e $data_directory) { make_path($data_directory); }

	my $qc_out_directory = join('/', $output_directory, 'QC');
	unless(-e $qc_out_directory) { make_path($qc_out_directory); }

	my $summary_out_directory = join('/', $output_directory, 'SUMMARY');
	unless(-e $summary_out_directory) { make_path($summary_out_directory); }

	my $dna_out_directory = join('/', $output_directory, 'Mutations');
	unless(-e $dna_out_directory) { make_path($dna_out_directory); }

	my $em_out_directory = join('/', $output_directory, 'Methylation');
	unless(-e $em_out_directory) { make_path($em_out_directory); }

	my $reports_directory = join('/', $output_directory, 'Reports');
	unless(-e $reports_directory) { make_path($reports_directory); }

	# check which tools have been requested
	my %tool_dirs = ();

	if ($run_wgs) {
		$tool_dirs{swgs_qc} = join('/', $wgs_directory, 'BAMQC', 'SequenceMetrics');
		$tool_dirs{swgs_ichor} = join('/', $wgs_directory, 'IchorCNA');
		}

	if ($run_dna) {
		$tool_dirs{dnaseq_qc} = join('/', $dna_directory, 'BAMQC', 'SequenceMetrics');
		$tool_dirs{dnaseq_germline} = join('/', $dna_directory, 'HaplotypeCaller');
		$tool_dirs{dnaseq_mutect2} = join('/', $dna_directory, 'MuTect2');
		$tool_dirs{dnaseq_pindel} = join('/', $dna_directory, 'Pindel');
		$tool_dirs{dnaseq_somaticsniper} = join('/', $dna_directory, 'SomaticSniper');
		$tool_dirs{dnaseq_strelka} = join('/', $dna_directory, 'Strelka');
		$tool_dirs{dnaseq_vardict} = join('/', $dna_directory, 'VarDict');
		$tool_dirs{dnaseq_mops} = join('/', $dna_directory, 'panelCNmops');
		$tool_dirs{dnaseq_msi} = join('/', $dna_directory, 'panelMSI');
		$tool_dirs{dnaseq_mavis} = join('/', $dna_directory, 'Mavis');
		}

	if ($run_em) {
		$tool_dirs{emseq_qc} = join('/', $em_directory, 'BAMQC');
		$tool_dirs{emseq_methyldackel} = join('/', $em_directory, 'MethylDackel');
		}

	# initiate objects
	my (@found_files, @job_ids);
	my ($base_file, $link_file, $run_script, $run_id);
	my ($ensemble_command, $ensemble_run_id, $ensemble_maf);

	# initiate command to combine/summarize dataset
	my $summarize_panel_command = join(' ',
		"Rscript $cwd/plot_panel_summary.R",
		'-p', $tool_data->{project_name},
		'-o', $summary_out_directory,
		'-r', $reports_directory,
		'-s', $args{dna_data_config},
		'-g', $tool_data->{ref_type}
		);

	my $summarize_mutations_command = join(' ',
		"Rscript $cwd/summarize_mutations.R",
		'-p', $tool_data->{project_name},
		'-o', $dna_out_directory,
		'-r', $reports_directory,
		'-y', $args{dna_data_config}
		);

	# find the required input files for 1XWGS
	if ($run_wgs) {

		# find sWGS QC data
		my $qc_dir = join('/', $data_directory, 'sWGS', 'QC');
		unless(-e $qc_dir) { make_path($qc_dir); }

		opendir(FILES, $tool_dirs{'swgs_qc'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		my @cov_files = grep { /WGSMetrics.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'swgs_qc'}, $cov_files[-1]);
		$link_file = join('/', $qc_dir, $cov_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# create some QC plots
		my $qc_command = "Rscript $cwd/plot_qc_metrics.R";
		$qc_command .= " " . join(' ',
			'-p', $tool_data->{project_name},
			'-o', $qc_out_directory,
			'-r', $reports_directory,
			'-s', $args{wgs_data_config},
			'-c', join('/', $qc_dir, $cov_files[-1]),
			'-t wgs'
			);

		# run command
		print $log "Submitting job to create QC plots...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'create_swgs_coverage_plots',
			cmd		=> $qc_command,
			modules		=> [$r_version],
			mem		=> '2G',
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'create_swgs_coverage_plots',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;


		# find sWGS ichorCNA estimates
		my $ichor_dir = join('/', $data_directory, 'sWGS', 'CNAs');
		unless(-e $ichor_dir) { make_path($ichor_dir); }

		opendir(FILES, $tool_dirs{'swgs_ichor'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		my @tf_files = grep { /ichorCNA_estimates.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'swgs_ichor'}, $tf_files[-1]);
		$link_file = join('/', $ichor_dir, $tf_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# add ichor file to summary command
		$summarize_panel_command .= " -i $base_file";
		$summarize_mutations_command .= " -i $base_file";

		# create some TF plots
		my $tf_command = "Rscript $cwd/plot_ichor_estimates.R";
		$tf_command .= " " . join(' ',
			'-p', $tool_data->{project_name},
			'-o', $qc_out_directory,
			'-i', $tool_dirs{'swgs_ichor'},
			'-r', $reports_directory,
			'-s', $args{wgs_data_config},
			'-t', join('/', $ichor_dir, $tf_files[-1]) 
			);

		# run command
		print $log "Submitting job to create TF plots...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'create_swgs_ichor_plots',
			cmd		=> $tf_command,
			modules		=> [$r_version],
			mem		=> '2G',
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'create_swgs_ichor_plots',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;

		}

	# find the required input files for DNA-Seq
	if ($run_dna) {

		# find DNA-Seq QC data
		my $qc_dir = join('/', $data_directory, 'DNASeq', 'Coverage');
		unless(-e $qc_dir) { make_path($qc_dir); }

		opendir(FILES, $tool_dirs{'dnaseq_qc'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		my @cov_files = grep { /HSMetrics.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_qc'}, $cov_files[-1]);
		$link_file = join('/', $qc_dir, $cov_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# find DNA-Seq contamination estimates
		my $qc2_dir = join('/', $data_directory, 'DNASeq', 'Contamination');
		unless(-e $qc2_dir) { make_path($qc2_dir); }

		opendir(FILES, join('/', $tool_dirs{'dnaseq_germline'}, 'cohort/sample_comparisons'));
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		my @contam_files = grep { /gtcheck.tsv$/ } @found_files;

		closedir(FILES);

		$base_file = join('/',
			$tool_dirs{'dnaseq_germline'},
			'cohort/sample_comparisons',
			$contam_files[-1]
			);

		$link_file = join('/', $qc2_dir, $contam_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# create some QC plots
		my $qc_command = "Rscript $cwd/plot_qc_metrics.R";
		$qc_command .= " " . join(' ',
			'-p', $tool_data->{project_name},
			'-o', $qc_out_directory,
			'-r', $reports_directory,
			'-s', $args{dna_data_config},
			'-c', join('/', $qc_dir, $cov_files[-1]),
			'-g', join('/', $qc2_dir, $contam_files[-1]),
			'-t dnaseq'
			);

		# run command
		print $log "Submitting job to create QC plots...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'create_dnaseq_coverage_plots',
			cmd		=> $qc_command,
			modules		=> [$r_version],
			mem		=> '2G',
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'create_dnaseq_coverage_plots',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;


		# add SNP directory to summary command
		$summarize_panel_command .= ' -v ' . join('/', $tool_dirs{'dnaseq_germline'}, 'cohort/VCF2MAF');

		# find DNA-Seq germline mutations
		my $germ_dir = join('/', $data_directory, 'DNASeq', 'Germline');
		unless(-e $germ_dir) { make_path($germ_dir); }

		opendir(FILES, join('/', $tool_dirs{'dnaseq_germline'}, 'CPSR'));
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		my @germ_files = grep { /mutations_for_cbioportal.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_germline'}, 'CPSR', $germ_files[-1]);
		$link_file = join('/', $germ_dir, $germ_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# add CPSR results to summary command
		$summarize_mutations_command .= " -g $base_file";

		# find DNA-Seq MSI estimates
		my $msi_dir = join('/', $data_directory, 'DNASeq', 'MSI');
		unless(-e $msi_dir) { make_path($msi_dir); }

		opendir(FILES, $tool_dirs{'dnaseq_msi'});
		@found_files = grep { /panelMSI_estimates.tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_msi'}, $found_files[-1]);
		$link_file = join('/', $msi_dir, $found_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# add MSI data to summary command
		$summarize_mutations_command .= " --msi $base_file";

		# find DNA-Seq somatic copy-number calls
		my $mops_dir = join('/', $data_directory, 'DNASeq', 'CNAs');
		unless(-e $mops_dir) { make_path($mops_dir); }

		opendir(FILES, $tool_dirs{'dnaseq_mops'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		my @cna_files = grep { /panelCN.mops_output.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_mops'}, $cna_files[-1]);
		$link_file = join('/', $mops_dir, $cna_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# add CN data to summary command
		$summarize_panel_command .= " -c $base_file";

		# find DNA-Seq SV calls
		my $sv_dir = join('/', $data_directory, 'DNASeq', 'SVs');
		unless(-e $sv_dir) { make_path($sv_dir); }

		opendir(FILES, $tool_dirs{'dnaseq_mavis'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		my @sv_files = grep { /sv_data_for_cbioportal.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_mavis'}, $sv_files[-1]);
		$link_file = join('/', $sv_dir, $sv_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# initiate command to find ensemble variants
		$ensemble_command = join(' ',
			"Rscript $cwd/../format_ensemble_mutations.R",
			'-p', $tool_data->{project_name},
			'-o', join('/', $data_directory, 'DNASeq', 'ENSEMBLE'),
			'-t targeted'
			);

		# find DNA-Seq somatic mutations (mutect2)
		my $mut_dir = join('/', $data_directory, 'DNASeq', 'Mutect2');
		unless(-e $mut_dir) { make_path($mut_dir); }

		opendir(FILES, $tool_dirs{'dnaseq_mutect2'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		my @mut_files = grep { /mutations_for_cbioportal.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_mutect2'}, $mut_files[-1]);
		$link_file = join('/', $mut_dir, $mut_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		$ensemble_command .= ' --mutect2 ' . $base_file;

		# find DNA-Seq somatic mutations (pindel)
		$mut_dir = join('/', $data_directory, 'DNASeq', 'Pindel');
		unless(-e $mut_dir) { make_path($mut_dir); }

		opendir(FILES, $tool_dirs{'dnaseq_pindel'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		@mut_files = grep { /mutations_for_cbioportal.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_pindel'}, $mut_files[-1]);
		$link_file = join('/', $mut_dir, $mut_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		$ensemble_command .= ' --pindel ' . $base_file;

		# find DNA-Seq somatic mutations (somaticsniper)
		$mut_dir = join('/', $data_directory, 'DNASeq', 'SomaticSniper');
		unless(-e $mut_dir) { make_path($mut_dir); }

		opendir(FILES, $tool_dirs{'dnaseq_somaticsniper'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		@mut_files = grep { /mutations_for_cbioportal.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_somaticsniper'}, $mut_files[-1]);
		$link_file = join('/', $mut_dir, $mut_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		$ensemble_command .= ' --somaticsniper ' . $base_file;

		# find DNA-Seq somatic mutations (strelka)
		$mut_dir = join('/', $data_directory, 'DNASeq', 'Strelka');
		unless(-e $mut_dir) { make_path($mut_dir); }

		opendir(FILES, $tool_dirs{'dnaseq_strelka'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		@mut_files = grep { /mutations_for_cbioportal.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_strelka'}, $mut_files[-1]);
		$link_file = join('/', $mut_dir, $mut_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		$ensemble_command .= ' --strelka ' . $base_file;

		# find DNA-Seq somatic mutations (vardict)
		$mut_dir = join('/', $data_directory, 'DNASeq', 'VarDict');
		unless(-e $mut_dir) { make_path($mut_dir); }

		opendir(FILES, $tool_dirs{'dnaseq_vardict'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		@mut_files = grep { /mutations_for_cbioportal.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'dnaseq_vardict'}, $mut_files[-1]);
		$link_file = join('/', $mut_dir, $mut_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		$ensemble_command .= ' --vardict ' . $base_file;
		}

	# find the required input files for EM-Seq
	if ($run_em) {

		# find EM-Seq QC data
		my $qc_dir = join('/', $data_directory, 'EMSeq', 'QC');
		unless(-e $qc_dir) { make_path($qc_dir); }

		opendir(FILES, $tool_dirs{'emseq_qc'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;

		my @cov_files = grep { /HSMetrics.tsv/ } @found_files;

		closedir(FILES);

		$base_file = join('/', $tool_dirs{'emseq_qc'}, $cov_files[-1]);
		$link_file = join('/', $qc_dir, $cov_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# create some QC plots
		my $qc_command = "Rscript $cwd/plot_qc_metrics.R";
		$qc_command .= " " . join(' ',
			'-p', $tool_data->{project_name},
			'-o', $qc_out_directory,
			'-r', $reports_directory,
			'-s', $args{em_data_config},
			'-c', join('/', $qc_dir, $cov_files[-1]),
			'-t emseq'
			);

		# run command
		print $log "Submitting job to create QC plots...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'create_emseq_coverage_plots',
			cmd		=> $qc_command,
			modules		=> [$r_version],
			mem		=> '2G',
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'create_emseq_coverage_plots',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;


		# find EM-Seq methylation data
		my $methyl_dir = join('/', $data_directory, 'EMSeq', 'Methylation');
		unless(-e $methyl_dir) { make_path($methyl_dir); }

		opendir(FILES, $tool_dirs{'emseq_methyldackel'});
		my @rdata_files = grep { /matrix.RData/ } readdir(FILES);
		@rdata_files = sort @rdata_files;
		closedir(FILES);

		opendir(FILES, $tool_dirs{'emseq_methyldackel'});
		@found_files = grep { /tsv/ } readdir(FILES);
		@found_files = sort @found_files;
		closedir(FILES);

		my @md_files = grep { /methylation_per_target_region.tsv/ } @found_files;

		$base_file = join('/', $tool_dirs{'emseq_methyldackel'}, $md_files[-1]);
		$link_file = join('/', $methyl_dir, $md_files[-1]);

		if ( -l $link_file) { unlink $link_file; }
		symlink($base_file, $link_file);

		# add methylation file to summary command
		$summarize_panel_command .= ' -m ' . join('/', $tool_dirs{'emseq_methyldackel'}, $rdata_files[-1]);

		# create some methylation plots
		my $em_command = "Rscript $cwd/plot_methylation.R";
		$em_command .= " " . join(' ',
			'-p', $tool_data->{project_name},
			'-o', $em_out_directory,
#			'-i', $tool_dirs{'swgs_ichor'},
			'-r', $reports_directory,
			'-s', $args{em_data_config},
			'-m', join('/', $tool_dirs{'emseq_methyldackel'}, $md_files[-1])
			);

		# run command
		print $log "Submitting job to create methylation plots...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'create_emseq_methylation_plots',
			cmd		=> $em_command,
			modules		=> [$r_version],
			mem		=> '2G',
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'create_emseq_methylation_plots',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;
		}

	# run analyses/plotting scripts
	if ($run_dna) {

		# find ENSEMBLE mutations
		$ensemble_maf = join('/',
			$data_directory,
			'DNASeq',
			'ENSEMBLE',
			'ensemble_mutation_data.tsv'
			);

		# add MAF to summary command
		$summarize_mutations_command .= " -m $ensemble_maf";

		# run command
		print $log "Submitting job to collect somatic variants...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'collect_somatic_variant_calls',
			cmd		=> $ensemble_command,
			modules		=> [$r_version],
			max_time	=> $max_time,
			mem		=> '8G',
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$ensemble_run_id = submit_job(
			jobname		=> 'collect_somatic_variant_calls',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $ensemble_run_id;

		}

	if ($run_dna && $run_em) {

		# add summary data to LOH command
		$summarize_mutations_command .= " -s " . join('/', $summary_out_directory, 'summarized_panel_data.RData');

		# run command
		print $log "Submitting job to summarize combined outputs...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'summarize_panel_outputs',
			cmd		=> $summarize_panel_command,
			modules		=> [$r_version],
			max_time	=> $max_time,
			mem		=> '8G',
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'summarize_panel_outputs',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;
		}


	# annotate/extract somatic mutations
	print $log "Submitting job to annotate/filter somatic variants...\n";
	$run_script = write_script(
		log_dir		=> $log_directory,
		name		=> 'summarize_mutation_data',
		cmd		=> $summarize_mutations_command,
		modules		=> [$r_version],
		dependencies	=> join(':', $ensemble_run_id, $run_id),
		max_time	=> $max_time,
		mem		=> '4G',
		hpc_driver	=> $args{hpc_driver},
		extra_args	=> [$hpc_group]
		);

	$run_id = submit_job(
		jobname		=> 'annotate_somatic_variants',
		shell_command	=> $run_script,
		hpc_driver	=> $args{hpc_driver},
		dry_run		=> $args{dry_run},
		log_file	=> $log
		);

	push @job_ids, $run_id;

	my $cohort_jobs = join(':', @job_ids);

	# generate report for each sample
	# we'll assume the dnaseq data config contains all tumour samples
	my $smp_data = LoadFile($args{dna_data_config});

	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient";

		# find tumour ids
		my @tumour_ids = keys %{$smp_data->{$patient}->{'tumour'}};

		next if (scalar(@tumour_ids) == 0);

		# create some directories
		my $patient_directory = join('/', $reports_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		# initiate command
		my $report_command = "cd $patient_directory\n\n";
		$report_command .= join(' ',
			"perl $cwd/generate_report.pl",
			'-o', $patient_directory
			);

		# run pdflatex (twice, to ensure proper page numbering)
		$report_command .= "\n\n" . join("\n",
			"pdflatex multiMMR_patient_report.tex",
			"pdflatex multiMMR_patient_report.tex"
			);

		$report_command .= "\n\n" . join("\n",
			'if [ 0 == $? ]; then',
			'  echo Generated and compiled LaTeX file with pdflatex',
			'else',
			'  echo Error running pdflatex: $?',
			'fi'
			);
	
		# run command
		print $log "\n >> Submitting job to generate report file...\n";
		$run_script = write_script(
			log_dir		=> $log_directory,
			name		=> 'create_report_' . $patient,
			cmd		=> $report_command,
			dependencies	=> $cohort_jobs,
			modules		=> [$perl,'texlive'],
			mem		=> '256M',
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'create_report_' . $patient,
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		push @job_ids, $run_id;
		}

	# if this is not a dry run OR there are jobs to assess (run or resumed with jobs submitted) then
	# collect job metrics (exit status, mem, run time)
	unless ( ($args{dry_run}) || (scalar(@job_ids) == 0) ) {

		# collect job stats
		my $collect_metrics = collect_job_stats(
			job_ids		=> join(',', @job_ids),
			outfile		=> $outfile,
			hpc_driver	=> $args{hpc_driver}
			);

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'output_job_metrics_' . $run_count,
			cmd	=> $collect_metrics,
			dependencies	=> join(':', @job_ids),
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
my ($wgs_tool_config, $wgs_data_config);
my ($dna_tool_config, $dna_data_config, $em_tool_config, $em_data_config);
my ($dry_run, $help, $no_wait, $output_directory);
my $hpc_driver = 'slurm';

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	'w|tool_wgs=s'	=> \$wgs_tool_config,
	't|tool_dna=s'	=> \$dna_tool_config,
	'f|tool_em=s'	=> \$em_tool_config,
	's|data_wgs=s'	=> \$wgs_data_config,
	'd|data_dna=s'	=> \$dna_data_config,
	'e|data_em=s'	=> \$em_data_config,
	'o|output_directory=s'	=> \$output_directory,
	'c|cluster=s'	=> \$hpc_driver,
	'dry-run'	=> \$dry_run,
	'no-wait'	=> \$no_wait
	);

if ($help) {
	my $help_msg = join("\n",
		"Options:",
		"\t--help|-h\tPrint this help message",
		"\t--tool_wgs|-w\t<string> tool config for sWGS (yaml format)",
		"\t--data_wgs|-s\t<string> data config for sWGS (yaml format)",
		"\t--tool_dna|-t\t<string> tool config for DNA-Seq (yaml format)",
		"\t--data_dna|-d\t<string> data config for DNA-Seq (yaml format)",
		"\t--tool_em|-f\t<string> tool config for EM-Seq (yaml format)",
		"\t--data_em|-e\t<string> data config for EM-Seq (yaml format)",
		"\t--output_directory|-o\t<string> path to output directory",
		"\t--cluster|-c\t<string> cluster scheduler (default: slurm)",
		"\t--dry-run\t<boolean> should jobs be submitted? (default: false)",
		"\t--no-wait\t<boolean> should we exit after job submission (true) or wait until all jobs have completed (false)? (default: false)"
		);

	print "$help_msg\n";
	exit;
	}

# do some quick error checks to confirm valid arguments	
if ( (!defined($wgs_tool_config)) | (!defined($wgs_data_config))) { 
	print("One or more config files for sWGS are missing; will skip this portion.\n");
	}

if ( (!defined($dna_tool_config)) | (!defined($dna_data_config))) { 
	print("One or more config files for DNA-Seq are missing; will skip this portion.\n");
	}

if ( (!defined($em_tool_config)) | (!defined($em_data_config))) { 
	print("One or more config files for EM-Seq are missing; will skip this portion.\n");
	}


main(
	wgs_tool_config		=> $wgs_tool_config,
	wgs_data_config		=> $wgs_data_config,
	dna_tool_config		=> $dna_tool_config,
	dna_data_config		=> $dna_data_config,
	em_tool_config		=> $em_tool_config,
	em_data_config		=> $em_data_config,
	output_directory	=> $output_directory,
	hpc_driver		=> $hpc_driver,
	dry_run			=> $dry_run,
	no_wait			=> $no_wait
	);

