#!/usr/bin/env perl
### utilities.pl ###################################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use File::Basename;
use File::Path qw(make_path);
use List::Util 'any';

### SHARED SUBROUTINES #############################################################################
# function to do standard error checking on input config files
sub error_checking {
	my %args = (
		tool_data	=> undef,
		data_type	=> 'dna',
		pipeline	=> 'main',
		@_
		);

	my $tool_data = $args{tool_data};
	my $pipeline  = $args{pipeline};
	my $data_type = $args{data_type};

	# check to see if tool-specific arguments are supplied and/or are correct

	# is ref_type either hg38 or hg19?
	# is this RNA?
	if ('rna' eq $data_type) {
		if ( ('GRCh38' ne $tool_data->{ref_type}) && ('hg38' ne $tool_data->{ref_type}) ) {
			die("RNA-Seq pipeline only configured for use with GRCh38/hg38.");
		}
	} else {
	# or DNA?
		if ( ('hg38' ne $tool_data->{ref_type}) && ('hg19' ne $tool_data->{ref_type}) && 
		('GRCh37' ne $tool_data->{ref_type}) && ('GRCh38' ne $tool_data->{ref_type})) {
			die("Unrecognized ref_type; must be one of hg19, GRCh37, hg38 or GRCh38.");
		}
	}


	my $is_ref_valid;

	# BWA (DNA only)
	if ('bwa' eq $pipeline) {

		# alignment method
		if ( ('emseq' eq $tool_data->{seq_type}) && ('bwameth' ne $tool_data->{bwa}->{aligner}) ) {
			die("EM-Seq is currently only compatible with BWA-METH alignment!");
		}

		if ( ('emseq' ne $tool_data->{seq_type}) && ('bwamem' ne $tool_data->{bwa}->{aligner}) ) {
			die("This pipeline is currently only compatible with BWA-MEM!");
		}

		# type of flowcell (for markduplicates)
		unless (defined($tool_data->{flowcell_type})) {
			print "tool_config: Option flowcell must be either random or patterned, defaulting to random\n";
			$tool_data->{flowcell_type} = 'random';
		}

		# fill in any missing parameters
		if (
			('Y' ne $tool_data->{bwa}->{parameters}->{merge}->{mark_dup}) && 
			('N' ne $tool_data->{bwa}->{parameters}->{merge}->{mark_dup})
			) {
			print "bwa_config: Option mark_dup must be either Y or N, defaulting to N\n";
			$tool_data->{bwa}->{parameters}->{merge}->{mark_dup} = 'N';
		}

		if (!defined($tool_data->{bwa}->{parameters}->{bwa}->{n_cpus}->{tumour})) {
			$tool_data->{bwa}->{parameters}->{bwa}->{n_cpus}->{tumour} = 4;
			print "bwa_config: Number of threads not defined for step bwa:tumour; setting to 4\n";
		}

		if (!defined($tool_data->{bwa}->{parameters}->{bwa}->{n_cpus}->{normal})) {
			$tool_data->{bwa}->{parameters}->{bwa}->{n_cpus}->{normal} = 4;
			print "bwa_config: Number of threads not defined for step bwa:normal; setting to 4\n";
		}

		if (!defined($tool_data->{reference}))  {
			die("Must supply path to reference genome!");
		}

		if (!defined($tool_data->{bwa}->{reference}))  {
			die("Must supply path to bwa index for indicated reference genome!");
		}

		$is_ref_valid = validate_ref(
			reference	=> $tool_data->{reference},
			pipeline	=> 'bwa',
			exts		=> [qw(.fa .fa.fai .dict)]
		);

		if ('emseq' eq $tool_data->{seq_type}) {
			$is_ref_valid = validate_ref(
				reference	=> $tool_data->{bwa}->{reference},
				pipeline	=> 'bwa',
				exts		=> [qw(.fa.bwameth.c2t .fa.bwameth.c2t.amb .fa.bwameth.c2t.ann .fa.bwameth.c2t.bwt .fa.bwameth.c2t.pac .fa.bwameth.c2t.sa)]
			);	
		} else {
			$is_ref_valid = validate_ref(
				reference	=> $tool_data->{bwa}->{reference},
				pipeline	=> 'bwa',
				exts		=> [qw(.fa .fa.amb .fa.ann .fa.bwt .fa.pac .fa.sa)]
			);
		}
	}

	# ensure compatibility with past versions
	if ( (defined($tool_data->{intervals_bed})) && !(defined($tool_data->{targets_bed})) ) {
		$tool_data->{targets_bed} = $tool_data->{intervals_bed};
		}

	# QC pipeline (gatk/picard functions)
	if ('qc' eq $pipeline) {

		if (!defined($tool_data->{reference})) { die("Must supply path to reference genome!"); }
		$is_ref_valid = validate_ref(
			reference	=> $tool_data->{reference},
			pipeline	=> 'gatk',
			exts		=> [qw(.fa .dict .fa.fai)]
		);

		if ('dna' eq $data_type) {
			if (!defined($tool_data->{targets_bed}) & 
				('wgs' ne $tool_data->{seq_type}) & ('emseq' ne $tool_data->{seq_type})) {
				print "gatk_config: WARNING: no target intervals provided.\n";
				print ">>If this is exome data, target regions are recommended!\n";
			}

			if ('emseq' eq $tool_data->{seq_type}) {

				if ( (defined($tool_data->{targets_bed})) & (!defined($tool_data->{baits_bed})) ) {
					$tool_data->{baits_bed} = $tool_data->{targets_bed};
					print "Only one of target/bait intervals files provided; using for both settings.\n";
				} elsif ( (!defined($tool_data->{targets_bed})) & (defined($tool_data->{baits_bed})) ) {
					$tool_data->{targets_bed} = $tool_data->{baits_bed};
					print "Only one of target/bait intervals files provided; using for both settings.\n";
				} elsif ( (!defined($tool_data->{targets_bed})) & (!defined($tool_data->{baits_bed})) ) {
					print "No target regions and/or baits bed file(s) provided; will run on whole-genome.\n";
				}
			}
		}
	}

	# GATK (indel realignment/recalibration + HaplotypeCaller + MuTect + MuTect2 + GATK-CNV)
	if ('gatk' eq $pipeline) {

		if (!defined($tool_data->{reference})) { die("Must supply path to reference genome!"); }
		$is_ref_valid = validate_ref(
			reference	=> $tool_data->{reference},
			pipeline	=> 'gatk',
			exts		=> [qw(.fa .dict .fa.fai)]
		);

		if ('dna' eq $data_type) {
			if (!defined($tool_data->{targets_bed}) & ('wgs' ne $tool_data->{seq_type})) {
				print "gatk_config: WARNING: no target intervals provided.\n";
				print ">>If this is exome data, target regions are recommended!\n";
			}
		}
	}

	# Novobreak
	if (('novobreak' eq $pipeline)) {
		if (!defined($tool_data->{reference})) { die("Must supply path to reference genome!"); }

		if ('wgs' eq $tool_data->{seq_type}) {
			print "NovoBreak is intended for WXS but this is WGS; will try to run but be prepared for failures!";
		}

		$is_ref_valid = validate_ref(
			reference	=> $tool_data->{bwa}->{reference},
			pipeline	=> 'bwa',
			exts		=> [qw(.fa .fa.amb .fa.ann .fa.bwt .fa.pac .fa.sa)]
		);
	}

	# SViCT
	if ('svict' eq $pipeline) {
		if (!defined($tool_data->{reference})) { die("Must supply path to reference genome!"); }
		unless ( ('ctdna' eq $tool_data->{sample_type}) && ('targeted' eq $tool_data->{seq_type}) ) {
			die("SViCT should only be run on targeted-panel ctDNA sequencing.");
		}
	}

	# ASCAT
	if ('ascat' eq $pipeline) {
		unless ( 'wgs' eq $tool_data->{seq_type} ) {
			die("ASCAT should only be run on WGS cohorts.");
		}
	}

	# panelCN.mops
	if ('mops' eq $pipeline) {
		unless ('targeted' eq $tool_data->{seq_type}) {
			die("panelCN.mops should only be run on targeted-panel data.");
		}

		if (!defined($tool_data->{targets_bed})) {
			die("Must provide path to target regions bed file.");
		}
	}

	# Strelka, VarScan, SomaticSniper and Delly
	my @pipeline_list = qw(strelka varscan delly somaticsniper vardict pindel);
	if ( any { /$pipeline/ } @pipeline_list ) {

		if (!defined($tool_data->{reference})) { die("Must supply path to reference genome!"); }

		my $intervals;
		if ('dna' eq $data_type) {

			# are intervals provided for exome/targeted seq?
			# are they properly formatted?
			if ( ('strelka' eq $pipeline) &&
				( any { /$tool_data->{seq_type}/ } qw(exome targeted) )
				) {

				$intervals = $tool_data->{targets_bed};
				$intervals =~ s/\.bed/_padding100bp.bed.gz/;

				if (!defined($tool_data->{targets_bed})) {
					die("Must supply path to target intervals!");
				} elsif ('Y' eq missing_file($intervals)) {
					die("Padded, bgzipped file: $intervals is missing. Please run format_intervals_bed.pl to ensure padding is added and file is bgzipped and tabix indexed.");
				}
			}

			if ( ( any { /$pipeline/ } qw(varscan somaticsniper vardict pindel) ) &&
				( any { /$tool_data->{seq_type}/ } qw(exome targeted) )
				) {

				$intervals = $tool_data->{targets_bed};
				$intervals =~ s/\.bed/_padding100bp.bed/;

				if (!defined($tool_data->{targets_bed})) {
					print "WARNING: no target intervals provided.\n";
					print ">>If this is exome data, target regions are recommend!\n";
				} elsif ('Y' eq missing_file($intervals)) {
					die("Padded file: $intervals is missing. Please run format_intervals_bed.pl to ensure padding is added.");
				}
			}

			if ( ('vardict' eq $pipeline) && ('wgs' eq $tool_data->{seq_type}) ) {
				$intervals = $tool_data->{vardict}->{intervals};
				if ( (!defined($intervals)) || ('Y' eq missing_file($intervals)) ) {
					die("For VarDict on WGS, intervals in the form of chr:start-end MUST be provided.");
				}
			}

		} elsif (('rna' eq $data_type) && ('strelka' eq $pipeline)) {

			if ('rna' ne $tool_data->{seq_type}) {
				$tool_data->{seq_type} = 'rna';
			}
		}
	}

	# STAR (RNA only)
	if ('star' eq $pipeline) {
		if (!defined($tool_data->{star_reference_dir}))  {
			die("star_config: Must supply path to reference genome directory!");
		}

		# type of flowcell (for markduplicates)
		unless (defined($tool_data->{flowcell_type})) {
			print "tool_config: Option flowcell must be either random or patterned, defaulting to random\n";
			$tool_data->{flowcell_type} = 'random';
		}
	}

	# STAR-Fusion (RNA only)
	if ('star-fusion' eq $pipeline) {
		if (!defined($tool_data->{star_fusion_reference_dir}))  {
			die("star_fusion_config: Must supply path to reference genome directory!");
		}
	}

	# RSEM (RNA only)
	if ('rsem' eq $pipeline) {
		if (!defined($tool_data->{rsem_reference}))  { die("rsem_config: Must supply path/to/reference/stem !"); }
		$is_ref_valid = validate_ref(
			reference	=> $tool_data->{rsem_reference},
			pipeline	=> 'rsem',
			exts		=> [qw(.idx.fa .grp .transcripts.fa .seq .chrlist)]
		);

		my @strand_options = qw(none forward reverse);
		if (!defined($tool_data->{rsem}->{strandedness})) {
			print "rsem_config: No option provided for 'strandedness'; setting to default: none.\n";
			$tool_data->{rsem}->{strandedness} = 'none';
		}

		if ( !any { /$tool_data->{rsem}->{strandedness}/ } @strand_options ) {
			print "rsem_config: Unrecognized 'strandedness' option: must be one of none, forward or reverse! Setting to default: none.\n";
			$tool_data->{rsem}->{strandedness} = 'none';
		}
	}

	return($tool_data);
}

# function to ensure all necessary reference files are available
sub validate_ref {
	my %args = (
		reference	=> undef,
		pipeline	=> undef,
		exts		=> undef,
		@_
		);

	my $ref_file_base;

	if ( ('bwa' eq $args{pipeline}) || ('gatk' eq $args{pipeline}) || ('variant_call' eq $args{pipeline}) ) {
		$ref_file_base = ($args{reference} =~ s/\.fa$//r);
	} else {
		$ref_file_base = $args{reference};
	}

	foreach my $ext (@{$args{exts}}) {
		unless (-e $ref_file_base . $ext) { die("Missing reference file: " . $ref_file_base . $ext); }
	}

	# if nothing fails, return Y
	return('Y');
}

# function to check for missing/empty files
sub missing_file {
	my @bad = grep { ! -e $_ || -s $_ <= 1}
		map { ref($_) eq 'HASH' ? values %{$_} : $_ }
		@_;

	if (@bad) {
		return('Y');
	}

	return('N');
}

# save commands to shell script in log directory
sub write_script {
	my %args = (
		log_dir		=> undef,
		name		=> undef,
		cmd		=> undef,
		modules 	=> [],
		dependencies	=> undef,
		max_time	=> '01:00:00',
		mem		=> '1G',
		cpus_per_task	=> 1,
		hpc_driver	=> 'slurm',
		extra_args	=> [],
		kill_on_error	=> 1,
		@_
		);

	if (!defined($args{max_time})) { $args{max_time} = '01:00:00'; }
	if (!defined($args{mem})) { $args{mem} = '1G'; }
	if (!defined($args{cpus_per_task})) { $args{cpus_per_task} = 1; }

	my $cmd_log_dir = join('/', $args{log_dir}, $args{name});
	unless(-e $cmd_log_dir) { make_path($cmd_log_dir); }

	# make a directory for error/log output
	my $job_log_dir = join('/', $cmd_log_dir, $args{hpc_driver});
	unless(-e $job_log_dir) { make_path($job_log_dir); }

	my @modules_list;
	for (my $i=0; $i < scalar (@{$args{modules}}); $i++) {
		next if (!defined($args{modules}->[$i]));
		next if ('' eq $args{modules}->[$i]);
		$modules_list[$i] = "module load $args{modules}->[$i]";
		}

	my $modules_to_load = join(";\n", @modules_list);

	my $script = join('/', $cmd_log_dir, 'script.sh');
	open (my $fh_script, '>', $script) or Carp::croak("Cannot open file $script: $!");

	my $job_params;

	# add in SBATCH parameters
	if ('slurm' eq $args{hpc_driver}) {

		print $fh_script "#!/bin/bash\n";

		$job_params = "#SBATCH " . join("\n#SBATCH ",
			'--job-name="' . $args{name} . '"',
			'-D ' . $job_log_dir,
			'-t ' . $args{max_time},
			'--mem ' . $args{mem},
			'-c ' . $args{cpus_per_task}
			);

		my ($size, $unit);
		if ($args{mem} =~ m/(\d+)([A-Z])/) {
			$size = $1;
			$unit = $2;
			if (($size > 28) && ($size < 61) && ('G' eq $unit)) {
				$job_params .= "\n#SBATCH -p himem";
			} elsif (($size >= 61) && ($size < 184) && ('G' eq $unit)) {
				$job_params .= "\n#SBATCH -p veryhimem";
			}
		}

		my ($days, $hours);
		if ($args{max_time} =~ m/(\d+)(-)(\d+)/) {
			$days = $1;
			$hours = $3;
			if (
				( ($days == 5) && ($hours ne '00') ) ||
				( ($days > 5) )
				) {
				$job_params .= "\n#SBATCH -p long";
			}
		}

		for (my $i=0; $i < scalar (@{$args{extra_args}}); $i++) {
			if (defined($args{extra_args}->[$i])) {
				$job_params .= "\n#SBATCH " . $args{extra_args}->[$i];
			}
		}

		if (defined($args{dependencies})) {

			# change any , to :
			$args{dependencies} =~ s/,/:/g;

			if ($args{dependencies} =~ m/:/) {
				my @parts = split(/:/, $args{dependencies});
				my @depends = grep { $_ ne '' } @parts;

				if (scalar(@depends) > 1) { $args{dependencies} = join(':', @depends); }
				elsif (scalar(@depends) == 1) { $args{dependencies} = $depends[0]; }
				else { $args{dependencies} = ''; }
			}

			if ('' ne $args{dependencies}) {

				if ($args{name} =~ m/job_metrics/) {
					$job_params .= "\n#SBATCH --dependency=afterany:" . $args{dependencies};
				} else {
					$job_params .= "\n#SBATCH --dependency=afterok:" . $args{dependencies};
					$job_params .= "\n#SBATCH --kill-on-invalid-dep=yes";
				}
			}
		}

	# add in PBS parameters
	} elsif ('pbs' eq $args{hpc_driver}) {

		$job_params = "#PBS " . join("\n#PBS ",
			'-S /bin/bash',
			'-N ' . $args{name},
			'-e ' . $job_log_dir . '/errorlog',
			'-k e',
			'-l nodes=1:ppn=' . $args{cpus_per_task}
			);

		my ($size, $unit);
		if ($args{mem} =~ m/(\d+)([A-Z])/) {
			$size = $1;
			$unit = $2;
			if ('G' eq $unit) { $unit = 'gb'; }
		#	if (($size > 28) && ($size < 61) && ('G' eq $unit)) {
		#		$job_params .= "\n#PBS -q himem";
		#	} elsif (($size >= 61) && ($size < 184) && ('G' eq $unit)) {
		#		$job_params .= "\n#PBS -q veryhimem";
		#	}

		$job_params .= "\n" . "#PBS -l mem=" . $size . $unit;
		}

		my ($days, $hours);
		if ($args{max_time} =~ m/(\d+)(-)(\d+)/) {
			$days = $1;
			$hours = $3;
		#	if (
		#		( ($days == 5) && ($hours ne '00') ) ||
		#		( ($days > 5) )
		#		) {
		#		$job_params .= "\n#PBS -q long";
		#	}

		$job_params .= "\n" . "#PBS -l walltime=" . $args{max_time};
		}

		for (my $i=0; $i < scalar (@{$args{extra_args}}); $i++) {
			if (defined($args{extra_args}->[$i])) {
				$job_params .= "\n#PBS " . $args{extra_args}->[$i];
			}
		}

		if (defined($args{dependencies})) {

			# change any , to :
			$args{dependencies} =~ s/,/:/g;

			if ($args{dependencies} =~ m/:/) {
				my @parts = split(/:/, $args{dependencies});
				my @depends = grep { $_ ne '' } @parts;

				if (scalar(@depends) > 1) { $args{dependencies} = join(':', @depends); }
				elsif (scalar(@depends) == 1) { $args{dependencies} = $depends[0]; }
				else { $args{dependencies} = ''; }
			}

			if ('' ne $args{dependencies}) {

				if ($args{name} =~ m/job_metrics/) {
					$job_params .= "\n#PBS -W depend=afterany:" . $args{dependencies};
				} else {
					$job_params .= "\n#PBS -W depend=afterok:" . $args{dependencies};
				}
			}
		}
	}

	print $fh_script $job_params . "\n\n";

	print $fh_script 'echo HOSTNAME: $HOSTNAME' . "\n";
	if ('slurm' eq $args{hpc_driver}) {
		print $fh_script 'echo SLURM NODE: $SLURMD_NODENAME' . "\n\n";
		}

	if ($args{kill_on_error}) {
		print $fh_script "set -e\n\n";
		}

	if (scalar(@modules_list) > 0) {
		print $fh_script $modules_to_load . "\n\n";
		}

	print $fh_script $args{cmd};
	close($fh_script);

	return($script);
	}

# Many java-based programs (ie, GATK) can fail without triggering a SLURM error
# To address this, we will add a final check to be run once these tasks finish,
# to verify that they completed successfully
sub check_java_output {
	my %args = (
		extra_cmd	=> undef,
		@_
		);

	my $java_check .= "\n" . join("\n",
		'if [ $? == 0 ]; then',
		'  echo "java task completed successfully"'
		);

	if (defined($args{extra_cmd})) {
		$java_check .= "\n  $args{extra_cmd}";
		}

	$java_check .= "\n" . join("\n",
		'else',
		'  echo "java task failed"',
		'  exit 1',
		'fi'
		);

	return($java_check);
	}

# execute commands/submit jobs
sub submit_job {
	my %args = (
		jobname		=> undef,
		shell_command	=> undef,
		dry_run		=> undef,
		hpc_driver	=> undef,
		log_file	=> undef,
		@_
		);

	my $log = $args{log_file};

	my $job_command;

	if ('slurm' eq $args{hpc_driver}) {
		$job_command = "sbatch " . $args{shell_command};
	} elsif ('pbs' eq $args{hpc_driver}) {
		$job_command = "qsub " . $args{shell_command};
	} else {
		die("Unrecognized HPC driver: currently only compatible with slurm or pbs");
	}

	print $log "\nCOMMAND IS: " . $job_command . "\n";

	my $job_id = $args{jobname};
	unless ($args{dry_run}) {
		$job_id = `$job_command`;
		chop $job_id;
		if ($job_id =~ m/Submitted batch job ([0-9]+)/) {
			$job_id = $1;
			} else {
			die("Failed to submit job with command $job_command");
			}

		print $log "Job number $job_id submitted.\n";
		} else {
		print $log "Job not submitted.\n";
		}

	# print 1 extra blank line in the log to separate steps
	print $log "\n";

	return($job_id);
	}

# command to extract job status / metrics
sub collect_job_stats {
	my %args = (
		job_ids		=> undef,
		outfile		=> undef,
		hpc_driver	=> undef,
		@_
		);

	my $sacct_command;

	if ('slurm' eq $args{hpc_driver}) {

		$sacct_command = join(' ',
			'sacct -P --delimiter=","',
			'--format="User,JobID,JobName,Start,End,AllocCPUS,CPUTime,MaxRSS,State,ExitCode"',
			'-j', $args{job_ids},
			'>', $args{outfile} . ';',
			'sed -i "s/,/\t/g"', $args{outfile}
			);

		$sacct_command .= "\n\n" . join(' ',
			'STATUS_COUNT=$(cut -f9', $args{outfile},
			"| awk '", '$1 != "COMPLETED" { print $0 }', "' | wc -l)",
			);

		$sacct_command .= "\n\n" . join("\n",
			'if (( $STATUS_COUNT > 1 )); then',
			'  exit 1;',
			'fi'
			);

		} elsif ('pbs' eq $args{hpc_driver}) {

		$sacct_command = join(' ',
			'qstat -f', $args{job_ids},
			'>', $args{outfile} . ';',
			);

		} else {
		$sacct_command = 'Job IDs to check: ' . $args{job_ids};
		}

	return($sacct_command);
	}

# command to check job status
sub check_final_status {
	my %args = (
		job_id	=> undef,
		@_
		);

	my $complete = 0;
	my $timeouts = 0;

	while (!$complete && $timeouts < 20 ) {
		sleep(30);
		my $status = `sacct --format='State' -j $args{job_id}`;

		# if we run into a server connection error (happens rarely with sacct)
		# increment timeouts (if we continue to repeatedly timeout, we will exit)
		if ($status =~ m/Connection timed out/) {
			$timeouts++;
			}
		else {
			my @status_lines = split("\n", $status);
			$status = $status_lines[2];

			# if final job has finished successfully:
			if ($status =~ m/COMPLETED/s) { $complete = 1; }
			# if the job is still pending or running, try again in a bit
			# but also reset timeouts, because we only care about consecutive timeouts
			elsif ($status =~ m/PENDING|RUNNING/) {
				$timeouts = 0;
				}
			# if none of the above, we will exit with an error
			else {
				die("Final accounting job: $args{job_id} finished with errors.");
				}
			}
		}

	print "Final accouting job: $args{job_id} finished successfully.";
	}

# format command to generate PON
sub generate_pon {
	my %args = (
		input		=> undef,
		output		=> undef,
		reference	=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		minN		=> 2,
		out_type	=> 'full',
		@_
		);

	my $pon_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $gatk_dir/GenomeAnalysisTK.jar -T CombineVariants',
		'-R', $args{reference},
		$args{input},
		'-o', $args{output},
		'--filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED',
		'--genotypemergeoption UNSORTED --filteredAreUncalled'
		);

	if ('trimmed' eq $args{out_type}) {
		$pon_command .= ' -minimalVCF -suppressCommandLineHeader --excludeNonVariants --sites_only';
		$pon_command .= " -minN $args{minN}";
		}

	return($pon_command);
	}

# find files recursively
sub _get_files {
	my ($dirs, $exten) = @_;

	my @files;
	my $want = sub {
		-e && /\Q$exten\E$/ && push @files, $File::Find::name
		};

	find($want, $dirs);
	return(@files);
	}

# format command to convert annotated VCF to MAF
sub get_vcf2maf_command {
	my $args = {
		input		=> undef,
		tumour_id	=> undef,
		tumour_vcf_id	=> undef,
		normal_id	=> undef,
		reference	=> undef,
		ref_type	=> undef,
		output		=> undef,
		tmp_dir		=> undef,
		parameters	=> undef,
		@_
		};

	my $ref_type;
	if ('hg19' eq $args->{ref_type}) {
		$ref_type = 'GRCh37';
		} elsif ('hg38' eq $args->{ref_type}) {
		$ref_type = 'GRCh38';
		} elsif ( ('GRCh38' eq $args->{ref_type}) || ('GRCh37' eq $args->{ref_type}) ) {
		$ref_type = $args->{ref_type};
		}

	my $maf_command = 'vcf2maf.pl';
	if (defined($args->{parameters}->{vcf2maf_path})) {
		$maf_command = "perl $args->{parameters}->{vcf2maf_path}"; 
		}

	$maf_command .= ' ' . join(' ',
		'--species homo_sapiens',
		'--ncbi-build', $ref_type,
		'--ref-fasta', $args->{reference},
		'--input-vcf', $args->{input},
		'--output-maf', $args->{output},
		'--tumor-id', $args->{tumour_id},
		'--vep-path', $args->{parameters}->{vep_path},
		'--vep-data', $args->{parameters}->{vep_data},
		'--vep-forks', $args->{parameters}->{n_cpus},
		'--filter-vcf', $args->{parameters}->{filter_vcf},
		'--buffer-size', $args->{parameters}->{buffer_size},
		'--tmp-dir', $args->{tmp_dir}
		);

	if (defined($args->{normal_id})) {
		$maf_command .= " --normal-id $args->{normal_id}";

		if ($args->{input} =~ m/Strelka|VarScan|MuTect2|SomaticSniper/) {
			$maf_command .= " --vcf-tumor-id TUMOR --vcf-normal-id NORMAL";
			}
	} else {
		if ($args->{input} =~ m/VarScan/) {
			$maf_command .= " --vcf-tumor-id Sample1";
		}
		if ( (defined($args->{tumour_vcf_id})) && ($args->{input} =~ m/MuTect2|Strelka/)) {
			$maf_command .= " --vcf-tumor-id $args->{tumour_vcf_id}";
		}
	}

	return($maf_command);
	}

1;
