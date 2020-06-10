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

	# check to see if common arguments are supplied and/or are correct
	if ('main' eq $pipeline) {

		# for del_intermediates, if not Y then default to N
		if (
			(!defined($tool_data->{del_intermediates})) ||
			('Y' ne $tool_data->{del_intermediates}) && ('N' ne $tool_data->{del_intermediates})
			) {
			print "Option del_intermediates is neither Y or N, defaulting to N\n";
			$tool_data->{del_intermediates} = 'N';
		}

		# check if this is a dry run or not
		if ((!defined($tool_data->{dry_run})) || ('Y' ne $tool_data->{dry_run})) {
			$tool_data->{dry_run} = 'N';
		}

		# check for compatible HPC driver; if not found, change dry_run to Y
		$tool_data->{HPC_driver} = lc $tool_data->{HPC_driver};
		my @compatible_drivers = qw(slurm);
		if ((!any { /$tool_data->{HPC_driver}/ } @compatible_drivers ) && ('N' eq $tool_data->{dry_run})) {
			print "Unrecognized HPC driver requested: setting dry_run to Y, jobs will not be submitted but commands will be written to file.\n";
			$tool_data->{dry_run} = 'Y';
		}

	# check to see if tool-specific arguments are supplied and/or are correct
	} else {

		# is ref_type either hg38 or hg19?
		if ( ('hg38' ne $tool_data->{ref_type}) && ('hg19' ne $tool_data->{ref_type}) ) {
			die("Unrecognized ref_type; must be one of hg19 or hg38.");
		}

		my $is_ref_valid;

		# BWA (DNA only)
		if ('bwa' eq $pipeline) {
			if ('bwamem' ne $tool_data->{aligner}) {
				die("This pipeline is currently only compatible with BWA-MEM!");
			}

			if (('Y' ne $tool_data->{mark_dup}) & ('N' ne $tool_data->{mark_dup})) {
				print "bwa_config: Option mark_dup must be either Y or N, defaulting to N\n";
				$tool_data->{mark_dup} = 'N';
			}

			if (!defined($tool_data->{reference}))  { die("Must supply path to reference genome!"); }
			$is_ref_valid = validate_ref(
				reference	=> $tool_data->{reference},
				pipeline	=> $pipeline,
				exts		=> [qw(.fa .fa.amb .fa.ann .fa.bwt .fa.fai .fa.pac .fa.sa)]
			);
		}

		# GATK (indel realignment/recalibration + HaplotypeCaller)
		if ('gatk' eq $pipeline) {

			if (!defined($tool_data->{reference})) { die("Must supply path to reference genome!"); }
			$is_ref_valid = validate_ref(
				reference	=> $tool_data->{reference},
				pipeline	=> $pipeline,
				exts		=> [qw(.fa .dict .fa.fai)]
			);

			if ('dna' eq $data_type) {
				if (!defined($tool_data->{intervals_bed})) {
					print "gatk_config: WARNING: no target intervals provided.\n";
					print ">>If this is exome data, please provide the target regions!\n";
				}
			}
		}

		# STAR (RNA only)
		if ('star' eq $pipeline) {
			if (!defined($tool_data->{reference_dir}))  {
				die("star_config: Must supply path to reference genome directory!");
			}

			if (('Y' ne $tool_data->{mark_dup}) && ('N' ne $tool_data->{mark_dup})) {
				print "star_config: Option mark_dup is neither Y or N, defaulting to N\n";
				$tool_data->{mark_dup} = 'N';
			}
		}

		# STAR-Fusion (RNA only)
		if ('star-fusion' eq $pipeline) {
			if (!defined($tool_data->{reference_dir}))  {
				die("star_fusion_config: Must supply path to reference genome directory!");
			}
		}

		# RSEM (RNA only)
		if ('rsem' eq $pipeline) {
			if (!defined($tool_data->{reference}))  { die("rsem_config: Must supply path/to/reference/stem !"); }
			$is_ref_valid = validate_ref(
				reference	=> $tool_data->{reference},
				pipeline	=> 'rsem',
				exts		=> [qw(.idx.fa .grp .transcripts.fa .seq .chrlist)]
			);

			my @strand_options = qw(none forward reverse);
			if (!defined($tool_data->{strandedness})) {
				print "rsem_config: No option provided for 'strandedness'; setting to default: none.\n";
				$tool_data->{strandedness} = 'none';
			}

			if ( !any { /$tool_data->{strandedness}/ } @strand_options ) {
				print "rsem_config: Unrecognized 'strandedness' option: must be one of none, forward or reverse! Setting to default: none.\n";
				$tool_data->{strandedness} = 'none';
			}
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
		dependencies	=> '',
		max_time	=> '01:00:00',
		mem		=> '1G',
		cpus_per_task	=> 1,
		hpc_driver	=> 'slurm',
		extra_args	=> undef,
		@_
		);

	my $cmd_log_dir = join('/', $args{log_dir}, $args{name});
	unless(-e $cmd_log_dir) { make_path($cmd_log_dir); }

	# make a directory for error/log output
	my $job_log_dir = join('/', $cmd_log_dir, $args{hpc_driver});
	unless(-e $job_log_dir) { make_path($job_log_dir); }

	my @modules_list;
	for (my $i=0; $i < scalar (@{$args{modules}}); $i++) {
		unless('' eq $args{modules}->[$i]) {
			$modules_list[$i] = "module load $args{modules}->[$i]";
			}
		}

	my $modules_to_load = join(";\n", @modules_list);

	my $script = join('/', $cmd_log_dir, 'script.sh');
	open (my $fh_script, '>', $script) or Carp::croak("Cannot open file $script: $!");
	print $fh_script "#!/bin/bash\n";

	# add in SBATCH parameters
	if ('slurm' eq $args{hpc_driver}) {

		my $sbatch_params = "#SBATCH " . join("\n#SBATCH ",
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
				$sbatch_params .= "\n#SBATCH -p himem";
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
				$sbatch_params .= "\n#SBATCH -p long";
			}
		}

		if (defined($args{extra_args})) {
			$sbatch_params .= "\n#SBATCH " . $args{extra_args};
			}

		if ('' ne $args{dependencies}) {

			if ($args{name} =~ m/job_metrics/) {
				$sbatch_params .= "\n#SBATCH --dependency=afterany:" . $args{dependencies};
				} else {
				$sbatch_params .= "\n#SBATCH --dependency=afterok:" . $args{dependencies};
				$sbatch_params .= "\n#SBATCH --kill-on-invalid-dep=yes";
				}
			}

		print $fh_script $sbatch_params . "\n\n";
		}

	print $fh_script $modules_to_load . "\n";
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
		'  echo "java task completed successfully!"'
		);

	if (defined($args{extra_cmd})) {
		$java_check .= "\n  $args{extra_cmd}";
		}

	$java_check .= "\n" . join("\n",
		'else',
		'  echo "java task failed!"',
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
		dry_run		=> 'N',
		hpc_driver	=> 'slurm',
		log_file	=> undef,
		@_
		);

	my $log = $args{log_file};

	my $job_command;
	unless('slurm' eq $args{hpc_driver}) {
		die("Unrecognized HPC driver: currently only compatible with slurm");
		}

	$job_command = "sbatch " . $args{shell_command};
	print $log "\nCOMMAND IS: " . $job_command . "\n";

	my $job_id = $args{jobname};
	if ('N' eq $args{dry_run}) {
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
		job_ids	=> undef,
		outfile	=> undef,
		@_
		);

	my $sacct_command = join(' ',
		'sacct -P --delimiter=","',
		'--format="User,JobID,Start,End,AllocCPUS,Elapsed,CPUTime,MaxRSS,ExitCode"',
		'-j', $args{job_ids},
		'>', $args{outfile} . ';',
		'sed -i "s/,/\t/g"', $args{outfile}
		);

	return($sacct_command);
	}

# format command to convert annotated VCF to MAF
sub get_vcf2maf_command {
	my %args = (
		input		=> undef,
		tumour_id	=> undef,
		normal_id	=> undef,
		reference	=> undef,
		ref_type	=> undef,
		output		=> undef,
		tmp_dir		=> undef,
		vcf2maf		=> undef,
		vep_path	=> undef,
		vep_data	=> undef,
		filter_vcf	=> undef,
		@_
		);

	my $ref_type;
	if ('hg19' eq $args{ref_type}) {
		$ref_type = 'GRCh37';
		} elsif ('hg38' eq $args{ref_type}) {
		$ref_type = 'GRCh38';
		}

	my $maf_command = join(' ',
		'perl', $args{vcf2maf},
		'--species homo_sapiens',
		'--ncbi-build', $ref_type,
		'--ref-fasta', $args{reference},
		'--input-vcf', $args{input},
		'--output-maf', $args{output},
		'--tumor-id', $args{tumour_id},
		'--vep-path', $args{vep_path},
		'--vep-data', $args{vep_data},
		'--vep-forks 1',
		'--filter-vcf', $args{filter_vcf},
		'--buffer-size 100',
		'--tmp-dir', $args{tmp_dir}
		);

	if (defined($args{normal_id})) {
		$maf_command .= " --normal-id $args{normal_id}";
		}

	return($maf_command);
	}

1;
