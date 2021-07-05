#!/usr/bin/env perl
### novobreak.pl ###################################################################################
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
use List::Util 'first';
use IO::Handle;

my $cwd = dirname(__FILE__);
require "$cwd/utilities.pl";

# define some global variables
our ($reference, $bwa_ref, $pon, $intervals_bed) = undef;

####################################################################################################
# version       author		comment
# 1.0		sprokopec       script to run NovoBreak

### USAGE ##########################################################################################
# novobreak.pl -t tool_config.yaml -d data_config.yaml -o /path/to/output/dir -c slurm --remove --dry_run
#
# where:
#	-t (tool.yaml) contains tool versions and parameters, reference information, etc.
#	-d (data.yaml) contains sample information (YAML file containing paths to BWA-aligned,
#		GATK-processed BAMs, generated by gatk.pl)
#	-o (/path/to/output/dir) indicates tool-specific output directory
#	-c indicates hpc driver (ie, slurm)
#	--remove indicates that intermediates will be removed
#	--dry_run indicates that this is a dry run

### DEFINE SUBROUTINES #############################################################################
# format command to split reference by chromosome pair
sub split_ref_command {
	my %args = (
		out_dir		=> undef,
		reference	=> undef,
		intervals	=> undef,
		@_
		);

	my $split_command = "cd $args{out_dir}\n\n";

	$split_command .= join("\n",
		"while IFS=' ' read -r CHRS || [ -n " . '"${CHRS}" ]; do',
		'  echo "Prepping fasta for $CHRS";',
		'  CHR_PAIR=$(echo $CHRS | ' . "sed s/' '/_/g );",
		'  if [ -s $CHR_PAIR.fa ]; then',
		'    echo ">> $CHR_PAIR.fa already exists";',
		'  else',
		"    samtools faidx $args{reference} " . '$CHRS > $CHR_PAIR.fa;',
		'    samtools faidx $CHR_PAIR.fa;',
		'    echo ">> $CHR_PAIR.fa complete!";',
		'  fi',
		"done < $args{intervals}"
		);

	$split_command .= "\n\necho 'Split REFERENCE completed successfully.' > split_reference.COMPLETE";

	return($split_command);
	}

# format command to run NovoBreak (parallel)
sub get_novobreak_command {
	my %args = (
		tumour		=> undef,
		normal		=> undef,
		reference	=> undef,
		output_stem	=> undef,
		@_
		);

	my $novobreak_command = join(' ',
		'novoBreak',
		'-i', $args{tumour},
		'-c', $args{normal},
		'-r', $args{reference},
		'-o', $args{output_stem}
		);

	return($novobreak_command);
	}

# format command to run novobreak (WXS)
sub novobreak_for_wxs_command {
	my %args = (
		tumour_id	=> undef,
		tumour_bam	=> undef,
		normal_bam	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $nb_command = "cd $args{tmp_dir}\n";

	$nb_command .= join("\n",
		"if [ -s $args{tumour_id}" . '_nb.out.md5]; then',
		"  echo Intermediate output file: $args{tumour_id}" . '_nb.out.md5 already exists.',
		"else",
		"  echo 'Running novoBreak...'"
		);

	$nb_command .= "\n  " . get_novobreak_command(
		tumour		=> $args{tumour_bam},
		normal		=> $args{normal_bam},
		reference	=> $reference,
		output_stem	=> $args{tumour_id} . '_nb.out'
		);

	$nb_command .= "\n  md5sum $args{tumour_id}\_nb.out > $args{tumour_id}\_nb.out.md5";

	$nb_command .="\nfi";

	return($nb_command);
	}

# format command to split bam by chromosome
sub novobreak_by_split_bam_command {
	my %args = (
		tumour_id	=> undef,
		normal_id	=> undef,
		tumour_bam	=> undef,
		normal_bam	=> undef,
		intervals	=> undef,
		ref_dir		=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $split_command = "cd $args{tmp_dir}\n" . join("\n",
		'CHRS=$(sed -n "$SLURM_ARRAY_TASK_ID"p ' . $args{intervals} . ');',
		"CHR_PAIR=\$(echo \$CHRS | sed s/' '/_/g);",
		"\n" . 'if [[ ! -e $CHR_PAIR ]]; then mkdir $CHR_PAIR; fi',
		'cd $CHR_PAIR',
		"\nif [ -s $args{tumour_id}" . '_${CHR_PAIR}.tsv.md5]; then',
		"  echo Final output file $args{tumour_id}" . '_${CHR_PAIR}.tsv.md5 already exists.',
		'  exit;',
		'fi',
		"\nif [ -s $args{tumour_id}" . '_${CHR_PAIR}.bam.bai ]; then',
		"  echo $args{tumour_id}" . '_${CHR_PAIR}.bam.bai already exists.',
		'else',
		"  samtools view -b $args{tumour_bam} " . '$CHRS' . " > $args{tumour_id}" . '_${CHR_PAIR}.bam;',
		"  samtools index $args{tumour_id}" . '_${CHR_PAIR}.bam;',
		'fi',
		"\nif [ -s $args{normal_id}" . '_${CHR_PAIR}.bam.bai ]; then',
		"  echo $args{normal_id}" . '_${CHR_PAIR}.bam.bai already exists.',
		'else',
		"  samtools view -b $args{normal_bam} " . '$CHRS' . " > $args{normal_id}" . '_${CHR_PAIR}.bam;',
		"  samtools index $args{normal_id}" . '_${CHR_PAIR}.bam;',
		'fi'
		);

	my $nb_command = $split_command . "\n\n" . get_novobreak_command(
		tumour		=> $args{tumour_id} . '_${CHR_PAIR}.bam',
		normal		=> $args{normal_id} . '_${CHR_PAIR}.bam',
		reference	=> join('/', $args{ref_dir}, '$CHR_PAIR.fa'),
		output_stem	=> $args{tumour_id} . '_${CHR_PAIR}.tsv'
		);

	$nb_command .= "\n\nmd5sum $args{tumour_id}\_\${CHR_PAIR}.tsv > $args{tumour_id}\_\${CHR_PAIR}.tsv.md5";

	$nb_command .="\n\n" . join("\n",
		"if [ -s $args{tumour_id}" . '_${CHR_PAIR}.tsv.md5]; then',
		"  rm $args{tumour_id}" . '_${CHR_PAIR}.bam*',
		"  rm $args{normal_id}" . '_${CHR_PAIR}.bam*',
		'else',
		"  echo " . "Final output: $args{tumour_id}" . '_${CHR_PAIR}.tsv.md5 is missing; something went wrong',
		'fi'
		);

	return($nb_command);
	}

# format command to check output of split
sub confirm_split_command {
	my %args = (
		tmp_dir		=> undef,
		intervals	=> undef,
		@_
		);

	my $confirm_command = "cd $args{tmp_dir}\n\n";

	$confirm_command .= join("\n",
		"while IFS=' ' read -r CHRS || [ -n " . '"${CHRS}" ]; do',
		'  echo "Checking for NovoBreak output for $CHRS";',
		'  CHR_PAIR=$(echo $CHRS | ' . "sed s/' '/_/g );",
		'  if [ -s $CHR_PAIR/*_${CHR_PAIR}.tsv ] && [ -s $CHR_PAIR/*_${CHR_PAIR}.tsv.md5 ]; then',
		'    echo ">> NovoBreak output for $CHR_PAIR already exists";',
		'  else',
		'    echo ">> NovoBreak output for $CHR_PAIR missing or incomplete";',
		'    exit 1;',
		'  fi',
		"done < $args{intervals}"
		);

	$confirm_command .= "\n\necho 'All NovoBreak runs completed successfully.' > novobreak.COMPLETE";

	return($confirm_command);
	}

# format command to prepare ssake 
sub get_prep_ssake_command {
	my %args = (
		output_file	=> undef,
		n_cpus		=> 1,
		@_
		);

	my $nb_output = $args{output_file};

	my $part1_command = join("\n",
		'  samtools bam2fq -1 read1.fq -2 read2.fq somaticreads.bam;',
		'  group_bp_reads.pl ' . $nb_output . ' read1.fq read2.fq > bp_reads.tsv;'
		);

	$part1_command .= "\n" . join("\n",
		'  cls=$(tail -1 bp_reads.tsv | cut -f1);',
		'  rec=$(echo $cls/' . $args{n_cpus} . ' | bc);',
		'  rec=$((rec+1));',
		'  awk -v rec=$rec ' . "'{print" . ' > int($1/rec)".txt"' . "}' bp_reads.tsv;",
		'  for file in *.txt; do',
		'    run_ssake.pl $file > /dev/null &',
		'  done',
		'  wait',
		"\n  awk 'length(\$0)>1' *.ssake.asm.out > ssake.fa;",
		'  bwa mem -t ' . $args{n_cpus} . " -M $bwa_ref ssake.fa > ssake.sam;",
		'  md5sum ssake.sam > ssake.sam.md5',
		'fi'
		);

	return($part1_command);
	}

# format command to infer breakpoints
sub get_infer_breakpoints_command {
	my %args = (
		output_file	=> undef,
		tumour_bam	=> undef,
		normal_bam	=> undef,
		n_cpus		=> 1,
		@_
		);

	my $nb_output = $args{output_file};

	my $part2_command = join("\n",
		'  infer_sv.pl ssake.sam > ssake.vcf',
		"  grep -v '^#' ssake.vcf | sed 's/|/\t/g' | sed 's/read//' |  awk '{ if (!x[\$1\$2]) { y[\$1\$2]=\$14; x[\$1\$2]=\$0 } else { if (\$14 > y[\$1\$2]) { y[\$1\$2]=\$14; x[\$1\$2]=\$0 }}} END { for (i in x) { print x[i]}}' | sort -k1,1 -k2,2n | perl -ne 'if (/TRA/) { print } elsif (/SVLEN=(\\d+)/) { if (\$1 > 100) { print \$_ }} elsif (/SVLEN=-(\\d+)/) { if (\$1 > 100 ) { print }}' > ssake.pass.vcf",
		"\n" . '  num=$(wc -l ssake.pass.vcf | cut -f1 -d' . "' ');",
		'  rec=$(echo $num/' . $args{n_cpus} . ' | bc);',
		'  rec=$((rec+1));',
		'  split -l $rec ssake.pass.vcf',
		'  for file in x??; do',
		'    infer_bp_v4.pl $file '. "$args{tumour_bam} $args{normal_bam}" . ' > $file.sp.vcf & ',
		'  done',
		'  wait',
		"\n" . "  grep '^#' ssake.vcf > header.txt;",
		'  filter_sv_icgc.pl split/*.sp.vcf | cat header.txt - > ' . $nb_output,
		"  md5sum $nb_output > $nb_output.md5",
		'fi'
		);

	return($part2_command);
	}

# format command to process novoBreak output (WXS)
sub get_process_novobreak_command {
	my %args = (
		tumour_id	=> undef,
		tumour_bam	=> undef,
		normal_bam	=> undef,
		tmp_dir		=> undef,
		n_cpus		=> 1,
		@_
		);

	my $nb_output = $args{tumour_id} . '_nb.out';

	my $part1_command = "cd $args{tmp_dir}\n";

	# part 1
	$part1_command .= "\n" . join("\n",
		'if [ -s ssake.sam ]; then',
		"  echo 'Intermediate output file: ssake.sam already exists';",
		'else',
		"  echo 'Extracting and aligning reads...';"
		);

	$part1_command .= "\n" . get_prep_ssake_command(
		output_file	=> $nb_output,
		n_cpus		=> $args{n_cpus}
		);

	# part 2
	my $final_output = join('/', $args{tmp_dir}, '..', $args{tumour_id} . '_novoBreak.pass.flt.vcf');

	my $part2_command .= "\n\n" . join("\n",
		"if [ -s $final_output.md5 ]; then",
		"  echo 'Final output file " . $args{tumour_id} . "_novoBreak.pass.flt.vcf already exists';",
		'else',
		"  echo 'Inferring variants and breakpoints...';"
		);

	$part2_command .= "\n" . get_infer_breakpoints_command(
		output_file	=> $final_output,
		tumour_bam	=> $args{tumour_bam},
		normal_bam	=> $args{normal_bam},
		n_cpus		=> $args{n_cpus}
		);

	my $final_command = $part1_command . "\n" . $part2_command;

	return($final_command);
	}

# format command to process novoBreak output (WGS)
sub get_process_novobreak_command_part1 {
	my %args = (
		tumour_id	=> undef,
		intervals	=> undef,
		tmp_dir		=> undef,
		n_cpus		=> 1,
		@_
		);

	my $split_command = join("\n",
		'CHRS=$(sed -n "$SLURM_ARRAY_TASK_ID"p ' . $args{intervals} . ');',
		"CHR_PAIR=\$(echo \$CHRS | sed s/' '/_/g);",
		"\n" . 'echo Processing output for: $CHR_PAIR',
		"cd $args{tmp_dir}/" . '$CHR_PAIR',
		"\n" . 'if [ -s ssake.sam.md5 ]; then',
		'  echo Final output file: ssake.sam already exists;',
		'else'
		);

	my $nb_output = $args{tumour_id} . '_${CHR_PAIR}.tsv';

	my $part1_command = get_prep_ssake_command(
		output_file	=> $nb_output,
		n_cpus		=> $args{n_cpus}
		);

	my $final_command = $split_command . "\n" . $part1_command;

	return($final_command);
	}

# format command to check output of ssake
sub confirm_ssake_command {
	my %args = (
		tmp_dir		=> undef,
		intervals	=> undef,
		@_
		);

	my $confirm_command = "cd $args{tmp_dir}\n\n";

	$confirm_command .= join("\n",
		"while IFS=' ' read -r CHRS || [ -n " . '"${CHRS}" ]; do',
		'  echo "Checking for NovoBreak output for $CHRS";',
		'  CHR_PAIR=$(echo $CHRS | ' . "sed s/' '/_/g );",
		'  if [ -s $CHR_PAIR/ssake.sam.md5 ]; then',
		'    echo ">> NovoBreak ssake output for $CHR_PAIR already exists";',
		'  else',
		'    echo ">> NovoBreak ssake output for $CHR_PAIR missing or incomplete";',
		'    exit 1;',
		'  fi',
		"done < $args{intervals}"
		);

	$confirm_command .= "\n\necho 'All SSAKE runs completed successfully.' > ssake.COMPLETE";

	return($confirm_command);
	}

# format command to process novoBreak output
sub get_process_novobreak_command_part2 {
	my %args = (
		tumour_id	=> undef,
		tumour_bam	=> undef,
		normal_bam	=> undef,
		intervals	=> undef,
		tmp_dir		=> undef,
		n_cpus		=> 1,
		@_
		);

	my $nb_output = $args{tumour_id} . '_${CHR_PAIR}_novoBreak.pass.flt.vcf';

	my $split_command = join("\n",
		'CHRS=$(sed -n "$SLURM_ARRAY_TASK_ID"p ' . $args{intervals} . ');',
		"CHR_PAIR=\$(echo \$CHRS | sed s/' '/_/g);",
		"\n" . 'echo Processing output for: $CHR_PAIR',
		"cd $args{tmp_dir}/" . '$CHR_PAIR',
		"\n" . "if [ -s $nb_output.md5 ]; then",
		"  echo Final output file: $nb_output already exists;",
		'else'
		);

	my $part2_command = get_infer_breakpoints_command(
		output_file	=> $nb_output,
		tumour_bam	=> $args{tumour_bam},
		normal_bam	=> $args{normal_bam},
		n_cpus		=> $args{n_cpus}
		);

	my $final_command = $split_command . "\n" . $part2_command;

	return($final_command);
	}

# format command to check output of infer breakpoints
sub confirm_breakpoint_command {
	my %args = (
		tmp_dir		=> undef,
		intervals	=> undef,
		@_
		);

	my $confirm_command = "cd $args{tmp_dir}\n\n";

	$confirm_command .= join("\n",
		"while IFS=' ' read -r CHRS || [ -n " . '"${CHRS}" ]; do',
		'  echo "Checking for NovoBreak output for $CHRS";',
		'  CHR_PAIR=$(echo $CHRS | ' . "sed s/' '/_/g );",
		'  if [ -s $CHR_PAIR/*_${CHR_PAIR}_novoBreak.pass.flt.vcf.md5 ]; then',
		'    echo ">> NovoBreak SV calls for $CHR_PAIR already exists";',
		'  else',
		'    echo ">> NovoBreak SV calls for $CHR_PAIR are missing or incomplete";',
		'    exit 1;',
		'  fi',
		"done < $args{intervals}"
		);

	$confirm_command .= "\n\necho 'All NovoBreak SV calling steps completed successfully.' > breakpoint.COMPLETE";

	return($confirm_command);
	}

# format command to check output of infer breakpoints
sub merge_and_filter_command {
	my %args = (
		tmp_dir		=> undef,
		output_file	=> undef,
		@_
		);

	my $merge_command = "cd $args{tmp_dir}\n" . join(' ',
		'vcf-concat */*${CHR_PAIR}_novoBreak.pass.flt.vcf',
		"| vcf-sort -c | perl $cwd/filter_novobreak_variants.pl",
		'>', $args{output_file}
		);
	
	return($merge_command);
	}

### MAIN ###########################################################################################
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

	# load tool config
	my $tool_data_orig = LoadFile($tool_config);
	my $tool_data = error_checking(tool_data => $tool_data_orig, pipeline => 'novobreak');
	my $date = strftime "%F", localtime;

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_NovoBreak_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_NovoBreak_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running NovoBreak SV calling pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n    Reference used: $tool_data->{reference}";

	$reference = $tool_data->{reference};
	$bwa_ref = $tool_data->{bwa}->{reference};

	print $log "\n    Output directory: $output_directory";
	print $log "\n  Sample config used: $data_config";
	print $log "\n---";

	# set tools and versions
	my $novobreak	= 'novoBreak/' . $tool_data->{novobreak_version};
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};
	my $vcftools	= 'vcftools/' . $tool_data->{vcftools_version};
	my $bwa		= 'bwa/' . $tool_data->{bwa_version};
	my $r_version	= 'R/' . $tool_data->{r_version};

	# get user-specified tool parameters
	my $parameters = $tool_data->{novobreak}->{parameters};

	# set chromosome list
	my $string;
	if ('exome' eq $tool_data->{seq_type}) {
		$string = 'exome';
		} elsif (defined($tool_data->{novobreak}->{chromosomes})) {
		$string = $tool_data->{novobreak}->{chromosomes};
		} elsif ( ('hg38' eq $tool_data->{ref_type}) || ('hg19' eq $tool_data->{ref_type})) {
		$string = 'chr' . join(',chr', 1..22) . ',chrX,chrY';
		} elsif ( ('GRCh37' eq $tool_data->{ref_type}) || ('GRCh37' eq $tool_data->{ref_type})) {
		$string = join(',', 1..22) . ',X,Y';
		}

	my @chroms = split(',', $string);

	# set some binaries
	my $is_multi_slurm = ((scalar(@chroms) > 1) && ('slurm' eq $args{hpc_driver}));
	my $is_wgs = ('wgs' eq $tool_data->{seq_type});

	### RUN ###########################################################################################
	my ($run_script, $run_id, $link, $prep_ref_run_id, $chr_file, $cleanup_run_id, $ref_dir);
	my (@all_jobs);

	# indicate chromosome pairs to check
	my $pair_count = 0;

	if ($is_wgs) {
		$chr_file = join('/', $output_directory, 'chromosome_list.txt');
		open (my $chr_list, '>', $chr_file) or die "Could not open $chr_file for writing.";

		foreach my $i (0..$#chroms) {
			my $chr1 = $chroms[$i];
			foreach my $j ($i..$#chroms) {
				my $chr2 = $chroms[$j];
				next if ( $chr1 eq $chr2 ); 
				print $chr_list "$chr1 $chr2\n";
				$pair_count++;
				}
			}

		close ($chr_list);

		# Prepare fasta file (chr pairs)
		$ref_dir = join('/', $output_directory, 'reference');
		unless(-e $ref_dir) { make_path($ref_dir); }

		my $split_ref_cmd = split_ref_command(
			reference	=> $reference,
			intervals	=> $chr_file,
			out_dir		=> $ref_dir
			);

		my $split_confirmation = join('/', $ref_dir, 'split_reference.COMPLETE');
		if ('Y' eq missing_file($split_confirmation)) {

			# record command (in log directory) and then run job
			print $log "Submitting job for SplitREF...\n";

			$run_script = write_script(
				log_dir	=> $log_directory,
				name	=> 'run_split_reference_by_chromosome_pair',
				cmd	=> $split_ref_cmd,
				modules	=> [$samtools],
				max_time	=> '04:00:00',
				mem		=> '1G',
				hpc_driver	=> $args{hpc_driver}
				);

			$prep_ref_run_id = submit_job(
				jobname		=> 'run_split_reference_by_chromosome_pair',
				shell_command	=> $run_script,
				hpc_driver	=> $args{hpc_driver},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			push @all_jobs, $prep_ref_run_id;
			} else {
			print $log "Skipping SplitREF because this has already been completed!\n";
			}
		}

	# get sample data
	my $smp_data = LoadFile($data_config);

	# process each sample in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\nInitiating process for PATIENT: $patient\n";

		# find bams
		my @normal_ids = keys %{$smp_data->{$patient}->{'normal'}};
		my @tumour_ids = keys %{$smp_data->{$patient}->{'tumour'}};

		if (scalar(@normal_ids) == 0) {
			print $log "\n>> No normal BAM provided, skipping patient.\n";
			next;
			}

		# create some directories
		my $patient_directory = join('/', $output_directory, $patient);
		unless(-e $patient_directory) { make_path($patient_directory); }

		my $link_directory = join('/', $patient_directory, 'bam_links');
		unless(-e $link_directory) { make_path($link_directory); }

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

		# for T/N pair
		foreach my $sample (@tumour_ids) {

			print $log "  SAMPLE: $sample\n\n";

			my $sample_directory = join('/', $patient_directory, $sample);
			unless(-e $sample_directory) { make_path($sample_directory); }

			my $tmp_directory = join('/', $sample_directory, 'TEMP');
			unless(-e $tmp_directory) { make_path($tmp_directory); }

			# indicate this should be removed at the end
			my $cleanup_cmd = "rm -rf $tmp_directory";

			$run_id = '';
			my $output_to_check;

			# run novoBreak
			if ($is_wgs && $is_multi_slurm) {

				# run novoBreak using split BAMs
				$output_to_check = join('/', $tmp_directory, 'novobreak.COMPLETE');

				if ('Y' eq missing_file($output_to_check)) {

					my $split_novo_command = novobreak_by_split_bam_command(
						tumour_id	=> $sample,
						normal_id	=> $normal_ids[0],
						tumour_bam	=> $smp_data->{$patient}->{tumour}->{$sample},
						normal_bam	=> $smp_data->{$patient}->{normal}->{$normal_ids[0]},
						intervals	=> $chr_file,
						ref_dir		=> $ref_dir,
						tmp_dir		=> $tmp_directory
						);

					# record command (in log directory) and then run job
					print $log "Submitting job for NovoBreak...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_split_novobreak_' . $sample,
						cmd	=> $split_novo_command,
						modules	=> [$samtools, $novobreak, 'perl'],
						dependencies	=> $prep_ref_run_id,
						max_time	=> $parameters->{novobreak}->{time},
						mem		=> $parameters->{novobreak}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> '--array=1-' . $pair_count . '%5'
						);

					$run_id = submit_job(
						jobname		=> 'run_split_novobreak_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;

					# check novoBreak output
					my $confirm_nb_cmd = confirm_split_command(
						tmp_dir		=> $tmp_directory,
						intervals	=> $chr_file,
						);

					# record command (in log directory) and then run job
					print $log "Submitting job for Check NovoBreak...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_check_novobreak_' . $sample,
						cmd	=> $confirm_nb_cmd,
						dependencies	=> $run_id,
						mem		=> '128M',
						hpc_driver	=> $args{hpc_driver}
						);

					$run_id = submit_job(
						jobname		=> 'run_check_novobreak_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log "Skipping NovoBreak because step is already complete!\n";
					}
				} elsif ($is_wgs && ! $is_multi_slurm) {
				# wip
				} else {
				# run novoBreak using full BAMs
				$output_to_check = join('/', $tmp_directory, $sample . '_nb.out.md5');

				if ('Y' eq missing_file($output_to_check)) {

					my $full_novo_command = novobreak_for_wxs_command(
						tumour_id	=> $sample,
						tumour_bam	=> $smp_data->{$patient}->{tumour}->{$sample},
						normal_bam	=> $smp_data->{$patient}->{normal}->{$normal_ids[0]},
						tmp_dir		=> $tmp_directory
						);

					# record command (in log directory) and then run job
					print $log "Submitting job for NovoBreak...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_novobreak_' . $sample,
						cmd	=> $full_novo_command,
						modules	=> [$samtools, $novobreak, 'perl'],
						max_time	=> $parameters->{novobreak}->{time},
						mem		=> $parameters->{novobreak}->{mem},
						hpc_driver	=> $args{hpc_driver}
						);

					$run_id = submit_job(
						jobname		=> 'run_novobreak_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log "Skipping NovoBreak because step is already complete!\n";
					}
				}

			# Run post-process steps
			if ($is_wgs && $is_multi_slurm) {
			
				$output_to_check = join('/', $tmp_directory, 'ssake.COMPLETE');

				if ('Y' eq missing_file($output_to_check)) {
				
					my $nb_part1_command = get_process_novobreak_command_part1(
						tumour_id	=> $sample,
						intervals	=> $chr_file,
						tmp_dir		=> $tmp_directory,
						n_cpus		=> $parameters->{process}->{n_cpu}
						);

					# record command (in log directory) and then run job
					print $log "Submitting job for NovoBreak SSAKE...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_novobreak_ssake_' . $sample,
						cmd	=> $nb_part1_command,
						modules	=> [$samtools, $novobreak, $bwa, 'perl'],
						dependencies	=> $run_id,
						max_time	=> $parameters->{ssake}->{time},
						mem		=> $parameters->{ssake}->{mem},
						cpus_per_task	=> $parameters->{ssake}->{n_cpu},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> '--array=1-' . $pair_count . '%5'
						);

					$run_id = submit_job(
						jobname		=> 'run_novobreak_ssake_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;

					# check process novoBreak (ssake) output
					my $confirm_ssake_cmd = confirm_ssake_command(
						tmp_dir		=> $tmp_directory,
						intervals	=> $chr_file,
						);

					# record command (in log directory) and then run job
					print $log "Submitting job for Check NovoBreak SSAKE...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_check_novobreak_ssake_' . $sample,
						cmd	=> $confirm_ssake_cmd,
						dependencies	=> $run_id,
						mem		=> '128M',
						hpc_driver	=> $args{hpc_driver}
						);

					$run_id = submit_job(
						jobname		=> 'run_check_novobreak_ssake' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log "Skipping Part 1 of post-process as this is already complete.\n";
					}

				# Run post-process part 2
				$output_to_check = join('/', $tmp_directory, 'breakpoint.COMPLETE');

				if ('Y' eq missing_file($output_to_check)) {

					my $nb_part2_command = get_process_novobreak_command_part2(
						tumour_id	=> $sample,
						tumour_bam	=> $smp_data->{$patient}->{tumour}->{$sample},
						normal_bam	=> $smp_data->{$patient}->{normal}->{$normal_ids[0]},
						intervals	=> $chr_file,
						tmp_dir		=> $tmp_directory,
						n_cpus		=> $parameters->{process}->{n_cpu} 
						);

					# record command (in log directory) and then run job
					print $log "Submitting job for NovoBreak process part 2...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_novobreak_find_breakpoints_' . $sample,
						cmd	=> $nb_part2_command,
						modules	=> [$samtools, $novobreak, $bwa, 'perl'],
						dependencies	=> $run_id,
						max_time	=> $parameters->{process}->{time},
						mem		=> $parameters->{process}->{mem},
						cpus_per_task	=> $parameters->{process}->{n_cpu},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> '--array=1-' . $pair_count . '%5'
						);

					$run_id = submit_job(
						jobname		=> 'run_novobreak_find_breakpoints_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;

					# check process novoBreak (infer breakpoints) output
					my $confirm_bp_cmd = confirm_breakpoint_command(
						tmp_dir		=> $tmp_directory,
						intervals	=> $chr_file,
						);

					# record command (in log directory) and then run job
					print $log "Submitting job for Check NovoBreak BREAKPOINTS...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_check_novobreak_breakpoints_' . $sample,
						cmd	=> $confirm_bp_cmd,
						dependencies	=> $run_id,
						mem		=> '128M',
						hpc_driver	=> $args{hpc_driver}
						);

					$run_id = submit_job(
						jobname		=> 'run_check_novobreak_breakpoints' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log "Skipping Part 2 of post-process as this is already complete.\n";
					}

				} elsif ($is_wgs && ! $is_multi_slurm) {
				# wip
				} else {
				# run post-process steps on full novoBreak output
				$output_to_check = join('/', $sample_directory, $sample . '_novoBreak.pass.flt.vcf');

				if ('Y' eq missing_file($output_to_check . '.md5')) {
				
					my $nb_process_command = get_process_novobreak_command(
						tumour_id	=> $sample,
						tumour_bam	=> $smp_data->{$patient}->{tumour}->{$sample},
						normal_bam	=> $smp_data->{$patient}->{normal}->{$normal_ids[0]},
						tmp_dir		=> $tmp_directory,
						n_cpus		=> $parameters->{process}->{n_cpu}
						);

					# record command (in log directory) and then run job
					print $log "Submitting job for NovoBreak post-process...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_novobreak_postprocess_' . $sample,
						cmd	=> $nb_process_command,
						modules	=> [$samtools, $novobreak, $bwa, 'perl'],
						dependencies	=> $run_id,
						max_time	=> $parameters->{process}->{time},
						mem		=> $parameters->{process}->{mem},
						cpus_per_task	=> $parameters->{process}->{n_cpu},
						hpc_driver	=> $args{hpc_driver}
						);

					$run_id = submit_job(
						jobname		=> 'run_novobreak_postprocess_' . $sample,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					push @final_outputs, $output_to_check;
					} else {
					print $log "Skipping post-process step as this is already complete!\n";
					}
				}

			# merge output (if split run)
			my $merged_output = join('/', $sample_directory, $sample . '_novoBreak.pass.flt.vcf');

			if ($is_wgs && ('Y' eq missing_file($merged_output . '.md5'))) {
				
				my $merge_command = merge_and_filter_command(
					tmp_dir		=> $tmp_directory,
					output_file	=> $merged_output
					);

				# record command (in log directory) and then run job
				print $log "Submitting job for MERGE step...\n";

				$run_script = write_script(
					log_dir	=> $log_directory,
					name	=> 'run_merge_and_filter_' . $sample,
					cmd	=> $merge_command,
					modules	=> ['perl', $vcftools],
					dependencies	=> $run_id,
					max_time	=> $parameters->{merge}->{time},
					mem		=> $parameters->{merge}->{mem},
					hpc_driver	=> $args{hpc_driver}
					);

				$run_id = submit_job(
					jobname		=> 'run_merge_and_filter_' . $sample,
					shell_command	=> $run_script,
					hpc_driver	=> $args{hpc_driver},
					dry_run		=> $args{dry_run},
					log_file	=> $log
					);

				push @patient_jobs, $run_id;
				push @all_jobs, $run_id;
				push @final_outputs, $merged_output;
				} else {
				print $log "Skipping final merge step as this is already completed!\n";
				}

			# should intermediate files be removed
			# run per patient
			if ($args{del_intermediates}) {

				if (scalar(@patient_jobs) == 0) {
					`rm -rf $tmp_directory`;

					} else {

					print $log "Submitting job to clean up temporary/intermediate files...\n";

					# make sure final output exists before removing intermediate files!
					my @files_to_check;
					foreach my $tmp ( @final_outputs ) {
						push @files_to_check, $tmp . '.md5';
						}

					$cleanup_cmd = join("\n",
						"if [ -s " . join(" ] && [ -s ", @files_to_check) . " ]; then",
						"  $cleanup_cmd",
						"else",
						'  echo "One or more FINAL OUTPUT FILES is missing; not removing intermediates"',
						"fi"
						);

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_cleanup_' . $patient,
						cmd	=> $cleanup_cmd,
						dependencies	=> join(':', @patient_jobs),
						mem		=> '256M',
						hpc_driver	=> $args{hpc_driver},
						kill_on_error	=> 0
						);

					$cleanup_run_id = submit_job(
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

			unless($args{dry_run}) { sleep(60); }
			}
		}

	# collate results
	my $collect_output = join(' ',
		"Rscript $cwd/collect_nb_output.R",
		'-d', $output_directory,
		'-p', $tool_data->{project_name}
		);

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'combine_variant_calls',
		cmd	=> $collect_output,
		modules	=> [$r_version],
		dependencies	=> join(':', @all_jobs),
		mem		=> '4G',
		max_time	=> '24:00:00',
		hpc_driver	=> $args{hpc_driver}
		);

	$run_id = submit_job(
		jobname		=> 'combine_variant_calls',
		shell_command	=> $run_script,
		hpc_driver	=> $args{hpc_driver},
		dry_run		=> $args{dry_run},
		log_file	=> $log
		);

	# if this is not a dry run OR there are jobs to assess (run or resumed with jobs submitted) then
	# collect job metrics (exit status, mem, run time)
	unless ( ($args{dry_run}) || (scalar(@all_jobs) == 0) ) {

		# collect job stats
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
			hpc_driver	=> $args{hpc_driver},
			kill_on_error	=> 0
			);

		$run_id = submit_job(
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
					die("Final VARSCAN accounting job: $run_id finished with errors.");
					}
				}
			}
		}

	# finish up
	print $log "\nProgramming terminated successfully.\n\n";
	close $log;
	}

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($tool_config, $data_config, $output_directory);
my $hpc_driver = 'slurm';
my ($remove_junk, $dry_run, $help, $no_wait);
my $panel_of_normals = undef;

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	'd|data=s'	=> \$data_config,
	't|tool=s'	=> \$tool_config,
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
		"\t--data|-d\t<string> data config (yaml format)",
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

# do some quick error checks to confirm valid arguments	
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