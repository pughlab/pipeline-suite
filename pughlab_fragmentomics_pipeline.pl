#!/usr/bin/env perl
### pughlab_fragmentomics_pipeline.pl ##############################################################
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
use File::Find;
use Data::Dumper;
use IO::Handle;

my $cwd = dirname(__FILE__);
require "$cwd/scripts/utilities.pl";

# define path to fragmentomics code location(s)
our ($fragmentomics_code_dir, $griffin_code_dir) = undef;

# define some global variables
our ($reference, $reference_peaks) = undef;
our ($delfi_healthy, $delfi_gaps, $delfi_vntrs, $delfi_tiles, $delfi_filters) = undef;
our ($tfbs_sites, $tcga_sites, $dhs_sites, $hk_sites) = undef;
our ($griffin_mappable, $griffin_mappable_name) = undef;

####################################################################################################
# version       author		comment
# 1.0		sprokopec	wrapper to run Derek's fragmentomics code

### DEFINE SUBROUTINES #############################################################################
# run downsample bam
sub create_downsample_command {
	my %args = (
		bam	=> undef,
		id	=> undef,
		n_reads	=> 50000000,
		outdir	=> undef,
		tmpdir	=> undef,
		@_
		);

	my $downsample_command = join("\n",
		"N_READS=" . $args{n_reads},
		"TOTAL_READS=\$(samtools view -c $args{bam})",
		"SCALE_FACTOR=\$(printf '%.4f\\n' \$(echo \"\$N_READS/\$TOTAL_READS\" | bc -l))",
		"echo 'downsampling to \$N_READS from \$TOTAL_READS (scale factor = \$SCALE_FACTOR)'",
		'if (( $(echo "$SCALE_FACTOR < 1" | bc -l) )); then'
		);

	$downsample_command .= "\n  " . join(' ',  
		'java -Xmx1g',
		'-Djava.io.tmpdir=' . $args{tmpdir},
		'-jar $picard_dir/picard.jar DownsampleSam',
		'I=' . $args{bam},
		'O=' . $args{outdir} . '/' . $args{id} . '_downsampled.bam',
		'P=$SCALE_FACTOR',
		'VALIDATION_STRINGENCY=LENIENT;'
		);

	$downsample_command .= "\nelse\n  " . join(' ',
		'ln -s',
		$args{bam},
		$args{outdir} . '/' . $args{id} . '_downsampled.bam;'
		);

	$downsample_command .= "\nfi";

	$downsample_command .= "\n" . join(' ',
		'samtools index',
		$args{outdir} . '/' . $args{id} . '_downsampled.bam;'
		);

	return($downsample_command);
	}

# run DELFI fragmentomics ratio
sub create_ratio_command {
	my %args = (
		bam	=> undef,
		id	=> undef,
		outdir	=> undef,
		@_
		);

	my $delfi_command = join(' ',
		"Rscript $fragmentomics_code_dir/ratio/runFrag.R",
		"--id", $args{id},
		"--bam", $args{bam},
		"--outdir", $args{outdir},
		"--libdir $fragmentomics_code_dir/ratio",
		"--filters $delfi_filters",
		"--gaps $delfi_gaps",
		"--VNTRs $delfi_vntrs",
		"--tiles $delfi_tiles",
		"--healthy $delfi_healthy"
		);

	return($delfi_command);
	}

# run VESSIES fragment score
sub create_fragmentscore_command {
	my %args = (
		id	=> undef,
		bam	=> undef,
		outdir	=> undef,
		@_
		);

	my $score_command = join(' ',
		"Rscript $fragmentomics_code_dir/fragment_score/scripts/04_patient_score.R",
		"--id", $args{id},
		"--bam", $args{bam},
		"--outdir", $args{outdir},
		"--ref $fragmentomics_code_dir/fragment_score/ref/vessies_reference_set.txt",
		"--libdir $fragmentomics_code_dir/fragment_score/scripts"
		);

	return($score_command);
	}

# format command to create de-duplicated BAM
sub get_deduplicate_command {
	my %args = (
		id		=> undef,
		bam		=> undef,
		outdir		=> undef,
		tmp_dir		=> undef,
		java_mem	=> undef,
		@_
		);

	my $output_stem = join('/', $args{outdir}, $args{id});

	# deduplicate BAM
	my $dedup_command = "# deduplicate bam";
	$dedup_command .= "\n" . join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $picard_dir/picard.jar MarkDuplicates',
		'I=' . $args{bam},
		'O=' . $output_stem . '_deduped.bam',
		'M=' . $output_stem . '_deduped.metrics',
		'REMOVE_DUPLICATES=true',
		'REMOVE_SEQUENCING_DUPLICATES=true'
		);

	# sort deduplicated BAM
	$dedup_command .= "\n\n# sort deduplicated bam";
	$dedup_command .= "\n" . join(' ',
		'samtools sort -n',
		$output_stem . '_deduped.bam',
		'-o', $output_stem . '_deduped_sorted.bam'
		);

	# convert to bedpe (all reads)
	$dedup_command .= "\n\n# convert deduplicated bam to bedpe";
	$dedup_command .= "\n" . join(' ',
		'samtools view -bf 0x2',
		$output_stem . '_deduped_sorted.bam',
		'| bedtools bamtobed -i stdin -bedpe',
		'>', $output_stem . '_all.bedpe'
		);

	# convert to bedpe (high quality reads)
	$dedup_command .= "\n\n# select high-quality reads and convert to bedpe";
	$dedup_command .= "\n" . join(' ',
		'samtools view -bf 0x2 -q 30',
		$output_stem . '_deduped_sorted.bam',
		'| bedtools bamtobed -i stdin -bedpe',
		'>', $output_stem . '_q30.bedpe'
		);

	# remove intermediates
	$dedup_command .= "\n\n# check for completeness and remove intermediates";
	$dedup_command .= "\n" . join("\n",
		'if [[ -f ' . $output_stem . '_q30.bedpe ]]; then',
		'  rm ' . $output_stem . '*deduped*',
		'fi'
		);

	return($dedup_command);
	}

# format command to create size-selected and de-duplicated BAM
sub get_sized_deduplicate_command {
	my %args = (
		id		=> undef,
		bam		=> undef,
		outdir		=> undef,
		size		=> 'genome',
		tmp_dir		=> undef,
		java_mem	=> undef,
		@_
		);

	my $output_stem = join('/', $args{outdir}, $args{id});

	my ($dedup_command, $new_bam);

	# select fragments of 167bp length (ie, 'normal' cfDNA)
	$dedup_command = "# extracting desired fragments: ($args{size})";
	if ('genome' eq $args{size}) {

		$new_bam = $output_stem . "_len167_sorted.bam";
		$output_stem = join('/', $args{outdir}, $args{id} . '_len167');

		$dedup_command .= "\n" . join(' ',
			"samtools view -h $args{bam}",
			"| awk 'substr(\$0,1,1)==\"@\" || (\$9==167) || (\$9==-167)'",
			"| samtools view -b",
			"> $output_stem.bam;"
			);

		# select fragments of 150bp length (ie, 'short' cfDNA)
		} elsif ('short' eq $args{size}) {

		$new_bam = $output_stem . "_len150_sorted.bam";
		$output_stem = join('/', $args{outdir}, $args{id} . '_len150');

		$dedup_command .= "\n" . join(' ',
			"samtools view -h $args{bam}",
			"| awk 'substr(\$0,1,1)==\"@\" || (\$9<=150 && \$9>=10) || (\$9>=-150 && \$9<=-10)'",
			"| samtools view -b",
			"> $output_stem.bam;"
			);
		}

	$dedup_command .= "\n\n# sort and index size-selected bam";
	$dedup_command .= "\nsamtools index $output_stem.bam";
	$dedup_command .= "\nsamtools sort $output_stem.bam -o $output_stem\_sorted.bam";
	$dedup_command .= "\nsamtools index $output_stem\_sorted.bam";	

	# deduplicate BAM
	$dedup_command .= "\n\n# deduplicate size-selected bam";
	$dedup_command .= "\n" . join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $picard_dir/picard.jar MarkDuplicates',
		'I=' . $new_bam,
		'O=' . $output_stem . '_deduped.bam',
		'M=' . $output_stem . '_deduped.metrics',
		'REMOVE_DUPLICATES=true',
		'REMOVE_SEQUENCING_DUPLICATES=true'
		);

	# sort deduplicated BAM
	$dedup_command .= "\n\n# sort deduplicated size-selected bam";
	$dedup_command .= "\n" . join(' ',
		'samtools sort -n',
		$output_stem . '_deduped.bam',
		'-o', $output_stem . '_deduped_sorted.bam'
		);

	# convert to bedpe (all reads)
	$dedup_command .= "\n\n# convert bam to bedpe";
	$dedup_command .= "\n" . join(' ',
		'samtools view -bf 0x2',
		$output_stem . '_deduped_sorted.bam',
		'| bedtools bamtobed -i stdin -bedpe',
		'>', $output_stem . '.bedpe'
		);

	# remove intermediates
	$dedup_command .= "\n\n# check for completeness and remove intermediates";
	$dedup_command .= "\n" . join("\n",
		"if [[ -f $output_stem.bedpe ]]; then",
		'  rm ' . $output_stem . '*bam*',
		'  rm ' . $output_stem . '*metrics',
		'fi'
		);

	return($dedup_command);
	}

# run nucleosome positioning
sub create_nucleosome_position_command {
	my %args = (
		id		=> undef,
		bedpe		=> undef,
		outdir		=> undef,
		@_
		);

	my $output_stem = join('/', $args{outdir}, $args{id});
	my $linked_bedpe = join('/', $args{outdir}, $args{id} . '.bedpe' );

	my $nuc_pos_command = "ln -s $args{bedpe} $linked_bedpe";

	$nuc_pos_command .= "\n\n# split bedpe by chromosome";
	$nuc_pos_command .= "\necho '>>> Running splitBedpe.sh <<<'";
	$nuc_pos_command .= "\n\n" . join(' ',
		"$fragmentomics_code_dir/nucleosome_peak/scripts/splitBedpe.sh",
		$linked_bedpe
		);

	$nuc_pos_command .= "\n\n# calculate distance to nucleosome peaks";
	$nuc_pos_command .= "\necho '>>> Running nucleosome_peaks_distance.R <<<'";
	$nuc_pos_command .= "\n\n" . join(' ',
		"Rscript $fragmentomics_code_dir/nucleosome_peak/scripts/nucleosome_peaks_distance.R",
		"--id", $args{id},
		"--path", $args{outdir},
		"--peaks", $reference_peaks, 
		"--outdir", $args{outdir}
		);

	$nuc_pos_command .= "\n\n# check for completeness and remove intermediates";
	$nuc_pos_command .= "\n" . join("\n",
		'if [[ -f ' . $output_stem . '_peak_distance.txt ]]; then',
		"  echo '>>> Nucleosome Position Profiling completed successfully! <<<'",
		'  rm ' . $output_stem . '*.bed*',
		'fi'
		);

	return($nuc_pos_command);
	}

# format command to extract insert size metrics
sub get_insert_sizes_command {
	my %args = (
		id		=> undef,
		bam		=> undef,
		outdir		=> undef,
		java_mem	=> undef,
		tmp_dir		=> undef,
		@_
		);

	my $output_stem = join('/', $args{outdir}, $args{id});

	my $fs_command = join(' ',
		'java -Xmx' . $args{java_mem},
		'-Djava.io.tmpdir=' . $args{tmp_dir},
		'-jar $picard_dir/picard.jar CollectInsertSizeMetrics',
		'I=' . $args{bam},
		'O=' . $output_stem . '_picard.txt',
		'H=' . $output_stem . '.pdf',
		'M=0 W=600'
		);

	return($fs_command);
	}

# format command to run Griffin GC functions
sub get_griffin_gc_command {
	my %args = (
		id	=> undef,
		bam	=> undef,
		outdir	=> undef,
		tmpdir	=> undef,
		n_cpu	=> 8,
		@_
		);

	my $gc_command;

	# using Griffin v0.1.0 (version provided with Pipeline-Suite)
	unless ($griffin_code_dir =~ m/v0.2.0/) {

		$gc_command = join(' ',
			"$griffin_code_dir/scripts/griffin_GC_counts.py",
			'--bam_file', $args{bam},
			'--bam_file_name', $args{id},
			'--mapable_regions', $griffin_mappable,
			'--ref_seq', $reference,
			'--chrom_sizes', "$griffin_code_dir/Ref/hg38.standard.chrom.sizes",
			'--out_dir', $args{outdir},
			'--map_q 20 --size_range 15 500 --CPU', $args{n_cpu}
			);

		$gc_command .= "\n\n" . join(' ',
			"$griffin_code_dir/scripts/griffin_GC_bias.py",
			'--bam_file_name', $args{id},
			'--mapable_name', $griffin_mappable_name,
			'--genome_GC_frequency', "$griffin_code_dir/Ref/genome_GC_frequency",
			'--out_dir', $args{outdir},
			'--size_range 15 500'
			);

		} else {

		# using Griffin v0.2.0 (currently doesn't work; testing in progress)
		$gc_command = "source activate base";
		$gc_command .= "\nconda activate griffin2";

		my $map_stem = join('/', $args{outdir}, 'mappability_bias', $args{id} . '.mappability_bias');

		$gc_command .= "\n\n" . join(' ',
			"$griffin_code_dir/scripts/griffin_mappability_correction.py",
			'--bam_file', $args{bam},
			'--bam_file_name', $args{id},
			'--output', $map_stem . '.txt',
			'--output_plot', $map_stem . '.pdf',
			'--mappability', "$griffin_code_dir/Ref/k50.Umap.MultiTrackMappability.hg38.bw",
			'--exclude_paths', "$griffin_code_dir/Ref/encode_unified_GRCh38_exclusion_list.bed",
			'--chrom_sizes', "$griffin_code_dir/Ref/hg38.standard.chrom.sizes",
			'--tmp_dir', $args{tmpdir},
			'--map_quality 20 --CPU', $args{n_cpu}
			);

		$gc_command .= "\n\n" . join(' ',
			"$griffin_code_dir/scripts/griffin_GC_counts.py",
			'--bam_file', $args{bam},
			'--bam_file_name', $args{id},
			'--mappable_regions_path', "$griffin_code_dir/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed",
			'--ref_seq', $reference,
			'--chrom_sizes', "$griffin_code_dir/Ref/hg38.standard.chrom.sizes",
			'--out_dir', $args{outdir},
			'--map_q 20 --size_range 15 500 --CPU', $args{n_cpu}
			);

		$gc_command .= "\n\n" . join(' ',
			"$griffin_code_dir/scripts/griffin_GC_bias.py",
			'--bam_file_name', $args{id},
			'--mappable_name k100_minus_exclusion_lists.mappable_regions.hg38',
			'--genome_GC_frequency', "$griffin_code_dir/Ref/genome_GC_frequency",
			'--out_dir', $args{outdir},
			'--size_range 15 500'
			);
		}

	return($gc_command);
	}

# format command to run Griffin nucleosome profiling
sub get_griffin_profiling_command {
	my %args = (
		id		=> undef,
		bam		=> undef,
		gc_bias		=> undef,
		map_bias	=> undef,
		n_cpu		=> 8,
		tmpdir		=> undef,
		outdir		=> undef,
		sites		=> undef,
		@_
		);

	my $griffin_command;

	# using Griffin v0.1.0 (version provided with Pipeline-Suite)
	unless ($griffin_code_dir =~ m/v0.2.0/) {

		$griffin_command = join(' ',
			"$griffin_code_dir/scripts/griffin_calc_coverage.py",
			'--sample_name', $args{id},
			'--bam', $args{bam},
			'--reference_genome', $reference,
			'--GC_bias', $args{gc_bias},
			'--background_normalization None',
			'--sites_yaml', $args{sites},
			'--results_dir', $args{outdir},
			'--chrom_column Chrom',
			'--strand_column Strand',
			'--chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22',
			'--norm_window -5000 5000',
			'--plot_window -1000 1000',
			'--fragment_length 165',
			'--step 15 --size_range 90 220 --map_quality 20',
			'--individual False --smoothing True --num_sites none',
			'--sort_by none --ascending none --cpu', $args{n_cpu},
			'--erase_intermediates False'
			);

		} else {

		# using Griffin v0.2.0 (currently doesn't work; testing in progress)
		$griffin_command = "source activate base";
		$griffin_command .= "\nconda activate griffin2";

		$griffin_command .= "\n\n" . join(' ',
			"$griffin_code_dir/scripts/griffin_coverage.py",
			'--sample_name', $args{id},
			'--bam', $args{bam},
			'--reference_genome', $reference,
			'--GC_bias', $args{gc_bias},
			'--mappability_bias', $args{map_bias},
			'--mappability_correction True',
			'--mappability_bw', "$griffin_code_dir/Ref/k50.Umap.MultiTrackMappability.hg38.bw",
			'--chrom_sizes_path', "$griffin_code_dir/Ref/hg38.standard.chrom.sizes",
			'--tmp_dir', $args{tmpdir},
			'--griffin_scripts', "$griffin_code_dir/scripts",
			'--sites_yaml', $args{sites},
			'--position_column position',
			'--size_range 90 220 --map_quality 20 --number_of_sites none',
			'--CPU', $args{n_cpu}
			);

		my $uncorrected_bw = join('/', $args{tmpdir}, $args{id}, 'tmp_bigWig', $args{id} . '.uncorrected.bw');
		my $corrected_bw = join('/', $args{tmpdir}, $args{id}, 'tmp_bigWig', $args{id} . '.GC_corrected.bw');
		my $map_corrected_bw = join('/', $args{tmpdir}, $args{id}, 'tmp_bigWig', $args{id} . '.GC_map_corrected.bw');

		my $exclude_encode = join('/', $griffin_code_dir, 'Ref', 'encode_unified_GRCh38_exclusion_list.bed');
		my $exclude_centromeres = join('/', $griffin_code_dir, 'Ref', 'hg38_centromeres.bed');
		my $exclude_gaps = join('/', $griffin_code_dir, 'Ref', 'hg38_gaps.bed');
		my $exclude_patches = join('/', $griffin_code_dir, 'Ref', 'hg38_fix_patches.bed');
		my $exclude_alts = join('/', $griffin_code_dir, 'Ref', 'hg38_alternative_haplotypes.bed');

		$griffin_command .= "\n\n" . join(' ',
			"$griffin_code_dir/scripts/griffin_merge_sites.py",
			'--sample_name', $args{id},
			'--uncorrected_bw_path', $uncorrected_bw,
			'--GC_corrected_bw_path', $corrected_bw,
			'--GC_map_corrected_bw_path', $map_corrected_bw,
			'--mappability_correction False',
			'--mappability_bw', "$griffin_code_dir/Ref/k50.Umap.MultiTrackMappability.hg38.bw",
			'--chrom_sizes_path', "$griffin_code_dir/Ref/hg38.standard.chrom.sizes",
			'--sites_yaml', $args{sites},
			'--griffin_scripts_dir', "$griffin_code_dir/scripts",
			'--tmp_dir', $args{tmpdir},
			'--results_dir', $args{outdir},
			'--position_column position',
			'--exclude_paths', $exclude_encode, $exclude_centromeres, $exclude_gaps, $exclude_patches, $exclude_alts,
			'--exclude_outliers True',
			'--exclude_zero_mappability True',
			'--step 15',
			'--CNA_normalization False',
			'--CPU', $args{n_cpu},
			'--individual False',
			'--smoothing True'
			);
		}

	$griffin_command .= "\n\n" . join(' ',
		"echo 'Griffin nucleosome profiling completed successfully' >",
		$args{outdir} . '/nucleosome_profiling.COMPLETE'
		);

	return($griffin_command);
	}

# format command to run End Motif profiling
sub get_end_motif_command {
	my %args = (
		id		=> undef,
		bedpe		=> undef,
		outdir		=> undef,
		@_
		);

	my $output_stem = join('/', $args{outdir}, $args{id});
	my $linked_bedpe = join('/', $args{outdir}, $args{id} . '.bedpe' );

	my $motif_command = "ln -s $args{bedpe} $linked_bedpe";

	# format bedpe to bed
	$motif_command .= "\n\n# convert bedpe to bed";
	$motif_command .= "\necho '>>> Running motif_format_bedpe.R <<<'";
	$motif_command .= "\n\n" . join(' ',
		"Rscript $fragmentomics_code_dir/end_motif/scripts/motif_format_bedpe.R",
		'--id', $args{id},
		'--bedpe', $linked_bedpe,
		'--outdir', $args{outdir}
		);

	# get fasta sequences
	$motif_command .= "\n\n# extract fasta sequences";
	$motif_command .= "\necho '>>> Running bedtools getfasta <<<'";
	$motif_command .= "\n\n" . join(' ',
		'bedtools getfasta',
		'-bedOut -fi', $reference,
		'-bed', $output_stem . '_5.bed',
		'>', $output_stem . '_fasta_5.bed'
		);

	$motif_command .= "\n\n" . join(' ',
		'bedtools getfasta',
		'-bedOut -fi', $reference,
		'-bed', $output_stem . '_3.bed',
		'>', $output_stem . '_fasta_3.bed'
		);

	# convert fasta to end motif context frequencies
	$motif_command .= "\n\n# run end motif context profiling";
	$motif_command .= "\necho '>>> Running motif_get_contexts.R <<<'";
	$motif_command .= "\n\n" . join(' ',
		"Rscript $fragmentomics_code_dir/end_motif/scripts/motif_get_contexts.R",
		"--id", $args{id},
		"--fasta_5", $output_stem . '_fasta_5.bed',
		"--fasta_3", $output_stem . '_fasta_3.bed',
		"--outdir", $args{outdir}
		);

	# remove intermediate files
	$motif_command .= "\n\n# check for completeness and remove intermediates";
	$motif_command .= "\n" . join("\n",
		'if [[ -f ' . $output_stem . '_motifs.txt ]]; then',
		"  echo '>>> End Motif Profiling completed successfully! <<<'",
		'  rm ' . $output_stem . '*.bed*',
		'fi'
		);

	return($motif_command);
	}

# format command to run breakpoint profiling
sub get_profile_breakpoints_command {
	my %args = (
		id		=> undef,
		bedpe		=> undef,
		outdir		=> undef,
		@_
		);

	my $output_stem = join('/', $args{outdir}, $args{id});
	my $linked_bedpe = join('/', $args{outdir}, $args{id} . '.bedpe' );

	my $bkpt_command = "ln -s $args{bedpe} $linked_bedpe";

	# format bedpe to bed
	$bkpt_command .= "\n\n# convert bedpe to bed";
	$bkpt_command .= "\necho '>>> Running breakpoint_format_bedpe.R <<<'";
	$bkpt_command .= "\n\n" . join(' ',
		"Rscript $fragmentomics_code_dir/breakpoint/scripts/breakpoint_format_bedpe.R",
		'--id', $args{id},
		'--bedpe', $linked_bedpe,
		'--outdir', $args{outdir}
		);

	# get fasta sequences
	$bkpt_command .= "\n\n# extract fasta sequences";
	$bkpt_command .= "\necho '>>> Running bedtools getfasta <<<'";
	$bkpt_command .= "\n\n" . join(' ',
		'bedtools getfasta',
		'-bedOut -fi', $reference,
		'-bed', $output_stem . '_5.bed',
		'>', $output_stem . '_fasta_5.bed'
		);

	$bkpt_command .= "\n\n" . join(' ',
		'bedtools getfasta',
		'-bedOut -fi', $reference,
		'-bed', $output_stem . '_3.bed',
		'>', $output_stem . '_fasta_3.bed'
		);

	# convert fasta to end motif context frequencies
	$bkpt_command .= "\n\n# run breakpoint context profiling";
	$bkpt_command .= "\necho '>>> Running breakpoint_get_contexts.R <<<'";
	$bkpt_command .= "\n\n" . join(' ',
		"Rscript $fragmentomics_code_dir/breakpoint/scripts/breakpoint_get_contexts.R",
		"--id", $args{id},
		"--fasta_5", $output_stem . '_fasta_5.bed',
		"--fasta_3", $output_stem . '_fasta_3.bed',
		"--outdir", $args{outdir}
		);

	# remove intermediate files
	$bkpt_command .= "\n\n# check for completeness and remove intermediates";
	$bkpt_command .= "\n" . join("\n",
		'if [[ -f ' . $output_stem . '_ratio.txt ]]; then',
		"  echo '>>> Breakpoint Profiling completed successfully! <<<'",
		'  rm ' . $output_stem . '*.bed*',
		'fi'
		);

	return($bkpt_command);
	}

# format command to run dinucleotide profiling
sub get_dinucleotide_profiling_command {
	my %args = (
		id		=> undef,
		bedpe		=> undef,
		outdir		=> undef,
		@_
		);

	my $output_stem = join('/', $args{outdir}, $args{id});

	# format bedpe to bed
	my $din_command = "# convert bedpe to bed";
	$din_command .= "\necho '>>> Running dinucleotide_format_bedpe.R <<<'";
	$din_command .= "\n\n" . join(' ',
		"Rscript $fragmentomics_code_dir/dinucleotide/scripts/dinucleotide_format_bedpe.R",
		'--id', $args{id},
		'--bedpe', $args{bedpe},
		'--outdir', $args{outdir}
		);

	# get fasta sequences
	$din_command .= "\n\n# extract fasta sequences";
	$din_command .= "\necho '>>> Running bedtools getfasta <<<'";
	$din_command .= "\n\n" . join(' ',
		'bedtools getfasta',
		'-bedOut -fi', $reference,
		'-bed', $output_stem . '.bed',
		'>', $output_stem . '_fasta.bed'
		);

	# convert fasta to end motif context frequencies
	$din_command .= "\n\n# run dinucleotide context profiling";
	$din_command .= "\necho '>>> Running dinucleotide_get_contexts.R <<<'";
	$din_command .= "\n\n" . join(' ',
		"Rscript $fragmentomics_code_dir/dinucleotide/scripts/dinucleotide_get_contexts.R",
		"--id", $args{id},
		"--fasta", $output_stem . '_fasta.bed',
		"--outdir", $args{outdir}
		);

	# remove intermediate files
	$din_command .= "\n\n# check for completeness and remove intermediates";
	$din_command .= "\n" . join("\n",
		'if [[ -f ' . $output_stem . '_contexts.txt ]]; then',
		"  echo '>>> Dinucleotide Profiling completed successfully! <<<'",
		#'  rm ' . $output_stem . '*.bed*',
		'fi'
		);

	return($din_command);
	}

### MAIN ###########################################################################################
sub main {
	my %args = (
		tool_config		=> undef,
		data_config		=> undef,
		output_directory	=> undef,
		hpc_driver		=> undef,
		dry_run			=> undef,
		no_wait			=> undef,
		@_
		);

	my $tool_config = $args{tool_config};
	my $data_config = $args{data_config};

	### PREAMBLE ######################################################################################
	unless($args{dry_run}) {
		print "Initiating Fragmenomics pipeline...\n";
		}

	# load tool config
	my $tool_data = LoadFile($tool_config);

	# organize output and log directories
	my $output_directory = $args{output_directory};
	$output_directory =~ s/\/$//;

	my $log_directory = join('/', $output_directory, 'logs');
	unless(-e $log_directory) { make_path($log_directory); }

	my $log_file = join('/', $log_directory, 'run_Fragmentomics_pipeline.log');

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

		$log_file = join('/', $log_directory, 'run_Fragmentomics_pipeline_' . $run_count . '.log');
		}

	# start logging
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";
	$log->autoflush;

	print $log "---\n";
	print $log "Running Fragmentomics pipeline.\n";
	print $log "\n  Tool config used: $tool_config";
	print $log "\n  >> Reference used: $tool_data->{reference}";
	print $log "\n  Sample config used: $data_config";
	print $log "\n  Output directory: $output_directory";
	print $log "\n---";

	# get user-specified tool parameters
	my $parameters = $tool_data->{parameters};
 
	# set tools and versions
	my $picard	= 'picard/' . $tool_data->{picard_version};
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};
	my $bedtools	= 'bedtools/' . $tool_data->{bedtools_version};
	my $python3	= 'python3/' . $tool_data->{python3_version};
	my $r_version	= 'R/' . $tool_data->{r_version};

	if (defined($tool_data->{fragmentomics_dir})) {
		$fragmentomics_code_dir = $tool_data->{fragmentomics_dir};
		} else {
		$fragmentomics_code_dir = "$cwd/scripts/fragmentomics";
		}

	if (defined($tool_data->{griffin_dir})) {
		$griffin_code_dir = $tool_data->{griffin_dir};
		} else {
		$griffin_code_dir = "$cwd/scripts/fragmentomics/griffin";
		}

	# set reference files
	$reference = $tool_data->{reference};

	$reference_peaks = defined($tool_data->{nucleosome_peaks}) ? $tool_data->{nucleosome_peaks} : "$fragmentomics_code_dir/nucleosome_peak/ref/hg38";

	$delfi_filters	= defined($tool_data->{delfi_filters}) ? $tool_data->{delfi_filters} : "$fragmentomics_code_dir/ratio/extdata/filters.hg38.rda";
	$delfi_gaps	= defined($tool_data->{delfi_gaps}) ? $tool_data->{delfi_gaps} : "$fragmentomics_code_dir/ratio/extdata/gaps.hg38.rda";
	$delfi_vntrs	= defined($tool_data->{delfi_vntrs}) ? $tool_data->{delfi_vntrs} : "$fragmentomics_code_dir/ratio/extdata/VNTRs.hg38.rda";
	$delfi_tiles	= defined($tool_data->{delfi_tiles}) ? $tool_data->{delfi_tiles} : "$fragmentomics_code_dir/ratio/extdata/hg38_tiles.bed";
	$delfi_healthy	= defined($tool_data->{delfi_pon}) ? $tool_data->{delfi_pon} : "$fragmentomics_code_dir/ratio/extdata/healthy.median.hg38.rda";

	if (defined($tool_data->{griffin_ref_mappable})) {
		$griffin_mappable = $tool_data->{griffin_ref_mappable};
		} else {
		$griffin_mappable = "$cwd/scripts/fragmentomics/griffin/Ref/repeat_masker.mapable.k50.Umap.hg38.bedGraph.gz";
		}

	$griffin_mappable_name = basename $griffin_mappable;
	$griffin_mappable_name =~ s/\.[^.]+$//;

	$tfbs_sites = defined($tool_data->{TFBS_sites}) ? $tool_data->{TFBS_sites} : "$griffin_code_dir/site_configs/TFBS_sites.yaml";
	$tcga_sites = defined($tool_data->{TCGA_sites}) ? $tool_data->{TCGA_sites} : "$griffin_code_dir/site_configs/TCGA_sites.yaml";
	$dhs_sites = defined($tool_data->{DHS_sites}) ? $tool_data->{DHS_sites} : "$griffin_code_dir/site_configs/DHS_sites.yaml";
	$hk_sites = defined($tool_data->{housekeeping_sites}) ? $tool_data->{housekeeping_sites} : "$griffin_code_dir/site_configs/housekeeping_sites.yaml";

	# get optional HPC group
	my $hpc_group = defined($tool_data->{hpc_group}) ? "-A $tool_data->{hpc_group}" : undef;
	
	### RUN ###########################################################################################
	# first, determine which steps/fragmentomics analyses to run
	my %tool_set = (
		'downsample'	=> 'Y',	#defined($parameters->{downsample}) ? 'Y' : 'N',
		'deduplicate'	=> 'N', #defined($parameters->{dedup}) ? 'Y' : 'N',
		'ratio'		=> defined($parameters->{ratio}) ? 'Y' : 'N',
		'score'		=> defined($parameters->{score}) ? 'Y' : 'N',
		'insertsize'	=> defined($parameters->{insertsize}) ? 'Y' : 'N',
		'griffin'	=> 'N',
		'griffin_gc'	=> defined($parameters->{griffin_gc}) ? 'Y' : 'N',
		'griffin_tfbs'	=> defined($parameters->{griffin_tfbs}) ? 'Y' : 'N',
		'griffin_tcga'	=> defined($parameters->{griffin_tcga}) ? 'Y' : 'N',
		'griffin_dhs'	=> defined($parameters->{griffin_dhs}) ? 'Y' : 'N',
		'griffin_hk'	=> defined($parameters->{griffin_housekeeping}) ? 'Y' : 'N',
		'nucleosome'	=> defined($parameters->{nucleosome_position}) ? 'Y' : 'N',
		'dinucleotide'	=> defined($parameters->{dinucleotide}) ? 'Y' : 'N',
		'end_motifs'	=> defined($parameters->{end_motifs}) ? 'Y' : 'N',
		'breakpoints'	=> defined($parameters->{breakpoints}) ? 'Y' : 'N'
		);

	if ( ('Y' eq $tool_set{'nucleosome'}) || 
		('Y' eq $tool_set{'end_motifs'}) || ('Y' eq $tool_set{'breakpoints'}) ) {
		$tool_set{'deduplicate'} = 'Y';
		}

	if ( ('Y' eq $tool_set{'griffin_gc'}) || ('Y' eq $tool_set{'griffin_tfbs'}) ||
		('Y' eq $tool_set{'griffin_tcga'}) || ('Y' eq $tool_set{'griffin_dhs'}) ||
		('Y' eq $tool_set{'griffin_hk'}) ) {
		$tool_set{'griffin'} = 'Y';
		}

	print $log "\nTools/steps to run:\n";
	print $log Dumper \%tool_set;

	# get sample data
	my $smp_data = LoadFile($data_config);

	unless($args{dry_run}) {
		print "Processing " . scalar(keys %{$smp_data}) . " patients.\n";
		}

	# initialize objects
	my ($run_script, $run_id);
	my @all_jobs;

	# process each patient in $smp_data
	foreach my $patient (sort keys %{$smp_data}) {

		print $log "\n---\nInitiating process for PATIENT: $patient";

		# find bams
		my @tumour_ids = sort keys %{$smp_data->{$patient}->{'tumour'}};

		if (scalar(@tumour_ids) == 0) {
			print $log "\n>> No tumour BAM provided, skipping patient.\n";
			next;
			}

		# create some directories
		my $patient_directory = join('/', $output_directory, $patient);
		my $tmp_directory = join('/', $patient_directory, 'TEMP');

		unless(-e $patient_directory) { make_path($patient_directory); }
		unless(-e $tmp_directory) { make_path($tmp_directory); }

		my (@final_outputs, @patient_jobs);

		my $cleanup_cmd = "  rm -rf $tmp_directory;";

		# process each sample provided for this patient
		foreach my $tumour (@tumour_ids) {

			# initiate some variables
			my $tumour_stem = $tumour;
			my $tumour_bam = $smp_data->{$patient}->{tumour}->{$tumour};

			my ($downsample_run_id, $dedup_run_id, $griffin_id) = '';

			# set up directory structure
			my $sample_directory = join('/', $patient_directory, $tumour);
			my $downsample_directory = join('/', $sample_directory, 'downsample');
			my $ratio_directory = join('/', $sample_directory, 'fragment_ratio');
			my $score_directory = join('/', $sample_directory, 'fragment_score');
			my $insertsize_directory = join('/', $sample_directory, 'insert_size');
			my $nucleosome_directory = join('/', $sample_directory, 'nucleosome_peaks');
			my $griffin_directory = join('/', $sample_directory, 'griffin');
			my $motif_directory = join('/', $sample_directory, 'end_motifs');
			my $breakpoint_directory = join('/', $sample_directory, 'breakpoints');
			my $di_directory = join('/', $sample_directory, 'dinucleotide');

			unless(-e $sample_directory) { make_path($sample_directory); }

			if ('Y' eq $tool_set{'downsample'}) {
				unless(-e $downsample_directory) { make_path($downsample_directory); }
				}
			if ('Y' eq $tool_set{'ratio'}) {
				unless(-e $ratio_directory) { make_path($ratio_directory); }
				}
			if ('Y' eq $tool_set{'score'}) {
				unless(-e $score_directory) { make_path($score_directory); }
				}
			if ('Y' eq $tool_set{'insertsize'}) {
				unless(-e $insertsize_directory) { make_path($insertsize_directory); }
				}
			if ('Y' eq $tool_set{'nucleosome'}) {
				unless(-e $nucleosome_directory) { make_path($nucleosome_directory); }
				}
			if ('Y' eq $tool_set{'griffin'}) {
				unless(-e $griffin_directory) { make_path($griffin_directory); }	
				}
			if ('Y' eq $tool_set{'end_motifs'}) {
				unless(-e $motif_directory) { make_path($motif_directory); }	
				}
			if ('Y' eq $tool_set{'breakpoints'}) {
				unless(-e $breakpoint_directory) { make_path($breakpoint_directory); }	
				}
			if ('Y' eq $tool_set{'dinucleotide'}) {
				unless(-e $di_directory) { make_path($di_directory); }	
				}

			# run downsample
			if ('Y' eq $tool_set{'downsample'}) {

				my $downsample_cmd = create_downsample_command(
					bam	=> $smp_data->{$patient}->{tumour}->{$tumour},
					id	=> $tumour,
					n_reads	=> $parameters->{downsample}->{n_reads},
					outdir	=> $downsample_directory,
					tmpdir	=> $tmp_directory
					);

				my $downsampled_bam = join('/',
					$downsample_directory,
					$tumour . '_downsampled.bam'
					);

				# note this stem to avoid repetitive code later
				$tumour_stem = $tumour . '_downsample';
				$tumour_bam = $downsampled_bam;

				# add to cleanup
				$cleanup_cmd .= "\n  rm $downsampled_bam;";

				# check if this should be run
				if ('Y' eq missing_file($downsampled_bam)) {

					# record command (in log directory) and then run job
					print $log "\n >> Submitting job for DownsampleBAM...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_downsample_bam_' . $tumour,
						cmd	=> $downsample_cmd,
						modules	=> [$picard, $samtools],
						max_time	=> $parameters->{downsample}->{time},
						mem		=> $parameters->{downsample}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$downsample_run_id = submit_job(
						jobname		=> 'run_downsample_bam_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $downsample_run_id;
					push @all_jobs, $downsample_run_id;
					} else {
					print $log "\n >> Skipping Downsample step because this has already been completed!\n";
					}
				}

			# run de-duplicate (and generate bedpe files)
			my ($deduped_bedpe, $deduped_bedpe_hc);

			if ('Y' eq $tool_set{'deduplicate'}) {

				my $dedup_cmd = get_deduplicate_command( 
					id		=> $tumour_stem,
					bam		=> $tumour_bam,
					outdir		=> $downsample_directory,
					tmp_dir		=> $tmp_directory,
					java_mem	=> $parameters->{dedup}->{java_mem}
					);

				$deduped_bedpe = join('/',
					$downsample_directory,
					$tumour_stem . '_all.bedpe'
					);

				$deduped_bedpe_hc = join('/',
					$downsample_directory,
					$tumour_stem . '_q30.bedpe'
					);

				# add to cleanup
				$cleanup_cmd .= "\n  rm $deduped_bedpe;";
				$cleanup_cmd .= "\n  rm $deduped_bedpe_hc;";

				# check if this should be run
				if ('Y' eq missing_file($deduped_bedpe_hc)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for DeduplicateBAM...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_deduplicate_and_make_bedpe_' . $tumour,
						cmd	=> $dedup_cmd,
						modules	=> [$picard, $samtools, $bedtools],
						dependencies	=> $downsample_run_id,
						max_time	=> $parameters->{dedup}->{time},
						mem		=> $parameters->{dedup}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$dedup_run_id = submit_job(
						jobname		=> 'run_deduplicate_and_make_bedpe_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $dedup_run_id;
					push @all_jobs, $dedup_run_id;
					} else {
					print $log " >> Skipping deduplicate/bedpe step because this has already been completed!\n";
					}
				}

			### run DELFI fragment ratio
			if ('Y' eq $tool_set{'ratio'}) {

				my $ratio_cmd = create_ratio_command(
					bam	=> $tumour_bam,
					id	=> $tumour_stem,
					outdir	=> $ratio_directory
					);

				my $ratio_output = join('/',
					$ratio_directory,
					$tumour_stem . '_combined.pdf'
					);

				# check if this should be run
				if ('Y' eq missing_file($ratio_output)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for DELFI...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_delfi_fragment_ratio_' . $tumour,
						cmd	=> $ratio_cmd,
						modules	=> [$r_version],
						dependencies	=> $downsample_run_id,
						max_time	=> $parameters->{ratio}->{time},
						mem		=> $parameters->{ratio}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_delfi_fragment_ratio_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping DELFI because this has already been completed!\n";
					}

				push @final_outputs, $ratio_output;
				}


			### run VESSIES fragment score
			if ('Y' eq $tool_set{'score'}) {

				my $score_cmd = create_fragmentscore_command(
					id	=> $tumour_stem,
					bam	=> $tumour_bam,
					outdir	=> $score_directory
					);

				my $score_output = join('/',
					$score_directory,
					$tumour_stem . '_score.txt'
					);

				# check if this should be run
				if ('Y' eq missing_file($score_output)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for VESSIES fragment score...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_vessies_fragment_score_' . $tumour,
						cmd	=> $score_cmd,
						modules	=> [$r_version],
						dependencies	=> $downsample_run_id,
						max_time	=> $parameters->{score}->{time},
						mem		=> $parameters->{score}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_vessies_fragment_score_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping VESSIES score because this has already been completed!\n";
					}

				push @final_outputs, $score_output;
				}


			### Extract insert size metrics
			if ('Y' eq $tool_set{'insertsize'}) {

				my $insertsize_cmd = get_insert_sizes_command(
					id		=> $tumour_stem,
					bam		=> $tumour_bam,
					outdir		=> $insertsize_directory,
					java_mem	=> $parameters->{insertsize}->{java_mem},
					tmp_dir		=> $tmp_directory
					);

				my $insert_sizes = join('/',
					$insertsize_directory,
					$tumour_stem . '_picard.txt'
					);

				# check if this should be run
				if ('Y' eq missing_file($insert_sizes)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for InsertSize...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_picard_insert_size_' . $tumour,
						cmd	=> $insertsize_cmd,
						modules	=> [$picard, $r_version],
						dependencies	=> $downsample_run_id,
						max_time	=> $parameters->{insertsize}->{time},
						mem		=> $parameters->{insertsize}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_picard_insert_size_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping InsertSizes because this has already been completed!\n";
					}

				push @final_outputs, $insert_sizes;
				}


			### Griffin Profiling
			my ($griffin_map_dir, $griffin_gc_dir, $map_bias_out, $gc_bias_out);

			if ('Y' eq $tool_set{'griffin'}) {

				# different versions require different bias files
				if ($griffin_code_dir =~ m/v0.2.0/) {

					$griffin_map_dir = join('/', $griffin_directory, 'mappability_bias');
					$griffin_gc_dir = join('/', $griffin_directory, 'GC_bias');

					unless(-e $griffin_map_dir) { make_path($griffin_map_dir); }
					unless(-e $griffin_gc_dir) { make_path($griffin_gc_dir); }

					$map_bias_out = join('/', $griffin_map_dir, $tumour_stem . '.mappability_bias.txt');
					$gc_bias_out = join('/', $griffin_gc_dir, $tumour_stem . '.GC_bias.txt');
					} else {

					$gc_bias_out = join('/',
						$griffin_directory,
						$griffin_mappable_name,
						'GC_bias',
						$tumour_stem . '.GC_bias.txt'
						);
					}

				# run GC correction functions
				my $griffin_cmd = get_griffin_gc_command(
					id	=> $tumour_stem,
					bam	=> $tumour_bam,
					outdir	=> $griffin_directory,
					tmpdir	=> $tmp_directory,
					n_cpu	=> $parameters->{griffin_gc}->{n_cpus}
					);

				# check if this should be run
				if ('Y' eq missing_file($gc_bias_out)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for GC correction...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_griffin_gc_correction_' . $tumour,
						cmd	=> $griffin_cmd,
						modules	=> [$python3, $bedtools],
						dependencies	=> $downsample_run_id,
						cpus_per_task	=> $parameters->{griffin_gc}->{n_cpus},
						max_time	=> $parameters->{griffin_gc}->{time},
						mem		=> $parameters->{griffin_gc}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$griffin_id = submit_job(
						jobname		=> 'run_griffin_gc_correction_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $griffin_id;
					push @all_jobs, $griffin_id;
					} else {
					print $log " >> Skipping Griffin GC Correction because this has already been completed!\n";
					}
				}

			# run Griffin TFBS profiling
			if ('Y' eq $tool_set{'griffin_tfbs'}) {

				my $griffin_tfbs_dir = join('/', $griffin_directory, 'TFBS');
				unless(-e $griffin_tfbs_dir) { make_path($griffin_tfbs_dir); }

				my $griffin_cmd = get_griffin_profiling_command(
					id	=> $tumour_stem,
					bam	=> $tumour_bam,
					sites	=> $tfbs_sites,
					outdir	=> $griffin_tfbs_dir,
					gc_bias		=> $gc_bias_out,
					map_bias	=> $map_bias_out,
					tmpdir	=> $tmp_directory,
					n_cpu	=> $parameters->{griffin_tfbs}->{n_cpus}
					);

				my $tfbs_output = join('/',
					$griffin_tfbs_dir,
					'nucleosome_profiling.COMPLETE'
					);

				# check if this should be run
				if ('Y' eq missing_file($tfbs_output)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for Griffin TFBS Coverage...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_griffin_tfbs_coverage_' . $tumour,
						cmd	=> $griffin_cmd,
						modules	=> [$python3, $bedtools],
						dependencies	=> $griffin_id,
						cpus_per_task	=> $parameters->{griffin_tfbs}->{n_cpus},
						max_time	=> $parameters->{griffin_tfbs}->{time},
						mem		=> $parameters->{griffin_tfbs}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_griffin_tfbs_coverage_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping Griffin TFBS Coverage because this has already been completed!\n";
					}

				push @final_outputs, $tfbs_output;
				}

			# run Griffin DHS profiling
			if ('Y' eq $tool_set{'griffin_dhs'}) {

				my $griffin_dhs_dir = join('/', $griffin_directory, 'DHS');
				unless(-e $griffin_dhs_dir) { make_path($griffin_dhs_dir); }

				my $griffin_cmd = get_griffin_profiling_command(
					id	=> $tumour_stem,
					bam	=> $tumour_bam,
					sites	=> $dhs_sites,
					outdir	=> $griffin_dhs_dir,
					gc_bias		=> $gc_bias_out,
					map_bias	=> $map_bias_out,
					tmpdir	=> $tmp_directory,
					n_cpu	=> $parameters->{griffin_dhs}->{n_cpus}
					);

				my $dhs_output = join('/',
					$griffin_dhs_dir,
					'nucleosome_profiling.COMPLETE'
					);

				# check if this should be run
				if ('Y' eq missing_file($dhs_output)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for Griffin DHS Coverage...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_griffin_dhs_coverage_' . $tumour,
						cmd	=> $griffin_cmd,
						modules	=> [$python3, $bedtools],
						dependencies	=> $griffin_id,
						cpus_per_task	=> $parameters->{griffin_dhs}->{n_cpus},
						max_time	=> $parameters->{griffin_dhs}->{time},
						mem		=> $parameters->{griffin_dhs}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_griffin_dhs_coverage_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping Griffin DHS Coverage because this has already been completed!\n";
					}

				push @final_outputs, $dhs_output;
				}

			# run Griffin TCGA tumour type profiling
			if ('Y' eq $tool_set{'griffin_tcga'}) {

				my $griffin_tcga_dir = join('/', $griffin_directory, 'TCGA');
				unless(-e $griffin_tcga_dir) { make_path($griffin_tcga_dir); }

				my $griffin_cmd = get_griffin_profiling_command(
					id	=> $tumour_stem,
					bam	=> $tumour_bam,
					sites	=> $tcga_sites,
					outdir	=> $griffin_tcga_dir,
					gc_bias		=> $gc_bias_out,
					map_bias	=> $map_bias_out,
					tmpdir	=> $tmp_directory,
					n_cpu	=> $parameters->{griffin_tcga}->{n_cpus}
					);

				my $tcga_output = join('/',
					$griffin_tcga_dir,
					'nucleosome_profiling.COMPLETE'
					);

				# check if this should be run
				if ('Y' eq missing_file($tcga_output)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for Griffin TCGA Coverage...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_griffin_tcga_coverage_' . $tumour,
						cmd	=> $griffin_cmd,
						modules	=> [$python3, $bedtools],
						dependencies	=> $griffin_id,
						cpus_per_task	=> $parameters->{griffin_tcga}->{n_cpus},
						max_time	=> $parameters->{griffin_tcga}->{time},
						mem		=> $parameters->{griffin_tcga}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_griffin_tcga_coverage_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping Griffin TCGA Coverage because this has already been completed!\n";
					}

				push @final_outputs, $tcga_output;
				}

			# run Griffin across housekeeping sites
			if ('Y' eq $tool_set{'griffin_hk'}) {

				my $griffin_hk_dir = join('/', $griffin_directory, 'housekeeping');
				unless(-e $griffin_hk_dir) { make_path($griffin_hk_dir); }

				my $griffin_cmd = get_griffin_profiling_command(
					id	=> $tumour_stem,
					bam	=> $tumour_bam,
					sites	=> $hk_sites,
					outdir	=> $griffin_hk_dir,
					gc_bias		=> $gc_bias_out,
					map_bias	=> $map_bias_out,
					tmpdir	=> $tmp_directory,
					n_cpu	=> $parameters->{griffin_housekeeping}->{n_cpus}
					);

				my $hk_output = join('/',
					$griffin_hk_dir,
					'nucleosome_profiling.COMPLETE'
					);

				# check if this should be run
				if ('Y' eq missing_file($hk_output)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for Griffin housekeeping Coverage...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_griffin_hk_coverage_' . $tumour,
						cmd	=> $griffin_cmd,
						modules	=> [$python3, $bedtools],
						dependencies	=> $griffin_id,
						cpus_per_task	=> $parameters->{griffin_housekeeping}->{n_cpus},
						max_time	=> $parameters->{griffin_housekeeping}->{time},
						mem		=> $parameters->{griffin_housekeeping}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_griffin_hk_coverage_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping Griffin housekeeping Coverage because this has already been completed!\n";
					}

				push @final_outputs, $hk_output;
				}


			### Nucleosome positioning
			if ('Y' eq $tool_set{'nucleosome'}) {

				my $position_cmd = create_nucleosome_position_command(
					id	=> $tumour_stem,
					bedpe	=> $deduped_bedpe,
					outdir	=> $nucleosome_directory
					);

				my $peaks_output = join('/',
					$nucleosome_directory,
					$tumour_stem . '_peak_distance.txt'
					);

				# check if this should be run
				if ('Y' eq missing_file($peaks_output)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for Nucleosome Positioning...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_nucleosome_positioning_' . $tumour,
						cmd	=> $position_cmd,
						modules	=> [$r_version],
						dependencies	=> $dedup_run_id,
						max_time	=> $parameters->{nucleosome_position}->{time},
						mem		=> $parameters->{nucleosome_position}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_nucleosome_positioning_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping NucleosomePositioning because this has already been completed!\n";
					}

				push @final_outputs, $peaks_output;
				}


			### End Motif profiling
			if ('Y' eq $tool_set{'end_motifs'}) {

				my $endmotif_cmd = get_end_motif_command(
					id	=> $tumour_stem,
					bedpe	=> $deduped_bedpe_hc,
					outdir	=> $motif_directory
					);

				my $motifs_output = join('/',
					$motif_directory,
					$tumour_stem . '_motifs.txt'
					);

				# check if this should be run
				if ('Y' eq missing_file($motifs_output)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for EndMotif Profiling...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_endmotif_profiling_' . $tumour,
						cmd	=> $endmotif_cmd,
						modules	=> [$samtools, $picard, $bedtools, $r_version],
						dependencies	=> $dedup_run_id,
						max_time	=> $parameters->{end_motifs}->{time},
						mem		=> $parameters->{end_motifs}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_endmotif_profiling_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping EndMotif profiling because this has already been completed!\n";
					}

				push @final_outputs, $motifs_output;
				}


			### Breakpoint profiling
			if ('Y' eq $tool_set{'breakpoints'}) {

				my $breakpoint_cmd = get_profile_breakpoints_command(
					id	=> $tumour_stem,
					bedpe	=> $deduped_bedpe,
					outdir	=> $breakpoint_directory
					);

				my $breakpoint_output = join('/',
					$breakpoint_directory,
					$tumour_stem . '_ratio.txt'
					);

				# check if this should be run
				if ('Y' eq missing_file($breakpoint_output)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for Breakpoint Profiling...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_breakpoint_profiling_' . $tumour,
						cmd	=> $breakpoint_cmd,
						modules	=> [$samtools, $picard, $bedtools, $r_version],
						dependencies	=> $dedup_run_id,
						max_time	=> $parameters->{breakpoints}->{time},
						mem		=> $parameters->{breakpoints}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_breakpoint_profiling_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping Breakpoint profiling because this has already been completed!\n";
					}

				push @final_outputs, $breakpoint_output;
				}


			### Dinucleotide profiling
			if ('Y' eq $tool_set{'dinucleotide'}) {

				# run size selection and de-duplicate (and generate bedpe files)
				my $dedup_cmd_part1 = get_sized_deduplicate_command( 
					id		=> $tumour_stem,
					bam		=> $tumour_bam,
					size 		=> 'genome',
					outdir		=> $di_directory,
					tmp_dir		=> $tmp_directory,
					java_mem	=> $parameters->{dedup}->{java_mem}
					);

				my $dedup_cmd_part2 = get_sized_deduplicate_command( 
					id		=> $tumour_stem,
					bam		=> $tumour_bam,
					size 		=> 'short',
					outdir		=> $di_directory,
					tmp_dir		=> $tmp_directory,
					java_mem	=> $parameters->{dedup}->{java_mem}
					);

				my $normal_bedpe = join('/',
					$di_directory,
					$tumour_stem . '_len167.bedpe'
					);

				my $short_bedpe = join('/',
					$di_directory,
					$tumour_stem . '_len150.bedpe'
					);

				my ($size_dedup_cmd, $required_bedpe);
				if ('genome' eq $parameters->{dinucleotide}->{size}) {
 					$size_dedup_cmd = $dedup_cmd_part1;
					$required_bedpe = $normal_bedpe;
					} elsif ('short' eq $parameters->{dinucleotide}->{size}) {
 					$size_dedup_cmd = $dedup_cmd_part2;
					$required_bedpe = $short_bedpe;
					} elsif ('both' eq $parameters->{dinucleotide}->{size}) {
 					$size_dedup_cmd = $dedup_cmd_part1 . "\n\n" . $dedup_cmd_part2;
					$required_bedpe = $short_bedpe;
					}

				# add to cleanup
				$cleanup_cmd .= "\n  rm $normal_bedpe;";
				$cleanup_cmd .= "\n  rm $short_bedpe;";

				# check if this should be run
				if ('Y' eq missing_file($required_bedpe)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for DeduplicateBAM...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_format_bedpe_for_dinucleotide_' . $tumour,
						cmd	=> $size_dedup_cmd,
						modules	=> [$picard, $samtools, $bedtools],
						dependencies	=> $downsample_run_id,
						max_time	=> $parameters->{dedup}->{time},
						mem		=> $parameters->{dedup}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$dedup_run_id = submit_job(
						jobname		=> 'run_format_bedpe_for_dinucleotide_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $dedup_run_id;
					push @all_jobs, $dedup_run_id;
					} else {
					print $log " >> Skipping deduplicate/bedpe step because this has already been completed!\n";
					}

				# format command for dinucleotide profilling
				my $dinucleotide_cmd_part1 = get_dinucleotide_profiling_command(
					id	=> $tumour_stem . '_len167',
					bedpe	=> $normal_bedpe,
					outdir	=> $di_directory
					);

				my $dinucleotide_cmd_part2 = get_dinucleotide_profiling_command(
					id	=> $tumour_stem . '_len150',
					bedpe	=> $short_bedpe,
					outdir	=> $di_directory
					);

				my $dinucleotide_normal = join('/',
					$di_directory,
					$tumour_stem . '_len167_contexts.txt'
					);

				my $dinucleotide_short = join('/',
					$di_directory,
					$tumour_stem . '_len150_contexts.txt'
					);

				my ($dinucleotide_cmd, $required_out);
				if ('genome' eq $parameters->{dinucleotide}->{size}) {
 					$dinucleotide_cmd = $dinucleotide_cmd_part1;
					$required_out = $dinucleotide_normal;
					} elsif ('short' eq $parameters->{dinucleotide}->{size}) {
 					$dinucleotide_cmd = $dinucleotide_cmd_part2;
					$required_out = $dinucleotide_short;
					} elsif ('both' eq $parameters->{dinucleotide}->{size}) {
 					$dinucleotide_cmd = $dinucleotide_cmd_part1 . "\n\n" . $dinucleotide_cmd_part2;
					$required_out = $dinucleotide_short;
					}

				# check if this should be run
				if ('Y' eq missing_file($required_out)) {

					# record command (in log directory) and then run job
					print $log " >> Submitting job for Dinucleotide Profiling...\n";

					$run_script = write_script(
						log_dir	=> $log_directory,
						name	=> 'run_dinucleotide_profiling_' . $tumour,
						cmd	=> $dinucleotide_cmd,
						modules	=> [$samtools, $picard, $bedtools, $r_version],
						dependencies	=> $dedup_run_id,
						max_time	=> $parameters->{dinucleotide}->{time},
						mem		=> $parameters->{dinucleotide}->{mem},
						hpc_driver	=> $args{hpc_driver},
						extra_args	=> [$hpc_group]
						);

					$run_id = submit_job(
						jobname		=> 'run_dinucleotide_profiling_' . $tumour,
						shell_command	=> $run_script,
						hpc_driver	=> $args{hpc_driver},
						dry_run		=> $args{dry_run},
						log_file	=> $log
						);

					push @patient_jobs, $run_id;
					push @all_jobs, $run_id;
					} else {
					print $log " >> Skipping Dinucleotide profiling because this has already been completed!\n";
					}

				push @final_outputs, $required_out;
				}
			}

		# should intermediate files be removed
		# run per patient
		if ($args{del_intermediates}) {

			print $log ">> Saving command to clean up temporary/intermediate files...\n";

			# make sure final output exists before removing intermediate files!
			$cleanup_cmd = join("\n",
				"if [ -s " . join(" ] && [ -s ", @final_outputs) . " ]; then",
				"$cleanup_cmd",
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

		print $log "\nFINAL OUTPUT:\n" . join("\n  ", @final_outputs) . "\n";
		}

	# collate results
	if (scalar(@all_jobs) > 0) {

		my $collect_output = join(' ',
			"Rscript $cwd/scripts/fragmentomics/collect_fragmentomics_output.R",
			'-d', $output_directory,
			'-p', $tool_data->{project_name}
			);

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'collect_fragmentomics_outputs',
			cmd	=> $collect_output,
			modules	=> ['R/4.1.0'],
			dependencies	=> join(':', @all_jobs),
			mem		=> '6G',
			max_time	=> '24:00:00',
			hpc_driver	=> $args{hpc_driver},
			extra_args	=> [$hpc_group]
			);

		$run_id = submit_job(
			jobname		=> 'collect_fragmentomics_outputs',
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
my ($tool_config, $data_config, $output_directory);
my $hpc_driver = 'slurm';
my ($dry_run, $help, $no_wait, $remove_intermediates);

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	'd|data=s'	=> \$data_config,
	't|tool=s'	=> \$tool_config,
	'o|out_dir=s'	=> \$output_directory,
	'c|cluster=s'	=> \$hpc_driver,
	'remove'	=> \$remove_intermediates,
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
	del_intermediates	=> $remove_intermediates,
	dry_run			=> $dry_run,
	no_wait			=> $no_wait
	);
