#!/usr/bin/env perl
### filter_novobreak_variants.pl ###################################################################
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
# 1.0		sprokopec       script to filter structural variants from NovoBreak

### USAGE ##########################################################################################
# filter_novobreak_variants.pl -v input.vcf -r REFERENCE -o /path/to/output/stem -t TUMOUR1,TUMOUR2 -n NORMAL1
#
# where:
#	-v (input.vcf) indicates input file to use (recalibrated gVCF from genotype_gvcfs.pl)
#	-o indicates the output stem (ie, /path/to/output/SAMPLE
#	-r REFERENCE indicates path to your reference fasta
#	-t indicates tumour ids to search for in header of input vcf
#	-n indicates normal ids to search for in header of input vcf

### DEFINE SUBROUTINES #############################################################################


### MAIN ###########################################################################################
sub main {
	my %args = (
		vcf	=> undef,
		output	=> undef,
		tumour	=> undef,
		normal	=> undef,
		ref	=> undef,
		@_
		);

	### RUN ###########################################################################################

	my ($variants, $stats);
	my $output_file = $args{output};

	# open the genotyped/recalibrated vcf file
	my $vcf_fh;
	if ($args{vcf} =~ /.gz$/) {
		open($vcf_fh, "gunzip -c $args{vcf} |") or die "Could not open pipe to $args{vcf}\n";
		}
	else {
		open($vcf_fh, $args{vcf}) or die "Could not open $args{vcf}\n";
		}

	# being checking each line
	my @metadata;
	my @positions;

	while (<$vcf_fh>) {
		my $line = $_;
		chomp($line);

		# if this line is a header line, store it and move on
		if ($line =~ /^##/) {
			push @metadata, $line;
			next;
			}

		my @info = split /\t/, $line;

		# if this line contains the field headings
		if ($line =~ /^#CHROM/) {
			push @metadata, join("\t", @info[0..9]);
			next;
			}

		# for all other lines:
		# apply quality threshold filtering, skip it
		my $qual = $info[5];
		if ($qual < 10) { next; }

		# find position details
		my $chr = $info[0];
		my $pos = $info[1];
		my $pos_to_check = $chr . '_' . $pos;

		# if the position has already been seen, skip it
		if (any { /$pos_to_check/ } @positions) {
			my @old_info = split /\t/, $variants->{$chr}->{$pos};
			next if ($qual <= $old_info[5]);
			} else {
			push @positions, $pos_to_check;
			}

		# if we make it this far, then it will be written to file
		$stats->{$chr}->{$pos} = join("\t", @info);
		$variants->{$chr}->{$pos} = join("\t", @info[0..9]);
		}

	close $vcf_fh;

	### VCF OUTPUT ###
	# now prepare to write these to file
	my $vcf_header = join("\n", @metadata);

	open (my $fh_output, '>', $output_file) or die "Could not open $output_file\n";

	print $fh_output $vcf_header . "\n";

	# to ensure output is sorted, retrieve the reference index
	my $index_file = $args{ref} . '.fai';
	open (my $fh_index, $index_file) or die "Could not open reference index $index_file\n";

	while (<$fh_index>) {

		my $line = $_;
		my @info_chr = split /\t/, $line;
		my $chr = $info_chr[0];

		foreach my $pos (sort {$a <=> $b} keys (%{$variants->{$chr}})) {
			print $fh_output $variants->{$chr}->{$pos} . "\n";
			}
		}	

	close $fh_output;

	### STATS OUTPUT
	# now prepare to write these to file
	my @nb_fields = qw(tumour normal left_chrom left_position right_chrom right_position sv_type sv_length connection_type consensus_sequence quality reads_used_for_assembly average_coverage tumor_bkpt1_depth tumor_bkpt1_sp_reads tumor_bkpt1_qual tumor_bkpt1_high_qual_sp_reads tumor_bkpt1_high_qual_qual normal_bkpt1_depth normal_bkpt1_sp_reads normal_bkpt1_qual normal_bkpt1_high_qual_sp_reads normal_bkpt1_high_qual_qual tumor_bkpt2_depth tumor_bkpt2_sp_reads tumor_bkpt2_qual tumor_bkpt2_high_qual_sp_reads tumor_bkpt2_high_qual_qual normal_bkpt2_depth normal_bkpt2_sp_reads normal_bkpt2_qual normal_bkpt2_high_qual_sp_reads normal_bkpt2_high_qual_qual tumor_bkpt1_discordant_reads normal_bkpt1_discordant_reads tumor_bkpt2_discordant_reads normal_bkpt2_discordant_reads);

	my $header = join("\t", @nb_fields);
	$output_file =~ s/.vcf/.tsv/;

	open (my $stats_output, '>', $output_file) or die "Could not open $output_file\n";

	print $stats_output $header . "\n";

	# to ensure output is sorted, retrieve the reference index
	open (my $stats_index, $index_file) or die "Could not open reference index $index_file\n";

	while (<$stats_index>) {

		my $line = $_;
		my @info_chr = split /\t/, $line;
		my $chr = $info_chr[0];

		foreach my $pos (sort {$a <=> $b} keys (%{$stats->{$chr}})) {

			my @stats_info = split /\t/, $stats->{$chr}->{$pos};
			my @info_fields = split /\;/, $stats_info[7];

			# split INFO field and extract key values
			my $chr2_idx   = first_index { /^CHR2=/ } @info_fields;
			my $pos2_idx   = first_index { /^END=/ } @info_fields;
			my $type_idx   = first_index { /^SVTYPE=/ } @info_fields;
			my $length_idx = first_index { /^SVLEN=/ } @info_fields;
			my $seq_idx    = first_index { /^CONSENSUS=/ } @info_fields;
			my $dir_idx    = first_index { /^CT=/ } @info_fields;

			my $chr2 = $info_fields[$chr2_idx];
			$chr2 =~ s/CHR2=//;

			my $pos2 = $info_fields[$pos2_idx];
			$pos2 =~ s/END=//;

			my $sv_type = $info_fields[$type_idx];
			$sv_type =~ s/SVTYPE=//;

			my $sv_length = $info_fields[$length_idx];
			$sv_length =~ s/SVLEN=//;

			my $sequence = $info_fields[$seq_idx];
			$sequence =~ s/CONSENSUS=//;

			my $direction = $info_fields[$dir_idx];
			$direction =~ s/CT=//;

			my $coverage = $stats_info[14];
			$coverage =~ s/cov//;

			# combine for output
			my $keep_info = join("\t",
				$args{tumour}, $args{normal},
				@stats_info[0..1],
				$chr2, $pos2, $sv_type, $sv_length, $direction, $sequence,
				$stats_info[5],
				$stats_info[13],
				$coverage,
				@stats_info[15..$#stats_info]
				);

			print $stats_output $keep_info . "\n";
			}
		}	

	close $stats_output;
	}

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($input_vcf, $tumour_id, $normal_id, $output, $reference, $help);

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	'v|vcf=s'	=> \$input_vcf,
	't|tumour=s'	=> \$tumour_id,
	'n|normal=s'	=> \$normal_id,
	'o|output=s'	=> \$output,
	'r|reference=s'	=> \$reference
	);

if ($help) {
	my $help_msg = join("\n",
		"Options:",
		"\t--help|-h\tPrint this help message",
		"\t--vcf|-v\t<string> input vcf file",
		"\t--tumour|-t\t<string> tumour ID",
		"\t--normal|-n\t<string> normal ID",
		"\t--output|-o\t<string> path to output file",
		"\t--reference|-r\t<string> path to reference used for alignments"
		);

	print $help_msg . "\n";
	exit;
	}

# do some quick error checks to confirm valid arguments	
main(
	vcf	=> $input_vcf,
	tumour	=> $tumour_id,
	normal	=> $normal_id,
	output	=> $output,
	ref	=> $reference
	);
