#!/usr/bin/env perl
### filter_germline_variants.pl ####################################################################
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
# 1.0		sprokopec       script to filter germline/somatic SNVs/INDELs from haplotypecaller
# 1.1		sprokopec	added help message

### USAGE ##########################################################################################
# filter_germline_variants.pl -v input.vcf -r REFERENCE -o /path/to/output/stem -t TUMOUR1,TUMOUR2 -n NORMAL1
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
		tumours	=> undef,
		normals	=> undef,
		@_
		);

	my (@tumour_samples, @normal_samples);
	if (defined($args{tumours})) {
		@tumour_samples = split /,/, $args{tumours};
		}
	if (defined($args{normals})) {
		@normal_samples = split /,/, $args{normals};
		}

	### RUN ###########################################################################################

	my ($variants);
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
	my (@normal_idx, @tumour_idx, $genotype);

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
			foreach my $norm (@normal_samples) {
				push @normal_idx, first_index { $_ eq $norm } @info;
				}
			foreach my $tum (@tumour_samples) {
				push @tumour_idx, first_index { $_ eq $tum } @info;
				}

			push @metadata, join("\t", @info[0..8], @info[@normal_idx], @info[@tumour_idx]);
			next;
			}

		# for all other lines:
		# if this variant didn't pass the filtering, skip it
		my $filter = $info[6];
		if ($filter ne 'PASS') { next; }

		# if the ALT allele is ambiguous (multiple options), skip it
		my $chr = $info[0];
		my $pos = $info[1];
		my $id  = $info[2];
		my $ref = $info[3];
		my $alt = $info[4];

		if ($alt =~ m/,/) { next; }

		# determine the variant type
		my $var_type = 'SNV';
		if ( (length($ref) > 1) || (length($alt) > 1)) { $var_type = 'INDEL'; }

		# extract key fields:
		my @format_elements = split /\:/, $info[8];
		my $gt_idx = first_index { $_ eq 'GT' } @format_elements;
		my $dp_idx = first_index { $_ eq 'DP' } @format_elements;

		my (@norm_gts, @tumour_gts, @norm_dp, @tumour_dp);
		my @entries;

		foreach my $idx (@normal_idx) {
			@entries = split /\:/, $info[($idx)];
			push @norm_gts, $entries[$gt_idx];
			push @norm_dp, $entries[$dp_idx];
			}
		foreach my $idx (@tumour_idx) {
			@entries = split /\:/, $info[($idx)];
			push @tumour_gts, $entries[$gt_idx];
			push @tumour_dp, $entries[$dp_idx];
			}

		# other filtering criteria:
		# for T/N pairs
		if ( (scalar(@normal_samples) > 0) && (scalar(@tumour_samples) > 0) ) {

			# skip, if ANY sample could not be called at this position
			if ( (any { $_ =~ m/\.\/\./ } @norm_gts) || (any { $_ =~ m/\.\/\./ } @tumour_gts) ) { next; }

			# or if depth is below a threshold (ie, 10 reads) for ANY sample
			if ( (any { $_ < 10 } @norm_dp) || (any { $_ < 10 } @tumour_dp) ) { next; }

			# or if clearly not a germline variant (ie, normal == 0/0) or multiple normals don't match
			if ( (any { $_ =~ m/0\/0/ } @norm_gts) ) { next; }
			if ( (any { $_ ne $norm_gts[0] } @norm_gts) ) { next; }

			# or if all genotypes are the reference allelle
			if ( (all { $_ =~ m/0\/0/ } @norm_gts) && (all { $_ =~ m/0\/0/ } @tumour_gts) ) { next; }

		# for tumour ONLY samples
		} elsif ( (scalar(@normal_samples) == 0) && (scalar(@tumour_samples) > 0) ) {

			# skip, if ANY sample could not be called at this position
			if ( (any { $_ =~ m/\.\/\./ } @tumour_gts) ) { next; }

			# or if depth is below a threshold (ie, 10 reads) for ANY sample
			if ( (any { $_ < 20 } @tumour_dp) ) { next; }

			# or if all gentoypes are the reference allele
			if ( (all { $_ =~ m/0\/0/ } @tumour_gts) ) { next; }

		# for normal ONLY samples (just in case)
		} elsif ( (scalar(@normal_samples) > 0) && (scalar(@tumour_samples) == 0) ) {

			# skip, if ANY sample could not be called at this position
			if ( (any { $_ =~ m/\.\/\./ } @norm_gts) ) { next; }

			# or if depth is below a threshold (ie, 10 reads) for ANY sample
			if ( (any { $_ < 20 } @norm_dp) ) { next; }

			# or if clearly not a germline variant (ie, normal == 0/0) or multiple normals don't match
			if ( (any { $_ ne $norm_gts[0] } @norm_gts) ) { next; }

			# or if all gentoypes are the reference allele
			if ( (all { $_ =~ m/0\/0/ } @norm_gts) ) { next; }

			}

		# if we make it this far, then it will be written as germline
		# or high confidence call (for tumour only)
		$variants->{$chr}->{$pos} = join("\t", @info[0..8], @info[@normal_idx], @info[@tumour_idx]);
		}

	close $vcf_fh;

	# now prepare to write these to file
	my $vcf_header = join("\n", @metadata);

	open (my $fh_output, '>', $output_file) or die "Could not open $output_file\n";

	print $fh_output $vcf_header . "\n";

	# to ensure output is sorted, retrieve the reference index
	my $index_file = $args{reference} . '.fai';
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
	}

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($input_vcf, $tumour_ids, $normal_ids, $output_stem, $reference, $help);

# get command line arguments
GetOptions(
	'h|help'		=> \$help,
	'v|vcf=s'		=> \$input_vcf,
	't|tumour=s'		=> \$tumour_ids,
	'n|normal=s'		=> \$normal_ids,
	'o|output-stem=s'	=> \$output_stem,
	'r|reference=s'		=> \$reference
	);

if ($help) {
	my $help_msg = join("\n",
		"Options:",
		"\t--help|-h\tPrint this help message",
		"\t--vcf|-v\t<string> input vcf file",
		"\t--tumour|-t\t<string> comma separated list of tumour IDs",
		"\t--normal|-n\t<string> comma separated list of normal IDs",
		"\t--output-stem|-o\t<string> filestem for output files",
		"\t--reference|-r\t<string> path to reference used for alignments"
		);

	print $help_msg;
	exit;
	}

# do some quick error checks to confirm valid arguments	
main(
	vcf		=> $input_vcf,
	tumours		=> $tumour_ids,
	normals		=> $normal_ids,
	output		=> $output_stem,
	reference	=> $reference
	);
