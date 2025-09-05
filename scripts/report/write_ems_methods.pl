#!/usr/bin/env perl
### write_ems_methods.pl ###########################################################################
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
use Data::Dumper;

my $cwd = dirname(__FILE__);
require "$cwd/../utilities.pl";

####################################################################################################
# version       author		comment
# 1.0		sprokopec       tool to automatically generate reports

### MAIN ###########################################################################################
sub main {
	my %args = (
		config		=> undef,
		directory	=> undef,
		@_
		);

	my $tool_data = LoadFile($args{config});

	### RUN ####################################################################################
	# for each tool (indicated in config file), read in and extract parameters
	my $methods = "\\section{Methods}\n";
	$methods .= "For all tools, default parameters were used unless otherwise indicated.\\newline\n";
	$methods .= "\\subsection{Alignment and Quality Checks}\n";

	my ($bwa, $gatk4, $methyldackel);
	my ($ref_type, $samtools, $picard, $bedtools, $vcftools, $bcftools, $sambamba);
	my ($k1000g, $mills, $kindels, $dbsnp, $hapmap, $omni, $cosmic, $pon, $gnomad);
	my ($vep, $vcf2maf);

	# check which tools have been requested
	my %tool_set = (
		'trim'	=> defined($tool_data->{trim_adapters}->{run}) ? $tool_data->{trim_adapters}->{run} : 'N',
		'bwa'	=> defined($tool_data->{bwa}->{run}) ? $tool_data->{bwa}->{run} : 'N',
		'bamqc'	=> defined($tool_data->{bamqc}->{run}) ? $tool_data->{bamqc}->{run} : 'N',
		'methyldackel'	=> defined($tool_data->{methyldackel}->{run}) ? $tool_data->{methyldackel}->{run} : 'N'
		);

	# how was BWA run?
	if ('Y' eq $tool_set{'bwa'}) {

		$bwa		= $tool_data->{bwa_meth_version};
		$samtools	= $tool_data->{samtools_version};
		$picard		= $tool_data->{picard_version};
		$sambamba	= $tool_data->{sambamba_version};
		$ref_type	= $tool_data->{ref_type};

		$methods .= "Fastq files were aligned to $ref_type using the BWA-METH algorithm (v$bwa). Resulting SAM files were coordinate sorted, converted to BAM format, indexed and filtered using using samtools (v$samtools).";

		if ('N' eq $tool_data->{bwa}->{parameters}->{merge}->{mark_dup}) {
			$methods .= "Lane- and library-level BAMs were not merged or duplicate marked.\n";
			} elsif ('sambamba' eq $tool_data->{bwa}->{parameters}->{merge}->{tool}) {
			$methods .= " Duplicate reads were marked and lane- and library-level BAMs were merged using Sambamba (v$sambamba).\\newline\n";
			} elsif ('picard' eq $tool_data->{bwa}->{parameters}->{merge}->{tool}) {
			$methods .= " Duplicate reads were marked and lane- and library-level BAMs were merged using Picard tools (v$picard).\\newline\n";
			}

		$methods .= "\\newline\n";
		} else {
		$methods .= "BWA not run.\\newline\n";
		}

	# how was QC run?
	if ('Y' eq $tool_set{'bamqc'}) {

		$picard = $tool_data->{picard_version};

		$methods .= "\\noindent\nPicard's ($picard) CollectAlignmentSummaryMetrics, CollectGcBiasMetrics and CollectWgsMetrics (or CollectHsMetrics) were run to obtain various alignment metrics.\\newline\n";
		} else {
		$methods .= "BAM quality checks (ContEst, coverage and callable bases) not performed.\\newline\n";
		}


	# Methyldackel
	if ('Y' eq $tool_set{'methyldackel'}) {

		$methyldackel = $tool_data->{methyldackel_version};

		my $sites = 'CpG sites';
		if ('Y' eq $tool_data->{methyldackel}->{parameters}->{extract}->{non_cpg}) {
			$sites = 'CpG|CHG|CHH sites';
			}

		# fill in methods
		$methods .= "\\subsubsection{MethylDackel (v$methyldackel)}\n";

		$methods .= "Methylation estimates were obtained for all $sites using MethylDackel extract with recommended paramters. MBIAS was evaluated but is not incorporated by default. Methylation levels (% of reads) for each CpG loci were determined using the bsseq package in R.\n";
		}

	# clean up special characters
	$methods =~ s/_/\\_/g;

	# write methods to file
	open(my $methods_file, '>', $args{directory} . "/methods.tex");
	print $methods_file <<EOI;
$methods
EOI
	close($methods_file);

	}
 
### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($help, $config, $directory);

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	't|tool=s'	=> \$config,
	'd|directory=s'	=> \$directory
	);

if ($help) {
	my $help_msg = join("\n",
		"Options:",
		"\t--help|-h\tPrint this help message",
		"\t--tool|-t\t<string> Master config file for the EM-Seq pipeline",
		"\t--directory|-d\t<string> path to output directory"
		) . "\n";

	print $help_msg;
	exit;
	}

# do some quick error checks to confirm valid arguments	
main(config => $config, directory => $directory);
