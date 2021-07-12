#!/usr/bin/env perl
### write_rna_methods.pl ###########################################################################
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
	$methods .= "\\subsection{Alignment and Quality Checks:}\n";

	my ($star, $rnaseqc, $gatk, $rsem, $star_fusion, $fusioncatcher);
	my ($samtools, $picard, $bedtools, $vcftools);
	my ($star_ref, $ref_type, $gtf);
	my ($k1000g, $mills, $kindels, $dbsnp, $hapmap, $omni, $cosmic);
	my ($vep, $vcf2maf, $bwa);

	$ref_type	= $tool_data->{ref_type};

	# how was STAR run?
	if ('Y' eq $tool_data->{star}->{run}) {

		$star		= $tool_data->{star_version};
		$rnaseqc	= $tool_data->{rna_seqc_version};
		$picard		= $tool_data->{picard_version};
		$bwa		= $tool_data->{bwa_path};

		my @parts	= split('\\/', $tool_data->{star_reference_dir});
		@parts		= split($star . '_', $parts[-1]);
		$star_ref	= $parts[-1];
		@parts		= split('\\/',$tool_data->{reference_gtf});
		$gtf		= $parts[-1];
		@parts		= split('\\/', $bwa);
		$bwa		= $parts[-1];

		$methods .= "Fastq files were aligned to the $ref_type transcriptome reference using STAR (v$star), with the following arguments: --twopassMode Basic, --outSAMtype BAM SortedByCoordinate, --outSAMunmapped Within, --outSAMprimaryFlag AllBestScore, --outFilterIntronMotifs RemoveNoncanonical, --alignSJDBoverhangMin 10, --alignMatesGapMax 100000, --alignIntronMax 100000, --alignSJstitchMismatchNmax 5 -1 5 5, --peOverlapNbasesMin 12, --peOverlapMMp 0.1, --chimOutType WithinBAM, --chimSegmentMin 10, --chimScoreJunctionNonGTAG -4, --chimMultimapNmax 20, --chimMultimapScoreRange 3, --chimNonchimScoreDropMin 10, --chimOutJunctionFormat 1, --chimJunctionOverhangMin 10, --quantMode GeneCounts TranscriptomeSAM. Where multiple fastq files were present (ie, multiple lanes), files were input together as a single alignment run.\\newline\n";
		$methods .= join("\n",
			"{\\scriptsize \\begin{itemize}",
			"  \\vspace{-0.2cm}\\item Reference: $star_ref",
			"\\end{itemize} }"
			) . "\n";

		$methods .= "\\noindent\nDuplicate reads were marked in the Aligned.sortedByCoord.out.bam file using Picard tools (v$picard). RNASeQC (v$rnaseqc) was run on each batch (cohort), using the duplicate-marked BAMs and $ref_type (genome) reference, with the following arguments: -bwa $bwa, -t $gtf, -singleEnd no.\\newline\n";
		} else {
		$methods .= "BWA not run.\\newline\n";
		}

	# how was GATK run?
	if ('Y' eq $tool_data->{gatk}->{run}) {

		$gatk = $tool_data->{gatk_version};

		# find reference files
		if ('hg38' eq $tool_data->{ref_type}) {

			$k1000g		= 'hg38bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz';
			$kindels	= 'hg38bundle/Homo_sapiens_assembly38.known_indels.vcf.gz';
			$mills		= 'hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz';
			$dbsnp		= 'hg38bundle/dbsnp_144.hg38.vcf.gz';

			} elsif ('hg19' eq $tool_data->{ref_type}) {

			$k1000g		= '1000G_phase1.snps.high_confidence.hg19.vcf';
			$kindels	= '1000G_phase1.indels.hg19.vcf';
			$mills		= 'Mills_and_1000G_gold_standard.indels.hg19.vcf';
			$dbsnp		= 'dbsnp_138.hg19.vcf';

			}

		if (defined($tool_data->{dbsnp})) {
			my @parts = split('\\/', $tool_data->{dbsnp});
			$dbsnp = $parts[-1];
			}

		$methods .= "\\noindent\nIndel realignment and base-quality recalibration were performed for each patient using GATK (v$gatk). First, SplitNCigarReads was run on duplicate-marked BAMs, using the $ref_type (genome) reference, ReassignOneMappingQuality (from 255 to 60) and UnmappedRead read filters, and -U ALLOW_N_CIGAR_READS. Known indels were provided for indel realignment target creation, and known snps provided for recalibration.\\newline\n";
		$methods .= join("\n",
			"{\\scriptsize \\begin{itemize}",
			"  \\vspace{-0.2cm}\\item Known INDELs: $mills",
			"  \\vspace{-0.2cm}\\item Known INDELs: $kindels",
			"  \\vspace{-0.2cm}\\item Known SNPs: $k1000g",
			"  \\vspace{-0.2cm}\\item Known SNPs: $dbsnp",
			"\\end{itemize} }"
			) . "\n";
		} else {
		$methods .= "GATK not run.\\newline\n";
		}

	# how was haplotypecaller run?
	if ('Y' eq $tool_data->{haplotype_caller}->{run}) {

		$methods .= "\\subsection{Variant calling:}\n";

		$gatk = $tool_data->{gatk_version};

		my @parts = split('\\/', $tool_data->{haplotype_caller}->{parameters}->{annotate}->{vep_path});
		$vep = $parts[-1];
		@parts = split('\\/', $tool_data->{haplotype_caller}->{parameters}->{annotate}->{vcf2maf_path});
		$vcf2maf = $parts[-2];

		$methods .= "Short variants (SNVs and Indels) were identified using GATK's (v$gatk) HaplotypeCaller as per GATK's RNA-Seq variant calling best practices. HaplotypeCaller was run using a minimum confidence threshold of 20 and -dontUseSoftClippedBases. VariantFiltration was used to remove low quality variants (QD \$<2.0\$ and FS \$>30\$) with a cluster size of 3 and cluster window size of 35.\\newline\n";
		$methods .= "\n\\noindent\nVariants were filtered and annotated using VEP (v$vep) and vcf2maf ($vcf2maf). Filters were applied to remove known common variants (ExAC nonTCGA version r1) and variants with coverage below 20x.\\newline\n";
		} else {
		$methods .= "No variant calling performed.\\newline\n";
		}

	# how were fusions called?
	$methods .= "\\subsection{Fusion Detection:}\n";

	if ('Y' eq $tool_data->{fusioncatcher}->{run}) {

		$fusioncatcher	= $tool_data->{fusioncatcher_version};

		# fill in methods
		$methods .= "\\subsubsection{Fusioncatcher (v$fusioncatcher):}\n";
		if ('hg38' eq $tool_data->{ref_type}) {
			$methods .= "Fusioncatcher was run on fastq files, using the GRCh38 reference, --skip-conversion-grch37 and --keep-viruses-alignments arguments.";
			} else {
			$methods .= "Fusioncatcher was run on fastq files, using the GRCh37 reference and --keep-viruses-alignments.";
			}
		$methods .= " Where multiple fastq files were present (ie, multiple lanes), files were input together as a single run.\\newline\n";
		} else {
		$methods .= "Fusioncatcher not run.\\newline\n";
		}

	if ('Y' eq $tool_data->{star_fusion}->{run}) {

		$star_fusion = $tool_data->{star_fusion_version};
		my $fusion_inspect = $tool_data->{star_fusion}->{parameters}->{FusionInspect};

		# fill in methods
		$methods .= "\\subsubsection{STAR-Fusion (v$star_fusion):}\n";
		$methods .= "STAR-fusion was run using the Chimeric.out.junction output by STAR (see above).";
		if (defined($fusion_inspect)) {
			$methods .= " FusionInspector ($fusion_inspect mode) was used to further validate detected fusions.";
			}
		$methods .= "\\newline\n";
		} else {
		$methods .= "STAR-Fusion was not run.\\newline\n";
		}

	# how was gene expression called?
	$methods .= "\\subsection{Gene Expression:}\n";

	if ('Y' eq $tool_data->{rsem}->{run}) {

		$rsem		= $tool_data->{rsem_version};
		my $strand	= $tool_data->{rsem}->{strandedness};
		my $rsem_ref	= $tool_data->{rsem_reference};
		my @parts	= split('\\/', $rsem_ref);
		$rsem_ref	= $parts[-1];

		# fill in methods
		$methods .= "\\textbf{RSEM (v$rsem):}\\newline\n";
		$methods .= "RSEM was run on the Aligned.toTranscriptome.out.bam file generated by STAR (see above) using rsem-calculate-expression with the RSEM $rsem_ref reference and the following arguments: --paired-end, --bam, --estimate-rspd, --output-genome-bam and --strandedness $strand.\\newline\n";
		} else {
		$methods .= "RSEM not run.\\newline\n";
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
		"\t--tool|-t\t<string> Master config file for the DNA pipeline",
		"\t--directory|-d\t<string> path to output directory"
		);

	print $help_msg . "\n";
	exit;
	}

# do some quick error checks to confirm valid arguments	
main(config => $config, directory => $directory);
