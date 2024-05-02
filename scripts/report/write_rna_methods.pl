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
	$methods .= "\\subsection{Alignment and Quality Checks}\n";

	my ($star, $rnaseqc, $gatk, $rsem, $star_fusion, $arriba, $fusioncatcher, $mavis);
	my ($samtools, $picard);
	my ($star_ref, $ref_type, $gtf);
	my ($k1000g, $mills, $kindels, $dbsnp, $hapmap, $omni, $cosmic);
	my ($vep, $vcf2maf, $bwa);

	# check which tools have been requested
	my %tool_set = (
		'star'		=> defined($tool_data->{star}->{run}) ? $tool_data->{star}->{run} : 'N',
		'gatk'		=> defined($tool_data->{gatk}->{run}) ? $tool_data->{gatk}->{run} : 'N',
		'haplotype_caller' => defined($tool_data->{haplotype_caller}->{run}) ? $tool_data->{haplotype_caller}->{run} : 'N',
		'rsem'		=> defined($tool_data->{rsem}->{run}) ? $tool_data->{rsem}->{run} : 'N',
		'star_fusion'	=> defined($tool_data->{star_fusion}->{run}) ? $tool_data->{star_fusion}->{run} : 'N',
		'arriba'	=> defined($tool_data->{arriba}->{run}) ? $tool_data->{arriba}->{run} : 'N',
		'fusioncatcher'	=> defined($tool_data->{fusioncatcher}->{run}) ? $tool_data->{fusioncatcher}->{run} : 'N',
		'mavis'		=> defined($tool_data->{mavis}->{run}) ? $tool_data->{mavis}->{run} : 'N'
		);

	# extract tool versions and common variables
	$ref_type = $tool_data->{ref_type};
	$picard = defined($tool_data->{picard_version}) ? $tool_data->{picard_version} : undef;
	$samtools = defined($tool_data->{samtools_version}) ? $tool_data->{samtools_version} : undef;

	$star = defined($tool_data->{star_version}) ? $tool_data->{star_version} : undef;
	$rnaseqc = defined($tool_data->{rna_seqc_version}) ? $tool_data->{rna_seqc_version} : undef;
	$gatk = defined($tool_data->{gatk_version}) ? $tool_data->{gatk_version} : undef;
	$rsem = defined($tool_data->{rsem_version}) ? $tool_data->{rsem_version} : undef;
	$star_fusion = defined($tool_data->{star_fusion_version}) ? $tool_data->{star_fusion_version} : undef;
	$arriba = defined($tool_data->{arriba_version}) ? $tool_data->{arriba_version} : undef;
	$fusioncatcher = defined($tool_data->{fusioncatcher_version}) ? $tool_data->{fusioncatcher_version} : undef;
	$mavis = defined($tool_data->{mavis_version}) ? $tool_data->{mavis_version} : undef;

	# how was STAR run?
	if ('Y' eq $tool_set{'star'}) {

		my @parts	= split('\\/', $tool_data->{star_reference_dir});
		@parts		= split($star . '_', $parts[-1]);
		$star_ref	= $parts[-1];
		$gtf		= basename($tool_data->{reference_gtf});
		@parts		= split('\\/', $tool_data->{bwa_path});
		$bwa		= $parts[-1];

		$methods .= "Fastq files were aligned to the $ref_type transcriptome reference using STAR (v$star), with the following arguments: --twopassMode Basic, --outSAMtype BAM SortedByCoordinate, --outSAMunmapped Within, --outSAMprimaryFlag AllBestScore, --outFilterIntronMotifs RemoveNoncanonical, --alignSJDBoverhangMin 10, --alignMatesGapMax 100000, --alignIntronMax 100000, --alignSJstitchMismatchNmax 5 -1 5 5, --peOverlapNbasesMin 12, --peOverlapMMp 0.1, --chimOutType WithinBAM, --chimSegmentMin 10, --chimScoreJunctionNonGTAG -4, --chimMultimapNmax 20, --chimMultimapScoreRange 3, --chimNonchimScoreDropMin 10, --chimOutJunctionFormat 1, --chimJunctionOverhangMin 10, --quantMode GeneCounts TranscriptomeSAM. Where multiple fastq files were present (ie, multiple lanes), files were input together as a single alignment run.\\newline\n";
		$methods .= join("\n",
			"{\\scriptsize \\begin{itemize}",
			"  \\vspace{-0.2cm}\\item Reference: $star_ref",
			"\\end{itemize} }"
			) . "\n";

		$methods .= "\\noindent\nDuplicate reads were marked in the Aligned.sortedByCoord.out.bam file using Picard tools (v$picard). RNASeQC (v$rnaseqc) was run on each batch (cohort), using the duplicate-marked BAMs and $ref_type (genome) reference, with the following arguments: -bwa $bwa, -t $gtf, -singleEnd no.\\newline\n";
		} else {
		$methods .= "Alignment step using STAR was not run.\\newline\n";
		}

	# how was GATK run?
	if ('Y' eq $tool_set{'gatk'}) {

		# find reference files
		if ('hg38' eq $ref_type) {

			$k1000g		= 'hg38bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz';
			$kindels	= 'hg38bundle/Homo_sapiens_assembly38.known_indels.vcf.gz';
			$mills		= 'hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz';
			$dbsnp		= 'hg38bundle/dbsnp_144.hg38.vcf.gz';

			} elsif ('hg19' eq $ref_type) {

			$k1000g		= '1000G_phase1.snps.high_confidence.hg19.vcf';
			$kindels	= '1000G_phase1.indels.hg19.vcf';
			$mills		= 'Mills_and_1000G_gold_standard.indels.hg19.vcf';
			$dbsnp		= 'dbsnp_138.hg19.vcf';

			}

		if (defined($tool_data->{dbsnp})) {
			$dbsnp = basename($tool_data->{dbsnp});
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
	if ('Y' eq $tool_set{'haplotype_caller'}) {

		$methods .= "\\subsection{Variant calling}\n";

		my @parts = split('\\/', $tool_data->{haplotype_caller}->{parameters}->{annotate}->{vep_path});
		$vep = $parts[-1];

		if (defined($tool_data->{vcf2maf_version})) {
			$vcf2maf = $tool_data->{vcf2maf_version};
			} elsif (defined($tool_data->{haplotype_caller}->{parameters}->{annotate}->{vcf2maf_path})) {
			@parts = split('\\/', $tool_data->{haplotype_caller}->{parameters}->{annotate}->{vcf2maf_path});
			$vcf2maf = $parts[-2];
			}

		$methods .= "Short variants (SNVs and Indels) were identified using GATK's (v$gatk) HaplotypeCaller as per GATK's RNA-Seq variant calling best practices. HaplotypeCaller was run using a minimum confidence threshold of 20 and -dontUseSoftClippedBases. VariantFiltration was used to remove low quality variants (QD \$<2.0\$ and FS \$>30\$) with a cluster size of 3 and cluster window size of 35.\\newline\n";
		$methods .= "\n\\noindent\nVariants were filtered and annotated using VEP (v$vep) and vcf2maf ($vcf2maf). Filters were applied to remove known common variants (ExAC nonTCGA version r1) and variants with coverage below 20x.\\newline\n";
		}

	# how were fusions called?
	$methods .= "\\subsection{Fusion Detection}\n";

	if ('Y' eq $tool_set{'arriba'}) {

		# fill in methods
		$methods .= "For Arriba, fastq files were first aligned to the $ref_type transcriptome reference using STAR (v$star), with the following arguments: --outSAMtype BAM Unsorted, --outSAMunmapped Within, --outFilterFilterMultimapNmax 50, --peOverlapNbasesMin 10, --alignSplicedMateMapLminOverLmate 0.5, --alignSJstitchMismatchNmax 5 -1 5 5, --chimSegmentMin 10, --chimOutType WithinBAM HardClip, --chimJunctionOverhangMin 10, --chimScoreDropMax 30, --chimScoreJunctionNonGTAG 0, --chimScoreSeparation 1, --chimSegmentReadGapMax 3, --chimMultimapNmax 50. Where multiple fastq files were present (ie, multiple lanes), files were input together as a single alignment run. Arriba (v$arriba) was then run using default parameters with lists of known fusions, false fusions (blacklist) and protein domains provided ($ref_type). Viral expression was also evaluated using Arriba's quantify\_virus\_expression.sh.\\newline\\newline\n";

		} else {
		$methods .= "Arriba was not run.\\newline\n";
		}

	if ('Y' eq $tool_set{'fusioncatcher'}) {

		# fill in methods
		if ('hg38' eq $ref_type) {
			$methods .= "Fusioncatcher (v$fusioncatcher) was run on fastq files, using the GRCh38 reference, --skip-conversion-grch37 and --keep-viruses-alignments arguments.";
			} else {
			$methods .= "Fusioncatcher (v$fusioncatcher) was run on fastq files, using the GRCh37 reference and --keep-viruses-alignments.";
			}
		$methods .= " Where multiple fastq files were present (ie, multiple lanes), files were input together as a single run.\\newline\n\\newline\n";
		} else {
		$methods .= "Fusioncatcher was not run.\\newline\n";
		}

	if ('Y' eq $tool_set{'star_fusion'}) {

		my $fusion_inspect = $tool_data->{star_fusion}->{parameters}->{FusionInspect};

		# fill in methods
		$methods .= "STAR-fusion (v$star_fusion) was run using the Chimeric.out.junction output by STAR (see above).";
		if (defined($fusion_inspect)) {
			$methods .= " FusionInspector ($fusion_inspect mode) was used to further validate detected fusions.";
			}
		$methods .= "\\newline\n\\newline\n";
		} else {
		$methods .= "STAR-Fusion was not run.\\newline\n";
		}

	if ('Y' eq $tool_set{'mavis'}) {
		$methods .= "Mavis (v$mavis) was run once for each patient, using available SV calls, with the $ref_type reference files provided by the developers and BWA indicated as the aligner to use.\\newline\n";
		} else {
		$methods .= "Mavis was not run.\\newline\n";
		}

	# how was gene expression called?
	$methods .= "\\subsection{Gene Expression}\n";

	if ('Y' eq $tool_set{'rsem'}) {
		my $strand	= $tool_data->{rsem}->{strandedness};
		my $rsem_ref	= basename($tool_data->{rsem_reference});

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
