#!/usr/bin/env perl
### write_targetseq_methods.pl #####################################################################
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

	my ($gatk, $gatk4);
	my ($mutect, $mutect2, $strelka, $manta, $varscan, $delly, $mavis, $somaticsniper, $vardict, $novobreak);
	my ($ref_type, $samtools, $picard, $intervals, $bedtools, $vcftools, $bcftools);
	my ($k1000g, $mills, $kindels, $dbsnp, $hapmap, $omni, $cosmic, $pon, $gnomad);
	my ($vep, $vcf2maf);

	# how was QC run?
	my $threshold = '3.0';
	my $t_depth = '20x';
	my $n_depth = '15x';

	if ('Y' eq $tool_data->{bamqc}->{run}) {

		$gatk = $tool_data->{gatk_version};
		$picard = $tool_data->{picard_version};
		$gatk4 = $tool_data->{gatk_cnv_version};
		$bedtools = $tool_data->{bedtools_version};

		my @parts = split('\\/', $tool_data->{gnomad});
		$gnomad = $parts[-1];

		if (defined($tool_data->{bamqc}->{parameters}->{contest}->{threshold})) {
			$threshold = $tool_data->{bamqc}->{parameters}->{contest}->{threshold};
			}

		if (defined($tool_data->{bamqc}->{parameters}->{callable_bases}->{min_depth}->{tumour})) {
			$t_depth = $tool_data->{bamqc}->{parameters}->{callable_bases}->{min_depth}->{tumour};
			$t_depth .= 'x';
			}
		if (defined($tool_data->{bamqc}->{parameters}->{callable_bases}->{min_depth}->{normal})) {
			$n_depth = $tool_data->{bamqc}->{parameters}->{callable_bases}->{min_depth}->{normal};
			$n_depth .= 'x';
			}

		$methods .= "\\noindent\nGATK's ($gatk4) GetPileupSummaries and CalculateContamination were used to estimate contamination for each sample, based on common germline SNPs from gnomAD. It is recommended that tumours with a contamination estimate \$>$threshold\\%\$ be excluded from downstream analyses. GATK's ($gatk) DepthOfCoverage was used to assess genome coverage, again on both sample and readgroup levels (with -omitBaseOutput -omitIntervals -omitLocusTable), and callable bases defined as those with a minimum of $t_depth (tumour) or $n_depth (normal) coverage using bedtools (v$bedtools). Picard's ($picard) CollectAlignmentSummaryMetrics and CollectSequencingArtifactMetrics were run to obtain various alignment metrics.\\newline\n";
		$methods .= join("\n",
			"{\\scriptsize \\begin{itemize}",
			"  \\vspace{-0.2cm}\\item $gnomad",
			"\\end{itemize} }"
			) . "\n";
		} else {
		$methods .= "BAM quality checks not performed.\\newline\n";
		}

	# how was haplotypecaller run?
	if ('Y' eq $tool_data->{haplotype_caller}->{run}) {

		$methods .= "\\subsection{Germline Variant calling}\n";

		$gatk = $tool_data->{gatk_version};
		my $cpsr = $tool_data->{cpsr_version};
		my $pcgr = $tool_data->{pcgr_version};
		my @parts = split('\\/', $tool_data->{intervals_bed});
		$intervals = $parts[-1];

		# find reference files
		if (defined($tool_data->{dbsnp})) {
			@parts = split('\\/', $tool_data->{dbsnp});
			$dbsnp = $parts[-1];
			} elsif ('hg38' eq $tool_data->{ref_type}) {
			$dbsnp = 'dbsnp_144.hg38.vcf.gz';
			} elsif ('hg19' eq $tool_data->{ref_type}) {
			$dbsnp = 'dbsnp_138.hg19.vcf';
			}

		# find reference files
		if ('hg38' eq $tool_data->{ref_type}) {

			$k1000g	= 'hg38bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz';
			$hapmap	= 'hg38bundle/hapmap_3.3.hg38.vcf.gz';
			$omni	= 'hg38bundle/1000G_omni2.5.hg38.vcf.gz';
			$mills	= 'hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz';

			} elsif ('hg19' eq $tool_data->{ref_type}) {

			$k1000g	= '1000G_phase1.snps.high_confidence.hg19.vcf';
			$hapmap	= 'hapmap_3.3.hg19.vcf';
			$omni	= '1000G_omni2.5.hg19.vcf';
			$mills	= 'Mills_and_1000G_gold_standard.indels.hg19.vcf';

			}

		$methods .= "Germline variants were identified using GATK's (v$gatk) HaplotypeCaller and GenotypeGVCFs as per GATK's germline variant calling best practices. HaplotypeCaller was run in GVCF mode, using a minimum confidence threshold of 30, variant index type and parameter of LINEAR and 128000. Known variants (dbSNP) were provided, as were target intervals with interval padding set to 100bp. Variants were combined across samples (using CombineGVCFs) and genotyped using GenotypeGVCFs. Variant score recalibration was performed for INDELs and SNPs separately. For INDELs, known indels were provided, and recalibrater run using tranche sensitivity thresholds of 90, 99, 99.9 and 100 percent and maximum 2 Gaussians for the positive model. For SNPs, known snps and dbSNP were provided, and recalibrater run using maximum 4 Gaussians for the positive model. Recalibration was applied using truth sensitivity filter of 99 using the above generated .tranches and .recal files.\\newline\n";
		$methods .= "Germline variant significance was assessed using CPSR (v$cpsr) and PCGR (v$pcgr). Variants with a CSPR score above 0 were carried forward, as were any variants previously identified as significant within the TCGA cohort (Huang et al.).\\newline\n"; 
		$methods .= join("\n",
			"{\\scriptsize \\begin{itemize}",
			"  \\vspace{-0.2cm}\\item dbSNP: $dbsnp",
			"  \\vspace{-0.2cm}\\item Target regions: $intervals",
			"  \\vspace{-0.2cm}\\item Known INDELs: $mills", # (known=true,training=true,truth=true,prior=12.0)",
			"  \\vspace{-0.2cm}\\item Known SNPs: $hapmap", # (known=false,training=true,truth=true,prior=15.0)",
			"  \\vspace{-0.2cm}\\item Known SNPs: $omni", # (known=false,training=true,truth=true,prior=12.0)",
			"  \\vspace{-0.2cm}\\item Known SNPs: $k1000g", # (known=false,training=true,truth=false,prior=10.0)",
			"  \\vspace{-0.2cm}\\item Known SNPs: $dbsnp", # (known=true,training=false,truth=false,prior=2.0)",
			"\\end{itemize} }"
			) . "\n";
		} else {
		$methods .= "Germline variants not called.\\newline\n";
		}

	# how were somatic SNVs called?
	$methods .= "\\subsection{Somatic Variant Calling}\n";

	if ('Y' eq $tool_data->{mutect}->{run}) {

		$mutect		= $tool_data->{mutect_version};
		$gatk		= $tool_data->{gatk_version};
		$vcftools	= $tool_data->{vcftools_version};

		if (defined($tool_data->{dbsnp})) {
			my @parts = split('\\/', $tool_data->{dbsnp});
			$dbsnp = $parts[-1];
			} elsif ('hg38' eq $tool_data->{ref_type}) {
			$dbsnp = 'dbsnp_144.hg38.vcf.gz';
			} elsif ('hg19' eq $tool_data->{ref_type}) {
			$dbsnp = 'dbsnp_138.hg19.vcf';
			}

		if (defined($tool_data->{cosmic})) {
			my @parts = split('\\/', $tool_data->{cosmic});
			$cosmic = $parts[-2];
			}

		if (defined($tool_data->{mutect}->{pon})) {
			$pon = basename($tool_data->{mutect}->{pon});
			}

		# annotation
		my @parts = split('\\/', $tool_data->{annotate}->{vep_path});
		$vep = $parts[-1];
		@parts = split('\\/', $tool_data->{annotate}->{vcf2maf_path});
		$vcf2maf = $parts[-2];

		# fill in methods
		$methods .= "\\subsubsection{MuTect (v$mutect)}\n";

		if (defined($tool_data->{mutect}->{pon})) {
			$methods .= "A pre-made panel of normals ($pon) was used for additional filtering.\\newline\n";
			} else {
			$methods .= "A panel of normals was produced using all available normal samples, with MuTect's artifact\_detection\_mode. COSMIC ($cosmic; likely somatic [keep]) and dbSNP ($dbsnp; likely germline [remove]) were used as known lists, and variant calling restricted to target regions (with interval padding set to 100bp). Variants were merged across samples using GATK's CombineVariants (v$gatk), removing any variant that did not pass MuTect and GATK's default quality criteria (FILTER field not equal to PASS), and keeping variants present in at least 2 samples.\\newline\n";
			}

		$methods .= "For somatic calls, MuTect was run on T/N pairs, or tumour-only samples using identical methods. COSMIC ($cosmic; likely somatic [keep]), dbSNP ($dbsnp; likely germline [remove]) and the panel of normals were provided as known lists, and variant calling restricted to target regions (with interval padding set to 100bp). Lastly, calls were filtered (using vcftools v$vcftools) to again remove any that did not pass MuTect's filter criteria (FILTER did not equal PASS).\\newline\n";
		} else {
		$methods .= "MuTect not run.\\newline\n";
		}

	if ('Y' eq $tool_data->{mutect2}->{run}) {

		$gatk = $tool_data->{gatk_version};
		$vcftools = $tool_data->{vcftools_version};

		if (defined($tool_data->{dbsnp})) {
			my @parts = split('\\/', $tool_data->{dbsnp});
			$dbsnp = $parts[-1];
			} elsif ('hg38' eq $tool_data->{ref_type}) {
			$dbsnp = 'dbsnp_144.hg38.vcf.gz';
			} elsif ('hg19' eq $tool_data->{ref_type}) {
			$dbsnp = 'dbsnp_138.hg19.vcf';
			}

		if (defined($tool_data->{cosmic})) {
			my @parts = split('\\/', $tool_data->{cosmic});
			$cosmic = $parts[-2];
			}

		if (defined($tool_data->{mutect2}->{pon})) {
			$pon = basename($tool_data->{mutect2}->{pon});
			}

		# annotation
		my @parts = split('\\/', $tool_data->{annotate}->{vep_path});
		$vep = $parts[-1];
		@parts = split('\\/', $tool_data->{annotate}->{vcf2maf_path});
		$vcf2maf = $parts[-2];

		# fill in methods
		$methods .= "\\subsubsection{MuTect2 (GATK v$gatk)}\n";

		if (defined($tool_data->{mutect2}->{pon})) {
			$methods .= "A pre-made panel of normals ($pon) was used for additional filtering.\\newline\n";
			} else {
			$methods .= "A panel of normals was produced using all available normal samples, with MuTect2's artifact\_detection\_mode. dbSNP ($dbsnp; likely germline [remove]) was provided and variant calling restricted to target regions (with interval padding set to 100bp). Variants were merged across samples using GATK's CombineVariants (v$gatk), removing any variant that did not pass MuTect and GATK's default quality criteria (FILTER field not equal to PASS), and keeping variants present in at least 2 samples.\\newline\n";
			}

		$methods .= "For somatic calls, MuTect2 was run on T/N pairs, or tumour-only samples using identical methods. COSMIC ($cosmic; likely somatic [keep]) , dbSNP ($dbsnp; likely germline [remove]) and the panel of normals were provided as known lists, and variant calling restricted to target regions (with interval padding set to 100bp). Lastly, calls were filtered (using vcftools v$vcftools) to again remove any that did not pass MuTect2's filter criteria (FILTER did not equal PASS).\\newline\n";
		} else {
		$methods .= "MuTect2 not run.\\newline\n";
		}

	if ('Y' eq $tool_data->{somaticsniper}->{run}) {

		$somaticsniper = $tool_data->{somaticsniper_version};
		my ($tool, $somaticsniper_version) = split('/', $somaticsniper);
		$somaticsniper = $somaticsniper_version;
		$samtools = $tool_data->{samtools_version};
		$vcftools = $tool_data->{vcftools_version};

		if (defined($tool_data->{somaticsniper}->{pon})) {
			$pon = $tool_data->{somaticsniper}->{pon};
			}

		# annotation
		my @parts = split('\\/', $tool_data->{annotate}->{vep_path});
		$vep = $parts[-1];
		@parts = split('\\/', $tool_data->{annotate}->{vcf2maf_path});
		$vcf2maf = $parts[-2];

		# fill in methods
		$methods .= "\\subsubsection{SomaticSniper (v$somaticsniper)}\n";
		$methods .= "SomaticSniper was run on each T/N pair using options -q 1 -Q 40 -G and -L . For filtering, bcftools (v$samtools) mpileup was run with results filtered for quality using bcftools vcfutils.pl varFilter -Q 20. SomaticSniper's filters were then applied as suggested (snpfilter.pl was first applied to the initial VCF using the normal indel pileup, with the results then filtered using the tumour pileup). The resulting positions were fed into bam-readcount (-b 15 -q 1) for the tumour BAM and SomaticSnipers fpfilter.pl and highconfidence.pl applied to remove probable false positives (min mapping quality = 40 and min somatic score = 40). ";

		if (defined($tool_data->{somaticsniper}->{pon})) {
			$methods .= "Remaining variants were filtered to target regions (\$\\pm\$ 50bp padding) and germline variants removed using the provided panel of normals ($pon) using bcftools (v$samtools) and vcftools (v$vcftools).\\newline\n";
			} else {
			$methods .= "Remaining variants were filtered to target regions only (\$\\pm\$ 50bp padding) using bcftools (v$samtools).\\newline\n";
			}

		} else {
		$methods .= "SomaticSniper not run.\\newline\n";
		}

	if ('Y' eq $tool_data->{vardict}->{run}) {

		$vardict	= $tool_data->{vardict_version};
		$bcftools	= $tool_data->{samtools_version};
		$vcftools	= $tool_data->{vcftools_version};
		$gatk		= $tool_data->{gatk_version};

		if (defined($tool_data->{vardict}->{pon})) {
			$pon = $tool_data->{vardict}->{pon};
			}

		# annotation
		my @parts = split('\\/', $tool_data->{annotate}->{vep_path});
		$vep = $parts[-1];
		@parts = split('\\/', $tool_data->{annotate}->{vcf2maf_path});
		$vcf2maf = $parts[-2];

		# fill in methods
		$methods .= "\\subsubsection{VarDict (v$vardict)}\n";
		$methods .= "VarDict was run using the recommended protocol for either paired (T/N) or tumour-only, with default parameters and target regions (\$\\pm\$ 50bp padding) provided. Variants were filtered by significance (p-value \$<\$ 0.05) using bcftools (v$bcftools).\\newline\n";

		if (defined($tool_data->{vardict}->{pon})) {
			$methods .= "For T/N, somatic (StrongSomatic and LikelySomatic) and germline variants were extracted.\\newline\n";
			} else {
			$methods .= "For T/N, somatic (StrongSomatic and LikelySomatic) and germline variants were extracted. Germline variants were used to generate a panel of normals (merged across samples using GATK's (v$gatk) CombineVariants, and keeping variants present in at least 2 samples).\\newline\n";
			}

		if (defined($tool_data->{vardict}->{pon})) {
			$methods .= "For tumour-only samples, somatic variants were extracted (StrongSomatic and LikelySomatic), and filtered further using a pre-made panel of normals ($pon) to remove probable germline variants (vcftools v$vcftools).\\newline\n";
			} else {
			$methods .= "For tumour-only samples, somatic variants were extracted (StrongSomatic and LikelySomatic), and filtered further using the PoN to remove probable germline variants (vcftools v$vcftools).\\newline\n";
			}

		} else {
		$methods .= "VarDict was not run.\\newline\n";
		}

	if ('Y' eq $tool_data->{strelka}->{run}) {

		$strelka	= $tool_data->{strelka_version};
		$manta		= $tool_data->{manta_version};
		$gatk		= $tool_data->{gatk_version};
		$vcftools	= $tool_data->{vcftools_version};

		if (defined($tool_data->{strelka}->{pon})) {
			$pon = $tool_data->{strelka}->{pon};
			}

		# annotation
		my @parts = split('\\/', $tool_data->{annotate}->{vep_path});
		$vep = $parts[-1];
		@parts = split('\\/', $tool_data->{annotate}->{vcf2maf_path});
		$vcf2maf = $parts[-2];

		# fill in methods
		$methods .= "\\subsubsection{Strelka (v$strelka)}\n";
		$methods .= "Strelka and Manta were run using default parameters in exome mode, with callable regions restricted to target regions (\$\\pm\$ 50bp padding).\\newline\n";

		if (defined($tool_data->{strelka}->{pon})) {
			$methods .= "A pre-made panel of normals ($pon) was used for additional filtering.\\newline\n";
			} else {
			$methods .= "A panel of normals was produced using all available normal samples, using Strelka's germline workflow. Resulting variants were merged across samples using GATK's CombineVariants (v$gatk), removing any variant that did not pass quality criteria (FILTER field not equal to PASS), and keeping variants present in at least 2 samples.\\newline\n";
			}

		$methods .= "For somatic variant detection, Strelka  was run on each sample following the developers recommended protocol. First, Manta (v$manta) was run on each T/N pair or tumour-only sample to identify a set of candidate small indels to be provided to Strelka for use in variant calling. Strelka's somatic workflow was run on T/N pairs, while the germline workflow was used for tumour-only samples. In both cases, resulting variant lists were filtered (using vcftools v$vcftools) to remove likely germline variants (found in the panel of normals) and poor quality calls (FILTER field did not equal PASS).\\newline\n";
		} else {
		$methods .= "Strelka not run.\\newline\n";
		}

	if ('Y' eq $tool_data->{varscan}->{run}) {

		$varscan	= $tool_data->{varscan_version};
		$gatk		= $tool_data->{gatk_version};
		$samtools	= $tool_data->{samtools_version};
		$vcftools	= $tool_data->{vcftools_version};

		if (defined($tool_data->{varscan}->{pon})) {
			$pon = $tool_data->{varscan}->{pon};
			}

		# annotation
		my @parts = split('\\/', $tool_data->{annotate}->{vep_path});
		$vep = $parts[-1];
		@parts = split('\\/', $tool_data->{annotate}->{vcf2maf_path});
		$vcf2maf = $parts[-2];

		# fill in methods
		$methods .= "\\subsubsection{VarScan (v$varscan)}\n";
		$methods .= "For variant calling in T/N pairs, samtools (v$samtools) mpileup was run on each T/N pair, using -B -q1 -d 10000 and restricting call regions to target regions (\$\\pm\$ 50bp padding). Positions with 0 coverage in both the tumour and normal were excluded and the resulting output provided to VarScan. VarScan somatic and processSomatic were used to generate lists of high-confidence germline and somatic variant positions. VarScan somatic was run again, using --output-vcf, and the resulting VCF was filtered (using bcftools v$samtools) to produce a high-confidence germline VCF file and a high-confidence somatic VCF file.\\newline\n";

		if (!defined($tool_data->{vardict}->{pon})) {
			$methods .= "To produce a panel of normals, high-confidence germline variants were merged across samples using GATK's CombineVariants (v$gatk), removing any variant that did not pass quality criteria (FILTER field not equal to PASS), and keeping variants present in at least 2 samples.\\newline\n";
			}

		$methods .= "For variant calling in tumour-only samples, samtools (v$samtools) mpileup was run, again using -B -q1 -d 10000 and restricting call regions to target regions (\$\\pm\$ 50bp padding). Positions with 0 coverage were excluded and the resulting output provided to VarScan mpileup2cns using --output-vcf and --variants. ";
		if (!defined($tool_data->{vardict}->{pon})) {
			$methods .= "Resulting variants were filtered (using vcftools v$vcftools) to remove germline variants (using the panel of normals).\\newline\n";
			} else {
			$methods .= "Resulting variants were filtered (using a pre-made panel of normals: $pon) to remove germline variants (using the panel of normals).\\newline\n";
			}

		} else {
		$methods .= "VarScan not run.\\newline\n";
		}

	if (defined($vep)) {
		$methods .= "\\subsubsection{Annotation}\n";
		$methods .= "Somatic short variants (SNVs and INDELs) were annotated using VEP (v$vep) and vcf2maf ($vcf2maf), with population frequencies for known common variants annotated from ExAC (nonTCGA version r1) and gnomAD.\\newline\n";
		$methods .= "\\noindent\nLastly, an ensemble approach was applied, such that variants meeting the following criteria were carried forward for downstream analyses:\\newline\n";
		# based on suggested criteria here:
		# 	https://www.nature.com/articles/s41598-020-69772-8
		$methods .= join("\n",
			"{\\scriptsize \\begin{itemize}",
			"  \\vspace{-0.2cm}\\item SNPs identified by 4 or more tools (of 6 total)",
			"  \\vspace{-0.2cm}\\item INDELs identified by 3 or more tools (of 5 total)",
			"  \\vspace{-0.2cm}\\item variants identified by MuTect2, with VAF \$<\$ 0.1",
			"  \\vspace{-0.2cm}\\item variants with intra-patient evidence (any of the above 3 criteria)",
			"\\end{itemize} }"
			) . "\n";
		$methods .= "\\noindent\nFlags were added to signify possible concerns, including:\\newline\n";
		$methods .= join("\n",
			"{\\scriptsize \\begin{itemize}",
			"  \\vspace{-0.2cm}\\item tumours without a matched normal", 
			"  \\vspace{-0.2cm}\\item low coverage (\$<\$ $t_depth in tumour; \$<\$ $n_depth in normal)", 
			"  \\vspace{-0.2cm}\\item low VAF \$<\$ 0.05",
			"  \\vspace{-0.2cm}\\item high population allele frequency (AF \$>\$ 0.001 in gnomAD, ExAc, other)",
			"\\end{itemize} }"
			) . "\n";
		} else {
		$methods .= "No somatic variant calling performed.\\newline\n";
		}

	# how were SVs called?
	$methods .= "\\subsection{Structural Variant Calling}\n";

	$delly = $tool_data->{delly_version};
	$manta = $tool_data->{manta_version};
	$novobreak = $tool_data->{novobreak_version};
	$samtools = $tool_data->{samtools_version};

	# fill in methods
	if ( ('Y' eq $tool_data->{delly}->{run}) && 
		('Y' eq $tool_data->{strelka}->{run}) && ('Y' eq $tool_data->{novobreak}->{run}) ) {
		$methods .= "Somatic structural variants (SVs) including large insertions/deletions, duplications, inversions and translocations were identified using Delly (v$delly), Manta (v$manta) and NovoBreak (v$novobreak).\\newline\n";
		} elsif ( ('Y' eq $tool_data->{delly}->{run}) && 
		('Y' eq $tool_data->{strelka}->{run}) && ('N' eq $tool_data->{novobreak}->{run}) ) {
		$methods .= "Somatic structural variants (SVs) including large insertions/deletions, duplications, inversions and translocations were identified using Delly (v$delly) and Manta (v$manta).\\newline\n";
		} elsif ( ('Y' eq $tool_data->{delly}->{run}) && 
		('N' eq $tool_data->{strelka}->{run}) && ('Y' eq $tool_data->{novobreak}->{run}) ) {
		$methods .= "Somatic structural variants (SVs) including large insertions/deletions, duplications, inversions and translocations were identified using Delly (v$delly) and NovoBreak (v$novobreak).\\newline\n";
		} elsif ( ('N' eq $tool_data->{delly}->{run}) && 
		('Y' eq $tool_data->{strelka}->{run}) && ('Y' eq $tool_data->{novobreak}->{run}) ) {
		$methods .= "Somatic structural variants (SVs) including large insertions/deletions, duplications, inversions and translocations were identified using Manta (v$manta) and NovoBreak (v$novobreak).\\newline\n";
		} elsif ( ('Y' eq $tool_data->{delly}->{run}) && 
		('N' eq $tool_data->{strelka}->{run}) && ('N' eq $tool_data->{novobreak}->{run}) ) {
		$methods .= "Somatic structural variants (SVs) including large insertions/deletions, duplications, inversions and translocations were identified using Delly (v$delly).\\newline\n";
		} elsif ( ('N' eq $tool_data->{delly}->{run}) && 
		('Y' eq $tool_data->{strelka}->{run}) && ('N' eq $tool_data->{novobreak}->{run}) ) {
		$methods .= "Somatic structural variants (SVs) including large insertions/deletions, duplications, inversions and translocations were identified using Manta (v$manta).\\newline\n";
		} elsif ( ('N' eq $tool_data->{delly}->{run}) && 
		('N' eq $tool_data->{strelka}->{run}) && ('Y' eq $tool_data->{novobreak}->{run}) ) {
		$methods .= "Somatic structural variants (SVs) including large insertions/deletions, duplications, inversions and translocations were identified using NovoBreak (v$novobreak).\\newline\n";
		} else {
		$methods .= "No structural variant calling performed or this was performed outside of the present pipeline.\\newline\n";
		}
	$methods .= "\\newline\n";

	if ('Y' eq $tool_data->{delly}->{run}) {
		$methods .= "Delly was run on each T/N pair or tumour-only sample, with variants filtered (per patient; -m 0 -a 0.1 -r 0.5 -v 10 -p) and merged (cohort; -m 0 -n 250000000 -b 0 -r 1.0) to identify a joint site list. Sites were genotyped in each sample, merged using bcftools (v$samtools) and finalized (per tumour, filtered against all available normals; -m 0 -a 0.1 -r 0.5 -v 10 -p), according to published best practices.\\newline\n\\newline\n";
		}

	if ('Y' eq $tool_data->{strelka}->{run}) {
		$methods .= "Manta was run on each T/N pair or tumour-only sample, using the --exome setting and --callRegions to restrict calls to target intervals.\\newline\n\\newline\n";
		}

	if ('Y' eq $tool_data->{novobreak}->{run}) {
		$methods .= "NovoBreak was run on each T/N pair (if available) using default parameters and recommended protocol; target intervals were not provided.\\newline\n\\newline\n";
		}

	if ('Y' eq $tool_data->{mavis}->{run}) {
		$ref_type	= $tool_data->{ref_type};
		$mavis		= $tool_data->{mavis_version};

		$methods .= "Mavis (v$mavis) was run once for each patient, using available SV calls, with the $ref_type reference files provided by the developers. BWA was indicated as the aligner, with the bwa-indexed reference file as above.\\newline\n";
		} else {
		$methods .= "Mavis not run.\\newline\n";
		}

	# for copy number variants
	if ('Y' eq $tool_data->{gatk_cnv}->{run}) {
		$gatk4 = $tool_data->{gatk_cnv_version};
		} else {
		$gatk4 = 0;
		}

	if (defined($varscan) & defined($gatk4)) {
		$methods .= "\\newline\n\\noindent Somatic copy-number aberrations (SCNAs) were identified using VarScan and Sequenza, as well as GATK's CNV pipeline.\\newline\n";
		} elsif (defined($varscan)) {
		$methods .= "\\newline\n\\noindent Somatic copy-number aberrations (SCNAs) were identified using VarScan and Sequenza.\\newline\n";
		}

	# for varscan/sequenza
	if (defined($varscan)) {
		$methods .= "VarScan was run as above on T/N pairs, using the copynumber and copyCaller tools. Input for Sequenza was formatted using VarScan2seqz (Sequenza v2.1, R v3.3.0). Sequenza functions (extract and fit) were run using Sequenza v3.0.0; R v3.6.1, with the modified copynumber package (v1.29.0.9000 - to ensure compatibility with hg38). Tuning of the gamma parameter within sequenza.extract was performed to obtain the optimal value, and sequenza.fit performed utilizing ploidy priors obtained from TGCA. Resulting copy number segments were mapped to gene IDs and ploidy-adjusted.\\newline\n";
		}

	if (defined($gatk4)) {
		$methods .= "GATK's CNV pipeline (v$gatk4) was run using best practice guidelines. Pipeline was run with default parameters except model step which was run using number of smoothing iterations = 1.\\newline";
		}

	if (!defined($varscan) & !defined($gatk4)) {
		$methods .= "No tools were run to detect copy-number alterations.\\newline\n";
		}

	# for MSI
	if ('Y' eq $tool_data->{other_tools}->{run_msi}) {

		my $msi	= 'msisensor-pro/1.2.0';

		# fill in methods
		$methods .= "\\subsection{MSI}\n";
		$methods .= "MSI-Sensor pro (v1.2.0) was run using recommended best practices with target regions provided to msi (T/N pairs) or pro (tumour only) commands and coverage threshold of 20x (as recommended).\\newline\n";
		} else {
		$methods .= "MSI-Sensor was not run.\\newline\n";
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
		"\t--directory|-d\t<string> path to output directory",
		) . "\n";

	print $help_msg;
	exit;
	}

# do some quick error checks to confirm valid arguments	
main(config => $config, directory => $directory);
