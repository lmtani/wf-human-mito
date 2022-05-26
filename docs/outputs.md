---
layout: default
title: Explaining the outputs
nav_order: 3
---

# 📦 Explaining the outputs

## BAM and VCF files

Each sample will have an output of:

- the alignment (BAM) between its sequencing reads and the mitochondrion genome
- the variants (VCF) found by using the Mutect2 program.

## The CSV file

This output has many data that can be useful in different ways.

> I've copy-paste some outputs definitions from [Picard documentation](https://broadinstitute.github.io/picard/picard-metric-definitions.html) and from [Haplockeck](https://mitoverse.readthedocs.io/haplocheck/haplocheck/#minor-heteroplasmy-level) here for convenience.

### From MarkDuplicates process

==SampleID==: sample name

==LIBRARY==: sample name (parsed from filename) and the barcode (from inside FASTQ file)

==UNPAIRED_READS_EXAMINED==: This is the number of mapped reads examined which did not have a mapped mate pair.

==READ_PAIRS_EXAMINED==: This is the number of mapped reads examined.

==SECONDARY_OR_SUPPLEMENTARY_RDS==: The number of reads that were either secondary or supplementary

==UNMAPPED_READS==: The total number of unmapped reads examined

==UNPAIRED_READ_DUPLICATES==: The number of fragments that were marked as duplicates.

==READ_PAIR_DUPLICATES==: The number of read pairs that were marked as duplicates.

==READ_PAIR_OPTICAL_DUPLICATES==: The number of read pairs duplicates that were caused by optical duplication. Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source.

==PERCENT_DUPLICATION==: The fraction of mapped sequence that is marked as duplicate.

==ESTIMATED_LIBRARY_SIZE==: The estimated number of unique molecules in the library based on PE duplication.

### From CollectWgsMetrics process

==GENOME_TERRITORY==: The number of non-N bases in the genome reference over which coverage will be evaluated.

==MEAN_COVERAGE==: The mean coverage in bases of the genome territory, after all filters are applied.

==SD_COVERAGE==: The standard deviation of coverage of the genome after all filters are applied.

==MEDIAN_COVERAGE==: The median coverage in bases of the genome territory, after all filters are applied.

==MAD_COVERAGE==: The median absolute deviation of coverage of the genome after all filters are applied.

==PCT_EXC_ADAPTER==:

==PCT_EXC_MAPQ==: The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20).

==PCT_EXC_DUPE==: The fraction of aligned bases that were filtered out because they were in reads marked as duplicates.

==PCT_EXC_UNPAIRED==: The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair.

==PCT_EXC_BASEQ==: The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20).

==PCT_EXC_OVERLAP==: The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.

==PCT_EXC_CAPPED==: The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x).

==PCT_EXC_TOTAL==: The total fraction of aligned bases excluded due to all filters.

==PCT_1X -> PCT_100X==: The fraction of bases that attained at least 1X (5X, 10X...) sequence coverage in post-filtering bases.

==FOLD_80_BASE_PENALTY==: The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets.

==FOLD_90_BASE_PENALTY==: Similar to above

==FOLD_95_BASE_PENALTY==:

==HET_SNP_SENSITIVITY==: The theoretical HET SNP sensitivity.

==HET_SNP_Q==: The Phred Scaled Q Score of the theoretical HET SNP sensitivity.

### From mtDnaServer process

Contamination:

==SampleHomoplasmies==:

==SampleHeteroplasmies==

==SampleMeanCoverage==: defines the mean coverage for the sample.

==HgMajor==

==HgQualityMajor==

==HgMinor==

==HgQualityMinor==

==HomoplasmiesMajor==

==HomoplasmiesMinor==

==HeteroplasmiesMajor==

==HeteroplasmiesMinor==

==MeanHetLevelMajor==: The major haplogroup is calculated by using Haplogrep. The input profile includes all homoplasmies and the major component of each heteroplasmy.

==MeanHetLevelMinor==: The minor heteroplasmy level is calculated by averaging the level of the minor component of each heteroplasmy. The sample HG00245 includes two minor components (3010A (0.011), 16356C (0.012)) to calculate the final heteroplasmy level of 0.011.

==HG_Distance==: This column defines the distance between the haplogroups of the major and minor profile using the graph structure of Phylotree 17


### From CollectAlignmentSummaryMetrics process

==CATEGORY==: Just being explicity that te reads metrics are only considering paired-end alignments

==TOTAL_READS==: Total of reads in the alignment

==PF_READS==: The total number of reads passing filter (PF), where the filter(s) can be platform/vendor quality controls

==PCT_PF_READS==: The fraction of reads passing filter, PF_READS/TOTAL_READS

==PF_NOISE_READS==: The number of PF reads that are marked as noise reads. A noise read is one which is composed entirely of A bases and/or N bases. These reads are marked as they are usually artifactual and are of no use in downstream analysis.

==PF_READS_ALIGNED==:

==PCT_PF_READS_ALIGNED==: The percentage of PF reads that aligned to the reference sequence. PF_READS_ALIGNED / PF_READS

==PF_ALIGNED_BASES==: The total number of aligned bases, in all mapped PF reads, that are aligned to the reference sequence.

==PF_HQ_ALIGNED_READS==: The number of PF reads that were aligned to the reference sequence with a mapping quality of Q20 or higher signifying that the aligner estimates a 1/100 (or smaller) chance that the alignment is wrong.

==PF_HQ_ALIGNED_BASES==: The number of bases aligned to the reference sequence in reads that were mapped at high quality. Will usually approximate PF_HQ_ALIGNED_READS * READ_LENGTH but may differ when either mixed read lengths are present or many reads are aligned with gaps.

==PF_HQ_ALIGNED_Q20_BASES==: The subset of PF_HQ_ALIGNED_BASES where the base call quality was Q20 or higher.

==PF_HQ_MEDIAN_MISMATCHES==: The median number of mismatches versus the reference sequence in reads that were aligned to the reference at high quality (i.e. PF_HQ_ALIGNED READS).

==PF_MISMATCH_RATE==: The rate of bases mismatching the reference for all bases aligned to the reference sequence.

==PF_HQ_ERROR_RATE==: The fraction of bases that mismatch the reference in PF HQ aligned reads.

==PF_INDEL_RATE==: The number of insertion and deletion events per 100 aligned bases. Uses the number of events as the numerator, not the number of inserted or deleted bases.

==MEAN_READ_LENGTH==: The mean read length of the set of reads examined. When looking at the data for a single lane with equal length reads this number is just the read length. When looking at data for merged lanes with differing read lengths this is the mean read length of all reads.

==SD_READ_LENGTH==: Standard deviation of the read length.

==MEDIAN_READ_LENGTH==: Median of the read length.

==MAD_READ_LENGTH==:

==MIN_READ_LENGTH==: Lenght of the shortest read

==MAX_READ_LENGTH==: Length of the longest read

==READS_ALIGNED_IN_PAIRS==: The number of aligned reads whose mate pair was also aligned to the reference.

==PCT_READS_ALIGNED_IN_PAIRS==: The fraction of reads whose mate pair was also aligned to the reference. READS_ALIGNED_IN_PAIRS / PF_READS_ALIGNED

==PF_READS_IMPROPER_PAIRS==: The number of (primary) aligned reads that are **not** "properly" aligned in pairs (as per SAM flag 0x2).

==PCT_PF_READS_IMPROPER_PAIRS==: The fraction of (primary) reads that are *not* "properly" aligned in pairs (as per SAM flag 0x2). PF_READS_IMPROPER_PAIRS / PF_READS_ALIGNED

==BAD_CYCLES==": The number of instrument cycles in which 80% or more of base calls were no-calls.

==STRAND_BALANCE==: The number of PF reads aligned to the positive strand of the genome divided by the number of PF reads aligned to the genome.

==PCT_CHIMERAS==: The fraction of reads that map outside of a maximum insert size (usually 100kb) or that have the two ends mapping to different chromosomes.

==PCT_ADAPTER==: The fraction of PF reads that are unaligned and match to a known adapter sequence right from the start of the read.

==PCT_SOFTCLIP==: the fraction of PF bases that are on (primary) aligned reads and are soft-clipped, as a fraction of the PF_ALIGNED_BASES (even though these are not aligned!)

==PCT_HARDCLIP==: The fraction of PF bases that are (on primary, aligned reads and) hard-clipped, as a fraction of the PF_ALIGNED_BASES (even though these are not aligned!)

==AVG_POS_3PRIME_SOFTCLIP_LENGTH==: The average length of the soft-clipped bases at the 3' end of reads.
