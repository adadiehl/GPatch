# GPatch
## Assemble contigs into a chromosome-scalse pseudo-assembly using alignments to a reference sequence.

Starting with alignments of contigs to a reference genome, produce a chromosome-scale pseudoassembly by patching gaps between mapped contigs with sequences from the reference.

## Dependencies
* Python >= v3.7
* samtools (https://github.com/samtools/samtools)
* biopython (https://biopython.org/)
* pysam (https://github.com/pysam-developers/pysam)
* minimap2 (https://github.com/lh3/minimap2)

We recommend using minimap2 for alignment, using the -a option to generate SAM output.

## Installation

We recommend installing with conda, into a new environment:
```
conda create -n GPatch -c conda-forge -c bioconda Bio pysam minimap2 samtools GPatch
```

Install with pip:
```
pip install GPatch
```

Installation from the github repository is not recommended. However, if you must, follow the steps below:
1) git clone https://github.com/adadiehl/GPatch
2) cd GPatch/
3) python3 -m pip install -e .


## Usage
```
GPatch.py [-h] -q SAM/BAM -r FASTA [-x STR] [-b FILENAME] [-m N] [-w PATH] [-d] [-t] [-s] [-l N] [-e] [-k]
          [--legacy_junctions] [--max_resid N] [--chain_mapq N] [--chain_minlen N] [--chain_maxgap N]
          [--chain_tol N] [--anchor N] [--min_overlap_identity F] [--indel_tol N] [--max_mm2_tests N]
          [--win_min N] [--win_max N] [--minimap2 PATH] [--threads N]
          [--no_extension] [--ext_min_overlap N] [--ext_region_pad N] [--ext_window N] [--ext_max_steps N]
```

Starting with alignments of contigs to a reference genome, produce a chromosome-scale pseudoassembly by patching gaps between mapped contigs with sequences from the reference. By default, reference chromosomes with no mapped contigs are printed to output unchanged. Use the --drop_missing option to disable this behavior. By default, patches are applied to the 5' and 3' telomere ends of pseudochromsomes if the first and last contig alignments do not extend to the start/end of the reference chromsome. In some cases, this may cause spurious duplications. Use the --no_extend option if this is a concern. Note that GPatch is designed to be run on single-haplotype or unphased genome assemblies. For phased assemblies, each haplotype should be separated into its own input FASTA file prior to alignment. GPatch can then be run separately on the BAM files for each haplotype to obtain phased pseudoassemblies, otherwise results will be unpredictable and likely incorrect.

#### Required Arguments
| Argument | Description |
|---|---|
| __-q SAM/BAM, --query_bam SAM/BAM__ | Path to SAM/BAM file containing non-overlapping contig mappings to the reference genome. |
| __-r FASTA, --reference_fasta FASTA__ | Path to reference genome fasta. |

#### Optional Arguments:
| Argument | Description |
|---|---|
| __-h, --help__ | Show this help message and exit. |
| __-x STR, --prefix STR__ | Prefix to add to output file names. Default=None |
| __-b FILENAME, --store_final_bam FILENAME__ | Store the final set of primary contig alignments to the given file name. Default: Do not store the final BAM. |
| __-m N, --min_qual_score N__ | Minimum mapping quality score to retain an alignment. Default=30 |
| __-w PATH, --whitelist PATH__ | Path to BED file containing whitelist regions: i.e., the inverse of blacklist regions. Supplying this will have the effect of excluding alignments that fall entirely within blacklist regions. Default=None |
| __-d, --drop_missing__ | Omit unpatched reference chromosome records from the output if no contigs map to them. Default: Unpatched chromosomes are printed to output unchanged. |
| __-t, --no_trim__ | Do not trim the 5' end of contigs that overlap the previously-placed contig; overlapping contigs are bookended untrimmed. Default: Trim the overlap, measured directly from the contig sequences (or, with --legacy_junctions, inferred from the soft-clip-extended alignment coordinates). |
| __-s, --scaffold_only__ | Pad gaps between placed contigs with strings of N characters instead of patching with sequence from the reference assembly. Effectively turns GPatch into a reference-guided scaffolding tool. Note that patches.bed will still be generated to document (inverse) mapped contig boundaries in reference frame. |
| __-l N, --gap_length N__ | Length of "N" gaps separating placed gontigs when using --scaffold_only. Has no effect when in default patching mode. Default=Estimate gap length from alignment. |
| __-e, --no_extend__ | Do not patch telomere ends of pseudochromosomes with reference sequence upstream of the first mapped contig and downstream of the last mapped contig. Default is to include 5' and 3' patches to extend telomeres to the ends implied by the alignment. |
| __-k, --keep_nested__ | Do not drop contigs with mapped positions nested entirely inside other mapped contigs. Instead, these will be bookended after the contig in which they are nested. Default is to drop contigs with mapped positions nested entirely within other mapped contigs. This option should be used with caution as these mappings cannot be placed unambigiously relative to other mapped contigs, thus including them is likely to lead to unpredictable and possibly incorrect results. Do not use this unless you are sure you know what you are doing! |

#### Junction-Handling Arguments:
| Argument | Description |
|---|---|
| __--legacy_junctions__ | Use the original (GPatch <= 0.4) junction logic: contig extents from the primary alignment extended by its full soft-clips, overlaps trimmed by the overlap of those extents, and gaps patched by their distance. No junctions.tsv is written. |
| __--max_resid N__ | Maximum unaligned sequence (bp) remaining at a contig end, after chaining, for the reference gap estimate at that end to be trusted. Default=10000 |
| __--chain_mapq N__ | Minimum mapping quality of supplementary alignments used to extend contig extents. Default=20 |
| __--chain_minlen N__ | Minimum aligned query length (bp) of supplementary alignments used to extend contig extents. Default=10000 |
| __--chain_maxgap N__ | Maximum gap (bp), in both query and reference, between chained alignments. Default=1000000 |
| __--chain_tol N__ | Maximum overlap (bp) allowed between chained alignments. Default=10000 |
| __--anchor N__ | Length (bp) of the exact anchors used to detect overlaps between adjacent contigs; also the minimum detectable overlap. Default=500 |
| __--min_overlap_identity F__ | Minimum identity over the whole overlap for a candidate overlap to be accepted. Default=0.999 |
| __--indel_tol N__ | Maximum length difference (bp) between the two anchor-implied overlap lengths for a candidate containing indels; such candidates are verified with minimap2. 0 disables. Default=50 |
| __--max_mm2_tests N__ | Maximum number of indel-containing candidate overlaps verified with minimap2 per junction. Default=20 |
| __--win_min N__ | Minimum size (bp) of the contig-end windows searched for overlaps. Default=1000000 |
| __--win_max N__ | Maximum size (bp) of the contig-end windows searched for overlaps. Default=4000000 |
| __--minimap2 PATH__ | Path to minimap2, used only to verify candidate overlaps containing indels. Default: minimap2 on PATH; if not found, such candidates are not considered. |
| __--threads N__ | Number of processes used for overlap detection. Default=1 |

#### Overlap-Extension Arguments:
| Argument | Description |
|---|---|
| __--no_extension__ | Do not insert unplaced contigs by overlap extension. |
| __--ext_min_overlap N__ | Minimum overlap (bp) for an extension step, or for a walk to be considered to have reached the opposite flanking contig. Short overlaps in tandem and palindromic repeats are often spurious. Default=5000 |
| __--ext_region_pad N__ | A candidate contig's primary alignment must lie within the gap region extended by this many bp beyond the flanking contigs' alignments. Default=1000000 |
| __--ext_window N__ | Size (bp) of the contig-end windows searched for extension overlaps. Default=2000000 |
| __--ext_max_steps N__ | Maximum number of contigs added by one walk. Default=50 |

## Junction Handling

At each junction between adjacent placed contigs, GPatch either trims the 5' end of the second contig (the contigs overlap), patches the gap between them with reference sequence, or bookends them. N characters are never introduced, except that patches are replaced by N runs in --scaffold_only mode.

1. __Contig extents.__ The primary alignment of each contig is chained to supplementary alignments of the same contig that are on the same chromosome and strand and colinear with it in both query and reference. Each contig end is then extrapolated onto the reference only by the unaligned sequence remaining beyond the outermost chained alignment (the "residual"). These extents are also used to identify nested contigs. With --legacy_junctions, the primary alignment's full soft-clips are used instead.
2. __Overlap measurement.__ Overlaps are measured directly from the contig sequences. Candidate overlaps come from exact matches of two anchors: the next contig's first bases within the end of the previous contig, and the previous contig's last bases within the start of the next contig. A candidate supported by both anchors is accepted if identity over the whole overlap is at least --min_overlap_identity.
3. __Decision.__ The extent-based gap estimate is *trusted* when the residuals at both contig ends of the junction are at most --max_resid, *one-sided* when only one exceeds it, and *untrusted* otherwise.

| Decision | When |
|---|---|
| trim_measured | Exactly one overlap was found; it is trimmed. |
| trim_ambiguous | Several overlap offsets were found (tandem repeats). The one closest to the estimated overlap is trimmed if the estimate is trusted or one-sided; otherwise the smallest. |
| patch | No overlap; the estimated gap is >= 0 and trusted. The gap is patched with reference sequence. |
| patch_one_sided | No overlap; the estimated gap is >= 0 and one-sided. The gap is patched with reference sequence. |
| bookend | No overlap; the estimate predicts an overlap or is untrusted. The contigs are joined directly. |

__Overlap extension.__ Contigs left out because of low mapping quality or nesting are often unambiguously connected to placed contigs by sequence overlap. After alignment-based placement, GPatch walks into each gap (and each chromosome end) from both flanking placed contigs. At each step, it adds an unplaced contig if it is the only one that overlaps the current contig end by at least --ext_min_overlap and its primary alignment is on the same reference sequence, within the gap region and in the orientation in which it would be placed. A walk stops at the first ambiguous or missing step, or when it reaches the contig on the other side of the gap; contigs claimed by more than one walk are rejected. Because the alignments of inserted contigs are unreliable, their reference extents are carried from the contig each was reached from through the measured overlap, and inherit that contig's residual for the decision rules. Every step and stop reason is written to extension.tsv. Disable with --no_extension.

Chromosome ends are patched when the residual at the terminal contig end is at most --max_resid (and --no_extend is not set). Contigs are never trimmed against a chromosome end. The decision and reasoning for every junction are written to junctions.tsv.

## Output

GPatch produces the following output files:
| File | Description |
|---|---|
| __patched.fasta__ | The final patched genome. |
| __contigs.bed__ | Location of contigs in the coordinate frame of the patched genome. |
| __patches.bed__ | Location of patches in the coordinate frame of the reference genome. |
| __unplaced.tsv__ | Contigs in the input BAM that are not in the pseudoassembly, with the reason (unmapped, low_mapq, whitelist, nested) and their primary alignment. |
| __unplaced.fasta__ | Sequences of the unplaced contigs, in their original orientation. |
| __extension.tsv__ | Every overlap-extension walk step (contig added, orientation, overlaps) and the reason each walk stopped, including contigs rejected because several walks claimed them. Not written with --legacy_junctions or --no_extension. |
| __junctions.tsv__ | The decision (trim, patch or bookend) made at each junction between adjacent contigs and at chromosome ends, with the contig extents, residual unaligned lengths, gap estimate, candidate overlaps and the reason for the decision. Not written with --legacy_junctions. |

## Helper Scripts

The scripts directory contains helper scripts for running GPatch and working with its output, including a shell script to automate a two-stage patching process, including initial alignment and patching steps, misjoin breakpoint prediction and contig-breaking, and subsequent realignment and patching of the split contigs, with dot-plots against the reference assembly created after both patching stages. In addition, scripts are provided to generate a chrom.sizes file for a patched pseudoassembly and a set of liftover chains that can be used to translate features mapped to the unpatched contigs to the patched pseudoassembly. Please see the README.md within the scripts for detailed usage information.

## Citing GPatch
Please use the following citation if you use this software in your work:

Fast and Accurate Draft Genome Patching with GPatch
Adam Diehl, Alan Boyle
bioRxiv 2025.05.22.655567; doi: https://doi.org/10.1101/2025.05.22.655567
