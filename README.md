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
usage: GPatch.py [-h] -q SAM/BAM -r FASTA [-x STR] [-b FILENAME] [-m N]
                 [-w PATH] [-d]
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
| __-s, --scaffold_only__ | Pad gaps between placed contigs with strings of N characters instead of patching with sequence from the reference assembly. Effectively turns GPatch into a reference-guided scaffolding tool. Note that patches.bed will still be generated to document (inverse) mapped contig boundaries in reference frame. |
| __-l N, --gap_length N__ | Length of "N" gaps separating placed gontigs when using --scaffold_only. Has no effect when in default patching mode. Default=Estimate gap length from alignment. |
| __-e, --no_extend__ | Do not patch telomere ends of pseudochromosomes with reference sequence upstream of the first mapped contig and downstream of the last mapped contig. Default is to include 5' and 3' patches to extend telomeres to the ends implied by the alignment. |
| __-k, --keep_nested__ | Do not drop contigs with mapped positions nested entirely inside other mapped contigs. Instead, these will be bookended after the contig in which they are nested. Default is to drop contigs with mapped positions nested entirely within other mapped contigs. This option should be used with caution as these mappings cannot be placed unambigiously relative to other mapped contigs, thus including them is likely to lead to unpredictable and possibly incorrect results. Do not use this unless you are sure you know what you are doing! |

## Output

GPatch produces three output files:
| File | Description |
|---|---|
| __patched.fasta__ | The final patched genome. |
| __contigs.bed__ | Location of contigs in the coordinate frame of the patched genome. |
| __patches.bed__ | Location of patches in the coordinate frame of the reference genome. |

## Helper Scripts

The scripts directory contains helper scripts for running GPatch and working with its output, including a shell script to automate a two-stage patching process, including initial alignment and patching steps, misjoin breakpoint prediction and contig-breaking, and subsequent realignment and patching of the split contigs, with dot-plots against the reference assembly created after both patching stages. In addition, scripts are provided to generate a chrom.sizes file for a patched pseudoassembly and a set of liftover chains that can be used to translate features mapped to the unpatched contigs to the patched pseudoassembly. Please see the README.md within the scripts for detailed usage information.

## Citing GPatch
Please use the following citation if you use this software in your work:

Fast and Accurate Draft Genome Patching with GPatch
Adam Diehl, Alan Boyle
bioRxiv 2025.05.22.655567; doi: https://doi.org/10.1101/2025.05.22.655567
