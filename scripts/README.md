# GPatch Helper Scripts
## Helper scripts

This directory contains helper scripts for running GPatch and working with its output, including a shell script to automate a two-stage patching process, including initial alignment and patching steps, misjoin breakpoint prediction and contig-breaking, and subsequent realignment and patching of the split contigs, with dot-plots against the reference assembly created after both patching stages. In addition, scripts are provided to generate a chrom.sizes file for a patched pseudoassembly and a set of liftover chains that can be used to translate features mapped to the unpatched contigs to the patched pseudoassembly.

## Dependencies
* Python >= v3.7
* samtools (https://github.com/samtools/samtools)
* biopython (https://biopython.org/)
* pysam (https://github.com/pysam-developers/pysam)
* minimap2 (https://github.com/lh3/minimap2)
* R (https://www.r-project.org/) and pafr (https://github.com/dwinter/pafr) for	dot-plot generation

## The scripts

### patch_genome.sh

Script to automate a GPatch workflow incuding the following analysis stages:

1) Initial alignment of contigs to the reference assembly
2) Initial GPatch patching of contigs into a patched pseuedoassembly
3) Realignment of the patched pseudoassembly to the reference assembly
4) Generation of dot-plots between the patched pseudoassembly and reference assembly
5) Automated identification of misjoin breakpoints and subsequent contig-breaking
6) Realignment of split contigs to the reference assembly
7) GPatch patching of split contigs into a final patched pseudoassembly
8) Realignment of the final patched pseudoassembly to the reference assembly
9) Generation of dot-plots between the final patched pseudoassembly and reference assembly

Stats on patched genome composition are calculated after both patching stages and collected in <PREFIX>.stats

#### Usage
```
bash patch_genome.sh <contigs.fa> <reference.fa> <REF_NAME> <PREFIX> <whitelist.bed> [GPATCH_ARGS]
```
Note that, before use, you must edit the PG_PATH and SCRIPTS_PATH variables within patch_genome.sh to reflect the local path of the src/GPatch and scripts subdirectories for your local copy of the GPatch github repository!

#### Required Arguments
| Argument | Description |
|---|---|
| __contigs.fa__ | Path to the contig-level assembly to be patched. |
| __reference.fa__ | Path to the reference assembly. |
| __REF_NAME__ | Name of the reference assembly, e.g., "hg38" |
| __PREFIX__ | Prefix to add to output file names, e.g., "NA12878.hap1" |
| __whitelist.bed__ | Path to file containing whitelist regions, i.e., the inverse of blacklist regions. Only contigs mapping within these regions will be retained in the patched pseudoassembly. This is usually only necessary for contig assemblies with many short, repetitive contigs. Supply an empty string (i.e., "") to disable. |

#### Optional Arguments
| Argument | Description |
|---|---|
| __GPATCH_ARGS__ | Additional optional arguments to GPatch. e.g., "-d -s" See GPatch usage information on the main README.md for available options. |

Also note that the number of mapping threads (defautl 24) and assembly mode (default asm20) for minimap2 may be specified by editing the THREADS and/or ASM_MODE variables within patch_genome.sh

### 

## Citing GPatch Helper Scripts
Please use the following citation if you use this software in your work:

Fast and Accurate Draft Genome Patching with GPatch
Adam Diehl, Alan Boyle
bioRxiv 2025.05.22.655567; doi: https://doi.org/10.1101/2025.05.22.655567
