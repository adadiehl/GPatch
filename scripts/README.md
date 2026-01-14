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

Stats on patched genome composition are calculated after both patching stages and collected in \<PREFIX\>.stats

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
| __GPATCH_ARGS__ | Additional optional arguments to GPatch. i.e., "-d -s" See GPatch usage information on the main README.md for available options. |

Also note that the number of mapping threads (defautl 24) and assembly mode (default asm20) for minimap2 may be specified by editing the THREADS and/or ASM_MODE variables within patch_genome.sh


### find_contig_breakpoints.py

Given an alignment of a patched pseudoassembly to a reference assembly, locate likely misjoins within the component contigs and generate a text file with their breakpoint coordinates.

#### Usage
```
find_contig_breakpoints.py [-h] -c BED -p PATH [-m INT] [-t INT] [-d INT] [-q INT] [-s INT] [-a]
```

#### Required Arguments
| Argument | Description |
|---|---|
| __-c BED, --contigs BED__ | Path to contigs.bed containing contig mappings within the current patched genome. |
| __-p PATH, --paf PATH__ | PAF alignment of the current patched genome to the reference genome. |

#### Optional Arguments
| Argument | Description |
|---|---|
| __-h, --help__ | show this help message and exit |
| __-m INT, --max_dist INT__ | Maximum distance, in bp, from the contig end to a potential breakpoint to consider it likely spurious and include this in breakpoints.txt. Default=1000 |
| __-t INT, --min_trans_dist INT__ | Minimum distance, in bp, from the query position to the mapped position in the reference to call a likely translocation and include this alignment in breakpoints.txt. Default=1000000 |
| __-d INT, --max_cluster_dist INT__ | Maximum distance between mapped positions to merge a partial alignment into a cluster. default = 0 |
| __-q INT, --max_query_dist INT__ | Maximum distance between query intervals to merge a partial alignment into a cluster. default = 10000 |
| __-s INT, --min_size INT__ | Minimum size, in bp, to call an event and include this alignment in breakpoints.txt. Default=1000000 |
| __-a, --all_chroms__  | Find breakpoints for rearrangements across all chromosomes, not just intrachromosomal rearrangements. Default=Find only intrachromosomal breakpoints. |

#### Output
| File | Description |
|---|---|
| __breakpoints.txt__ | A text file describing the breakpoint coordinates of likely misjoins. See below for column definitions. |

| Column | Description |
|---|---|
| 1-3 | Coordinates of misjoin breakpoints within the patched pseudoassembly. |
| 4 | Name of the contig harboring the misjoin. |
| 5 | Score. Not used. |
| 6 | Orientation of the contig mapping within the patched pseudoassembly. |
| 7-8 | Contig-level start and end positions of the placed contig fragment. |
| 9-10 | Coordinates of misjoin breakpoints within the contig. |
| 11 | Misjoin type. TRS=Translocation, INT=Inversion, INT=Inverted Translocation, TRD=Translocation from different chromosome. |

Note that this file may be edited to remove individual breakpoints if it is believed they represent legitimate structural variants in the contig-level assembly.


### break_contigs.py

Given a contig-level assembly and breakpoints.txt file from find_contig_breakpoints.py, break contigs at the boundaries of likely misjoins and write results to a new fasta file.

#### Usage
```
break_contigs.py [-h] -f FASTA -b PATH [-m FASTA]
```

#### Required Arguments
| Argument | Description |
|---|---|
| __-f FASTA, --fasta FASTA__ | Path to genome assembly fasta containing contigs to break. |
| __-b PATH, --breakpoints PATH__ | Path to text file containing breakpoints for each contig to break. Only supports one breakpoint per contig, at present. Contigs without an entry in this file will be printed to output unchanged. |

#### Optional Arguments
| Argument | Description |
|---|---|
| __-h, --help__ | show this help message and exit |
| __-m FASTA, --min_frag_len FASTA__ | Minimum fragment length when breaking a contig, in bp. Default=1000 |

#### Output
Output in FASTA format are written to stdout.


### create_dotplots.Rscript

Generate dot plots between all patched pseudochromosomes and corresponding reference chromosomes.

#### Usage
```
create_dotplots.Rscript <patched_to_reference.paf> <CTG_NAME> <REF_NAME>
```

#### Required Arguments
| Argument | Description |
|---|---|
| __patched_to_reference.paf__ | PAF format alignment of the patched pseudoassembly to the reference assembly. |
| __CTG_NAME__ | Name of the contig-level draft assembly, e.g., NA12878.1. Used as axis labels in the dot-plot and in ouput file naming. |
| __REF_NAME__ | Name of the reference assembly, e.g., hg38. Used as axis labels in the dot-plot and in ouput file naming. |

#### Output
Output is written as individual image files, in PDF format, for each matching chromosome in the patched pseudoassenbly and reference assembly, with naming convention <CTG_NAME>.dotplot.<REF_NAME>.<CHR>.pdf.


## Citing GPatch Helper Scripts
Please use the following citation if you use this software in your work:

Fast and Accurate Draft Genome Patching with GPatch
Adam Diehl, Alan Boyle
bioRxiv 2025.05.22.655567; doi: https://doi.org/10.1101/2025.05.22.655567
