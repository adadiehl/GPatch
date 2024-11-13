#!/bin/bash

# Shell script to automate two-stage genome patching, including initial alignment
# and patching, automated suspicious-rearrangement-breakpoint location,
# contig breaking, realignment, and final patching. Dot plots are produced for
# all chromsomes before and after contig-breaking.

ASSEMBLY_FASTA=$1
REFERENCE_FASTA=$2
REFERENCE_NAME=$3
PREFIX=$4
WHITELIST=$5

PG_PATH=/data/projects/adadiehl/genome_patching/patch_genome/src/patch_genome
SCRIPTS_PATH=/data/projects/adadiehl/genome_patching/patch_genome/scripts

# Initial mapping of the assembly genome and first-round of patching.
echo "Initial mapping to reference..."
#minimap2 -x asm20 -t 24 -a $REFERENCE_FASTA $ASSEMBLY_FASTA | samtools view -b - > $PREFIX.bam 2> $PREFIX.patch_genome.err

echo "Initial genome patching..."
#time $PG_PATH/patch_genome.py -q $PREFIX.bam -r $REFERENCE_FASTA -x $PREFIX -w $WHITELIST 2>> $PREFIX.patch_genome.err

# Stats and dot plots
#printf "pre-break contig base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $PREFIX.contigs.bed) > $PREFIX.stats
#printf "pre-break patch base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $PREFIX.patches.bed) >> $PREFIX.stats
#printf "Initial unpatched assembly genome size (bp): %d\n" $(grep -v ">" $ASSEMBLY_FASTA | awk 'BEGIN{N=0}{N+=length($0)}END{print N}') >> $PREFIX.stats
#printf "pre-break total patched genome size (bp): %d\n" $(grep -v ">" $PREFIX.patched.fasta | awk 'BEGIN{N=0}{N+=length($0)}END{print N}') >> $PREFIX.stats

# Alignment to reference and dot-plots.
echo "Aligning patched genome to reference..."
minimap2 -x asm20 -t 24 $REFERENCE_FASTA $PREFIX.patched.fasta > $PREFIX.patched.$REFERENCE_NAME.paf 2>> $PREFIX.patch_genome.err
echo "Generating initial dot-plots..."
$SCRIPTS_PATH/create_dotplots.Rscript $PREFIX.patched.$REFERENCE_NAME.paf $PREFIX $REFERENCE_NAME

# Automated location of suspicious rearrangement breakpoints that
# are likely misassemblies. Note we overwrite original results here.
echo "Finding suspicious rearrangement breakpoints..."
$SCRIPTS_PATH/find_contig_breakpoints.py -c $PREFIX.contigs.bed -p $PREFIX.patched.$REFERENCE_NAME.paf -d 25000 -q 25000 > $PREFIX.breakpoints.txt 2>> $PREFIX.patch_genome.err
echo "Breaking contigs at suspicious rearrangement breakpoints..."
$SCRIPTS_PATH/break_contigs.py -f $ASSEMBLY_FASTA -b $PREFIX.breakpoints.txt > $PREFIX.break_contigs.fa 2>> $PREFIX.patch_genome.err

echo "Remapping to reference..."
minimap2 -x asm20 -t 24 -a $REFERENCE_FASTA $PREFIX.break_contigs.fa | samtools view -b - > $PREFIX.break_contigs.bam
time $PG_PATH/patch_genome.py -q $prefix.break_contigs.bam -r $REFERENCE_FASTA -x $PREFIX -w $WHITELIST 2>> $PREFIX.patch_genome.err

# Stats and dot plots
printf "post-break contig base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $PREFIX.contigs.bed) >> $PREFIX.stats
printf "post-break patch base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $PREFIX.patches.bed) >> $PREFIX.stats
printf "post-break total patched genome size (bp): %d\n" $(grep -v ">" $PREFIX.patched.fasta | awk 'BEGIN{N=0}{N+=length($0)}END{print N}') >> $PREFIX.stats

# Alignment to reference and dot-plots.
echo "Mapping the final patched genome against the reference..."
minimap2 -x asm20 -t 24 $REFERENCE_FASTA $PREFIX.patched.fasta > $PREFIX.patched.$REFERENCE_NAME.paf 2>> $PREFIX.patch_genome.err
echo "Generating dot plots"
$SCRIPTS_PATH/create_dotplots.Rscript $PREFIX.patched.$REFERENCE_NAME.paf $PREFIX.break_contigs $REFERENCE_NAME
