#!/bin/bash

# Shell script to automate two-stage genome patching, including initial alignment
# and patching, automated suspicious-rearrangement-breakpoint location,
# contig breaking, realignment, and final patching. Dot plots are produced for
# all chromsomes before and after contig-breaking.

ASSEMBLY_FASTA=$1
REFERENCE_FASTA=$2
REFERENCE_NAME=$3
PREFIX=$4
WHITELIST=$5  # Supply emtpy string ("") to disable whitelist
GPATCH_ARGS=$6  # For supplying extra args to GPatch: "-d -t", e.g.

# Edit the following two lines to indicate the local path to GPatch source and scripts directories
PG_PATH=/data/projects/adadiehl/genome_patching/GPatch/src/GPatch
SCRIPTS_PATH=/data/projects/adadiehl/genome_patching/GPatch/scripts

# Set this to the number of compute threads to use for mapping steps
THREADS=24

CBREAK_MAXDIST=20000000
CBREAK_MAXCDIST=100000
CBREAK_MAXQDIST=100000
CBREAK_MINSIZE=1000000

# Initial mapping of the assembly genome and first-round of patching.
echo "Initial mapping to reference..."
minimap2 -x asm20 -t $THREADS -a $REFERENCE_FASTA $ASSEMBLY_FASTA | samtools view -b - > $PREFIX.bam 2> $PREFIX.GPatch.err

echo "Initial genome patching..."
if [ "${WHITELIST}" != "" ]; then
    /usr/bin/time -v $PG_PATH/GPatch.py -q $PREFIX.bam -r $REFERENCE_FASTA -x $PREFIX -w $WHITELIST $GPATCH_ARGS 2>> $PREFIX.GPatch.err
else
    /usr/bin/time -v $PG_PATH/GPatch.py -q $PREFIX.bam -r $REFERENCE_FASTA -x $PREFIX $GPATCH_ARGS 2>> $PREFIX.GPatch.err
fi

# Stats and dot plots
printf "pre-break contig base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $PREFIX.contigs.bed) > $PREFIX.stats 2>> $PREFIX.GPatch.err
printf "pre-break patch base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $PREFIX.patches.bed) >> $PREFIX.stats 2>> $PREFIX.GPatch.err
printf "Initial unpatched assembly genome size (bp): %d\n" $(grep -v ">" $ASSEMBLY_FASTA | awk 'BEGIN{N=0}{N+=length($0)}END{print N}') >> $PREFIX.stats 2>> $PREFIX.GPatch.err
printf "pre-break total patched genome size (bp): %d\n" $(grep -v ">" $PREFIX.patched.fasta | awk 'BEGIN{N=0}{N+=length($0)}END{print N}') >> $PREFIX.stats 2>> $PREFIX.GPatch.err

# Alignment to reference and dot-plots.
echo "Aligning patched genome to reference..."
minimap2 -x asm20 -t $THREADS $REFERENCE_FASTA $PREFIX.patched.fasta > $PREFIX.patched.$REFERENCE_NAME.paf 2>> $PREFIX.GPatch.err
echo "Generating initial dot-plots..."
$SCRIPTS_PATH/create_dotplots.Rscript $PREFIX.patched.$REFERENCE_NAME.paf $PREFIX $REFERENCE_NAME 2>> $PREFIX.GPatch.err

# Automated location of suspicious rearrangement breakpoints that
# are likely misassemblies. Note we overwrite original results here.
echo "Finding suspicious rearrangement breakpoints..."
$SCRIPTS_PATH/find_contig_breakpoints.py -c $PREFIX.contigs.bed -p $PREFIX.patched.$REFERENCE_NAME.paf -d $CBREAK_MAXCDIST -q $CBREAK_MAXQDIST -m $CBREAK_MAXDIST -s $CBREAK_MINSIZE -a > $PREFIX.breakpoints.txt 2>> $PREFIX.GPatch.err
echo "Breaking contigs at suspicious rearrangement breakpoints..."
$SCRIPTS_PATH/break_contigs.py -f $ASSEMBLY_FASTA -b $PREFIX.breakpoints.txt > $PREFIX.cbreak.fa 2>> $PREFIX.GPatch.err

echo "Remapping to reference..."
minimap2 -x asm20 -t $THREADS -a $REFERENCE_FASTA $PREFIX.cbreak.fa | samtools view -b - > $PREFIX.cbreak.bam 2>> $PREFIX.GPatch.err

echo "Stage two genome patching..."
if [ "${WHITELIST}" != "" ]; then
    time $PG_PATH/GPatch.py -q $PREFIX.cbreak.bam -r $REFERENCE_FASTA -x $PREFIX.cbreak -w $WHITELIST $GPATCH_ARGS 2>> $PREFIX.GPatch.err
else
    time $PG_PATH/GPatch.py -q $PREFIX.cbreak.bam -r $REFERENCE_FASTA -x $PREFIX.cbreak $GPATCH_ARGS 2>> $PREFIX.GPatch.err
fi

# Stats and dot plots
printf "post-break contig base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $PREFIX.cbreak.contigs.bed) >> $PREFIX.stats 2>> $PREFIX.GPatch.err
printf "post-break patch base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $PREFIX.cbreak.patches.bed) >> $PREFIX.stats 2>> $PREFIX.GPatch.err
printf "post-break total patched genome size (bp): %d\n" $(grep -v ">" $PREFIX.cbreak.patched.fasta | awk 'BEGIN{N=0}{N+=length($0)}END{print N}') >> $PREFIX.stats 2>> $PREFIX.GPatch.err

# Alignment to reference and dot-plots.
echo "Mapping the final patched genome against the reference..."
minimap2 -x asm20 -t $THREADS $REFERENCE_FASTA $PREFIX.cbreak.patched.fasta > $PREFIX.cbreak.patched.$REFERENCE_NAME.paf 2>> $PREFIX.GPatch.err
echo "Generating dot plots"
$SCRIPTS_PATH/create_dotplots.Rscript $PREFIX.cbreak.patched.$REFERENCE_NAME.paf $PREFIX.cbreak $REFERENCE_NAME 2>> $PREFIX.GPatch.err
