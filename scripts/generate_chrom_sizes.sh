#!/bin/bash

# Generate a chrom.sizes file for the given genome fasta. Assumes
# uncompressed FASTA input.

FASTA_IN=$1
CS_OUT=$2

awk 'BEGIN {RS = ">"} {if (NF > 0) {chrom_name = $1; sub(/ .*/, "", chrom_name); seq_length = 0; for (i = 2; i <= NF; i++) {seq_length += length($i)}; printf "%s\t%d\n", chrom_name, seq_length}}' $FASTA_IN > $CS_OUT
