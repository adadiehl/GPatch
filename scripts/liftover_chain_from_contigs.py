#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser

__version__ = "0.0.1"

"""
Build a liftover chain for lifting annotations from the draft genome (target)
to the patched genome (query).
"""

def parse_chrom_sizes(cs_file):
    """
    Parse chrom.sizes file into a dict.
    """
    ret = {}
    with open(cs_file, "r") as csf:
        for line in csf:
            rec = line.strip().split()
            ret[rec[0]] = int(rec[1])
    return ret
    

def main():
    parser = ArgumentParser(description="Build a liftover chain for lifting annotations from the draft genome (target) to the patched genome (query). ")
    parser.add_argument('-c', '--contigs_bed', metavar='BED', type=str,
                        required=True, help='Path to the contigs.bed file created by GPatch for the final patched genome.')
    parser.add_argument('-s', '--chrom_sizes', metavar='CHROM.SIZES', type=str,
                        required=True, help='Path to chrom.sizes file for the final patched genome. Can be generated with build_chrom_sizes.sh.')
    
    args = parser.parse_args()

    # Parse the chrom_sizes file into a dict
    chrom_sizes = parse_chrom_sizes(args.chrom_sizes)

    # Loop over records in the contigs.bed, creating a chain for each record.
    # For all records, "query" is the genome we are mapping to (patched) and "target"
    # is the genome we are mapping from (contigs).
    with open(args.contigs_bed, "r") as contigs:
        chain_id = 1
        for line in contigs:
            rec = line.strip().split()
            # Target name is the name of the contig in the draft genome.
            target_name = rec[3]
            # Target size is the length of the contig
            target_size = int(rec[7])
            # Target strand is always +
            target_strand = '+'
            target_start = int(rec[6])
            target_end = int(rec[7])
            # Query name is the name of the patched chromosome
            query_name = rec[0]
            # Get the size of the query chromosome in the patched genome
            query_size = chrom_sizes[rec[0]]
            # Query strand is the orientation of the contig within the patched genome (?)
            query_strand = rec[5]
            query_start = int(rec[1])
            query_end =	int(rec[2])
            # Chain score is an arbitrary int (meant to be the alignment score for a block).
            # These are used to select between overlapping chains for a region. Since we only
            # have one chain per region, we can use an arbitrary values. We will go with the
            # target_size.
            chain_score = target_size

            # Write a chain record to stdout.
            sys.stdout.write("chain %d %s %d %s %d %d %s %d %s %d %d %d\n" % (chain_score,
                                                                              target_name,
                                                                              target_size,
                                                                              target_strand,
                                                                              target_start,
                                                                              target_end,
                                                                              query_name,
                                                                              query_size,
                                                                              query_strand,
                                                                              query_start,
                                                                              query_end,
                                                                              chain_id))
            sys.stdout.write("%d\t%d\t%d\n" % (target_size, 0, 0))
            sys.stdout.write("%d\n" % (0))

            chain_id += 1
    
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
