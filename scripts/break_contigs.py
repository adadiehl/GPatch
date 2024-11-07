#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "0.0.1"

"""
Break contigs at breakpoints indicated in the given text file.
"""

def parse_breakpoints(breakpoints_f):
    """
    Parse the breakpoints file into a dict.
    """
    ret = {}
    with open(breakpoints_f, "r") as breakpoints:
        for line in breakpoints:
            rec = line.strip().split()
            #sys.stderr.write("%s\n" % (rec))
            ret[rec[3]] = rec
    return ret
    


def main():
    parser = ArgumentParser(description="Break contigs at breakpoints indicated in the given text file.")
    parser.add_argument('-f', '--fasta', metavar='FASTA', type=str,
                        required=True, help='Path to genome assembly fasta containing contigs to break.')
    parser.add_argument('-b', '--breakpoints', metavar='PATH', type=str,
                        required=True, help='Path to text file containing breakpoints for each contig to break. Only supports one breakpoint per contig, at present. Contigs without an entry in this file will be printed to output unchanged.')
    parser.add_argument('-m', '--min_frag_len', metavar='FASTA', type=str,
                        required=False, default=1000,
                        help='Minimum fragment length when breaking a contig, in bp. Default=1000')

    
    args = parser.parse_args()

    breakpoints = parse_breakpoints(args.breakpoints)

    for fasta in SeqIO.parse(args.fasta, "fasta"):
        if fasta.id in breakpoints:
            # Break the contig at the indicated breakpoints and
            # print a separate record for each fragment.
            seq = fasta.seq
            bp = breakpoints[fasta.id]
            dist3 = int(bp[2]) - int(bp[8])  # Distance from left breakpoint to contig mapped end
            dist1 = int(bp[7]) - dist3   # Distance from contig mapped start to left breakpoint
            dist2 = dist1 + (int(bp[9]) - int(bp[8])) # Distance from contig start to right breakpoint

            # Check for negative values in dist1.
            if dist1 < 0:
                dist1 = 0
            # Check for illegal values for dist2.
            if dist2 > len(seq):
                dist2 = len(seq)

            frag = ["","",""]
            if bp[5] == "+":
                frag[0] = seq[0:dist1]
                frag[1] = seq[dist1:dist2]
                frag[2] = seq[dist2:len(seq)]
            else:
                frag[0] = seq[len(seq)-dist1:len(seq)]
                frag[1] = seq[len(seq)-dist2:len(seq)-dist1]
                frag[2] = seq[0:len(seq)-dist2]

            # Just swap out names and sequences in the current seq object to print results.
            root_id = fasta.id
            for i in range(3):
                if len(frag[i]) < args.min_frag_len:
                    continue
                fasta.seq = frag[i]
                fasta.id = root_id + '_' + str(i)
                sys.stdout.write("%s\n" % fasta.format("fasta"))
            
        else:
            # Print the record unchanged.
            sys.stdout.write("%s\n" % (fasta.format("fasta")))
    
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
