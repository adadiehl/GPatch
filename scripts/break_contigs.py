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
            if rec[3] in ret:
                ret[rec[3]].append(rec)
            else:
                ret[rec[3]] = [rec]
    # Breakpoints must be sorted in ascending order or
    # not all breakpoints will be used!
    for contig in ret:
        ret[contig] = sorted(ret[contig], key = lambda x: (x[8], x[9]))
    return ret


def extract_oriented_fragment(seq, start, end, is_reverse=False):
    """
    Extract the substring from the given sequence, returning
    the raw sequence if is_reverse == False and the reverse-
    complement if is_reverse == True.
    """
    subseq = seq[start:end]
    if is_reverse:
        subseq = subseq.reverse_complement()
    return subseq


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
            #sys.stderr.write("%d\n" % (len(seq)))
            root_id = fasta.id
            pos = 0    # Start position in fasta sequence
            frag = 0   # Fragment number
            contig_bp = breakpoints[fasta.id]
            # Check for minus strand mappings. These need special-handling.
            # These should be looped over in reverse, but for now we will
            # just reverse-complement the whole seq and then do the same to
            # each subseq before writing fasta output to restore the same
            # strandedness.
            is_reverse = False
            if contig_bp[0][5] == "-":
                is_reverse = True
                seq = seq.reverse_complement()
            for bp in contig_bp:
                #sys.stderr.write("%s\n" % (bp))
                # Need to calculate start and end of fragment from the start
                # of the contig from genomic coordinates.
                fstart = int(bp[8]) - int(bp[1])
                fend = int(bp[9]) - int(bp[1])
                #sys.stderr.write("%d\t%d\t%d\n" % (pos, fstart, fend))
                if fstart < 0:
                    # If the mapped fragment start is upstream of the contig
                    # start, set fstart to zero to avoid errors.
                    fstart = 0
                if fstart < pos:
                    # Fragment overlaps previous fragment. Truncate at pos.
                    fstart = pos
                    if fstart >= fend:
                        # Nested fragment. Skip.
                        continue
                if fstart > pos:
                    # Handle the fragment upstream of the current breakpoint.
                    subseq = extract_oriented_fragment(seq, pos, fstart-1, is_reverse=is_reverse)
                    fasta.seq = subseq
                    fasta.id = root_id + '_' + str(frag)
                    sys.stdout.write("%s\n" % fasta.format("fasta"))
                    frag += 1
                if fend > len(seq):
                    # Don't try to grab sequence past the end of the
                    # contig!
                    fend = len(seq)

                subseq = extract_oriented_fragment(seq, fstart, fend, is_reverse=is_reverse)
                fasta.seq = subseq
                fasta.id = root_id + '_' + str(frag)
                sys.stdout.write("%s\n" % fasta.format("fasta"))
                frag += 1
                pos = fend

            # Handle the last fragment in the contig
            if pos != len(seq):
                #sys.stderr.write("%d\t%d\n" % (pos, len(seq)))
                subseq = extract_oriented_fragment(seq, pos, len(seq), is_reverse=is_reverse)
                fasta.seq = subseq
                fasta.id = root_id + '_' + str(frag)
                sys.stdout.write("%s\n" % fasta.format("fasta"))
            
        else:
            # Print the record unchanged.
            sys.stdout.write("%s\n" % (fasta.format("fasta")))
    
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
