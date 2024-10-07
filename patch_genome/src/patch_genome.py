#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
import random

__version__ = "0.0.1"

"""
Given a name-sorted BAM file, cluster mapped reads by name and
use the maximal interval defined by the primary and supplementary
reads on the same chromosome and strand to define the contig
break-points in the reference genome.
"""

def create_alignment_dict(query_bam, min_qual_score):
    """
    Create a dict of alignments, keyed on read name. This gives us
    the clusters of primary and supplementary alignments we need to
    define the BED intervals that bound contigs in the refernece
    genome.
    """
    ret = {}    
    for aln in query_bam:
        if aln.is_unmapped:
            # Skip unmapped contigs
            continue
        if aln.mapping_quality < min_qual_score:
            # Skip low-quality mappings
            continue
        if aln.query_name in ret:
            ret[aln.query_name].append(aln)
        else:
            ret[aln.query_name] = [aln]        
    return(ret)


def find_breakpoints(aligned_contigs, query_bam, max_merge_dist, min_mapped_fraction, max_expansion):
    """
    Deterimine maximal BED intervals for contig mapping chrom, start, and stop
    based on the reference start and end coordinates of primary and supplementary
    contig mappings.
    """
    contig_breakpoints = {}
    primary_alignments = []
    for cname in aligned_contigs:
        #sys.stderr.write("%s\n" % (cname))
        # First identify primary and supplementary alignments
        supplementary_aln = []
        to_keep = []
        for aln in aligned_contigs[cname]:
            if aln.is_supplementary or aln.is_secondary:
                supplementary_aln.append(aln)
            else:
                # Primary alignment. Store as first index.
                to_keep.append(aln)

        # Skip contigs where the primary alignment is absent (due to
        # not passing the quality filter, for example).
        if len(to_keep) == 0:
            continue

        # Next, find all supplementary alignments on the same chromosome
        # and strand as the primary alignment.
        for aln in supplementary_aln:
            if aln.reference_name == to_keep[0].reference_name:
                if (aln.is_reverse and to_keep[0].is_reverse) or (not aln.is_reverse and not to_keep[0].is_reverse):
                    to_keep.append(aln)

        # Next, order the alignents by mapped position. This requires
        # a BAM file intermediate.
        sorted_bam_name = sorted_bam_from_aln_list(to_keep, query_bam.header)
        pysam.index(sorted_bam_name)

        # Read the position-sorted alignments into a list for
        # further processing.
        sorted_aln = []
        sorted_bam = pysam.AlignmentFile(sorted_bam_name, "rb")
        for aln in sorted_bam:
            sorted_aln.append(aln)
        sorted_bam.close()
        os.remove(sorted_bam_name)
        os.remove(sorted_bam_name + '.bai')

        # Locate the primary alignment in the sorted list and check distances
        # between supplementary alignments to exclude any that are excessively
        # far away from the primary alignment site.
        primary_idx = 0
        for idx, aln in enumerate(sorted_aln):
            if not aln.is_secondary and not aln.is_supplementary:
                # Primary alignment
                primary_idx = idx
        final_aln = [sorted_aln[primary_idx]]
        # Starting at the indeces one-off the primary alignment, check that the
        # distance between upstream supplementary alignments is within the
        # allowed range.
        i = primary_idx - 1
        prev_start = 0
        while i >= 0:
            dist = abs(sorted_aln[i+1].reference_start - sorted_aln[i].reference_end)
            if dist <= max_merge_dist:
                final_aln.append(sorted_aln[i])
            else:
                # Since these are sorted by position, we can (and should!) stop looking after
                # we find any distances greater than the threshold.
                break
            i -= 1
        # Next check downstream alignments.
        i = primary_idx + 1
        while i < len(sorted_aln):
            dist = abs(sorted_aln[i].reference_start - sorted_aln[i-1].reference_end)
            if dist <= max_merge_dist:
                final_aln.append(sorted_aln[i])
            else:
                break
            i += 1

        # Collect the reference starts and ends and sort to locate the 5'-most
        # start and 3'-most end to define our insertion interval.
        ref_starts = []
        ref_ends = []
        for aln in final_aln:
            ref_starts.append(aln.reference_start)
            ref_ends.append(aln.reference_end)
        ref_starts.sort()
        ref_ends.sort()

        # Check that the aligned length represents at least the
        # minimum fraction required (by default, at least 50% of the
        # contig length must be encompassed by the interval.)
        mapped_len = abs(ref_ends[-1]-ref_starts[0])
        if (mapped_len / to_keep[0].query_length) <= min_mapped_fraction:
            continue

        # Also check that the mapped interval has not expanded in length
        # by more than the max_expansion param.
        if (mapped_len / to_keep[0].query_length) > max_expansion:
            continue

        # We now have what we need. Store the result.
        """
        sys.stderr.write("%s\t%d\t%d\t%s\t%d\t%d\n" % (to_keep[0].reference_name, ref_starts[0],
                                                       ref_ends[-1], to_keep[0].query_name,
                                                       to_keep[0].query_length, ref_ends[-1]-ref_starts[0]))
        """
        # Store the primary alignment for later processing.
        primary_alignments.append(to_keep[0])

        # Store the BED interval for the maximal contig mapping.
        contig_breakpoints[to_keep[0].query_name]=[to_keep[0].reference_name, ref_starts[0], ref_ends[-1], to_keep[0].query_name]

    return primary_alignments, contig_breakpoints


def sorted_bam_from_aln_list(aln_list, bam_header, sorted_bam_name=None):
    """
    Given a list of pysam AlignedSegment objects, sort by position with
    a pysam system call. Returns the name of the sorted file.
    """
    tmp_root = ''.join(random.sample('0123456789', 10))
    tmp_bam_name = tmp_root + '.tmp.bam'
    if sorted_bam_name is None:
        sorted_bam_name = tmp_root + '.tmp.sorted.bam'
    tmp_bam = pysam.AlignmentFile(tmp_bam_name, "wb", header=bam_header)
    for aln in aln_list:
        tmp_bam.write(aln)
    tmp_bam.close()
    pysam.sort("-o", sorted_bam_name, tmp_bam_name)
    os.remove(tmp_bam_name)
    return sorted_bam_name


def main():
    parser = ArgumentParser(description="Given a name-sorted BAM file, cluster mapped reads by name and use the maximal interval defined by the primary and supplementary reads on the same chromosome and strand to define the contig break-points in the reference genome. ")
    parser.add_argument('-q', '--query_bam', metavar='SAM/BAM', type=str,
                        required=True, help='Path to SAM/BAM file containing non-overlapping contig mappings to the reference genome.')
    parser.add_argument('-r', '--reference_fasta', metavar='FASTA', type=str,
                        required=True, help='Path to reference genome fasta.')
    parser.add_argument('-x', '--prefix', metavar='BED', type=str,
                        required=False, default="",
                        help='Prefix to add to output file names. Default=None')
    parser.add_argument('-b', '--store_final_bam', metavar='FILENAME', type=str,
                        required=False, default=None,
                        help='Store the final set of primary contig alignments to the fiven file name. Default: Do not store the final BAM.')
    parser.add_argument('-m', '--min_qual_score', metavar='N', type=int,
                        required=False, default=30,
                        help='Minimum mapping quality score to retain an alignment. Default=30')
    parser.add_argument('-d', '--max_merge_dist', metavar='N', type=int,
                        required=False, default=100000,
                        help='Maximum distance to merge adjacent alignments within a cluster. Default=100000')
    parser.add_argument('-f', '--min_mapped_fraction', metavar='FLOAT', type=float,
                        required=False, default=0.5,
                        help='Minimum fraction of the contig length that must be encompassed by the mapped interval. Default=0.5')
    parser.add_argument('-e', '--max_expansion', metavar='FLOAT', type=float,
                        required=False, default=2.0,
                        help='Maximum factor by which the mapped interval is allowed to expand relative to the contig length. I.e., 2.0 allows the mapped interval to be up to double the length of the contig. Default=2.0')
    
    args = parser.parse_args()

    query_bam = pysam.AlignmentFile(args.query_bam)
    
    # Load up the BAM into a dict, keyed on contig name.
    aligned_contigs = create_alignment_dict(query_bam, args.min_qual_score)

    # Determine insertion break-points for each contig based
    # on primary and supplementary mappings.
    primary_alignments, contig_breakpoints = find_breakpoints(aligned_contigs, query_bam, args.max_merge_dist, args.min_mapped_fraction, args.max_expansion)

    # Sort the useful primary alignments by position via pysam sort.
    sorted_primary_alignments_bam = sorted_bam_from_aln_list(primary_alignments, query_bam.header, sorted_bam_name=args.store_final_bam)
    pysam.index(sorted_primary_alignments_bam)

    # Set up output streams for fasta, patches.bed, and contigs.bed
    pf_fname = "patched.fasta"
    pb_fname = "patches.bed"
    cb_fname = "contigs.bed"
    if args.prefix != "":
        pf_fname = args.prefix + '.' + pf_fname
        pb_fname = args.prefix + '.' + pb_fname
        cb_fname = args.prefix + '.' + cb_fname
    patched_fasta = open(pf_fname, "w")    # Patched fasta
    patches_bed = open(pb_fname, "w")      # Reference-frame patch coordinates
    contigs_bed = open(cb_fname, "w")      # Patched genome frame contig coordinates

    # Open up the sorted primary alignments bam for reading.
    sorted_primary_alignments = pysam.AlignmentFile(sorted_primary_alignments_bam, "rb")
    
    for ref_seq in SeqIO.parse(args.reference_fasta, "fasta"):
        # Select all contigs mapped to this sequence.
        contigs = sorted_primary_alignments.fetch(ref_seq.id)
        
        pos = 0 # Tracks position on the reference chromosome        
        patched_seq = ""
        for contig in contigs:
            rstart = contig_breakpoints[contig.query_name][1]
            if rstart > 0:
                rstart -= 1

            # Check for overlapping/bookended contig mappings. These will not
            # need any intevening patch sequence.
            if rstart > pos:
                patched_seq = patched_seq + ref_seq.seq[pos:rstart]

            # For now, we won't worry about tyring to merge overlapping contig
            # ends. We will just bookend them into the sequence. This might be
            # an area to revisit in the future.
            contig_start = len(patched_seq)
            patched_seq = patched_seq + Seq(contig.query_sequence)

            # Gather what we need to write BED coordinates for patch and contig
            qstrand = "+"
            if contig.is_reverse:
                # Since BAM sequence is already reverse-complemented, we don't have to!
                qstrand = "-"
                
            # Write patch coordinates in reference frame to patches_bed
            patches_bed.write("%s\t%d\t%d\n" % (ref_seq.id, pos, rstart))
            # Write contig coordinates in patched-genome frame to contigs_bed
            contigs_bed.write("%s\t%d\t%d\t%s\t%s\n" % (ref_seq.id, contig_start, len(patched_seq), contig.query_name, qstrand))

            # Update the current position in the reference sequence.
            pos = contig_breakpoints[contig.query_name][2]
            
        # Once the above loop finishes, we need to add the terminal segment
        # from the reference genome.
        # It is important to note that, if no contigs were placed on the current
        # reference chromosome, the value of 'pos' will remain zero, so patched_seq
        # ends up being identical to the reference sequence in the output.
        if pos == 0:
            sys.stderr.write("No contigs map to %s. Printing the reference sequence to output unchanged.\n" % (ref_seq.id))
        patched_seq = patched_seq + ref_seq.seq[pos:len(ref_seq.seq)]
        patches_bed.write("%s\t%d\t%d\n" % (ref_seq.id, pos, len(ref_seq.seq)))

        # Swap in the patched sequence for the reference sequence and print the result 
        ref_seq.seq = patched_seq
        patched_fasta.write("%s\n" % (ref_seq.format("fasta")))

    # Close open file handles, etc.
    sorted_primary_alignments.close()
    patched_fasta.close()
    patches_bed.close()
    contigs_bed.close()

    # Clean up any leftover files as needed.
    if args.store_final_bam is None:
        os.remove(sorted_primary_alignments_bam)
        os.remove(sorted_primary_alignments_bam + '.bai')
    
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
