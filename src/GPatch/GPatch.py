#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
import random

__version__ = "0.3.4"

"""
Given a name-sorted BAM file, cluster mapped reads by name and
use the maximal interval defined by the primary alignment's
cigar string to infer the contig break-points in the reference 
genome.
"""

def create_primary_alignments_list(query_bam, min_qual_score):
    """
    Extract all primary alignments from the BAM file and return
    as a list.
    """
    ret = []
    for aln in query_bam:
        if aln.is_unmapped:
            # Skip unmapped contigs
            continue
        if aln.mapping_quality < min_qual_score:
            # Skip low-quality mappings
            continue
        if not (aln.is_supplementary or aln.is_secondary):
            ret.append(aln)
    return ret


def contig_breakpoints_from_cigar(primary_alignments):
    """
    Infer the contig breakpoints using soft-clips on the primary
    alignments.
    """
    contig_breakpoints = {}
    for contig in primary_alignments:
        # Default to mapped coordinates
        chrom = contig.reference_name
        start = contig.reference_start
        end = contig.reference_end
        # Check for left and right soft-clips
        if contig.cigartuples[0][0] == 4:
            # left soft-clip
            start -= contig.cigartuples[0][1]
        if contig.cigartuples[-1][0] == 4:
            # right soft-clip
            end += contig.cigartuples[-1][1]
        contig_breakpoints[contig.query_name] = [chrom, start, end, contig.query_name]

    return contig_breakpoints
    

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


def whitelist_filter(sorted_alignments_bam, whitelist, sorted_bam_name=None):
    """
    Given a sorted BAM file, use a system call to SAMTools to
    extract only alignments overlapping intervals within the
    whitelist BED.
    """    
    tmp_root = ''.join(random.sample('0123456789', 10))
    tmp_bam_name = tmp_root + '.tmp.bam'
    if sorted_bam_name is None:
        sorted_bam_name = tmp_root + '.tmp.sorted.bam'
    # Because pysam.view is buggy about creating output files with -o, we need
    # to touch the output file before calling pysam view.
    pb = open(tmp_bam_name, 'w')
    pb.close()
    pysam.view("-L", whitelist, "-h", "-b", "-o", tmp_bam_name, sorted_alignments_bam, catch_stdout=False)
    # Sort and index the result. Not sure if this is really necessary,
    # since input is already sorted. In principle, sort order may not be
    # maintainted, though -- not sure this is guaranteed.
    pysam.sort("-o", sorted_bam_name, tmp_bam_name)
    pysam.index(sorted_bam_name)
    os.remove(tmp_bam_name)
    return sorted_bam_name


def contigs_to_sorted_interval_list(contigs, contig_breakpoints):
    """
    Given an iterator over contig alignments from pysam.fetch,
    generate a list of intervals, sorted on position. Alignments are
    assumed to be on the same chromosome.
    """
    # Push contig alignment space start, end, and name to a list
    ints = []
    for contig in contigs:
        start = contig_breakpoints[contig.query_name][1]
        end = contig_breakpoints[contig.query_name][2]
        name = contig.query_name
        ints.append((start, end, name))
    # Sort the list by position
    ret = sorted(ints, key=lambda x: (x[0], x[1]))
    return ret


def cluster_contig_positions(sorted_contig_intervals, drop_nested=False):
    """
    Given a list of contig intervals, sorted on start then end
    position, cluster intervals based on linear overlap and
    return a list of clusters, optionally leaving out any nested
    intervals.
    """
    clusters = []
    for contig in sorted_contig_intervals:
        #sys.stderr.write("{0}\n".format(contig))
        merged = False
        for cluster in clusters:
            #sys.stderr.write("\t%s\n" % (cluster))
            if contig[0] <= cluster[-1][1]:
                if drop_nested:
                    if not contig[1] <= cluster[-1][1]:
                        cluster.append(contig)
                    else:
                        #sys.stderr.write("\tnested\n")
                        pass
                else:
                    cluster.append(contig)
                #sys.stderr.write("\tmerged\n")
                merged = True
                break
        if not merged:
            clusters.append([contig])
    return clusters


def main():
    parser = ArgumentParser(description="Starting with alignments of contigs to a reference genome, produce a chromosome-scale pseudoassembly by patching gaps between mapped contigs with sequences from the reference.")
    parser.add_argument('-q', '--query_bam', metavar='SAM/BAM', type=str,
                        required=True, help='Path to SAM/BAM file containing non-overlapping contig mappings to the reference genome.')
    parser.add_argument('-r', '--reference_fasta', metavar='FASTA', type=str,
                        required=True, help='Path to reference genome fasta.')
    parser.add_argument('-x', '--prefix', metavar='STR', type=str,
                        required=False, default="",
                        help='Prefix to add to output file names. Default=None')
    parser.add_argument('-b', '--store_final_bam', metavar='FILENAME', type=str,
                        required=False, default=None,
                        help='Store the final set of primary contig alignments to the given file name. Default: Do not store the final BAM.')
    parser.add_argument('-m', '--min_qual_score', metavar='N', type=int,
                        required=False, default=30,
                        help='Minimum mapping quality score to retain an alignment. Default=30')
    parser.add_argument('-w', '--whitelist', metavar='PATH', type=str,
                        required=False, default=None,
                        help='Path to BED file containing whitelist regions: i.e., the inverse of blacklist regions. Supplying this will have the effect of excluding alignments that fall entirely within blacklist regions. Default=None')
    parser.add_argument('-d', '--drop_missing', type=bool,
                        required=False, default=False,
                        help='Omit unpatched reference chromosome records from the output if no contigs map to them. Default: Unpatched chromosomes are printed to output unchanged.')
    
    args = parser.parse_args()

    query_bam = pysam.AlignmentFile(args.query_bam)
    
    # Load up the BAM into a dict, keyed on contig name.
    sys.stderr.write("Loading alignments...\n")
    primary_alignments = create_primary_alignments_list(query_bam, args.min_qual_score)

    # Determine insertion break-points for each contig based
    # on the cigar string for the primary alignment.
    sys.stderr.write("Locating contig breakpoints...\n")
    contig_breakpoints = contig_breakpoints_from_cigar(primary_alignments)

    # Sort the useful primary alignments by position via pysam sort.
    sys.stderr.write("Sorting primary alignments...\n")
    sorted_primary_alignments_bam = sorted_bam_from_aln_list(primary_alignments, query_bam.header, sorted_bam_name=args.store_final_bam)
    pysam.index(sorted_primary_alignments_bam)

    # Exclude alignments entirely within blacklist regions if requested.
    if args.whitelist is not None:
        sys.stderr.write("Filtering against the whitelist...\n")
        filtered_bam = whitelist_filter(sorted_primary_alignments_bam, args.whitelist)
        # Swap the filtered bam in for the sorted primary alignments bam
        os.remove(sorted_primary_alignments_bam)
        os.remove(sorted_primary_alignments_bam + ".bai")
        os.rename(filtered_bam, sorted_primary_alignments_bam)
        os.rename(filtered_bam + '.bai', sorted_primary_alignments_bam + ".bai")

    # Check for nested/overlapping alignment spaces among the remaining contigs.
    sys.stderr.write("Checking for nested or overlapping contigs...\n")
    sorted_primary_alignments = pysam.AlignmentFile(sorted_primary_alignments_bam, "rb")
    final_alignments = []
    for ref_seq in SeqIO.parse(args.reference_fasta, "fasta"):
        # Select all contigs mapped to this sequence.
        contigs = list(sorted_primary_alignments.fetch(ref_seq.id))

        # Convert the list of contigs into a list of tuples: (start, end, name),
        # representing the inferred contig breakpoints, sorted by start, then end position.
        sorted_contig_intervals = contigs_to_sorted_interval_list(contigs, contig_breakpoints)
        
        # Cluster the linear intervals on overlap, leaving out nested contigs
        clusters = cluster_contig_positions(sorted_contig_intervals, drop_nested=True)

        # Next thing to do is to push all retained contigs across clusters into
        # a final list of useful contigs. Since pysam.fetch does not allow
        # retrieval by alignment name, the easiest way to do this is through a
        # dict of names that will allow us to make a single pass through the
        # contig list to select what we need.
        retained_alignments = {}
        for cluster in clusters:
            for interval in cluster:
                retained_alignments[interval[2]] = interval
        for contig in contigs:
            if contig.query_name in retained_alignments:
                final_alignments.append(contig)

    sorted_primary_alignments.close()
                    
    # Write the final set of alignments to a sorted BAM file
    sorted_primary_alignments_bam = sorted_bam_from_aln_list(final_alignments, query_bam.header, sorted_bam_name=args.store_final_bam)                
    pysam.index(sorted_primary_alignments_bam)
        
    # Set up output streams for fasta, patches.bed, and contigs.bed
    sys.stderr.write("Patching the genome...\n")
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
        contigs = list(sorted_primary_alignments.fetch(ref_seq.id))

        # If no contigs map to this sequence, either print it unchanged
        # or omit it, depending on args.
        if len(contigs) == 0:
            if args.drop_missing:
                sys.stderr.write("No contigs map to %s. Omitting from output.\n" % (ref_seq.id))
            else:
                sys.stderr.write("No contigs map to %s. Printing the reference sequence to output unchanged.\n" % (ref_seq.id))
                patched_fasta.write("%s\n" % (ref_seq.format("fasta")))
            continue
        
        pos = 0 # Tracks position on the reference chromosome        
        patched_seq = ""
        for contig in contigs:
            rstart = contig_breakpoints[contig.query_name][1]
            if rstart > 0:
                rstart -= 1

            # Check for overlapping/bookended contig mappings.
            qstart = 0  # Start of contig interval to append
            if rstart > pos:
                # No overlap. Append the patch, in all lower-case
                # for easy identification.
                patched_seq = patched_seq + ref_seq.seq[pos:rstart].lower()
                # Write patch coordinates in reference frame to patches_bed
                patches_bed.write("%s\t%d\t%d\n" % (ref_seq.id, pos, rstart))
            else:
                # Handle overlapping contig ends by trimming the 5' end of
                # this contig sequence by the length of the overlap.
                qstart = pos - rstart
                
            # Append the contig sequence to the patched sequence string
            contig_start = len(patched_seq)
            # All patch sequences are printed in upper case.
            # Overlaps with the previous contig are handled by truncating
            # the 3' end of this contig sequence.
            patched_seq = patched_seq + Seq(contig.query_sequence)[qstart:len(Seq(contig.query_sequence))].upper()

            # Gather what we need to write BED coordinates for this contig
            qstrand = "+"
            if contig.is_reverse:
                # Since BAM sequence is already reverse-complemented, we don't have to!
                qstrand = "-"
                
            # Write contig coordinates in patched-genome frame to contigs_bed
            contigs_bed.write("%s\t%d\t%d\t%s\t.\t%s\t%d\t%d\n" % (ref_seq.id, contig_start, len(patched_seq), contig.query_name, qstrand, qstart, len(Seq(contig.query_sequence))))

            # Update the current position in the reference sequence if the
            # end position of the current contig is 3' of the current pos.
            if contig_breakpoints[contig.query_name][2] > pos:
                # This should always evaluate true in the absence of nested contigs.
                pos = contig_breakpoints[contig.query_name][2]
            
        # Once the above loop finishes, we need to add the terminal segment
        # from the reference genome.
        patched_seq = patched_seq + ref_seq.seq[pos:len(ref_seq.seq)]
        patches_bed.write("%s\t%d\t%d\n" % (ref_seq.id, pos, len(ref_seq.seq)))

        # Swap in the patched sequence for the reference sequence and print the result 
        ref_seq.seq = patched_seq
        patched_fasta.write("%s\n" % (ref_seq.format("fasta")))

    # Close open file handles, etc.
    query_bam.close()
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
