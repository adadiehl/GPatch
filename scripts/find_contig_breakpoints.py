#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "0.0.1"

"""
Examine contigs.bed and a paf alignment to the reference genome from the given 
patched genome to identify potentially spurious large-scale rearrangement
breakpoints for contig breaking purposes.
"""

def parse_paf(paf_f, matched_only = False):
    """
    Parse the PAF alignment file into a dict keyed on query chromosome.
    If matched_only == True, will return only alignents between the same
    query and target chromosomes.
    """
    ret = {}
    with open(paf_f, "r") as paf:
        for line in paf:
            rec = line.strip().split()
            #sys.stderr.write("%s\n" % (rec))
            if matched_only and rec[0] != rec[5]:
                continue
            if rec[0] in ret:
                ret[rec[0]].append(rec)
            else:
                ret[rec[0]] = [rec]
    return ret
    

def find_overlapping_alignments(contig, paf):
    """
    Find all partial alignments overlapping the given contig and
    those on different chromosomes.
    """
    ret = []
    if contig[0] not in paf:
        return ret
    for aln in paf[contig[0]]:
        # loop over alignments from the same query chromosome.
        #sys.stderr.write("%s\n" % (aln))
        if overlaps(contig, aln):
            # Query interval overlap the contig's mapped interval
            ret.append(aln)
    return ret


def get_sorted_starts_ends(cluster):
    """
    Get sorted start and end positions for a cluster
    """
    starts = []
    ends = []
    for aln in cluster:
        starts.append(int(aln[2]))
        ends.append(int(aln[3]))
    starts.sort()
    ends.sort()
    return starts, ends
    

def overlaps(contig, aln):
    """
    Check for simple overlap between contig ends alignment query interval.
    """
    if (int(aln[2]) >= int(contig[1]) and int(aln[2]) <= int(contig[2])) or (int(aln[3]) >= int(contig[1]) and int(aln[3]) <= int(contig[2])):
        return True
    return False


def cluster_on_mapped_then_query_pos(alignments, max_cluster_dist = 0, max_query_dist = 10000):
    """
    Given a set of alignments, create clusters based on mapped
    position, query position, and strand, such that all alignments
    in the cluster likely come from the same event. Assumes all
    alignments are between the same query and target chromosome.
    """
    sorted_alignments = sorted(alignments, key = lambda x: (x[5], int(x[7]), int(x[8])))
    #sys.stderr.write("%s\n\n" % (sorted_alignments))
    clusters = []
    for aln in sorted_alignments:
        #sys.stderr.write("%s\n" % (aln))
        merged = False
        for idx, cluster in enumerate(clusters):
            #sys.stderr.write("\n%d Before: %s\n" % (idx,cluster))
            #sys.stderr.write("%s\n" % (aln))
            #sys.stderr.write("%d\t%d\t%d\n" % (int(aln[7]), int(cluster[-1][8]), (int(cluster[-1][8]) + max_cluster_dist)))
            #sys.stderr.write("%s\n" % (cluster[-1]))
            #sys.stderr.write("%s\t%s\n" % (aln[4], cluster[-1][4]))
            if aln[5] != cluster[-1][5]:
                # Alignment to a different chromosome.
                merged = False
            elif int(aln[7]) <= (int(cluster[-1][8]) + max_cluster_dist) and aln[5] == cluster[-1][5]:
                # Mapped intervals overlap. Check for nesting, then strand.
                if int(aln[8]) <= int(cluster[-1][8]) and aln[4] == cluster[-1][4]:
                    # Nested. Need to do some more checks to see if we want to merge
                    # as rearrangments can sometimes "hide" inside nested mappings.
                    #sys.stderr.write("nested\n\n")
                    # Check to see if the alignment ends within the max_cluster_dist
                    if abs(int(cluster[-1][8]) - int(aln[8])) >= max_cluster_dist:
                        merged = False
                elif aln[4] == cluster[-1][4]:
                    # Same strand. Check distance between query positions.
                    if abs(int(aln[2]) - int(cluster[-1][3])) <= max_query_dist:
                        cluster.append(aln)
                        merged = True
                        #sys.stderr.write("After: %s\n\n" % (cluster))
                        break
            
        if not merged:
            clusters.append([aln])
            #sys.stderr.write("\nINIT: %s\n\n" % ([aln]))
    return clusters


def main():
    parser = ArgumentParser(description="Break contigs at breakpoints indicated in the given text file.")
    parser.add_argument('-c', '--contigs', metavar='BED', type=str,
                        required=True, help='Path to contigs.bed containing contig mappings within the current patched genome.')
    parser.add_argument('-p', '--paf', metavar='PATH', type=str,
                        required=True, help='PAF alignment of the current patched genome to the reference genome.')
    parser.add_argument('-m', '--max_dist', metavar='INT', type=int,
                        required=False, default=1000,
                        help='Maximum distance, in bp, from the contig end to a potential breakpoint to consider it likely spurious and include this in breakpoints.txt. Default=1000')
    parser.add_argument('-t', '--min_trans_dist', metavar='INT', type=int,
                        required=False, default=1000000,
                        help='Minimum distance, in bp, from the query position to the mapped position in the reference to call a likely translocation and include this alignment in breakpoints.txt. Default=1000000')
    parser.add_argument('-d', '--max_cluster_dist', metavar='INT', type=int,
                        required=False, default=0,
                        help='Maximum distance between mapped positions to merge a partial alignment into a cluster. default = 0')
    parser.add_argument('-q', '--max_query_dist', metavar='INT', type=int,
                        required=False, default=10000,
                        help='Maximum distance between query intervals to merge a partial alignment into a cluster. default = 10000')
    parser.add_argument('-s', '--min_size', metavar='INT', type=int,
                        required=False, default=1000000,
                        help='Minimum size, in bp, to call an event and include this alignment in breakpoints.txt. Default=1000000')
    parser.add_argument('-a', '--all_chroms', default=False, action="store_true",
                        required=False, help="Find breakpoints for rearrangements across all chromosomes, not just intrachromosomal rearrangements. Default=Find only intrachromosomal breakpoints.")
    
    args = parser.parse_args()

    # Load up all partial alignments overlapping the contig
    matched_only = True
    if args.all_chroms == True:
        matched_only = False
    paf = parse_paf(args.paf, matched_only = matched_only)

    with open(args.contigs, "r") as contigs:
        for line in contigs:
            contig = line.strip().split()

            # Gather overlapping alignments on the same query and target chromosome
            alignments = find_overlapping_alignments(contig, paf)
            if len(alignments) == 0:
                continue
            #sys.stderr.write("%s\n" % (contig))
            #sys.stderr.write("%s\n\n" % (alignments))

            # Cluster alignments based on mapped position, then query strand and position
            # in an attempt to group alignments by event.
            clustered_alignments = cluster_on_mapped_then_query_pos(alignments, max_cluster_dist = args.max_cluster_dist, max_query_dist = args.max_query_dist)
            #for cluster in clustered_alignments:
            #    sys.stderr.write("%s\n" % (cluster))

            # Go through clusters to isolate those with partial overlap of the right and/or
            # left end of the contig (terminating within max_dist of the contig end)
            suspicious_clusters = []
            for cluster in clustered_alignments:
                #sys.stderr.write("%s\n" % (cluster))
                # Sort query starts and ends independently to account for possible nested intervals
                starts, ends = get_sorted_starts_ends(cluster)
                ldist = abs(starts[0] - int(contig[1])) 
                rdist = abs(ends[-1] - int(contig[2]))
                length = ends[-1] - starts[0]
                #sys.stderr.write("%d\t%d\n" % (ldist, rdist))
                if (ldist <= args.max_dist or rdist <= args.max_dist) and length >= args.min_size:
                    #sys.stderr.write("%s\n" % (cluster))
                    suspicious_clusters.append(cluster)

            """
            This turns out to cause problems with missed rearrangements in more-complex
            rearrangement scenarios, so skip it!            
            # Check for overlapping/duplicated query intervals. Sometimes this happens when
            # things multimap or in more-complex scenarios (e.g., breakpoints based on multiple
            # mapped clusters. When this happens, we will adjust the start position of the
            # overlapping cluster.
            final_clusters = []
            for cluster in suspicious_clusters:
                sys.stderr.write("%s\n" % (cluster))
                keep = True
                cstarts, cends = get_sorted_starts_ends(cluster)
                sys.stderr.write("%s\t%s\n" % (cstarts, cends))
                for fcluster in final_clusters:                    
                    fstarts, fends = get_sorted_starts_ends(fcluster)
                    sys.stderr.write("%s\t%s\n" % (fstarts, fends))
                    if cstarts[0] <= fends[-1]:
                        keep = False
                        break
                if keep == True:
                    sys.stderr.write("%s\n" % (cluster))
                    final_clusters.append(cluster)
            """
            final_clusters = suspicious_clusters
                    
            # Next check to be sure this is a rearrangement breakpoint based on mapped strand
            # (inversions) and distance between query start and mapped start (translocations)
            for cluster in final_clusters:
                #sys.stderr.write("%s\n" % (cluster))
                is_bp = False
                code = ""
                if cluster[0][4] == "-" and abs(int(cluster[0][2]) - int(cluster[0][7])) >= args.min_trans_dist:
                    # Inverted translocation
                    is_bp = True
                    code = "INT"
                elif cluster[0][4] == "-":
                    # Inversion
                    is_bp = True
                    code = "INV"
                elif abs(int(cluster[0][2]) - int(cluster[0][7])) >= args.min_trans_dist:
                    # Transversion
                    is_bp = True
                    code = "TRS"
                if is_bp:
                    starts, ends = get_sorted_starts_ends(cluster)
                    sys.stdout.write("%s\t%d\t%d\t%s\n" % ("\t".join(contig), starts[0], ends[-1], code))
                
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
