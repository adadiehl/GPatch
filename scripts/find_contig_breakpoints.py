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
    for aln in paf[contig[0]]:
        # loop over alignments from the same query chromosome.
        #sys.stderr.write("%s\n" % (aln))
        if overlaps(contig, aln):
            # Query interval overlap the contig's mapped interval
            ret.append(aln)
    return ret
    

def overlaps(contig, aln):
    """
    Check for simple overlap between contig ends alignment query interval.
    """
    if (int(aln[2]) >= int(contig[1]) and int(aln[2]) <= int(contig[2])) or (int(aln[3]) >= int(contig[1]) and int(aln[3]) <= int(contig[2])):
        return True
    return False


def cluster_alignments(alignments, max_cluster_dist = 0):
    """
    Cluster linear intervals closer together than max_cluster_dist.
    Reports a pseudo-alignment that spans the maximal query interval
    defined by each cluster.
    """
    sorted_alignments = sorted(alignments, key = lambda x: (int(x[2]), int(x[3])))
    clusters = []
    for aln in sorted_alignments:
        merged = False
        for cluster in clusters:
            if int(aln[2]) <= int(cluster[-1][3]) + max_cluster_dist:
                cluster.append(aln)
                merged = True
                break
        if not merged:
            clusters.append([aln])
    #sys.stderr.write("%s\n" % (clusters))
    return clusters
            

def cluster_on_mapped_then_query_pos(alignments, max_cluster_dist = 0, max_query_dist = 10000):
    """
    Given a set of alignments, create clusters based on mapped
    position, query position, and strand, such that all alignments
    in the cluster likely come from the same event. Assumes all
    alignments are between the same query and target chromosome.
    """
    sorted_alignments = sorted(alignments, key = lambda x: (int(x[7]), int(x[8])))
    clusters = []
    for aln in sorted_alignments:
        #sys.stderr.write("%s\n" % (aln))
        merged = False
        for idx, cluster in enumerate(clusters):
            #sys.stderr.write("\n%d Before: %s\n" % (idx,cluster))
            #sys.stderr.write("%s\n" % (aln))
            #sys.stderr.write("%d\t%d\t%d\n" % (int(aln[7]), int(cluster[-1][8]), (int(cluster[-1][8]) + max_cluster_dist)))
            #sys.stderr.write("%s\t%s\n" % (aln[4], cluster[-1][4]))
            if int(aln[7]) <= (int(cluster[-1][8]) + max_cluster_dist):
                # Intervals overlap. Check for nesting, then strand.
                #sys.stderr.write("%s\n" % (aln))
                if int(aln[8]) <= int(cluster[-1][8]):
                    # Nested. Don't use!
                    #sys.stderr.write("nested\n\n")
                    merged = True
                    break
                elif aln[4] == cluster[-1][4]:
                    # Same strand. Check distance between query positions.
                    if int(aln[2]) - int(cluster[-1][3]) <= max_query_dist:
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
    
    args = parser.parse_args()

    # Load up all partial alignments between the same query and target chromosome
    paf = parse_paf(args.paf, matched_only = True)

    with open(args.contigs, "r") as contigs:
        for line in contigs:
            contig = line.strip().split()

            # Gather overlapping alignments on the same query and target chromosome
            alignments = find_overlapping_alignments(contig, paf)
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
                starts = []
                ends = []
                for aln in cluster:
                    starts.append(int(aln[2]))
                    ends.append(int(aln[3]))
                starts.sort()
                ends.sort()
                ldist = abs(starts[0] - int(contig[1])) 
                rdist = abs(ends[-1] - int(contig[2]))
                if ldist <= args.max_dist or rdist <= args.max_dist:
                    suspicious_clusters.append(cluster)

            # Check for overlapping/duplicated query intervals. Sometimes this happens when
            # things multipmap. Note this keeps the first cluster, regardless of if it's the
            # maximal interval for the given query window.
            final_clusters = []
            for cluster in suspicious_clusters:
                keep = True
                for fcluster in final_clusters:
                    if cluster[0][2] <= fcluster[-1][3]:
                        keep = False
                        break
                if keep == True:
                    final_clusters.append(cluster)
                    
                    
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
                    starts = []
                    ends = []
                    for aln in cluster:
                        starts.append(int(aln[2]))
                        ends.append(int(aln[3]))
                    starts.sort()
                    ends.sort()
                    sys.stdout.write("%s\t%d\t%d\t%s\n" % ("\t".join(contig), starts[0], ends[-1], code))
                
                
            """
            suspicious_alignments = []
            for aln in alignments:
                #sys.stderr.write("%s\n" % (aln))
                ldist = abs(int(aln[2]) - int(contig[1]))   # Distance from left end of 3'-most partial alignment to left end of contig
                rdist = abs(int(aln[3]) - int(contig[2]))   # Distance from right end of 5'-most partial alignment to right end of contig
                if ldist <= args.max_dist or rdist <= args.max_dist:
                    if aln[0] != aln[5]:
                        # Translocation on different chromosome
                        code = "TRD"
                        suspicious_alignments.append(aln)
                    else:
                        if aln[4] == "-":
                            # Inversion
                            code = "INV"
                            suspicious_alignments.append(aln)
                        elif abs(int(aln[2]) - int(aln[7])) >= args.min_trans_dist:
                            code = "TRS"
                            suspicious_alignments.append(aln)

            if len(suspicious_alignments) > 0:
                starts = []
                ends = []
                for aln in suspicious_alignments:
                    starts.append(int(aln[2]))
                    ends.append(int(aln[3]))
                starts.sort()
                ends.sort()
                sys.stdout.write("%s\t%d\t%d\t%s\n" % ("\t".join(contig), starts[0], ends[-1], code))
            """

            """
            # Group overlapping alignments into linear clusters
            clustered_alignments = cluster_alignments(alignments, max_cluster_dist = args.max_cluster_dist)
            
            
            for cluster in clustered_alignments:
                #sys.stderr.write("%s\n" % (cluster))
                starts = []
                ends = []
                for aln in cluster:
                    #sys.stderr.write("%s\n" % (aln))
                    starts.append(int(aln[2]))
                    ends.append(int(aln[3]))
                ldist = abs(starts[0] - int(contig[1]))   # Distance from left end of 3'-most partial alignment to left end of contig
                rdist = abs(ends[-1] - int(contig[2]))   # Distance from right end of 5'-most partial alignment to right end of contig
                #sys.stderr.write("%d\t%d\n" % (ldist, rdist))
                if ldist <= args.max_dist or rdist <= args.max_dist:
                    # Breakpoint near the end of a conti             
                    sys.stderr.write("%s\n" % (cluster))
                    sys.stdout.write("%s\t%d\t%d\n" % ("\t".join(contig), starts[0], ends[-1]))                
            """
            
    
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
