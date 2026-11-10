#!/usr/bin/env python3

import sys, os, re
import bisect, shutil, subprocess, tempfile
from concurrent.futures import ProcessPoolExecutor
from argparse import ArgumentParser, SUPPRESS
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
import random

__version__ = "0.5.0"

"""
Given a name-sorted BAM file, cluster mapped reads by name and
use the maximal interval defined by the primary alignment's
cigar string to infer the contig break-points in the reference 
genome.
"""

def create_primary_alignments_list(query_bam, min_qual_score, dropped=None):
    """
    Extract all primary alignments from the BAM file and return
    as a list. If a dict is given as dropped, contigs excluded here
    are recorded in it as name -> (reason, detail).
    """
    #i = 0
    #j = 0
    #k = 0
    #l = 0
    ret = []
    for aln in query_bam:
        #i += 1
        if aln.is_unmapped:
            # Skip unmapped contigs
            #j += 1
            if dropped is not None:
                dropped[aln.query_name] = ("unmapped", "")
            continue
        if aln.mapping_quality < min_qual_score:
            # Skip low-quality mappings
            #k += 1
            if dropped is not None and not (aln.is_supplementary or aln.is_secondary):
                dropped[aln.query_name] = ("low_mapq", "MAPQ %d < %d" % (aln.mapping_quality, min_qual_score))
            continue
        if not (aln.is_supplementary or aln.is_secondary):
            #l += 1
            ret.append(aln)
    #sys.stderr.write("%d\t%d\t%d\t%d\t%d\n" % (i, j, k, l, len(ret)))
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
        #sys.stderr.write("%s\t%s\t%s\n" % (chrom, start, end))
        # Check for left and right soft-clips
        #sys.stderr.write("%s\n%s\n" % (contig.cigartuples[0], contig.cigartuples[-1]))
        if contig.cigartuples[0][0] == 4:
            # left soft-clip
            start -= contig.cigartuples[0][1]
        if contig.cigartuples[-1][0] == 4:
            # right soft-clip
            end += contig.cigartuples[-1][1]
        contig_breakpoints[contig.query_name] = [chrom, start, end, contig.query_name]
        #sys.stderr.write("%s\n" % (contig_breakpoints))

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


def cluster_contig_positions(sorted_contig_intervals, drop_nested=False, nested=None):
    """
    Given a list of contig intervals, sorted on start then end
    position, cluster intervals based on linear overlap and
    return a list of clusters, optionally leaving out any nested
    intervals. If a dict is given as nested, dropped nested contigs
    are recorded in it as name -> name of the containing contig.
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
                        #sys.stderr.write("\t%s\n\t%s\n" % (cluster[-1], contig))
                        if nested is not None:
                            nested[contig[2]] = cluster[-1][2]
                else:
                    cluster.append(contig)
                #sys.stderr.write("\tmerged\n")
                merged = True
                break
        if not merged:
            clusters.append([contig])
    return clusters


"""
Junction handling.

For each pair of adjacent placed contigs, GPatch decides whether to trim the
5' end of the second contig (the contigs overlap), patch the gap between them
with reference sequence, or bookend them. N characters are never introduced
(except in --scaffold_only mode, where patches are replaced by N runs).

1) Contig extents on the reference. The primary alignment is chained to
   supplementary alignments of the same contig on the same chromosome and
   strand that are colinear with it in both query and reference. Each contig
   end is extrapolated 1:1 only by the residual unaligned query sequence
   beyond the outermost chained alignment.

2) Direct overlap measurement. Candidate overlaps between the end of the
   previous contig and the start of the next (both as placed) come from exact
   matches of two anchors: the next contig's first bases in the previous
   contig's tail, and the previous contig's last bases in the next contig's
   head. Candidates supported by both anchors at the same offset (or within
   --indel_tol, verified with minimap2) are accepted when identity over the
   whole overlap is >= --min_overlap_identity. A general-purpose aligner is
   not used to search for overlaps because, in palindromes and tandem arrays,
   a short true overlap is reported as a low-scoring secondary of a longer
   paralogous alignment and filtered out.

The extent-based gap estimate is "trusted" when the residual unaligned
sequence at both contig ends of a junction is <= --max_resid, "one-sided" when
only one end exceeds it, and "untrusted" otherwise. Decisions:
  trim_measured    one overlap found -> trim it
  trim_ambiguous   several overlap offsets (tandem repeats) -> the one closest
                   to the estimated overlap if the estimate is trusted or
                   one-sided, otherwise the smallest
  patch            no overlap; estimate >= 0 and trusted -> reference patch
  patch_one_sided  no overlap; estimate >= 0 and one-sided -> reference patch
  bookend          no overlap; estimate negative or untrusted
Chromosome ends are patched when the residual at that contig end is
<= --max_resid; contigs are never trimmed against a chromosome end.
"""


class AlignedSpan(object):
    """Reference and query (placement-orientation) coordinates of one alignment."""
    __slots__ = ("chrom", "rs", "re", "qs", "qe", "qlen", "is_reverse", "mapq")

    def __init__(self, aln):
        lead = 0
        for op, n in aln.cigartuples:
            if op in (4, 5):   # soft or hard clip
                lead += n
            else:
                break
        self.chrom = aln.reference_name
        self.rs = aln.reference_start
        self.re = aln.reference_end
        self.qlen = aln.infer_read_length()
        self.qs = lead
        self.qe = lead + aln.query_alignment_length
        self.is_reverse = aln.is_reverse
        self.mapq = aln.mapping_quality


def collect_supplementary_spans(bam_path, names, min_mapq, min_len):
    """Supplementary alignments usable for chaining, keyed on contig name."""
    supps = {}
    with pysam.AlignmentFile(bam_path) as bam:
        for aln in bam:
            if not aln.is_supplementary or aln.is_unmapped or aln.query_name not in names:
                continue
            if aln.mapping_quality < min_mapq or aln.query_alignment_length < min_len:
                continue
            supps.setdefault(aln.query_name, []).append(AlignedSpan(aln))
    return supps


def chained_extent(primary, supps, maxgap, tol):
    """
    Chain colinear supplementary alignments onto the primary alignment and
    return the extrapolated reference extent of the contig together with the
    residual unaligned query length at each end.
    """
    cands = [s for s in supps if s.chrom == primary.chrom and s.is_reverse == primary.is_reverse]
    left = right = primary
    n = 1
    while True:
        best = None
        for s in cands:
            dq, dr = left.qs - s.qe, left.rs - s.re
            if s.qs < left.qs and s.rs < left.rs and -tol <= dq <= maxgap and -tol <= dr <= maxgap:
                if best is None or s.qs < best.qs:
                    best = s
        if best is None:
            break
        left, n = best, n + 1
    while True:
        best = None
        for s in cands:
            dq, dr = s.qs - right.qe, s.rs - right.re
            if s.qe > right.qe and s.re > right.re and -tol <= dq <= maxgap and -tol <= dr <= maxgap:
                if best is None or s.qe > best.qe:
                    best = s
        if best is None:
            break
        right, n = best, n + 1
    return {"start": left.rs - left.qs, "end": right.re + (primary.qlen - right.qe),
            "resid5": left.qs, "resid3": primary.qlen - right.qe, "nchain": n}


def find_all(hay, needle):
    """All start positions of needle in hay."""
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def identity_at_least(x, y, min_id, chunk=1024):
    """True if equal-length strings x and y are >= min_id identical."""
    max_mm = int((1.0 - min_id) * len(x))
    mm = 0
    for i in range(0, len(x), chunk):
        a, b = x[i:i + chunk], y[i:i + chunk]
        if a != b:
            mm += sum(1 for c, d in zip(a, b) if c != d)
            if mm > max_mm:
                return False
    return True


def identity_minimap2(x, y, minimap2):
    """
    Identity of a near-equal-length sequence pair (an overlap containing
    indels) via minimap2; the alignment must cover >= 99% of both.
    """
    tmp = tempfile.mkdtemp(prefix="gpatch_")
    try:
        q, t = os.path.join(tmp, "q.fa"), os.path.join(tmp, "t.fa")
        with open(q, "w") as o:
            o.write(">q\n%s\n" % x)
        with open(t, "w") as o:
            o.write(">t\n%s\n" % y)
        out = subprocess.run([minimap2, "-x", "asm5", "-c", "-t", "1", t, q],
                             capture_output=True, text=True, check=True).stdout
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    best = 0.0
    for line in out.splitlines():
        f = line.split("\t")
        qlen, qs, qe, strand = int(f[1]), int(f[2]), int(f[3]), f[4]
        tlen, ts, te, nmatch, blen = int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
        if strand == "+" and qe - qs >= 0.99 * qlen and te - ts >= 0.99 * tlen:
            best = max(best, nmatch / blen)
    return best


def measure_overlaps(tail, head, prev_len, cur_len, L, p):
    """
    Anchored suffix/prefix overlap detection between the end of the previous
    contig (tail: its last <= win_max bases, as placed) and the start of the
    next contig (head: its first <= win_max bases). Returns
    (sorted candidate trims, final window size, note).
    """
    k = p["anchor"]
    note = ""
    while True:
        lp, lc = min(L, prev_len), min(L, cur_len)
        ps, cs = tail[len(tail) - lp:], head[:lc]
        ov1 = sorted(set(lp - i for i in find_all(ps, cs[:k])))
        ov2 = sorted(set(j + k for j in find_all(cs, ps[-k:])))
        s2 = set(ov2)
        trims, n_tested = [], 0
        for o in ov1:
            if o in s2:
                n_tested += 1
                if identity_at_least(ps[lp - o:], cs[:o], p["min_id"]):
                    trims.append(o)
        if not trims and p["indel_tol"] > 0 and p["minimap2"]:
            near = []
            for o1 in ov1:
                lo = bisect.bisect_left(ov2, o1 - p["indel_tol"])
                hi = bisect.bisect_right(ov2, o1 + p["indel_tol"])
                near.extend((abs(o1 - o2), o1, o2) for o2 in ov2[lo:hi] if o2 != o1)
            for _, o1, o2 in sorted(near)[:p["max_mm2_tests"]]:
                n_tested += 1
                if identity_minimap2(ps[lp - o1:], cs[:o2], p["minimap2"]) >= p["min_id"]:
                    trims.append(o2)   # trim = position in the next contig of the previous contig's end
        # The overlap may extend beyond the window if an anchor lies at the
        # far edge of it.
        if not trims and L < p["win_max"] and ((ov1 and ov1[-1] > lp - k and lp < prev_len) or
                                               (ov2 and ov2[-1] > lc - k and lc < cur_len)):
            L = min(2 * L, p["win_max"])
            note = "window doubled"
            continue
        trims = [t for t in sorted(set(trims)) if t < cur_len]
        note = ("%s; " % note if note else "") + "anchors %d/%d, tested %d" % (len(ov1), len(ov2), n_tested)
        return trims, L, note


def _measure_overlaps_job(job):
    return measure_overlaps(*job)


def decide_junction(prev_ext, cur_ext, pos, trims, max_resid):
    """
    Decide how to join two adjacent contigs. pos is the reference position of
    the end of the placed sequence so far. Returns (decision, trim, reason).
    """
    est = cur_ext["start"] - pos   # > 0 gap, < 0 overlap
    big = []
    if prev_ext["resid3"] > max_resid:
        big.append("prev residual 3' unaligned %d" % prev_ext["resid3"])
    if cur_ext["resid5"] > max_resid:
        big.append("residual 5' unaligned %d" % cur_ext["resid5"])
    trust = ("trusted", "one-sided", "untrusted")[len(big)]
    if len(trims) == 1:
        return "trim_measured", trims[0], "single overlap"
    if len(trims) > 1 and trust != "untrusted":
        t = min(trims, key=lambda x: (abs(x + est), x))
        return "trim_ambiguous", t, ("%d overlap offsets; closest to %s estimated overlap %d used"
                                     % (len(trims), trust, -est))
    if len(trims) > 1:
        return "trim_ambiguous", trims[0], "%d overlap offsets; estimate untrusted, smallest used" % len(trims)
    if est >= 0 and trust == "trusted":
        return "patch", 0, "no overlap; gap estimate trusted"
    if est >= 0 and trust == "one-sided":
        return "patch_one_sided", 0, "no overlap; one-sided gap estimate used (%s)" % big[0]
    rs = (["estimate predicts overlap of %d but none found" % -est] if est < 0 else []) + big
    return "bookend", 0, "no overlap; unusable gap estimate: " + ", ".join(rs)


def main():
    parser = ArgumentParser(description="Starting with alignments of contigs to a reference genome, produce a chromosome-scale pseudoassembly by patching gaps between mapped contigs with sequences from the reference. By default, reference chromosomes with no mapped contigs are printed to output unchanged. Use the --drop_missing option to disable this behavior. By default, patches are applied to the 5' and 3' telomere ends of pseudochromsomes if the first and last contig alignments do not extend to the start/end of the reference chromsome. In some cases, this may cause spurious duplications. Use the --no_extend option if this is a concern. Note that GPatch is designed to be run on single-haplotype or unphased genome assemblies. For phased assemblies, each haplotype should be separated into its own input FASTA file prior to alignment. GPatch can then be run separately on the BAM files for each haplotype to obtain phased pseudoassemblies, otherwise results will be unpredictable and likely incorrect.")
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
    parser.add_argument('-d', '--drop_missing',
                        required=False, default=False, action="store_true",
                        help='Omit unpatched reference chromosome records from the output if no contigs map to them. Default: Unpatched chromosomes are printed to output unchanged.')
    parser.add_argument('-t', '--no_trim',
                        required=False, default=False, action="store_true",
                        help='Do not trim the 5-prime end of contigs that overlap the previously-placed contig; overlapping contigs are bookended untrimmed. Default: Trim the overlap, measured directly from the contig sequences (or, with --legacy_junctions, inferred from the soft-clip-extended alignment coordinates).')
    parser.add_argument('-s', '--scaffold_only',
                        required=False, default=False, action="store_true",
                        help='Pad gaps between placed contigs with strings of N characters instead of patching with sequence from the reference assembly. Effectively turns GPatch into a reference-guided scaffolding tool. Note that patches.bed will still be generated to document (inverse) mapped contig boundaries in reference frame.')
    parser.add_argument('-l', '--gap_length', metavar='N', type=int,
                        required=False, default=-1,
                        help='Length of "N" gaps separating placed gontigs when using --scaffold_only. Has no effect when in default patching mode. Default=Estimate gap length from alignment.')
    parser.add_argument('-e', '--no_extend',
                        required=False, default=False, action="store_true",
                        help='Do not patch telomere ends of pseudochromosomes with reference sequence upstream of the first mapped contig and downstream of the last mapped contig. Default is to include 5\' and 3\' patches to extend telomeres to the ends implied by the alignment.')
    parser.add_argument('-k', '--keep_nested',
                        required=False, default=False, action="store_true",
                        help='Do not drop contigs with mapped positions nested entirely inside other mapped contigs. Instead, these will be bookended after the contig in which they are nested. Default is to drop contigs with mapped positions nested entirely within other mapped contigs. This option should be used with caution as these mappings cannot be placed unambigiously relative to other mapped contigs, thus including them is likely to lead to unpredictable and possibly incorrect results. Do not use this unless you are sure you know what you are doing!')
    junc = parser.add_argument_group("junction handling", "Options controlling how adjacent placed contigs are joined (trim, patch or bookend). See junctions.tsv for the decision made at each junction.")
    junc.add_argument('--legacy_junctions',
                      required=False, default=False, action="store_true",
                      help='Use the original (GPatch <= 0.4) junction logic: contig extents from the primary alignment extended by its full soft-clips, overlaps trimmed by the overlap of those extents, and gaps patched by their distance. No junctions.tsv is written. Default: chained extents and directly measured overlaps.')
    junc.add_argument('--max_resid', metavar='N', type=int, default=10000,
                      help='Maximum unaligned sequence (bp) remaining at a contig end, after chaining, for the reference gap estimate at that end to be trusted. Default=10000')
    junc.add_argument('--chain_mapq', metavar='N', type=int, default=20,
                      help='Minimum mapping quality of supplementary alignments used to extend contig extents. Default=20')
    junc.add_argument('--chain_minlen', metavar='N', type=int, default=10000,
                      help='Minimum aligned query length (bp) of supplementary alignments used to extend contig extents. Default=10000')
    junc.add_argument('--chain_maxgap', metavar='N', type=int, default=1000000,
                      help='Maximum gap (bp), in both query and reference, between chained alignments. Default=1000000')
    junc.add_argument('--chain_tol', metavar='N', type=int, default=10000,
                      help='Maximum overlap (bp) allowed between chained alignments. Default=10000')
    junc.add_argument('--anchor', metavar='N', type=int, default=500,
                      help='Length (bp) of the exact anchors used to detect overlaps between adjacent contigs; also the minimum detectable overlap. Default=500')
    junc.add_argument('--min_overlap_identity', metavar='F', type=float, default=0.999,
                      help='Minimum identity over the whole overlap for a candidate overlap to be accepted. Default=0.999')
    junc.add_argument('--indel_tol', metavar='N', type=int, default=50,
                      help='Maximum length difference (bp) between the two anchor-implied overlap lengths for a candidate containing indels; such candidates are verified with minimap2. 0 disables. Default=50')
    junc.add_argument('--max_mm2_tests', metavar='N', type=int, default=20,
                      help='Maximum number of indel-containing candidate overlaps verified with minimap2 per junction. Default=20')
    junc.add_argument('--win_min', metavar='N', type=int, default=1000000,
                      help='Minimum size (bp) of the contig-end windows searched for overlaps. Default=1000000')
    junc.add_argument('--win_max', metavar='N', type=int, default=4000000,
                      help='Maximum size (bp) of the contig-end windows searched for overlaps. Default=4000000')
    junc.add_argument('--minimap2', metavar='PATH', type=str, default=None,
                      help='Path to minimap2, used only to verify candidate overlaps containing indels. Default: minimap2 on PATH; if not found, such candidates are not considered.')
    junc.add_argument('--threads', metavar='N', type=int, default=1,
                      help='Number of processes used for overlap detection. Default=1')
    ext = parser.add_argument_group("overlap extension", "After alignment-based placement, unplaced contigs (low mapping quality or nested) are inserted into gaps by extending from the flanking placed contigs through unambiguous sequence overlaps. Not used with --legacy_junctions. See extension.tsv for every walk step and stop reason.")
    ext.add_argument('--no_extension',
                     required=False, default=False, action="store_true",
                     help='Do not insert unplaced contigs by overlap extension.')
    ext.add_argument('--ext_min_overlap', metavar='N', type=int, default=5000,
                     help='Minimum overlap (bp) for an extension step, or for a walk to be considered to have reached the opposite flanking contig. Short overlaps in tandem and palindromic repeats are often spurious. Default=5000')
    ext.add_argument('--ext_region_pad', metavar='N', type=int, default=1000000,
                     help='A candidate contig\'s primary alignment must lie within the gap region extended by this many bp beyond the flanking contigs\' alignments. Default=1000000')
    ext.add_argument('--ext_window', metavar='N', type=int, default=2000000,
                     help='Size (bp) of the contig-end windows searched for extension overlaps. Default=2000000')
    ext.add_argument('--ext_max_steps', metavar='N', type=int, default=50,
                     help='Maximum number of contigs added by one walk. Default=50')

    args = parser.parse_args()

    query_bam = pysam.AlignmentFile(args.query_bam)
    
    # Load up the BAM into a dict, keyed on contig name.
    sys.stderr.write("Loading alignments...\n")
    unplaced = {}   # contig name -> (reason, detail) for contigs left out of the pseudoassembly
    primary_alignments = create_primary_alignments_list(query_bam, args.min_qual_score, dropped=unplaced)

    # Determine insertion break-points for each contig: in legacy mode from
    # the soft-clips of the primary alignment; otherwise from the primary
    # alignment chained with colinear supplementary alignments, extrapolated
    # only by the remaining unaligned sequence. These extents are used for
    # the nested-contig check and for junction handling.
    sys.stderr.write("Locating contig breakpoints...\n")
    extents = None
    if args.legacy_junctions:
        contig_breakpoints = contig_breakpoints_from_cigar(primary_alignments)
    else:
        extents = compute_extents(args, primary_alignments)
        contig_breakpoints = {}
        for aln in primary_alignments:
            e = extents[aln.query_name]
            contig_breakpoints[aln.query_name] = [aln.reference_name, e["start"], e["end"], aln.query_name]

    # Sort the useful primary alignments by position via pysam sort.
    sys.stderr.write("Sorting primary alignments...\n")
    sorted_primary_alignments_bam = sorted_bam_from_aln_list(primary_alignments, query_bam.header, sorted_bam_name=args.store_final_bam)
    pysam.index(sorted_primary_alignments_bam)

    # Exclude alignments entirely within blacklist regions if requested.
    if args.whitelist is not None:
        sys.stderr.write("Filtering against the whitelist...\n")
        filtered_bam = whitelist_filter(sorted_primary_alignments_bam, args.whitelist)
        with pysam.AlignmentFile(filtered_bam) as fb:
            kept = set(aln.query_name for aln in fb)
        for aln in primary_alignments:
            if aln.query_name not in kept:
                unplaced[aln.query_name] = ("whitelist", "alignment entirely outside whitelist regions")
        # Swap the filtered bam in for the sorted primary alignments bam
        os.remove(sorted_primary_alignments_bam)
        os.remove(sorted_primary_alignments_bam + ".bai")
        os.rename(filtered_bam, sorted_primary_alignments_bam)
        os.rename(filtered_bam + '.bai', sorted_primary_alignments_bam + ".bai")

    # Check for nested/overlapping alignment spaces among the remaining contigs.
    if not args.keep_nested:
        sys.stderr.write("Checking for nested contigs...\n")
    sorted_primary_alignments = pysam.AlignmentFile(sorted_primary_alignments_bam, "rb")
    final_alignments = []
    for ref_seq in SeqIO.parse(args.reference_fasta, "fasta"):
        # Select all contigs mapped to this sequence.
        contigs = list(sorted_primary_alignments.fetch(ref_seq.id))
        #sys.stderr.write("%s\t%s\n" % (ref_seq.id, contigs))

        # Convert the list of contigs into a list of tuples: (start, end, name),
        # representing the inferred contig breakpoints, sorted by start, then end position.
        sorted_contig_intervals = contigs_to_sorted_interval_list(contigs, contig_breakpoints)
        
        # Cluster the linear intervals on overlap, leaving out nested contigs
        drop_nested = True
        if args.keep_nested:
            drop_nested = False
        nested = {}
        clusters = cluster_contig_positions(sorted_contig_intervals, drop_nested=drop_nested, nested=nested)
        for name, container in nested.items():
            unplaced[name] = ("nested", "extent within that of %s" % container)

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

    # Clean up the initial sorted BAM and index files.
    os.remove(sorted_primary_alignments_bam)
    os.remove(sorted_primary_alignments_bam + '.bai')
    
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

    # Placement order on each reference sequence: by alignment position, plus
    # (default mode) contigs inserted by overlap extension.
    placement = {}
    for chrom in sorted_primary_alignments.references:
        contigs = list(sorted_primary_alignments.fetch(chrom))
        if contigs:
            placement[chrom] = contigs
    if not args.legacy_junctions and not args.no_extension:
        et_fname = "extension.tsv"
        if args.prefix != "":
            et_fname = args.prefix + '.' + et_fname
        ref_lengths = dict(zip(sorted_primary_alignments.references, sorted_primary_alignments.lengths))
        candidates = collect_extension_candidates(args.query_bam, unplaced)
        placement, inserted = extend_placement(args, placement, candidates, ref_lengths, et_fname)
        propagate_extents(extents, inserted, placement)
        if inserted and args.store_final_bam is not None:
            # Include inserted contigs in the stored final BAM.
            sorted_primary_alignments.close()
            os.remove(sorted_primary_alignments_bam + '.bai')
            sorted_bam_from_aln_list([c for cs in placement.values() for c in cs], query_bam.header,
                                     sorted_bam_name=args.store_final_bam)
            pysam.index(sorted_primary_alignments_bam)
            sorted_primary_alignments = pysam.AlignmentFile(sorted_primary_alignments_bam, "rb")

    junctions_tsv = None
    if not args.legacy_junctions:
        jt_fname = "junctions.tsv"
        if args.prefix != "":
            jt_fname = args.prefix + '.' + jt_fname
        junctions_tsv = open(jt_fname, "w")  # Per-junction decisions
        junctions_tsv.write("# GPatch %s %s\n" % (__version__, " ".join(sys.argv[1:])))
        junctions_tsv.write("chrom\tprev_contig\tcontig\tpos\tcontig_ext_start\tprev_resid3\tcontig_resid5"
                            "\tprev_chain_n\tcontig_chain_n\test_dist\twindow_bp\toverlap_candidates"
                            "\tdecision\ttrim_bp\tpatch_ref_start\tpatch_ref_end\tpatch_bp\treason\n")
        overlaps = prepare_junctions(args, placement, extents)

    for ref_seq in SeqIO.parse(args.reference_fasta, "fasta"):
        # Select all contigs placed on this sequence.
        contigs = placement.get(ref_seq.id, [])

        # If no contigs map to this sequence, either print it unchanged
        # or omit it, depending on args.
        if len(contigs) == 0:
            if args.drop_missing:
                sys.stderr.write("No contigs map to %s. Omitting from output.\n" % (ref_seq.id))
            else:
                sys.stderr.write("No contigs map to %s. Printing the reference sequence to output unchanged.\n" % (ref_seq.id))
                patched_fasta.write("%s\n" % (ref_seq.format("fasta")))
            continue

        if args.legacy_junctions:
            patched_seq = legacy_patch_sequence(ref_seq, contigs, contig_breakpoints, args, patches_bed, contigs_bed,
                                                unplaced=unplaced)
        else:
            patched_seq = patch_sequence(ref_seq, contigs, extents, overlaps, args,
                                         patches_bed, contigs_bed, junctions_tsv)

        # Swap in the patched sequence for the reference sequence and print the result
        ref_seq.seq = patched_seq
        patched_fasta.write("%s\n" % (ref_seq.format("fasta")))

    # Close open file handles, etc.
    query_bam.close()
    sorted_primary_alignments.close()
    patched_fasta.close()
    patches_bed.close()
    contigs_bed.close()
    if junctions_tsv is not None:
        junctions_tsv.close()

    # Report contigs left out of the pseudoassembly, with their sequences.
    with open(cb_fname) as f:
        placed = set(line.split("\t")[3] for line in f)
    write_unplaced(args, placed, unplaced)

    # Clean up any leftover files as needed.
    if args.store_final_bam is None:
        os.remove(sorted_primary_alignments_bam)
        os.remove(sorted_primary_alignments_bam + '.bai')

    sys.stderr.write("Done!\n")


def compute_extents(args, primary_alignments):
    """
    Chained reference extents (see chained_extent) for each primary
    alignment, keyed on contig name.
    """
    spans = {}
    for aln in primary_alignments:
        if aln.query_sequence is None or len(aln.query_sequence) != aln.infer_read_length():
            sys.stderr.write("ERROR: the primary alignment of %s does not contain the full contig sequence (hard-clipped?). Re-run the alignment with soft-clipping for primary alignments, or use --legacy_junctions.\n" % aln.query_name)
            sys.exit(1)
        spans[aln.query_name] = AlignedSpan(aln)
    supps = collect_supplementary_spans(args.query_bam, set(spans), args.chain_mapq, args.chain_minlen)
    return dict((name, chained_extent(span, supps.get(name, []), args.chain_maxgap, args.chain_tol))
                for name, span in spans.items())


def suffix_prefix_overlaps(tail, head, len_a, len_b, k, min_id, min_ovl):
    """
    Lengths (>= min_ovl) of exact-anchored overlaps between the end of
    contig a (tail: its last bases, as placed) and the start of contig b
    (head: its first bases), verified at >= min_id identity.
    """
    lp, lc = len(tail), len(head)
    if lp < k or lc < k:
        return []
    ov1 = set(lp - i for i in find_all(tail, head[:k]))
    if not ov1:
        return []
    ov2 = set(j + k for j in find_all(head, tail[-k:]))
    return [o for o in sorted(ov1 & ov2)
            if min_ovl <= o < min(len_a, len_b) and identity_at_least(tail[lp - o:], head[:o], min_id)]


def collect_extension_candidates(bam_path, unplaced):
    """
    Primary alignments of unplaced contigs eligible for overlap extension
    (left out for low mapping quality or nesting), keyed on reference name.
    """
    names = set(n for n, (reason, _) in unplaced.items() if reason in ("low_mapq", "nested"))
    cands = {}
    with pysam.AlignmentFile(bam_path) as bam:
        for aln in bam:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.query_name not in names:
                continue
            if aln.query_sequence is None or len(aln.query_sequence) != aln.infer_read_length():
                continue
            cands.setdefault(aln.reference_name, []).append(aln)
    return cands


def extend_placement(args, placement, candidates, ref_lengths, log_fname):
    """
    Insert unplaced contigs into the gaps between placed contigs (and at
    chromosome ends) by walking through unambiguous sequence overlaps.

    For each gap between adjacent placed contigs A and B, a forward walk
    starts at A and a backward walk at B (aiming for the end of the forward
    walk). At each step, the candidates are unplaced contigs whose primary
    alignment is on the same reference sequence, within the gap region
    (padded by --ext_region_pad) and used in the orientation of that
    alignment. A candidate is added if it is the only one whose 5' end
    (forward) or 3' end (backward) overlaps the current contig end by
    >= --ext_min_overlap. A walk stops when no candidate or more than one
    qualifies, when it reaches the opposite flanking contig, or after
    --ext_max_steps contigs. Contigs claimed by more than one walk are
    rejected and walks are truncated before them.

    Returns (new placement, inserted), where inserted maps each inserted
    contig to (direction, anchor contig name, overlap).
    """
    sys.stderr.write("Extending placements through contig overlaps...\n")
    W = args.ext_window
    seq_ends = {}

    def ends(aln):
        if aln.query_name not in seq_ends:
            s = aln.query_sequence.upper()
            seq_ends[aln.query_name] = (s[:W], s[-W:], len(s))
        return seq_ends[aln.query_name]

    def ovl(a, b):
        _, ta, la = ends(a)
        hb, _, lb = ends(b)
        return suffix_prefix_overlaps(ta, hb, la, lb, args.anchor, args.min_overlap_identity, args.ext_min_overlap)

    log = []
    walks = {}
    steps = {}   # (chrom, gap, direction, name) -> (anchor name, overlap)
    for chrom, placed in placement.items():
        L = ref_lengths[chrom]
        for gi in range(len(placed) + 1):
            A = placed[gi - 1] if gi > 0 else None
            B = placed[gi] if gi < len(placed) else None
            lo = (A.reference_start - args.ext_region_pad) if A is not None else 0
            hi = (B.reference_end + args.ext_region_pad) if B is not None else L
            pool = [c for c in candidates.get(chrom, []) if c.reference_start >= lo and c.reference_end <= hi]
            res = {"fwd": [], "bwd": []}
            fwd_closed = False
            for direction in ("fwd", "bwd"):
                if direction == "fwd":
                    start, target = A, B
                else:
                    start, target = B, (res["fwd"][-1] if res["fwd"] else A)
                if start is None or (direction == "bwd" and fwd_closed):
                    continue
                fwd_names = set(c.query_name for c in res["fwd"])
                cur, used, why = start, set(), "max steps"
                for step in range(args.ext_max_steps):
                    if target is not None:
                        closed = ovl(cur, target) if direction == "fwd" else ovl(target, cur)
                        if closed:
                            why = "reached %s" % target.query_name
                            fwd_closed = fwd_closed or direction == "fwd"
                            break
                    hits = []
                    for x in pool:
                        if x.query_name in used or x.query_name in fwd_names:
                            continue
                        o = ovl(cur, x) if direction == "fwd" else ovl(x, cur)
                        if o:
                            hits.append((x, o))
                    if len(hits) != 1:
                        why = "no candidate" if not hits else "ambiguous: %s" % ",".join(h[0].query_name for h in hits)
                        break
                    x, o = hits[0]
                    res[direction].append(x)
                    used.add(x.query_name)
                    steps[(chrom, gi, direction, x.query_name)] = (cur.query_name, min(o))
                    log.append([chrom, gi, A.query_name if A is not None else "START",
                                B.query_name if B is not None else "END", direction, step, cur.query_name,
                                x.query_name, "-" if x.is_reverse else "+", ",".join(str(v) for v in o)])
                    cur = x
                log.append([chrom, gi, A.query_name if A is not None else "START",
                            B.query_name if B is not None else "END", direction, "stop", cur.query_name,
                            ".", ".", why])
            walks[(chrom, gi)] = res

    claims = {}
    for r in walks.values():
        for d in ("fwd", "bwd"):
            for x in r[d]:
                claims[x.query_name] = claims.get(x.query_name, 0) + 1
    conflicted = set(n for n, c in claims.items() if c > 1)
    for (chrom, gi), r in walks.items():
        for d in ("fwd", "bwd"):
            for i, x in enumerate(r[d]):
                if x.query_name in conflicted:
                    for y in r[d][i:]:
                        log.append([chrom, gi, ".", ".", d, "rejected", ".", y.query_name, ".",
                                    "claimed by several walks" if y.query_name in conflicted
                                    else "follows a rejected contig"])
                    r[d] = r[d][:i]
                    break

    new_placement, inserted = {}, {}
    for chrom, placed in placement.items():
        order = []
        for gi in range(len(placed) + 1):
            r = walks[(chrom, gi)]
            for d, lst in (("fwd", r["fwd"]), ("bwd", list(reversed(r["bwd"])))):
                for x in lst:
                    order.append(x)
                    anchor, o = steps[(chrom, gi, d, x.query_name)]
                    inserted[x.query_name] = (d, anchor, o)
            if gi < len(placed):
                order.append(placed[gi])
        new_placement[chrom] = order

    with open(log_fname, "w") as f:
        f.write("# GPatch %s %s\n" % (__version__, " ".join(sys.argv[1:])))
        f.write("chrom\tgap\tleft\tright\tdirection\tstep\tfrom\tadded\tstrand\toverlaps_or_reason\n")
        for row in log:
            f.write("\t".join(str(v) for v in row) + "\n")
    sys.stderr.write("Overlap extension inserted %d contigs (%d rejected as claimed by several walks).\n"
                     % (len(inserted), len(conflicted)))
    return new_placement, inserted


def propagate_extents(extents, inserted, placement):
    """
    Reference extents of contigs inserted by overlap extension, carried from
    the contig each was reached from through the measured overlap (their own
    alignments have low mapping quality and are not used). Forward: start =
    anchor end - overlap; backward: end = anchor start + overlap. The
    residual unaligned length of the anchor end is inherited at both ends.
    """
    lengths = dict((c.query_name, len(c.query_sequence)) for cs in placement.values() for c in cs)
    pending = dict(inserted)
    while pending:
        progress = False
        for name, (d, anchor, o) in list(pending.items()):
            if anchor in pending:
                continue
            A = extents[anchor]
            if d == "fwd":
                s0 = A["end"] - o
                extents[name] = {"start": s0, "end": s0 + lengths[name], "resid5": A["resid3"],
                                 "resid3": A["resid3"], "nchain": 0}
            else:
                e0 = A["start"] + o
                extents[name] = {"start": e0 - lengths[name], "end": e0, "resid5": A["resid5"],
                                 "resid3": A["resid5"], "nchain": 0}
            del pending[name]
            progress = True
        if not progress:
            raise RuntimeError("cyclic overlap-extension anchors: %s" % ",".join(pending))


def write_unplaced(args, placed, unplaced):
    """
    Write contigs present in the input BAM but absent from the
    pseudoassembly to unplaced.fasta (original contig orientation) and
    unplaced.tsv (reason and primary alignment).
    """
    uf_fname, ut_fname = "unplaced.fasta", "unplaced.tsv"
    if args.prefix != "":
        uf_fname = args.prefix + '.' + uf_fname
        ut_fname = args.prefix + '.' + ut_fname
    n = bp = 0
    with pysam.AlignmentFile(args.query_bam) as bam, open(uf_fname, "w") as uf, open(ut_fname, "w") as ut:
        ut.write("contig\tlength\treason\tdetail\tmapq\tref\tref_start\tref_end\tstrand\n")
        for aln in bam:
            if aln.is_secondary or aln.is_supplementary or aln.query_name in placed:
                continue
            seq = aln.query_sequence or ""
            if aln.is_reverse:
                seq = str(Seq(seq).reverse_complement())
            reason, detail = unplaced.get(aln.query_name, ("not_placed", "reference sequence not in FASTA or contig skipped"))
            if aln.is_unmapped:
                ut.write("%s\t%d\t%s\t%s\t.\t.\t.\t.\t.\n" % (aln.query_name, len(seq), reason, detail))
            else:
                ut.write("%s\t%d\t%s\t%s\t%d\t%s\t%d\t%d\t%s\n" % (aln.query_name, len(seq), reason, detail,
                         aln.mapping_quality, aln.reference_name, aln.reference_start, aln.reference_end,
                         "-" if aln.is_reverse else "+"))
            uf.write(">%s\n" % aln.query_name)
            for i in range(0, len(seq), 80):
                uf.write(seq[i:i + 80] + "\n")
            n += 1
            bp += len(seq)
    sys.stderr.write("%d contigs (%d bp) were not placed; see %s and %s.\n" % (n, bp, ut_fname, uf_fname))


def prepare_junctions(args, placement, extents):
    """
    Measure overlaps at every junction between adjacent placed contigs,
    using the contig extents to size the search windows. placement maps each
    reference sequence to its placed contig alignments, in order. Returns the
    overlaps keyed on (prev_name, name).
    """
    sys.stderr.write("Measuring junction overlaps...\n")
    order = placement

    minimap2 = args.minimap2 if args.minimap2 is not None else shutil.which("minimap2")
    if args.indel_tol > 0 and not minimap2:
        sys.stderr.write("WARNING: minimap2 not found; candidate overlaps containing indels will not be considered.\n")
    params = {"anchor": args.anchor, "min_id": args.min_overlap_identity, "indel_tol": args.indel_tol,
              "max_mm2_tests": args.max_mm2_tests, "win_max": args.win_max, "minimap2": minimap2}
    keys, jobs = [], []
    for chrom, contigs in order.items():
        for prev, cur in zip(contigs, contigs[1:]):
            P, C = extents[prev.query_name], extents[cur.query_name]
            est_overlap = P["end"] - C["start"]
            L0 = min(args.win_max, max(args.win_min, 2 * max(est_overlap, 0) + P["resid3"] + C["resid5"]))
            pseq, cseq = prev.query_sequence.upper(), cur.query_sequence.upper()
            keys.append((prev.query_name, cur.query_name))
            jobs.append((pseq[-args.win_max:], cseq[:args.win_max], len(pseq), len(cseq), L0, params))
    if args.threads > 1:
        with ProcessPoolExecutor(args.threads) as ex:
            results = list(ex.map(_measure_overlaps_job, jobs))
    else:
        results = [_measure_overlaps_job(j) for j in jobs]
    return dict(zip(keys, results))


def patch_sequence(ref_seq, contigs, extents, overlaps, args, patches_bed, contigs_bed, junctions_tsv):
    """
    Build the patched sequence for one reference chromosome, deciding at each
    junction whether to trim, patch or bookend (see the junction-handling
    notes above). Returns the patched sequence as a Seq.
    """
    chrom = ref_seq.id
    L = len(ref_seq.seq)
    pieces = []
    plen = 0
    pos = None    # reference position of the end of the sequence placed so far
    prev = None

    def row(*v):
        junctions_tsv.write("\t".join(str(x) for x in v) + "\n")

    def patch_seq(a, b):
        if args.scaffold_only:
            return "n" * (args.gap_length if args.gap_length >= 0 else b - a)
        return str(ref_seq.seq[a:b]).lower()

    for contig in contigs:
        name = contig.query_name
        C = extents[name]
        trim, patch = 0, None
        if prev is None:
            s = C["start"]
            if args.no_extend:
                dec, why = "bookend", "chromosome start; --no_extend"
            elif s > 0 and C["resid5"] <= args.max_resid:
                dec, why, patch = "patch", "chromosome start; extent trusted", (0, min(s, L))
            elif s > 0:
                dec, why = "bookend", "chromosome start; residual 5' unaligned %d > --max_resid" % C["resid5"]
            else:
                dec, why = "bookend", "chromosome start; extrapolated start %d <= 0 (no trim at chromosome ends)" % s
            row(chrom, "START", name, 0, s, ".", C["resid5"], ".", C["nchain"], s, ".", ".",
                dec, 0, *(patch if patch else (".", ".")), (patch[1] - patch[0]) if patch else 0, why)
        else:
            P = extents[prev]
            trims, win, note = overlaps[(prev, name)]
            est = C["start"] - pos
            dec, trim, why = decide_junction(P, C, pos, trims, args.max_resid)
            if dec.startswith("trim") and args.no_trim:
                dec, why, trim = "bookend", "overlap of %d not trimmed (--no_trim); %s" % (trim, why), 0
            if dec.startswith("patch"):
                patch = (min(pos, L), min(C["start"], L))
            row(chrom, prev, name, pos, C["start"], P["resid3"], C["resid5"], P["nchain"], C["nchain"],
                est, win, ",".join(str(t) for t in trims) or ".", dec, trim,
                *(patch if patch else (".", ".")), (patch[1] - patch[0]) if patch else 0,
                why + "; " + note)
        if patch and patch[1] > patch[0]:
            # Patches are printed in lower case for easy identification.
            ps = patch_seq(patch[0], patch[1])
            pieces.append(ps)
            plen += len(ps)
            patches_bed.write("%s\t%d\t%d\n" % (chrom, patch[0], patch[1]))
        # Contig sequence (already in placement orientation in the BAM), in upper case.
        seq = contig.query_sequence
        start = plen
        pieces.append(seq[trim:].upper())
        plen += len(seq) - trim
        qstrand = "-" if contig.is_reverse else "+"
        contigs_bed.write("%s\t%d\t%d\t%s\t.\t%s\t%d\t%d\n" % (chrom, start, plen, name, qstrand, trim, len(seq)))
        pos = C["end"] if pos is None else max(pos, C["end"])
        prev = name

    # Terminal patch from the end of the last contig to the chromosome end.
    P = extents[prev]
    patch = None
    if args.no_extend:
        dec, why = "bookend", "chromosome end; --no_extend"
    elif pos < L and P["resid3"] <= args.max_resid:
        dec, why, patch = "patch", "chromosome end; extent trusted", (max(pos, 0), L)
    elif pos < L:
        dec, why = "bookend", "chromosome end; residual 3' unaligned %d > --max_resid" % P["resid3"]
    else:
        dec, why = "bookend", "chromosome end; extrapolated end %d >= chromosome length" % pos
    row(chrom, prev, "END", pos, L, P["resid3"], ".", P["nchain"], ".", L - pos, ".", ".",
        dec, 0, *(patch if patch else (".", ".")), (patch[1] - patch[0]) if patch else 0, why)
    if patch and patch[1] > patch[0]:
        pieces.append(patch_seq(patch[0], patch[1]))
        patches_bed.write("%s\t%d\t%d\n" % (chrom, patch[0], patch[1]))
    return Seq("".join(pieces))


def legacy_patch_sequence(ref_seq, contigs, contig_breakpoints, args, patches_bed, contigs_bed, unplaced=None):
    """
    Original (GPatch <= 0.4) junction handling (--legacy_junctions): contig extents from the
    primary alignment extended by its soft-clips; overlapping contigs trimmed
    by the overlap of those extents (unless --no_trim); gaps patched by their
    distance. Returns the patched sequence.
    """
    pos = 0 # Tracks position on the reference chromosome
    patched_seq = ""
    for contig in contigs:
        rstart = contig_breakpoints[contig.query_name][1]
        """
        # This introduced a one-off error! We don't need to do this
        # because all coordinates are already in proper zero-based form.
        if rstart > 0:
            rstart -= 1
        """

        # Check for overlapping/bookended contig mappings.
        qstart = 0  # Start of contig interval to append
        #sys.stderr.write("%s\t%s\t%s\n" % (rstart, pos, pos-rstart))
        if rstart > pos:
            # No overlap. Append the patch, in all lower-case
            # for easy identification.
            patch = ""
            if args.scaffold_only:
                # Use a string of N characters the length of the patch instead
                # of an actual patch.
                if args.gap_length >= 0:
                    patch = "n" * args.gap_length
                else:
                    patch = "n" * (rstart-pos)
            else:
                patch = ref_seq.seq[pos:rstart].lower()
            if pos == 0 and args.no_extend:
                pass
            else:
                patched_seq = patched_seq + patch
                # Write patch coordinates in reference frame to patches_bed
                patches_bed.write("%s\t%d\t%d\n" % (ref_seq.id, pos, rstart))
        else:
            # Handle overlapping contig ends by trimming the 5' end of
            # this contig sequence by the length of the overlap. Note this
            # evaluates to zero if contigs are bookended (i.e., pos and
            # rstart are equal.)
            # Trimming by alignment-coordinate overlap is unreliable when
            # contig ends are soft-clipped; the default (non-legacy) mode
            # measures overlaps from the contig sequences instead.
            if not args.no_trim:
                qstart = pos - rstart
                if qstart > len(Seq(contig.query_sequence)):
                    # This should not happen, as it is indicative of a nested
                    # mapping, which should be filtered out prior to sequence-
                    # building. However, I have seen instances of nested mappings
                    # that somehow slip through the cracks. This is a hack to fix
                    # those instances until I can figure out why it happens.
                    if unplaced is not None:
                        unplaced[contig.query_name] = ("nested", "trim would exceed contig length")
                    continue
            
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
        #sys.stderr.write("%s\t%s\t%s\n" % (pos, contig_breakpoints[contig.query_name][2], len(ref_seq.seq)))
        if contig_breakpoints[contig.query_name][2] > pos:
            # This should always evaluate true in the absence of nested contigs.
            pos = contig_breakpoints[contig.query_name][2]
            
    # Once the above loop finishes, we need to add the terminal segment
    # from the reference genome.
    if pos < len(ref_seq.seq) and not args.no_extend:
        #sys.stderr.write("%s\n" % (ref_seq.seq[pos:len(ref_seq.seq)]))
        patched_seq = patched_seq + ref_seq.seq[pos:len(ref_seq.seq)]
        patches_bed.write("%s\t%d\t%d\n" % (ref_seq.id, pos, len(ref_seq.seq)))

    return patched_seq


if __name__ == '__main__':
    main()
    exit(0)
