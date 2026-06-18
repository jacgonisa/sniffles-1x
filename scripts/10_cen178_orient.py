#!/usr/bin/env python3
"""Step 10 — CEN178 array orientation track for the genome maps.
Maps the 178-bp CEN178 consensus to each centromere reference (both strands, all
repeat hits) with minimap2, then majority-votes strand per 1-kb bin:
  + strand hits -> FORWARD array (red),  - strand hits -> REVERSE array (blue).
-> results/cen178_orientation.tsv  (hap chrom start end strand)
Run with nextflow_env python (pysam) + minimap2."""
import os, subprocess, tempfile, sys
from collections import defaultdict
import pysam
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts")
from common import HAPS, CEN, REF, OUT

MM = "/home/jg2070/minimap2/minimap2"
CEN178 = ("AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAA"
          "GACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAG"
          "AAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG")
BIN = 1000


def orient_chrom(refseq, subseq, td):
    """Return list of minimap2 (start,end,strand) hits of CEN178 in subseq."""
    sf = os.path.join(td, "sub.fa"); qf = os.path.join(td, "q.fa")
    open(sf, "w").write(f">sub\n{subseq}\n")
    open(qf, "w").write(f">CEN178\n{CEN178}\n")
    out = subprocess.run([MM, "-k11", "-w3", "-N", "1000000", "-p", "0",
                          "--secondary=yes", sf, qf],
                         capture_output=True, text=True).stdout
    hits = []
    for ln in out.splitlines():
        c = ln.split("\t")
        if len(c) < 9:
            continue
        hits.append((int(c[7]), int(c[8]), c[4]))   # tstart, tend, strand
    return hits


def main():
    fa = {h: pysam.FastaFile(REF[h][0]) for h in HAPS}
    rows = []
    with tempfile.TemporaryDirectory() as td:
        for hap in HAPS:
            for chrom, (a, b) in CEN[hap].items():
                sub = fa[hap].fetch(chrom, a, b)
                hits = orient_chrom(None, sub, td)
                # majority strand per BIN
                vote = defaultdict(lambda: [0, 0])  # bin -> [plus, minus]
                for s, e, strand in hits:
                    bi = ((s + e) // 2) // BIN
                    vote[bi][0 if strand == "+" else 1] += 1
                # collapse consecutive same-strand bins into blocks
                bins = sorted(vote)
                blocks = []
                for bi in bins:
                    pl, mi = vote[bi]
                    st = "+" if pl >= mi else "-"
                    s0 = a + bi * BIN; e0 = a + (bi + 1) * BIN
                    if blocks and blocks[-1][2] == st and s0 - blocks[-1][1] <= BIN:
                        blocks[-1][1] = e0
                    else:
                        blocks.append([s0, e0, st])
                npl = sum(1 for bi in bins if vote[bi][0] >= vote[bi][1])
                print(f"{hap} {chrom}: {len(hits)} hits, {len(bins)} bins ({npl}+ / {len(bins)-npl}-), {len(blocks)} blocks")
                for s0, e0, st in blocks:
                    rows.append((hap, chrom, s0, e0, st))
    with open(f"{OUT}/cen178_orientation.tsv", "w") as f:
        f.write("hap\tchrom\tstart\tend\tstrand\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print("DONE_ORIENT")


if __name__ == "__main__":
    main()
