#!/usr/bin/env python3
"""Step 2 — Sniffles2's own SV-signal algorithm, run at the single-molecule level.

Sniffles needs >=2 reads only in its *clustering* step. Its per-read signal
extraction (leadprov) and topology classifier (sv.classify_splits) are inherently
single-molecule. We import those directly and emit every per-read "lead" with NO
clustering / min-support -> one row = one single-molecule SV.

  - INLINE leads: CIGAR I/D operations >= 50 bp           (read_iterindels loop)
  - SPLIT  leads: primary + SA supplementaries -> classify_splits -> DEL/DUP/INV/INS/BND

-> results/leadprov_sm.tsv
Run with nextflow_env python (has pysam + sniffles).  Self-test: --selftest
"""
import sys, os, pysam
from types import SimpleNamespace
sys.path.insert(0, "/home/jg2070/miniforge3/envs/nextflow_env/lib/python3.13/site-packages")
from sniffles import sv
from sniffles.leadprov import Lead, CIGAR_analyze
from common import SAMPLES, HAPS, CEN, bam_path, OUT

# minimal config mirroring Sniffles defaults; only attrs classify_splits reads.
# minsvlen_screen=50 -> the colleague's ">=50 bp" threshold (vs sniffles' 45).
CFG = SimpleNamespace(minsvlen_screen=50, long_ins_length=2500,
                      bnd_min_split_length=1000, dev_seq_cache_maxlen=0)
MINSV = 50
MAPQ_MIN = 10


def cigar_leads(read):
    """read_iterindels: yield (svtype, ref_pos, svlen) for CIGAR I/D >= MINSV."""
    pos_ref = read.reference_start
    for op, oplen in read.cigartuples:
        if op == pysam.CINS and oplen >= MINSV:
            yield ("INS", pos_ref, oplen)
        elif op == pysam.CDEL and oplen >= MINSV:
            yield ("DEL", pos_ref + oplen, -oplen)
        # advance ref for M/D/N/=/X
        if op in (pysam.CMATCH, pysam.CDEL, pysam.CREF_SKIP, pysam.CEQUAL, pysam.CDIFF):
            pos_ref += oplen


def split_leads(read, contig):
    """read_itersplits: build leads from primary + SA, run sv.classify_splits."""
    if not read.has_tag("SA"):
        return
    supps = [p.split(",") for p in read.get_tag("SA").split(";") if p]
    qry_start = (read.query_length - read.query_alignment_end) if read.is_reverse else read.query_alignment_start
    leads = [Lead(0, read.query_name, contig, read.reference_start,
                  read.reference_start + read.reference_length, qry_start,
                  qry_start + read.query_alignment_length,
                  "-" if read.is_reverse else "+", read.mapping_quality, 0, "SPLIT_PRIM", "?")]
    for refname, pos, strand, cigar, mapq, nm in supps:
        try:
            rs_fwd, rs_rev, refspan, readspan = CIGAR_analyze(cigar)
        except Exception:
            return
        pos0 = int(pos) - 1
        qs = rs_rev if strand == "-" else rs_fwd
        leads.append(Lead(0, read.query_name, refname, pos0, pos0 + refspan,
                          qs, qs + readspan, strand, int(mapq), 0, "SPLIT_SUP", "?"))
    leads = sv.classify_splits(read, leads, CFG, contig)
    for ld in leads:
        for svtype, svstart, arg in (ld.svtypes_starts_lens or []):
            if svtype == "NOSV":
                continue
            svlen = None if svtype == "BND" else arg
            yield (svtype, svstart, svlen, ld.mapq)


def run():
    out = open(f"{OUT}/leadprov_sm.tsv", "w")
    out.write("sample\thap\ttissue\tchrom\tpos\tsvtype\tsvlen\tsource\tmapq\tread\n")
    for sample, tis in SAMPLES:
        for hap in HAPS:
            bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
            n = 0
            for chrom, (a, b) in CEN[hap].items():
                for r in bam.fetch(chrom, a, b):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary:
                        continue
                    if r.mapping_quality < MAPQ_MIN or not (a <= r.reference_start < b):
                        continue
                    for svtype, pos, svlen in cigar_leads(r):
                        if a <= pos < b:
                            out.write(f"{sample}\t{hap}\t{tis}\t{chrom}\t{pos}\t{svtype}\t{svlen}\tINLINE\t{r.mapping_quality}\t{r.query_name}\n"); n += 1
                    for svtype, pos, svlen, mapq in split_leads(r, chrom):
                        if mapq < MAPQ_MIN:
                            continue
                        if a <= pos < b:
                            sl = "" if svlen is None else svlen
                            out.write(f"{sample}\t{hap}\t{tis}\t{chrom}\t{pos}\t{svtype}\t{sl}\tSPLIT\t{mapq}\t{r.query_name}\n"); n += 1
            print(f"{sample} {hap}: {n} single-molecule leads")
            bam.close()
    out.close(); print("DONE_LEADPROV")


def selftest():
    """Assert sv.classify_splits yields DEL, DUP, INV on hand-built split topologies."""
    class R:  # minimal read stand-in for classify_splits (only query_sequence used, and only if seq cached)
        query_sequence = "A" * 100000
    def lead(contig, rs, re_, qs, qe, strand):
        return Lead(0, "r", contig, rs, re_, qs, qe, strand, 60, 0, "SPLIT", "?")
    def types(leads):
        got = sv.classify_splits(R(), leads, CFG, "Chr1")
        return [t for ld in got for (t, _s, _a) in (ld.svtypes_starts_lens or [])]
    # DEL: two fwd fragments, ref gap >> query gap
    assert "DEL" in types([lead("Chr1", 1000, 2000, 0, 1000, "+"),
                            lead("Chr1", 5000, 6000, 1000, 2000, "+")]), "DEL fail"
    # DUP: second fragment starts before first ends on ref (fwd), overlap
    assert "DUP" in types([lead("Chr1", 5000, 6000, 0, 1000, "+"),
                            lead("Chr1", 4000, 5000, 1000, 2000, "+")]), "DUP fail"
    # INV: opposite strands
    assert "INV" in types([lead("Chr1", 1000, 2000, 0, 1000, "+"),
                           lead("Chr1", 3000, 4000, 1000, 2000, "-")]), "INV fail"
    print("SELFTEST_OK")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        run()
