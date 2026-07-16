#!/usr/bin/env python3
"""Step 18 — split-and-map on the chromosome ARMS (unique-sequence baseline for §12).

Extends the CEN-vs-ARM control to the split-and-map route: runs the same
contrast-frontier split + winnowmap re-map + classify_splits on candidate reads anchored
in the distal arm windows (5 Mb past the centromere → transition excluded). Gives the arm
(no-satellite) rate for DUP/INV/large-DEL/BND so their centromere enrichment gets a proper
unique-sequence baseline. Compares to the CEN split-and-map rate.
-> results/arm_splitmap_control.tsv  (+ printed).
Run with nextflow_env python (pysam + winnowmap)."""
import importlib, subprocess, tempfile, os, csv
from collections import Counter, defaultdict
from types import SimpleNamespace
import pysam
sm = importlib.import_module("03_split_and_map")
from sniffles import sv
from sniffles.leadprov import Lead
from common import SAMPLES, HAPS, CEN, REF, CHRLEN, GROUPS, bam_path, OUT, refkey

ARM_BUFFER = 5_000_000; ARM_WIN = 3_000_000
DE_MIN = 0.005; NM_MIN = 50
TYPES = ["DEL", "INS", "DUP", "INV", "BND"]
_DUMMY = sm._Dummy()


def arm_windows(rk):
    w = {}
    for chrom, (a, b) in CEN[rk].items():
        s = b + ARM_BUFFER; e = min(s + ARM_WIN, CHRLEN[rk][chrom] - 100_000)
        if e - s > 200_000:
            w[chrom] = (s, e)
    return w


def run_one(sample, hap, td):
    rk = refkey(sample, hap)
    wins = arm_windows(rk)
    bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
    fa = os.path.join(td, f"{sample}_{hap}.fa"); meta = {}; idx = 0; arm_reads = 0
    with open(fa, "w") as out:
        for chrom, (a, b) in wins.items():
            for r in bam.fetch(chrom, a, b):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or not (a <= r.reference_start < b):
                    continue
                arm_reads += 1
                de = r.get_tag("de") if r.has_tag("de") else 0.0
                nm = r.get_tag("NM") if r.has_tag("NM") else 0
                if not (de >= DE_MIN or nm >= NM_MIN or r.has_tag("SA")):
                    continue
                seq = r.query_sequence
                if seq is None or len(seq) < 2 * sm.MIN_FRAG + 200:
                    continue
                bs = sm.best_split(r)
                if bs is None:
                    continue
                qsplit, _c = bs
                if qsplit < sm.MIN_FRAG or len(seq) - qsplit < sm.MIN_FRAG:
                    continue
                meta[idx] = (r.query_name, chrom, qsplit, (a, b))
                out.write(f">{idx}_A\n{seq[:qsplit]}\n>{idx}_B\n{seq[qsplit:]}\n"); idx += 1
    bam.close()
    if idx == 0:
        return arm_reads, Counter()
    ref, rep = REF[rk]; frag_bam = fa.replace(".fa", ".bam")
    subprocess.run(f"{sm.WIN} -W {rep} -ax map-pb -t {sm.THREADS} {ref} {fa} 2>/dev/null | "
                   f"{sm.ST} sort -@4 -o {frag_bam} - && {sm.ST} index {frag_bam}",
                   shell=True, check=True, executable="/bin/bash")
    fb = pysam.AlignmentFile(frag_bam, "rb"); halves = {}
    for r in fb.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        if r.query_alignment_length < sm.MIN_FRAG or r.mapping_quality < sm.MAPQ_MIN:
            continue
        i, half = r.query_name.rsplit("_", 1); halves.setdefault(int(i), {})[half] = r
    fb.close()
    cnt = Counter()
    for i, hv in halves.items():
        if "A" not in hv or "B" not in hv:
            continue
        qname, main_contig, qsplit, (a, b) = meta[i]; A, B = hv["A"], hv["B"]
        leads = [Lead(0, qname, A.reference_name, A.reference_start, A.reference_end,
                      A.query_alignment_start, A.query_alignment_end, "-" if A.is_reverse else "+", A.mapping_quality, 0, "x", "?"),
                 Lead(0, qname, B.reference_name, B.reference_start, B.reference_end,
                      qsplit + B.query_alignment_start, qsplit + B.query_alignment_end, "-" if B.is_reverse else "+", B.mapping_quality, 0, "x", "?")]
        for ld in sv.classify_splits(_DUMMY, leads, sm.CFG, main_contig):
            for svtype, svstart, arg in (ld.svtypes_starts_lens or []):
                if svtype != "NOSV" and a <= svstart < b:
                    cnt[svtype] += 1
    return arm_reads, cnt


def main():
    arm_cnt = defaultdict(Counter); arm_reads = defaultdict(int)
    with tempfile.TemporaryDirectory() as td:
        for sample, tis in SAMPLES:
            for hap in HAPS:
                nr, cnt = run_one(sample, hap, td)
                arm_reads[sample] += nr
                for k, v in cnt.items():
                    arm_cnt[sample][k] += v
                print(f"{sample} {hap}: arm reads {nr}, split-and-map {dict(cnt)}")

    # CEN split-and-map rate from split_and_map.tsv + cen_read_counts
    cen_cnt = defaultdict(Counter); cen_reads = defaultdict(int)
    for r in csv.DictReader(open(f"{OUT}/split_and_map.tsv"), delimiter="\t"):
        cen_cnt[r["sample"]][r["svtype"]] += 1
    with open(f"{OUT}/cen_read_counts.tsv") as f:
        next(f)
        for ln in f:
            s, h, n = ln.split(); cen_reads[s] += int(n)

    cols = ["group", "svtype", "cen_per_Mreads", "arm_per_Mreads", "enrichment_CEN_over_ARM"]
    with open(f"{OUT}/arm_splitmap_control.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        print(f"{'group':16}{'type':5}{'CEN/Mr':>9}{'ARM/Mr':>9}{'CEN÷ARM':>9}  (split-and-map route)")
        for samp in GROUPS:
            cm = cen_reads[samp] / 1e6; am = arm_reads[samp] / 1e6
            for t in TYPES:
                cr = cen_cnt[samp][t] / cm if cm else 0; ar = arm_cnt[samp][t] / am if am else 0
                enr = round(cr / ar, 1) if ar else float("inf")
                f.write(f"{samp}\t{t}\t{cr:.1f}\t{ar:.1f}\t{enr}\n")
                print(f"{samp:16}{t:5}{cr:9.1f}{ar:9.1f}{str(enr):>9}")
    print(f"arm reads: {dict(arm_reads)}")
    print("DONE_ARM_SPLITMAP")


if __name__ == "__main__":
    main()
