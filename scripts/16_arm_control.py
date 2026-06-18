#!/usr/bin/env python3
"""Step 16 — chromosome-ARM control (unique sequence, no satellite).

Runs the same per-read leadprov detection (CIGAR + SA split, via 02) on reads anchored in
the chromosome ARMS, then compares the per-million-read SV rate to the centromere. Arms are
unique sequence, so this is the no-satellite background: real arm variants + sequencing/mapping
artifacts, but ~no CEN178 cross-mapping. The gap CEN − ARM = the centromere-specific signal.

IMPORTANT: arm windows start 5 Mb PAST the centromere end, so the CEN↔ARM transition
(pericentromere) is excluded — we are not measuring junction/transition reads.

-> results/arm_control.tsv  (type · CEN /Mreads · ARM /Mreads · enrichment)
Run with nextflow_env python."""
import importlib, csv
from collections import Counter, defaultdict
import pysam
lp = importlib.import_module("02_leadprov_sm")
from common import SAMPLES, HAPS, CEN, bam_path, OUT

CHRLEN = {"col": {"Chr1": 32640075, "Chr2": 23012915, "Chr3": 26150667, "Chr4": 22582341, "Chr5": 30170985},
          "ler": {"Chr1": 32485061, "Chr2": 21328600, "Chr3": 27335240, "Chr4": 22700724, "Chr5": 30661135}}
ARM_BUFFER = 5_000_000   # exclude pericentromere / CEN-ARM transition
ARM_WIN = 3_000_000      # arm window size per chromosome (distal to CEN)
TYPES = ["DEL", "INS", "DUP", "INV", "BND"]


def arm_windows(hap):
    """Distal-arm window per chrom, 5 Mb past CEN end (transition excluded)."""
    w = {}
    for chrom, (a, b) in CEN[hap].items():
        s = b + ARM_BUFFER
        e = min(s + ARM_WIN, CHRLEN[hap][chrom] - 100_000)
        if e - s > 200_000:
            w[chrom] = (s, e)
    return w


def main():
    arm_cnt = defaultdict(lambda: Counter())   # tissue -> svtype -> n
    arm_reads = defaultdict(int)               # tissue -> reads
    for sample, tis in SAMPLES:
        for hap in HAPS:
            bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
            for chrom, (a, b) in arm_windows(hap).items():
                for r in bam.fetch(chrom, a, b):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary:
                        continue
                    if r.mapping_quality < lp.MAPQ_MIN or not (a <= r.reference_start < b):
                        continue
                    arm_reads[tis] += 1
                    for svtype, pos, svlen in lp.cigar_leads(r):
                        if a <= pos < b:
                            arm_cnt[tis][svtype] += 1
                    for svtype, pos, svlen, mapq, mate in lp.split_leads(r, chrom):
                        if mapq >= lp.MAPQ_MIN and a <= pos < b:
                            arm_cnt[tis][svtype] += 1
            bam.close()
            print(f"{sample} {hap}: arm done")

    # CEN leadprov rate (from leadprov_sm.tsv) for comparison
    cen_cnt = defaultdict(lambda: Counter()); cen_reads = defaultdict(int)
    for r in csv.DictReader(open(f"{OUT}/leadprov_sm.tsv"), delimiter="\t"):
        cen_cnt[r["tissue"]][r["svtype"]] += 1
    with open(f"{OUT}/cen_read_counts.tsv") as f:
        next(f)
        for ln in f:
            s, h, nr = ln.split()
            cen_reads["leaf" if s == "wt_leaf" else "pollen"] += int(nr)

    rows = []
    for tis in ("leaf", "pollen"):
        cm = cen_reads[tis] / 1e6; am = arm_reads[tis] / 1e6
        for t in TYPES:
            cr = cen_cnt[tis][t] / cm if cm else 0
            ar = arm_cnt[tis][t] / am if am else 0
            rows.append({"tissue": tis, "svtype": t, "cen_per_Mreads": round(cr, 1),
                         "arm_per_Mreads": round(ar, 1), "enrichment_CEN_over_ARM": round(cr / ar, 1) if ar else float("inf")})

    with open(f"{OUT}/arm_control.tsv", "w") as f:
        cols = ["tissue", "svtype", "cen_per_Mreads", "arm_per_Mreads", "enrichment_CEN_over_ARM"]
        f.write("\t".join(cols) + "\n")
        for d in rows:
            f.write("\t".join(str(d[c]) for c in cols) + "\n")

    print(f"\narm reads: leaf {arm_reads['leaf']}, pollen {arm_reads['pollen']}  (windows: CEN_end+{ARM_BUFFER//10**6}Mb, {ARM_WIN//10**6}Mb wide)")
    print(f"{'tissue':8}{'type':5}{'CEN/Mr':>9}{'ARM/Mr':>9}{'CEN÷ARM':>9}")
    for d in rows:
        print(f"{d['tissue']:8}{d['svtype']:5}{d['cen_per_Mreads']:9.1f}{d['arm_per_Mreads']:9.1f}{str(d['enrichment_CEN_over_ARM']):>9}")
    print("DONE_ARM")


if __name__ == "__main__":
    main()
