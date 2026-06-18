#!/usr/bin/env python3
"""Step 11 — are the vertical lines (recurrent positions) hotspots or artifacts?

A vertical line = many reads with the same SV type at ~the same coordinate. For a
SINGLE-molecule caller this is the key discriminator:
  - high VAF (support / spanning coverage) + tight size  => the read population agrees
    with each other but DISAGREES with the reference => a FIXED difference between the
    sample and the reference assembly (germline/assembly discrepancy), NOT a somatic
    hotspot. It is "recurrent" only because every spanning read shows the same fixed event.
  - low VAF + variable breakpoints/size => a minority, recurrent signal => genuine
    somatic-hotspot candidate (or satellite mapping ambiguity).

Bins calls per (sample,hap,chrom,200bp,svtype), counts supporting reads, fetches
spanning coverage from the BAM, computes VAF, and classifies.
-> results/recurrent_loci.tsv  (+ printed verdict).
Run with nextflow_env python."""
import csv, statistics as st, sys
from collections import defaultdict
import pysam
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts")
from common import bam_path, OUT, in_phase

BINW = 200
MINSUP = 10          # a "vertical line"
VAF_FIXED = 0.30     # >= this = fixed difference vs reference


def main():
    rows = []
    with open(f"{OUT}/sm_sv_calls.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["pos"] = int(r["pos"]); r["svlen"] = int(r["svlen"]) if r["svlen"] not in ("", "None") else 0
            rows.append(r)
    # bin
    binned = defaultdict(list)
    for r in rows:
        if r["svtype"] == "BND":
            continue
        key = (r["sample"], r["hap"], r["chrom"], (r["pos"] // BINW) * BINW, r["svtype"])
        binned[key].append(r)

    bams = {}
    def cov_at(sample, hap, chrom, pos):
        k = (sample, hap)
        if k not in bams:
            bams[k] = pysam.AlignmentFile(bam_path(sample, hap), "rb")
        return bams[k].count(chrom, pos, pos + 1,
                             read_callback=lambda r: not (r.is_secondary or r.is_supplementary or r.is_unmapped))

    out = []
    for (sample, hap, chrom, pos, svt), g in binned.items():
        sup = len(g)
        if sup < MINSUP:
            continue
        center = pos + BINW // 2
        cov = cov_at(sample, hap, chrom, center)
        vaf = sup / cov if cov else 0
        sizes = [abs(x["svlen"]) for x in g]
        med = int(st.median(sizes))
        spread = max(sizes) - min(sizes)
        ph, rem = in_phase(med)
        cls = "FIXED_vs_ref" if vaf >= VAF_FIXED else "hotspot_candidate"
        out.append({"sample": sample, "hap": hap, "chrom": chrom, "pos": center, "svtype": svt,
                    "support": sup, "coverage": cov, "vaf": round(vaf, 3),
                    "median_size": med, "size_spread": spread,
                    "in_phase": int(ph), "class": cls})

    out.sort(key=lambda d: -d["support"])
    cols = ["sample", "hap", "chrom", "pos", "svtype", "support", "coverage", "vaf",
            "median_size", "size_spread", "in_phase", "class"]
    with open(f"{OUT}/recurrent_loci.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for d in out:
            f.write("\t".join(str(d[c]) for c in cols) + "\n")

    nfixed = sum(1 for d in out if d["class"] == "FIXED_vs_ref")
    ncand = len(out) - nfixed
    print(f"recurrent loci (>= {MINSUP} reads): {len(out)}  ->  FIXED_vs_ref {nfixed}, hotspot_candidate {ncand}")
    print(f"{'sample':9}{'hap':4}{'chrom':6}{'pos':>10}{'type':5}{'sup':>5}{'cov':>5}{'vaf':>6}{'medsz':>7}{'class':>18}")
    for d in out[:15]:
        print(f"{d['sample']:9}{d['hap']:4}{d['chrom']:6}{d['pos']:10d}{d['svtype']:5}{d['support']:5d}{d['coverage']:5d}{d['vaf']:6.2f}{d['median_size']:7d}{d['class']:>18}")
    print("DONE_RECUR")


if __name__ == "__main__":
    main()
