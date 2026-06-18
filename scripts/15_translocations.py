#!/usr/bin/env python3
"""Step 15 — summarize the translocation (BND) calls.

BND = a read whose two fragments map to DIFFERENT contigs (sv.classify_splits). The partner
locus is in the `mate` column (mate_contig:mate_ref_start). We categorize each mate:
  other_CEN          mate in another chromosome's centromere  -> likely CEN178 cross-mapping
  other_chrom_arm    mate on another chromosome, outside its CEN -> possible real junction
  unplaced_organellar mate on a ptg*/organellar contig          -> mapping artifact
-> results/translocations.tsv (+ printed summary).
Run with nextflow_env python."""
import csv
from collections import Counter
from common import CEN, OUT


def cat(hap, mate):
    if not mate or ":" not in mate:
        return "none"
    mc, mp = mate.rsplit(":", 1)
    try:
        mp = int(mp)
    except ValueError:
        return "none"
    if mc in CEN[hap]:
        a, b = CEN[hap][mc]
        return "other_CEN" if a <= mp < b else "other_chrom_arm"
    return "unplaced_organellar"


def main():
    out = []
    for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t"):
        if r["svtype"] != "BND":
            continue
        c = cat(r["hap"], r.get("mate", ""))
        out.append({"sample": r["sample"], "hap": r["hap"], "tissue": r["tissue"],
                    "chrom": r["chrom"], "pos": r["pos"], "mate": r.get("mate", ""),
                    "category": c, "methods": r["methods"], "mapq": r["mapq"], "read": r["read"]})
    cols = ["sample", "hap", "tissue", "chrom", "pos", "mate", "category", "methods", "mapq", "read"]
    with open(f"{OUT}/translocations.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for d in sorted(out, key=lambda d: (d["category"], d["sample"], d["chrom"])):
            f.write("\t".join(str(d[c]) for c in cols) + "\n")

    print(f"BND / translocation calls: {len(out)}")
    print("by category:", dict(Counter(d["category"] for d in out)))
    print("by tissue  :", dict(Counter(d["tissue"] for d in out)))
    print("by method  :", dict(Counter(d["methods"] for d in out)))

    # inter-CEN BND = satellite cross-mapping noise background (per million CEN reads).
    # These are fragments demonstrably landing in the WRONG centromere (shared CEN178),
    # so their rate is an empirical noise floor for single-molecule satellite split calls.
    denom = {}
    try:
        with open(f"{OUT}/cen_read_counts.tsv") as f:
            next(f)
            for ln in f:
                s, h, nr = ln.split(); denom[(s, h)] = int(nr)
    except FileNotFoundError:
        denom = None
    with open(f"{OUT}/crossmap_background.tsv", "w") as f:
        f.write("tissue\tinterCEN_BND\tcen_reads\tper_million_reads\n")
        print("--- satellite cross-mapping background (inter-CEN BND / million reads) ---")
        for tis, samp in (("leaf", "wt_leaf"), ("pollen", "wt_pollen")):
            n = sum(1 for d in out if d["tissue"] == tis and d["category"] == "other_CEN")
            reads = sum(denom.get((samp, h), 0) for h in ("col", "ler")) if denom else 0
            rate = n / reads * 1e6 if reads else 0
            f.write(f"{tis}\t{n}\t{reads}\t{rate:.1f}\n")
            print(f"  {tis}: {n} inter-CEN BND / {reads} reads = {rate:.1f} per million")
    print("DONE_TRANSLOC")


if __name__ == "__main__":
    main()
