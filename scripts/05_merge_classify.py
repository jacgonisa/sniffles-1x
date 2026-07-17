#!/usr/bin/env python3
"""Step 5 — union of the two single-molecule detectors, dedup, annotate, concordance.

Inputs : results/leadprov_sm.tsv  (step 2, sniffles topology per read)
         results/split_and_map.tsv (step 3, re-mapped fragments)
         results/stock/*.vcf        (step 4, stock sniffles --minsupport 1)
Dedup  : same read + svtype + pos within 100 bp = one event (keep method flags).
Annotate: CEN178 register (in_phase / monomer remainder) — the unequal-HR signature.
Concord : flag whether a stock-sniffles call of same type sits within 200 bp.
-> results/sm_sv_calls.tsv  + prints a concordance summary.
"""
import os, glob, csv
from collections import defaultdict
from common import OUT, in_phase, in_cen, refkey

POS_DEDUP = 100
POS_STOCK = 200


def load(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            r["pos"] = int(r["pos"])
            r["svlen"] = int(r["svlen"]) if r["svlen"] not in ("", "None") else 0
            rows.append(r)
    return rows


def load_stock():
    """{(sample,hap): [(chrom,pos,svtype), ...]}"""
    st = defaultdict(list)
    for vcf in glob.glob(f"{OUT}/stock/*.vcf"):
        base = os.path.basename(vcf)[:-4]            # wt_leaf_col
        sample, hap = base.rsplit("_", 1)
        for ln in open(vcf):
            if ln.startswith("#"):
                continue
            c = ln.split("\t")
            chrom, pos, info = c[0], int(c[1]), c[7]
            svt = next((x[7:] for x in info.split(";") if x.startswith("SVTYPE=")), "NA")
            st[(sample, hap)].append((chrom, pos, svt))
    return st


def stock_hit(stock, key, chrom, pos, svtype):
    for c, p, t in stock.get(key, []):
        if c == chrom and t == svtype and abs(p - pos) <= POS_STOCK:
            return True
    return False


def main():
    rows = load(f"{OUT}/leadprov_sm.tsv") + load(f"{OUT}/split_and_map.tsv")
    stock = load_stock()

    # dedup within (sample,hap,read,svtype) by position
    groups = defaultdict(list)
    for r in rows:
        groups[(r["sample"], r["hap"], r["read"], r["svtype"])].append(r)

    merged = []
    for (sample, hap, read, svtype), g in groups.items():
        g.sort(key=lambda r: r["pos"])
        clusters = []
        for r in g:
            if clusters and abs(r["pos"] - clusters[-1][0]["pos"]) <= POS_DEDUP:
                clusters[-1].append(r)
            else:
                clusters.append([r])
        for cl in clusters:
            _lbl = {"INLINE": "CIGAR", "SPLIT": "SPLITREAD", "SPLITMAP": "SPLITANDMAP"}
            methods = sorted({_lbl.get(r["source"], r["source"]) for r in cl})
            rep = max(cl, key=lambda r: abs(r["svlen"]))     # representative = largest
            ph, rem = in_phase(rep["svlen"]) if rep["svlen"] else (False, "")
            merged.append({
                "sample": sample, "hap": hap, "tissue": rep["tissue"], "chrom": rep["chrom"],
                "pos": rep["pos"], "svtype": svtype, "svlen": rep["svlen"],
                "methods": "+".join(methods),
                "in_phase": int(ph) if rep["svlen"] else "", "monomer_rem": rem,
                "in_cen": int(in_cen(refkey(sample, hap), rep["chrom"], rep["pos"])),
                "mapq": rep["mapq"], "read": read, "mate": rep.get("mate", ""),
                "stock_match": int(stock_hit(stock, (sample, hap), rep["chrom"], rep["pos"], svtype)),
            })

    cols = ["sample", "hap", "tissue", "chrom", "pos", "svtype", "svlen",
            "methods", "in_phase", "monomer_rem", "in_cen", "mapq", "read", "mate", "stock_match"]
    with open(f"{OUT}/sm_sv_calls.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for m in sorted(merged, key=lambda r: (r["sample"], r["hap"], r["chrom"], r["pos"])):
            out.write("\t".join(str(m[c]) for c in cols) + "\n")

    # summary
    n = len(merged)
    print(f"merged single-molecule SV calls: {n}")
    by_method = defaultdict(int); by_type = defaultdict(int); stockm = 0; inph = 0; inph_den = 0
    for m in merged:
        by_method[m["methods"]] += 1
        by_type[m["svtype"]] += 1
        stockm += m["stock_match"]
        if m["in_phase"] != "":
            inph_den += 1; inph += m["in_phase"]
    print("by method:", dict(by_method))
    print("by type  :", dict(by_type))
    print(f"stock-sniffles concordant: {stockm}/{n} ({100*stockm/max(n,1):.1f}%)")
    if inph_den:
        print(f"CEN178 in-register (whole-monomer): {inph}/{inph_den} ({100*inph/inph_den:.1f}%)")
    print("DONE_MERGE")


if __name__ == "__main__":
    main()
