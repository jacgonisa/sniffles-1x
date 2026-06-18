#!/usr/bin/env python3
"""Step 17 — where do the SV calls come from, and how in-register is each source?

Three detection routes (one read each):
  CIGAR        inline I/D in a single alignment
  SPLITREAD    the ALIGNER already split the read (SA tag)            — Sniffles-compatible
  SPLITANDMAP  WE split the read at the contrast frontier and re-mapped — not in stock Sniffles
A call's `methods` is the union of routes that found it. in_register uses the whole-CEN178-monomer
heuristic (|svlen| mod 178; the `in_phase` column); the rigorous TRASH flanking-monomer version is
in step 13 / report §9 (singletons).

-> results/source_breakdown.tsv (+ printed).
Run with nextflow_env python."""
import csv
from collections import defaultdict, Counter
from common import OUT

ORDER = ["CIGAR", "SPLITREAD", "SPLITANDMAP", "CIGAR+SPLITANDMAP"]


def main():
    by = defaultdict(lambda: {"n": 0, "reg": 0, "regden": 0, "types": Counter()})
    bytis = defaultdict(lambda: defaultdict(lambda: {"n": 0, "reg": 0, "regden": 0}))
    for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t"):
        m = r["methods"]; ip = r["in_phase"]
        d = by[m]; d["n"] += 1; d["types"][r["svtype"]] += 1
        if ip in ("0", "1"):     # in_phase computed (svlen present, not BND)
            d["regden"] += 1; d["reg"] += int(ip == "1")
        t = bytis[r["tissue"]][m]; t["n"] += 1
        if ip in ("0", "1"):
            t["regden"] += 1; t["reg"] += int(ip == "1")

    cols = ["methods", "calls", "in_register", "pct_in_register", "DEL", "INS", "DUP", "INV", "BND"]
    with open(f"{OUT}/source_breakdown.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        keys = [k for k in ORDER if k in by] + [k for k in by if k not in ORDER]
        for k in keys:
            d = by[k]; pct = 100 * d["reg"] / d["regden"] if d["regden"] else 0
            f.write(f"{k}\t{d['n']}\t{d['reg']}/{d['regden']}\t{pct:.1f}\t"
                    + "\t".join(str(d['types'][t]) for t in ('DEL', 'INS', 'DUP', 'INV', 'BND')) + "\n")

    print(f"{'route':22}{'calls':>7}{'in-register':>14}{'%':>7}")
    for k in keys:
        d = by[k]; pct = 100 * d["reg"] / d["regden"] if d["regden"] else 0
        print(f"{k:22}{d['n']:7d}{(str(d['reg'])+'/'+str(d['regden'])):>14}{pct:7.1f}")
    print("--- by tissue ---")
    for tis in ("leaf", "pollen"):
        for k in keys:
            t = bytis[tis][k]
            if t["n"]:
                pct = 100 * t["reg"] / t["regden"] if t["regden"] else 0
                print(f"  {tis:7} {k:20} calls={t['n']:5d} in_register={pct:.1f}%")
    print("DONE_SOURCE")


if __name__ == "__main__":
    main()
