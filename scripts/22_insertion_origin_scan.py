#!/usr/bin/env python3
"""Batch scan (no plots): trace the origin of EVERY INS >= MINSZ, tabulate donor distance,
and flag any whose inserted fragment maps FAR from the insertion site (distal / other-chrom).
-> results/insertion_origin_scan.tsv  (+ printed distance distribution & distal hits)."""
import os, csv, subprocess, tempfile
from collections import defaultdict, Counter
import pysam
from common import OUT, REF, CEN, bam_path, refkey

MINIMAP2 = "/home/jg2070/minimap2/minimap2"
MINSZ = 1000


def find_ins(cigar, target):
    q = 0
    for op, ln in cigar:
        if op in (0, 7, 8): q += ln
        elif op == 1:
            if abs(ln - target) <= 3: return q, q + ln
            q += ln
        elif op == 4: q += ln
    return None, None


def main():
    rows = [r for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t")
            if r["svtype"] == "INS" and "CIGAR" in r["methods"]
            and r["svlen"] not in ("", "None") and abs(int(r["svlen"])) >= MINSZ]
    # dedup nearby recurrent (same sample/hap/chrom within 2 kb) -> keep one
    seen, evs = [], []
    for r in sorted(rows, key=lambda r: -abs(int(r["svlen"]))):
        k = (r["sample"], r["hap"], r["chrom"]); pos = int(r["pos"])
        if any(kk == k and abs(pos - p) < 2000 for kk, p in seen): continue
        seen.append((k, pos)); evs.append(r)
    print(f"{len(evs)} distinct INS >= {MINSZ} bp to trace")

    # extract inserted seq per event, batched per refkey
    recs = defaultdict(list); meta = {}
    for i, r in enumerate(evs):
        s, h, chrom, pos, sz = r["sample"], r["hap"], r["chrom"], int(r["pos"]), abs(int(r["svlen"]))
        rk = refkey(s, h); sid = f"e{i}"
        bam = pysam.AlignmentFile(bam_path(s, h), "rb")
        for rd in bam.fetch(chrom, max(0, pos - 60000), pos + 60000):
            if rd.query_name != r["read"] or rd.is_unmapped or not rd.cigartuples: continue
            qs, qe = find_ins(rd.cigartuples, sz)
            if qs is None: continue
            recs[rk].append((sid, rd.query_sequence[qs:qe]))
            meta[sid] = {"sample": s, "hap": h, "chrom": chrom, "pos": pos, "sz": sz, "rk": rk}
            break
        bam.close()

    hits = defaultdict(list)
    with tempfile.TemporaryDirectory() as td:
        for rk, rr in recs.items():
            if not rr: continue
            fa = f"{td}/{rk}.fa"; paf = f"{td}/{rk}.paf"
            with open(fa, "w") as f:
                for sid, seq in rr: f.write(f">{sid}\n{seq}\n")
            subprocess.run([MINIMAP2, "-x", "map-hifi", "-c", "--secondary=yes", "-N", "8",
                            "-t", "12", REF[rk][0], fa], stdout=open(paf, "w"),
                           stderr=subprocess.DEVNULL, check=True)
            for ln in open(paf):
                c = ln.split("\t")
                if len(c) < 12: continue
                hits[c[0]].append({"chrom": c[5], "start": int(c[7]), "strand": c[4],
                                   "ident": int(c[9]) / max(int(c[10]), 1), "alen": int(c[10])})
    for sid in hits: hits[sid].sort(key=lambda d: -d["alen"] * d["ident"])

    cats = Counter(); dist_bins = Counter(); distal = []
    with open(f"{OUT}/insertion_origin_scan.tsv", "w") as f:
        f.write("sample\thap\tchrom\tpos\tins_bp\tn_hits\tbest_chrom\tbest_pos\tstrand\tident\tdist_bp\tcategory\n")
        for sid, m in meta.items():
            hs = hits.get(sid, [])
            if not hs:
                cats["no_hit"] += 1; continue
            best = hs[0]; same = best["chrom"] == m["chrom"]
            d = abs(best["start"] - m["pos"]) if same else None
            cen = CEN[m["rk"]].get(m["chrom"])
            if not same: cat = "other_chrom"
            elif d < 10_000: cat = "local_<10kb"
            elif d < 100_000: cat = "local_10-100kb"
            elif cen and cen[0] <= best["start"] < cen[1]: cat = "same_CEN_>100kb"
            else: cat = "distal_>100kb"
            cats[cat] += 1
            dist_bins["other_chrom" if not same else
                      ("<1kb" if d < 1000 else "1-10kb" if d < 10000 else
                       "10-100kb" if d < 100000 else ">100kb")] += 1
            if not same or (d is not None and d >= 100_000):
                distal.append((m, best, d, cat))
            f.write(f"{m['sample']}\t{m['hap']}\t{m['chrom']}\t{m['pos']}\t{m['sz']}\t{len(hs)}\t"
                    f"{best['chrom']}\t{best['start']}\t{best['strand']}\t{best['ident']:.3f}\t"
                    f"{d if d is not None else 'NA'}\t{cat}\n")

    print("\n=== origin category ==="); [print(f"  {k:18} {v}") for k, v in cats.most_common()]
    print("\n=== donor distance (best hit) ==="); [print(f"  {k:12} {v}") for k, v in
          sorted(dist_bins.items(), key=lambda x: ["<1kb","1-10kb","10-100kb",">100kb","other_chrom"].index(x[0]))]
    print(f"\n=== {len(distal)} events with DISTAL / other-chrom origin ===")
    for m, best, d, cat in sorted(distal, key=lambda x: (x[2] is None, -(x[2] or 0)))[:25]:
        dtxt = f"{d/1000:.0f} kb" if d is not None else f"-> {best['chrom']}"
        print(f"  {m['sample']:14} {m['hap']} {m['chrom']}:{m['pos']:>9,} {m['sz']:>5}bp  "
              f"origin {best['chrom']}:{best['start']:>9,} ({best['ident']*100:.0f}% id)  {dtxt}  [{cat}]")
    print("DONE_INSSCAN")


if __name__ == "__main__":
    main()
