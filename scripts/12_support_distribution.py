#!/usr/bin/env python3
"""Step 12 — distribution of READ SUPPORT per SV locus (not Sniffles VAF).

For a single-molecule caller the meaningful quantity is: how many independent reads
carry the same event. We cluster per-read calls into loci (same sample+hap+chrom+svtype,
positions within 100 bp) and count supporting reads per locus:
  support = 1  -> exclusively single-molecule (singleton; the cleanest somatic candidate)
  support high -> recurrent (fixed-vs-reference or hotspot; see step 11)

-> results/support_distribution.tsv  (tissue, support, n_loci)
   results/singleton_events.tsv       (every support==1 locus = the 1x events)
   results/figures/support_distribution.png
Run with nextflow_env python."""
import os, csv, base64, io, zlib
from collections import defaultdict, Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import OUT, HAPS, in_phase

SAMPLE_ROWS = ["wt_leaf", "wt_pollen"]

FIGDIR = f"{OUT}/figures"; os.makedirs(FIGDIR, exist_ok=True)
GAP = 100   # bp; calls within this on the same chrom/type = one locus
TISSUE = {"wt_leaf": "leaf", "wt_pollen": "pollen"}
TCOL = {"leaf": "#4C9A2A", "pollen": "#E8820C"}
# support-bin layout for the x axis
BINS = [(1, 1, "1"), (2, 2, "2"), (3, 3, "3"), (4, 4, "4"), (5, 5, "5"),
        (6, 10, "6-10"), (11, 20, "11-20"), (21, 50, "21-50"), (51, 10**9, ">50")]


def load():
    rows = []
    with open(f"{OUT}/sm_sv_calls.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["pos"] = int(r["pos"]); r["svlen"] = int(r["svlen"]) if r["svlen"] not in ("", "None") else 0
            rows.append(r)
    return rows


def load_denom():
    d = {}
    with open(f"{OUT}/cen_read_counts.tsv") as f:
        next(f)
        for ln in f:
            s, h, n = ln.split(); d[(s, h)] = int(n)
    return d


def match_reads(rows, denom):
    """Downsample leaf to the same CEN-read count as pollen, per haplotype (deterministic
    per-read keep). Equalizes read budget so the support tail is comparable."""
    N = {h: min(denom[("wt_leaf", h)], denom[("wt_pollen", h)]) for h in HAPS}
    p = {(s, h): N[h] / denom[(s, h)] for s in SAMPLE_ROWS for h in HAPS}
    return [r for r in rows
            if p[(r["sample"], r["hap"])] >= 1
            or (zlib.crc32(r["read"].encode()) & 0xffffffff) / 2**32 < p[(r["sample"], r["hap"])]], N


def cluster_loci(rows):
    """Yield loci dicts: support (distinct reads) + a representative call."""
    groups = defaultdict(list)
    for r in rows:
        if r["svtype"] == "BND":
            continue
        groups[(r["sample"], r["hap"], r["chrom"], r["svtype"])].append(r)
    for (sample, hap, chrom, svt), g in groups.items():
        g.sort(key=lambda r: r["pos"])
        cluster = [g[0]]
        for r in g[1:]:
            if r["pos"] - cluster[-1]["pos"] <= GAP:
                cluster.append(r)
            else:
                yield _locus(sample, hap, chrom, svt, cluster); cluster = [r]
        yield _locus(sample, hap, chrom, svt, cluster)


def _locus(sample, hap, chrom, svt, cluster):
    reads = {r["read"] for r in cluster}
    rep = max(cluster, key=lambda r: abs(r["svlen"]))
    return {"sample": sample, "hap": hap, "tissue": TISSUE[sample], "chrom": chrom,
            "pos": rep["pos"], "svtype": svt, "svlen": rep["svlen"], "support": len(reads),
            "read": rep["read"]}


def binidx(s):
    for i, (lo, hi, _) in enumerate(BINS):
        if lo <= s <= hi:
            return i
    return len(BINS) - 1


def distribution(loci):
    dist = defaultdict(lambda: Counter())
    for L in loci:
        dist[L["tissue"]][L["support"]] += 1
    return dist


def panel(ax, dist, title):
    x = np.arange(len(BINS)); w = 0.4
    for k, tis in enumerate(("leaf", "pollen")):
        h = [sum(dist[tis][s] for s in dist[tis] if lo <= s <= hi) for lo, hi, _ in BINS]
        ax.bar(x + (k - 0.5) * w, h, w, color=TCOL[tis], label=tis)
        for xi, hv in zip(x + (k - 0.5) * w, h):
            if hv:
                ax.text(xi, hv, str(hv), ha="center", va="bottom", fontsize=7)
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels([b[2] for b in BINS], fontsize=8)
    ax.set_xlabel("supporting reads per SV locus  (1 = single-molecule)")
    ax.axvspan(-0.5, 0.5, color="#FDEBD0", zorder=0)
    ax.set_title(title)


def main():
    rows = load()
    denom = load_denom()
    loci_full = list(cluster_loci(rows))
    matched_rows, N = match_reads(rows, denom)
    loci_match = list(cluster_loci(matched_rows))

    # full distribution table
    dist_full = distribution(loci_full)
    dist_match = distribution(loci_match)
    with open(f"{OUT}/support_distribution.tsv", "w") as f:
        f.write("set\ttissue\tsupport\tn_loci\n")
        for tag, dist in (("all_reads", dist_full), ("read_budget_matched", dist_match)):
            for tis in ("leaf", "pollen"):
                for s in sorted(dist[tis]):
                    f.write(f"{tag}\t{tis}\t{s}\t{dist[tis][s]}\n")

    # singleton (1x) events FROM FULL DATA (discovery set for annotation)
    singles = [L for L in loci_full if L["support"] == 1]
    with open(f"{OUT}/singleton_events.tsv", "w") as f:
        cols = ["sample", "hap", "tissue", "chrom", "pos", "svtype", "svlen", "in_phase", "read"]
        f.write("\t".join(cols) + "\n")
        for L in sorted(singles, key=lambda d: (d["sample"], d["hap"], d["chrom"], d["pos"])):
            ph = int(in_phase(L["svlen"])[0]) if L["svlen"] else ""
            f.write("\t".join(str(x) for x in
                    [L["sample"], L["hap"], L["tissue"], L["chrom"], L["pos"], L["svtype"], L["svlen"], ph, L["read"]]) + "\n")

    # two-panel figure: all reads vs read-budget-matched
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.4), sharey=True)
    panel(axes[0], dist_full, "All reads (leaf ~14× deeper → longer tail)")
    panel(axes[1], dist_match, f"Read-budget matched (leaf downsampled to pollen ≈{N['col']//1000}k/hap)")
    axes[0].set_ylabel("number of SV loci (log)")
    axes[1].legend(title="tissue")
    fig.suptitle("Read-support distribution per centromere SV locus — leaf vs pollen")
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=130, bbox_inches="tight"); plt.close(fig)
    open(f"{FIGDIR}/support_distribution.png", "wb").write(b.getvalue())

    # summary
    print(f"read-budget matched to: {N}")
    for tag, loci in (("ALL READS", loci_full), ("MATCHED", loci_match)):
        print(f"--- {tag} ---")
        print(f"{'tissue':8}{'loci':>8}{'1x':>8}{'1x%':>7}{'>=2x':>7}{'maxsup':>8}")
        for tis in ("leaf", "pollen"):
            tl = [L for L in loci if L["tissue"] == tis]
            n = len(tl); one = sum(1 for L in tl if L["support"] == 1)
            mx = max((L["support"] for L in tl), default=0)
            print(f"{tis:8}{n:8d}{one:8d}{100*one/max(n,1):7.1f}{n-one:7d}{mx:8d}")
    print(f"singleton (1x) events (full data) written: {len(singles)}")
    print("DONE_SUPPORT")


if __name__ == "__main__":
    main()
