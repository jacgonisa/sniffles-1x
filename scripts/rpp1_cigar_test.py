#!/usr/bin/env python3
"""Is the RPP1/DM2 cluster CIGAR-indel (DEL/INS) rate different from surrounding arm regions, and
between leaf and pollen? CIGAR-only (the trustworthy inline-indel channel; INV/BND are split-based
and excluded here). Col haplotype, WT leaf + WT pollen.

Design: count CIGAR DEL/INS (>=50bp) per read in the RPP1 cluster window and in many matched-size
(310 kb) arm control windows tiled across all 5 chromosome arms (pericentromere + telomeres excluded).
Rate = (DEL+INS) events per read. Tests:
  (1) cluster vs arm background, per tissue — empirical percentile among arm windows + Poisson exact
      (observed cluster events vs expected from the pooled-arm rate at the cluster's read depth).
  (2) leaf vs pollen — Poisson rate-ratio (conditional-binomial) test, in the cluster and in pooled arms,
      plus a paired Wilcoxon of per-window rates.
-> results/rpp1/rpp1_cigar_test.tsv + printed + rpp1_cigar_test.png
Run: /home/jg2070/miniforge3/envs/nextflow_env/bin/python rpp1_cigar_test.py"""
import os, csv, random, importlib
import pysam
from scipy import stats
lp = importlib.import_module("02_leadprov_sm")
from common import CEN, CHRLEN

CEN_WINDOWS = [(c, a, b) for c, (a, b) in CEN["col"].items()]  # the 5 centromeres (satellite hotspot ref)

ROOT = "/mnt/ssd-4tb/HIFI_NAMIL"
OUT = f"{ROOT}/single_molecule_sv/results/rpp1"
WINSZ = 310_000
CLUSTER = ("Chr3", 19_150_000, 19_460_000)
CEN_BUF, TEL_BUF = 3_000_000, 300_000
N_CTRL = 60
DATASETS = [("wt_leaf", "leaf"), ("wt_pollen", "pollen")]


def bam_path(sample):
    return f"{ROOT}/sv_calling/aligned/{sample}/strict90/col_all.bam"


def arm_windows():
    """Matched-size arm windows on both arms of every chrom, excluding pericentromere + telomeres."""
    w = []
    for c, (ca, cb) in CEN["col"].items():
        L = CHRLEN["col"][c]
        for s in range(TEL_BUF, ca - CEN_BUF - WINSZ, WINSZ):          # left arm
            w.append((c, s, s + WINSZ))
        for s in range(cb + CEN_BUF, L - TEL_BUF - WINSZ, WINSZ):      # right arm
            w.append((c, s, s + WINSZ))
    # drop any window overlapping the cluster
    cc, ca_, cb_ = CLUSTER
    w = [(c, a, b) for (c, a, b) in w if not (c == cc and a < cb_ and b > ca_)]
    random.Random(7).shuffle(w)
    return w[:N_CTRL]


def count(bam, chrom, a, b):
    """(reads, DEL, INS) for CIGAR indels >=50bp with junction in [a,b), reads anchored in [a,b)."""
    reads = dele = ins = 0
    for r in bam.fetch(chrom, a, b):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < lp.MAPQ_MIN:
            continue
        if not (a <= r.reference_start < b):
            continue
        reads += 1
        for svtype, pos, svlen in lp.cigar_leads(r):
            if a <= pos < b:
                if svtype == "DEL":
                    dele += 1
                elif svtype == "INS":
                    ins += 1
    return reads, dele, ins


def rate_ratio_test(n1, e1, n2, e2):
    """Poisson rate-ratio test: rate1 (n1 events / e1 exposure) vs rate2. Conditional-binomial p."""
    N = n1 + n2
    if N == 0 or e1 == 0 or e2 == 0:
        return float("nan"), float("nan")
    p = e1 / (e1 + e2)
    rr = (n1 / e1) / (n2 / e2) if n2 and e2 and (n2 / e2) else float("inf")
    return rr, stats.binomtest(n1, N, p).pvalue


def poisson_vs_background(obs, reads, bg_rate):
    """Two-sided Poisson exact: obs cluster events vs expected = bg_rate*reads."""
    exp = bg_rate * reads
    if exp <= 0:
        return exp, float("nan")
    p = min(1.0, 2 * min(stats.poisson.cdf(obs, exp), stats.poisson.sf(obs - 1, exp)))
    return exp, p


def main():
    os.makedirs(OUT, exist_ok=True)
    ctrl = arm_windows()
    data = {}          # tissue -> {"cluster":(reads,del,ins), "ctrl":[(reads,del,ins),...]}
    rows = []
    for sample, tis in DATASETS:
        bam = pysam.AlignmentFile(bam_path(sample), "rb")
        cl = count(bam, *CLUSTER)
        cw = [count(bam, c, a, b) for (c, a, b) in ctrl]
        cen = [count(bam, c, a, b) for (c, a, b) in CEN_WINDOWS]      # pooled over the 5 CENs
        cen_pool = (sum(x[0] for x in cen), sum(x[1] for x in cen), sum(x[2] for x in cen))
        bam.close()
        data[tis] = {"cluster": cl, "ctrl": cw, "cen": cen_pool}
        rows.append(("RPP1_cluster", tis, *cl))
        for (c, a, b), cc in zip(ctrl, cw):
            rows.append((f"{c}:{a}", tis, *cc))

    with open(f"{OUT}/rpp1_cigar_test.tsv", "w") as f:
        f.write("region\ttissue\treads\tDEL\tINS\tDELINS_per_1k_reads\n")
        for reg, tis, rd, de, ins in rows:
            rate = 1000 * (de + ins) / rd if rd else 0
            f.write(f"{reg}\t{tis}\t{rd}\t{de}\t{ins}\t{rate:.3f}\n")

    print(f"=== CIGAR DEL/INS rate — RPP1 cluster vs {len(ctrl)} matched arm windows (310kb), Col hap ===\n")
    per1k = {}
    for tis in ("leaf", "pollen"):
        cl = data[tis]["cluster"]; cw = data[tis]["ctrl"]
        cl_rate = 1000 * (cl[1] + cl[2]) / cl[0]
        ctrl_rates = sorted(1000 * (d + i) / rd for rd, d, i in cw if rd)
        # pooled arm background rate (events per read)
        bg = sum(d + i for _, d, i in cw) / max(sum(rd for rd, _, _ in cw), 1)
        exp, ppois = poisson_vs_background(cl[1] + cl[2], cl[0], bg)
        pct = 100 * sum(1 for r in ctrl_rates if r <= cl_rate) / len(ctrl_rates)
        med = ctrl_rates[len(ctrl_rates) // 2]
        per1k[tis] = (cl_rate, ctrl_rates, cl)
        print(f"[{tis}] cluster: {cl[0]} reads, {cl[1]} DEL + {cl[2]} INS = {cl[1]+cl[2]} events, "
              f"rate={cl_rate:.2f}/1k reads")
        print(f"       arm windows: median={med:.2f}/1k  (pooled arm rate {1000*bg:.2f}/1k)")
        print(f"       cluster vs arm background: expected {exp:.1f} events, observed {cl[1]+cl[2]} "
              f"-> Poisson p={ppois:.3f};  cluster is at the {pct:.0f}th percentile of arm windows\n")

    # compartment comparison: CEN (satellite hotspot) vs arm background vs RPP1 cluster
    print("=== compartment comparison — CIGAR DEL/INS per 1000 reads (Col hap) ===")
    print(f"{'compartment':16}{'leaf /1k':>10}{'pollen /1k':>12}{'pollen/leaf':>13}")
    def rate(t):
        rd, d, i = t; return 1000 * (d + i) / rd if rd else 0
    comps = [("centromere", data["leaf"]["cen"], data["pollen"]["cen"]),
             ("arm background", (sum(rd for rd, _, _ in data["leaf"]["ctrl"]),
                                 sum(d for _, d, _ in data["leaf"]["ctrl"]),
                                 sum(i for _, _, i in data["leaf"]["ctrl"])),
                                (sum(rd for rd, _, _ in data["pollen"]["ctrl"]),
                                 sum(d for _, d, _ in data["pollen"]["ctrl"]),
                                 sum(i for _, _, i in data["pollen"]["ctrl"]))),
             ("RPP1 cluster", data["leaf"]["cluster"], data["pollen"]["cluster"])]
    comp_rows = []
    for name, lt, pt in comps:
        lr, pr = rate(lt), rate(pt)
        print(f"{name:16}{lr:10.2f}{pr:12.2f}{(pr/lr if lr else float('nan')):13.2f}")
        comp_rows.append((name, lt, pt, lr, pr))
    cenL, cenP = data["leaf"]["cen"], data["pollen"]["cen"]
    print(f"  -> CEN is {rate(cenL)/max(rate(comps[1][1]),1e-9):.0f}x the arm background in leaf, "
          f"{rate(cenP)/max(rate(comps[1][2]),1e-9):.0f}x in pollen\n")

    # leaf vs pollen
    print("=== leaf vs pollen (CIGAR DEL/INS rate ratio) ===")
    clL, clP = data["leaf"]["cluster"], data["pollen"]["cluster"]
    rr, p = rate_ratio_test(clL[1] + clL[2], clL[0], clP[1] + clP[2], clP[0])
    print(f"  RPP1 cluster: leaf {1000*(clL[1]+clL[2])/clL[0]:.2f} vs pollen "
          f"{1000*(clP[1]+clP[2])/clP[0]:.2f} /1k  -> rate ratio(L/P)={rr:.2f}, p={p:.3f}")
    poolL = (sum(d + i for _, d, i in data["leaf"]["ctrl"]), sum(rd for rd, _, _ in data["leaf"]["ctrl"]))
    poolP = (sum(d + i for _, d, i in data["pollen"]["ctrl"]), sum(rd for rd, _, _ in data["pollen"]["ctrl"]))
    rr2, p2 = rate_ratio_test(poolL[0], poolL[1], poolP[0], poolP[1])
    print(f"  pooled arms : leaf {1000*poolL[0]/poolL[1]:.2f} vs pollen {1000*poolP[0]/poolP[1]:.2f} /1k  "
          f"-> rate ratio(L/P)={rr2:.2f}, p={p2:.4g}  (n_L={poolL[0]}, n_P={poolP[0]} events)")
    rr3, p3 = rate_ratio_test(cenL[1] + cenL[2], cenL[0], cenP[1] + cenP[2], cenP[0])
    print(f"  centromere  : leaf {rate(cenL):.2f} vs pollen {rate(cenP):.2f} /1k  "
          f"-> rate ratio(L/P)={rr3:.2f}, p={p3:.4g}  (n_L={cenL[1]+cenL[2]}, n_P={cenP[1]+cenP[2]} events)")
    # paired per-window Wilcoxon
    Lr = [1000 * (d + i) / rd if rd else 0 for rd, d, i in data["leaf"]["ctrl"]]
    Pr = [1000 * (d + i) / rd if rd else 0 for rd, d, i in data["pollen"]["ctrl"]]
    try:
        wp = stats.wilcoxon(Lr, Pr).pvalue
    except ValueError:
        wp = float("nan")
    print(f"  paired per-window leaf vs pollen rate: Wilcoxon p={wp:.4g}")

    plot(per1k, comp_rows)
    print("\nDONE_CIGAR_TEST")


def plot(per1k, comp_rows):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, numpy as np
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.3))
    for ax, tis, col in [(axes[0], "leaf", "#2c7fb8"), (axes[1], "pollen", "#d95f0e")]:
        cl_rate, ctrl_rates, cl = per1k[tis]
        ax.hist(ctrl_rates, bins=15, color=col, alpha=0.6, label="arm windows")
        ax.axvline(cl_rate, color="red", lw=2.5, label=f"RPP1 cluster ({cl_rate:.2f})")
        ax.axvline(np.median(ctrl_rates), color="black", ls="--", lw=1, label=f"arm median ({np.median(ctrl_rates):.2f})")
        ax.set_xlabel("CIGAR DEL+INS per 1000 reads"); ax.set_ylabel("# arm windows")
        ax.set_title(f"{tis}  (cluster {cl[1]}+{cl[2]} events / {cl[0]} reads)"); ax.legend(fontsize=8)
    # compartment bar: CEN vs arm vs cluster, log scale (CEN dwarfs the arms)
    ax = axes[2]
    names = [c[0] for c in comp_rows]; x = np.arange(len(names)); w = 0.38
    ax.bar(x - w / 2, [c[3] for c in comp_rows], w, color="#2c7fb8", label="leaf")
    ax.bar(x + w / 2, [c[4] for c in comp_rows], w, color="#d95f0e", label="pollen")
    ax.set_yscale("log"); ax.set_ylabel("CIGAR DEL+INS per 1000 reads (log)")
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right")
    ax.set_title("compartment comparison"); ax.legend(fontsize=8)
    for xi, c in zip(x, comp_rows):
        ax.text(xi - w / 2, c[3], f"{c[3]:.2f}", ha="center", va="bottom", fontsize=7)
        ax.text(xi + w / 2, c[4], f"{c[4]:.2f}", ha="center", va="bottom", fontsize=7)
    fig.suptitle("RPP1/DM2 CIGAR indel rate vs arm windows, and vs the centromere (Col hap)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(f"{OUT}/rpp1_cigar_test.png", dpi=130); plt.close(fig)
    print(f"figure -> {OUT}/rpp1_cigar_test.png")


if __name__ == "__main__":
    main()
