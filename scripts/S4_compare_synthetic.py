#!/usr/bin/env python3
"""Synthetic-control step 4 — real vs synthetic SV rate = the mapping/sequencing artefact floor.
Any SV on synthetic (assembly-derived) reads is an artefact; real − synthetic = biological signal.
CEN rates per Mb of CEN-mapped reads (from sm_sv_calls + cen_mapped_mb); ARM rates per million arm reads
(from arm_control.tsv). Also the in-register whole-monomer DEL/INS rate (the unequal-HR signature) — its
floor should be ~0.
-> results/synthetic_comparison.tsv  + printed. (figure lives in the combined report.)"""
import os, csv
from collections import defaultdict, Counter

ROOT = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv"
REAL, SYN = f"{ROOT}/results", f"{ROOT}/results_synthetic"
GROUPS = ["wt_leaf", "cenh3ox_leaf", "wt_pollen", "cenh3ox_pollen"]
TYPES = ["DEL", "INS", "DUP", "INV", "BND"]


def cen_mb(out):
    d = defaultdict(float)
    p = f"{out}/cen_mapped_mb.tsv"
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            d[r["sample"]] += float(r["cen_mapped_mb"])
    return d


def cen_rates(out):
    """per group: {type: calls/Mb}, plus INREG = in-register DEL/INS per Mb."""
    mb = cen_mb(out)
    cnt = defaultdict(Counter)
    p = f"{out}/sm_sv_calls.tsv"
    if not os.path.exists(p):
        return {}, mb
    for r in csv.DictReader(open(p), delimiter="\t"):
        g = r["sample"]; cnt[g][r["svtype"]] += 1
        if r["svtype"] in ("DEL", "INS") and r.get("in_phase") == "1":
            cnt[g]["INREG"] += 1
    rate = {g: {t: (cnt[g][t] / mb[g] if mb.get(g) else 0) for t in TYPES + ["INREG"]} for g in GROUPS}
    return rate, mb


def arm_rates(out):
    d = defaultdict(dict)
    p = f"{out}/arm_control.tsv"
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            d[r["group"]][r["svtype"]] = float(r["arm_per_Mreads"])
    return d


def main():
    rc, _ = cen_rates(REAL); sc, _ = cen_rates(SYN)
    ra, sa = arm_rates(REAL), arm_rates(SYN)
    rows = []
    for g in GROUPS:
        for t in ["INREG"] + TYPES:
            rr, ss = rc.get(g, {}).get(t, 0), sc.get(g, {}).get(t, 0)
            rows.append(("CEN", g, t, rr, ss))
        for t in TYPES:
            rr, ss = ra.get(g, {}).get(t, 0), sa.get(g, {}).get(t, 0)
            rows.append(("ARM", g, t, rr, ss))
    with open(f"{REAL}/synthetic_comparison.tsv", "w") as f:
        f.write("compartment\tgroup\ttype\treal_rate\tsynthetic_floor\tbiological_excess\tpct_floor\n")
        for comp, g, t, rr, ss in rows:
            exc = rr - ss; pct = 100 * ss / rr if rr else 0
            f.write(f"{comp}\t{g}\t{t}\t{rr:.3f}\t{ss:.3f}\t{exc:.3f}\t{pct:.0f}\n")

    print(f"{'comp':5}{'group':16}{'type':7}{'real':>9}{'synth(floor)':>14}{'%floor':>8}")
    for comp, g, t, rr, ss in rows:
        if comp == "CEN" and t in ("INREG", "DEL", "INS", "DUP", "BND"):
            print(f"{comp:5}{g:16}{t:7}{rr:9.3f}{ss:14.3f}{(100*ss/rr if rr else 0):7.0f}%")
    print("... full table -> results/synthetic_comparison.tsv")

    figure(rc, sc, ra, sa)
    print("DONE_SYNTH_COMPARE")


def figure(rc, sc, ra, sa):
    """CEN + ARM real rate (bars) vs synthetic floor (line). Floor ~0 => all signal is real.
    Two panels because CEN is per Mb of CEN-mapped seq and ARM is per million arm reads."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, numpy as np
    col = {"wt_leaf": "#2c7fb8", "cenh3ox_leaf": "#7fcdbb",
           "wt_pollen": "#d95f0e", "cenh3ox_pollen": "#fec44f"}
    lab = {"INREG": "in-register\nDEL/INS", "DEL": "DEL", "INS": "INS",
           "DUP": "DUP", "INV": "INV", "BND": "BND"}
    panels = [("CEN", ["INREG", "DEL", "INS", "DUP", "BND"], rc, sc,
               "real events per Mb of CEN-mapped read seq"),
              ("ARM", ["DEL", "INS", "DUP", "INV", "BND"], ra, sa,
               "real events per million arm reads")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    for ax, (comp, cats, real_d, syn_d, ylab) in zip(axes, panels):
        x = np.arange(len(cats)); w = 0.2; floor_max = 0.0
        for i, g in enumerate(GROUPS):
            real = [real_d.get(g, {}).get(t, 0) for t in cats]
            floor_max = max(floor_max, max((syn_d.get(g, {}).get(t, 0) for t in cats), default=0))
            ax.bar(x + (i - 1.5) * w, real, w, color=col[g], label=g.replace("_", " "))
        ax.axhline(floor_max, ls="--", lw=1.5, color="red")
        ax.text(len(cats) - 0.4, floor_max, f"synthetic floor = {floor_max:.3f}  ",
                color="red", va="bottom", ha="right", fontsize=9, fontweight="bold")
        ax.set_xticks(x); ax.set_xticklabels([lab[c] for c in cats])
        ax.set_ylabel(ylab); ax.set_title(f"{comp}")
    axes[0].legend(fontsize=8, ncol=2)
    fig.suptitle("Single-molecule SV rate vs synthetic mapping-artefact floor — CEN and ARMS\n"
                 "(reads simulated from variant-free assemblies, identical pipeline; floor = 0 in both)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(f"{REAL}/synthetic_control.png", dpi=130); plt.close(fig)
    print(f"figure -> {REAL}/synthetic_control.png")


if __name__ == "__main__":
    main()
