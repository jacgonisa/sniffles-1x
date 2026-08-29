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
SYND = f"{ROOT}/results_synthetic_dirty"   # dirty floor (mapping + sequencing), error tail matched to real
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
    rc, _ = cen_rates(REAL); sc, _ = cen_rates(SYN); dc, _ = cen_rates(SYND)
    ra, sa, da = arm_rates(REAL), arm_rates(SYN), arm_rates(SYND)
    rows = []
    for g in GROUPS:
        for t in ["INREG"] + TYPES:
            rows.append(("CEN", g, t, rc.get(g, {}).get(t, 0),
                         sc.get(g, {}).get(t, 0), dc.get(g, {}).get(t, 0)))
        for t in TYPES:
            rows.append(("ARM", g, t, ra.get(g, {}).get(t, 0),
                         sa.get(g, {}).get(t, 0), da.get(g, {}).get(t, 0)))
    # biological_excess is computed over the MAPPING floor (clean); the dirty run is a sequencing
    # stress-test whose MAGNITUDE overshoots real (glitch model harsher than real CCS failure) — its
    # value is the TYPE pattern (only INS/DEL are ever generated). seq_stress kept for that column.
    with open(f"{REAL}/synthetic_comparison.tsv", "w") as f:
        f.write("compartment\tgroup\ttype\treal_rate\tmapping_floor\tseq_stress\t"
                "biological_excess\tpct_mapping_floor\n")
        for comp, g, t, rr, cs, ds in rows:
            exc = rr - cs; pct = 100 * cs / rr if rr else 0
            f.write(f"{comp}\t{g}\t{t}\t{rr:.3f}\t{cs:.3f}\t{ds:.3f}\t{exc:.3f}\t{pct:.0f}\n")

    print(f"{'comp':5}{'group':16}{'type':7}{'real':>9}{'map.floor':>10}{'seq.stress':>11}")
    for comp, g, t, rr, cs, ds in rows:
        if comp == "ARM" or (comp == "CEN" and t in ("INREG", "DEL", "INS", "DUP", "BND")):
            print(f"{comp:5}{g:16}{t:7}{rr:9.3f}{cs:10.3f}{ds:11.1f}")
    print("... full table -> results/synthetic_comparison.tsv")

    detector_breakdown()
    periodicity()
    figure(rc, sc, dc, ra, sa, da)
    print("DONE_SYNTH_COMPARE")


def periodicity():
    """Do CIGAR DEL/INS sizes show 178-bp CEN178 periodicity? Real CEN = 90% whole-monomer (unequal-HR
    signature); the dirty/artefact CIGAR indels are random-sized. This periodicity is a biological
    fingerprint sequencing error does NOT reproduce."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, numpy as np

    def sizes(path):
        if not os.path.exists(path):
            return np.array([])
        return np.array([abs(int(r["svlen"])) for r in csv.DictReader(open(path), delimiter="\t")
                         if "CIGAR" in r["methods"] and r["svtype"] in ("DEL", "INS")])
    real = sizes(f"{REAL}/sm_sv_calls.tsv"); dirty = sizes(f"{SYND}/sm_sv_calls.tsv")
    if not len(real) or not len(dirty):
        return
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.3))
    # panel 1: phase = |svlen| mod 178, density
    for d, c, lab in [(real, "#138D75", f"real CEN (n={len(real)})"),
                      (dirty, "#b30000", f"dirty artefact (n={len(dirty)})")]:
        a1.hist(d % 178, bins=np.arange(0, 179, 8), density=True, histtype="step", lw=2, color=c, label=lab)
    a1.axvline(0, ls=":", c="grey"); a1.axvline(178, ls=":", c="grey")
    a1.set_xlabel("|svlen| mod 178 bp  (0 / 178 = whole-monomer)"); a1.set_ylabel("density")
    a1.set_title("CEN178 phase of CIGAR DEL/INS"); a1.legend(fontsize=8)
    # panel 2: size comb 0..1500, density; mark 178 multiples
    for d, c, lab in [(real, "#138D75", "real CEN"), (dirty, "#b30000", "dirty artefact")]:
        a2.hist(d[d <= 1500], bins=np.arange(0, 1501, 18), density=True, histtype="step", lw=1.8, color=c, label=lab)
    for k in range(1, 9):
        a2.axvline(178 * k, ls=":", c="grey", lw=0.8)
    a2.set_xlabel("|svlen| (bp) — dotted lines = multiples of 178"); a2.set_ylabel("density")
    a2.set_title("CIGAR indel size spectrum"); a2.legend(fontsize=8)
    fig.suptitle("178-bp periodicity: real CEN CIGAR indels are whole-monomer (90%); "
                 "sequencing-error CIGAR indels are random-sized", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(f"{REAL}/cigar_periodicity.png", dpi=130); plt.close(fig)
    print(f"periodicity figure -> {REAL}/cigar_periodicity.png")


def detector_breakdown():
    """Which detector (CIGAR inline-indel vs split-based) produced the calls, per run. Shows the dirty
    floor is ~100% CIGAR and the split detectors (source of DUP/INV/BND) are error-immune."""
    from collections import Counter

    def one(path):
        c = Counter()
        if not os.path.exists(path):
            return c, 0
        rows = list(csv.DictReader(open(path), delimiter="\t"))
        for r in rows:
            m = r["methods"]
            det = "CIGAR" if "CIGAR" in m else ("SPLITANDMAP" if "SPLITANDMAP" in m else "SPLITREAD")
            c[det] += 1
        return c, len(rows)

    runs = [("real", f"{REAL}/sm_sv_calls.tsv"), ("mapping_floor", f"{SYN}/sm_sv_calls.tsv"),
            ("seq_stress", f"{SYND}/sm_sv_calls.tsv")]
    with open(f"{REAL}/synthetic_detector_source.tsv", "w") as f:
        f.write("run\tCIGAR\tSPLITREAD\tSPLITANDMAP\ttotal\n")
        print(f"\n{'run':14}{'CIGAR':>8}{'SPLITREAD':>11}{'SPLITANDMAP':>13}{'total':>8}")
        for name, path in runs:
            c, n = one(path)
            f.write(f"{name}\t{c['CIGAR']}\t{c['SPLITREAD']}\t{c['SPLITANDMAP']}\t{n}\n")
            print(f"{name:14}{c['CIGAR']:8}{c['SPLITREAD']:11}{c['SPLITANDMAP']:13}{n:8}")


def figure(rc, sc, dc, ra, sa, da):
    """CEN + ARM real rate (bars) over the mapping floor (green line, =0). The sequencing stress-test
    magnitude overshoots real ~200x (glitch model harsher than real CCS), so it is NOT plotted to scale;
    instead each type is annotated by whether sequencing error can generate it at all
    ('e' = INS/DEL, error-generable; '0' = DUP/INV/BND, never produced under any error level)."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, numpy as np
    col = {"wt_leaf": "#2c7fb8", "cenh3ox_leaf": "#7fcdbb",
           "wt_pollen": "#d95f0e", "cenh3ox_pollen": "#fec44f"}
    lab = {"INREG": "in-register\nDEL/INS", "DEL": "DEL", "INS": "INS",
           "DUP": "DUP", "INV": "INV", "BND": "BND"}
    panels = [("CEN", ["INREG", "DEL", "INS", "DUP", "BND"], rc, sc, dc,
               "real events per Mb of CEN-mapped read seq"),
              ("ARM", ["DEL", "INS", "DUP", "INV", "BND"], ra, sa, da,
               "real events per million arm reads")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, (comp, cats, real_d, clean_d, dirty_d, ylab) in zip(axes, panels):
        x = np.arange(len(cats)); w = 0.2
        cflmax = max((clean_d.get(g, {}).get(t, 0) for g in GROUPS for t in cats), default=0)
        ymax = max((real_d.get(g, {}).get(t, 0) for g in GROUPS for t in cats), default=1)
        for i, g in enumerate(GROUPS):
            real = [real_d.get(g, {}).get(t, 0) for t in cats]
            ax.bar(x + (i - 1.5) * w, real, w, color=col[g],
                   label=g.replace("_", " ") if comp == "CEN" else None)
        ax.axhline(cflmax, ls="--", lw=1.3, color="green")
        ax.text(len(cats) - 0.4, cflmax + 0.02 * ymax, f"mapping floor = {cflmax:.3f}  ",
                color="green", va="bottom", ha="right", fontsize=8.5, fontweight="bold")
        # seq-stress generability tag under each type
        for xi, t in zip(x, cats):
            got = sum(dirty_d.get(g, {}).get(t, 0) for g in GROUPS) > 0
            ax.annotate("seq-err: yes" if got else "seq-err: 0",
                        (xi, 0), (xi, -0.13 * ymax), ha="center", fontsize=7.5,
                        color="#b30000" if got else "#0a7d0a", annotation_clip=False)
        ax.set_ylim(-0.18 * ymax, ymax * 1.12)
        ax.set_xticks(x); ax.set_xticklabels([lab[c] for c in cats])
        ax.set_ylabel(ylab); ax.set_title(f"{comp}")
    axes[0].legend(fontsize=7.5, ncol=2)
    fig.suptitle("Single-molecule SV rate vs synthetic controls — CEN and ARMS\n"
                 "green line = mapping floor (reads at real divergence, = 0 for all types); "
                 "'seq-err' tag = whether inflated HiFi error can generate that type at all",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(f"{REAL}/synthetic_control.png", dpi=130); plt.close(fig)
    print(f"figure -> {REAL}/synthetic_control.png")


if __name__ == "__main__":
    main()
