#!/usr/bin/env python3
"""Overlay CHARLA pollen crossover breakpoints on the RPP1/DM2 gene annotation (Col coords).
Gene track (all liftoff genes in the window; the 4 RPP1-cluster TNLs highlighted + labelled) +
crossover breakpoint markers (pollen, coloured by category; leaf shown if any). Small n — this is a
map, not a test.
-> results/rpp1/rpp1_crossover_overlay.png
Run: /home/jg2070/miniforge3/envs/nextflow_env/bin/python rpp1_crossover_overlay.py"""
import os, csv, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

ROOT = "/mnt/ssd-4tb/HIFI_NAMIL"
OUT = f"{ROOT}/single_molecule_sv/results/rpp1"
GFF = f"{ROOT}/01_genomes/Col-HiFi/Col-0.ragtag_scaffolds.fa_liftoff.edta.gff3"
WIN = ("Chr3", 19_100_000, 19_650_000)
CLUSTER_GENES = {"AT3G44400", "AT3G44480", "AT3G44630", "AT3G44670"}
LABEL = {"AT3G44480": "RPP1\n(AT3G44480)", "AT3G44400": "AT3G44400",
         "AT3G44630": "AT3G44630", "AT3G44670": "AT3G44670"}
CO = {"leaf": f"{ROOT}/01_f1leaf-wt/01-mask_0/output/10-call_recombination_sites/wt_leaf_threshold_50.crossovers.tsv",
      "pollen": f"{ROOT}/03_f1pollen-wt/01-mask_0/output/10-call_recombination_sites/wt_pollen_threshold_50.crossovers.tsv"}
CATCOL = {"4b-2": "#d62728", "4b-5": "#1f77b4"}


def genes():
    chrom, a, b = WIN
    g = []
    for ln in open(GFF):
        if ln.startswith("#"):
            continue
        c = ln.rstrip("\n").split("\t")
        if len(c) < 9 or c[0] != chrom or c[2] != "gene":
            continue
        s, e = int(c[3]), int(c[4])
        if e < a or s > b:
            continue
        m = re.search(r"Name=([^;]+)", c[8])
        g.append((s, e, c[6], m.group(1) if m else "?"))
    return g


def crossovers(tis):
    chrom, a, b = WIN
    out = []
    for r in csv.DictReader(open(CO[tis]), delimiter="\t"):
        if r["decision"] != "crossover" or r["col_crossover_chr"] != chrom:
            continue
        s, e = int(r["col_crossover_start"]), int(r["col_crossover_end"])
        if a <= (s + e) / 2 <= b:
            out.append((s, e, r["category"]))
    return out


def main():
    chrom, a, b = WIN
    gs = genes()
    pol, leaf = crossovers("pollen"), crossovers("leaf")
    fig, ax = plt.subplots(figsize=(14, 4.6))

    # gene track at y=0
    ci = 0
    for s, e, strand, name in sorted(gs):
        cluster = name in CLUSTER_GENES
        y = 0
        col = "#B8860B" if cluster else "#bbbbbb"
        dx = (e - s) if strand == "+" else -(e - s)
        x0 = s if strand == "+" else e
        ax.add_patch(FancyArrow(x0, y, dx, 0, width=0.18 if cluster else 0.06,
                                length_includes_head=True, head_width=0.32 if cluster else 0.10,
                                head_length=min(2500, abs(dx) * 0.4), color=col,
                                alpha=0.95 if cluster else 0.5, zorder=3 if cluster else 1))
        if cluster:
            ly = -0.5 - 0.42 * (ci % 2)          # stagger labels so neighbours don't collide
            ax.plot([(s + e) / 2, (s + e) / 2], [-0.12, ly + 0.12], color="#7a5c00", lw=0.5, zorder=2)
            ax.text((s + e) / 2, ly, LABEL[name], ha="center", va="top", fontsize=8.5,
                    fontweight="bold", color="#7a5c00")
            ci += 1

    # crossover markers
    def draw(cos, y, marker, lab):
        for i, (s, e, cat) in enumerate(sorted(cos)):
            mid = (s + e) / 2
            ax.plot([s, e], [y, y], color=CATCOL.get(cat, "grey"), lw=4, solid_capstyle="butt", zorder=5)
            ax.plot(mid, y, marker, color=CATCOL.get(cat, "grey"), ms=9, zorder=6)
            ax.vlines(mid, 0.2, y, color=CATCOL.get(cat, "grey"), lw=0.8, ls=":", zorder=2)
    draw(pol, 1.6, "v", "pollen")
    if leaf:
        draw(leaf, 2.2, "s", "leaf")

    ax.axhline(0.2, color="black", lw=0.5)
    ax.set_ylim(-1.5, 2.6); ax.set_yticks([0, 1.6])
    ax.set_yticklabels(["genes", f"pollen COs (n={len(pol)})"])
    ax.set_xlim(a, b); ax.set_xlabel(f"{chrom} position (Col coords, bp)")
    ax.set_title(f"CHARLA pollen crossover breakpoints over the RPP1/DM2 cluster  "
                 f"{chrom}:{a:,}-{b:,}\n(gold = RPP1-like TNL array; leaf had 0 crossovers here; "
                 f"red=4b-2, blue=4b-5)")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], color=CATCOL["4b-2"], lw=4, label="crossover 4b-2"),
                       Line2D([0], [0], color=CATCOL["4b-5"], lw=4, label="crossover 4b-5")],
              fontsize=8, loc="upper left")
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(f"{OUT}/rpp1_crossover_overlay.png", dpi=140); plt.close(fig)
    print(f"{len(gs)} genes, {len(pol)} pollen + {len(leaf)} leaf crossovers in window")
    print(f"figure -> {OUT}/rpp1_crossover_overlay.png")


if __name__ == "__main__":
    main()
