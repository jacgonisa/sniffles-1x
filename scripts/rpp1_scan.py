#!/usr/bin/env python3
"""RPP1 / DM2 cluster de-novo SV scan (chr3 TNL array; Ian's request).

Single-molecule SV detection (same detectors as the pipeline: CIGAR I/D>=50 + SA split-read
via 02_leadprov_sm) restricted to the RPP1/DM2 cluster window on Chr3, comparing POLLEN (probe
for de-novo events) vs LEAF (deep somatic baseline / paralog-artefact control), Col haplotype
(exact gene coords + RPP1=At3g44480 is a Col gene).

De-novo candidate = a large cluster SV carried by pollen reads that is NOT recurrent in the deep
leaf baseline (low VAF, not a fixed assembly/paralog difference). Repetitive TNL cluster => expect
paralog-mismap artefacts (sharp gene-edge boundaries, recurrent in BOTH tissues) — the leaf control
is what separates those from real de-novo change.

-> results/rpp1/rpp1_events.tsv (all SV leads) + rpp1_summary printed + rpp1_cluster.png
Run: /home/jg2070/miniforge3/envs/nextflow_env/bin/python rpp1_scan.py"""
import os, csv, importlib
from collections import defaultdict, Counter
import pysam
lp = importlib.import_module("02_leadprov_sm")

ROOT = "/mnt/ssd-4tb/HIFI_NAMIL"
OUT = f"{ROOT}/single_molecule_sv/results/rpp1"
BIG = 1000                       # "structural" event size for the de-novo catalogue
WIN = ("Chr3", 19_150_000, 19_460_000)
GENES = {"AT3G44400": (19_179_993, 19_185_327), "AT3G44630": (19_378_961, 19_384_613),
         "AT3G44670": (19_403_010, 19_407_575), "RPP1/AT3G44480": (19_432_477, 19_438_182)}
# col haplotype only (Col-HiFi carries the annotation); pollen = probe, leaf = baseline
DATASETS = [("wt_pollen", "col", "pollen"), ("wt_leaf", "col", "leaf")]


def bam_path(sample, hap):
    return f"{ROOT}/sv_calling/aligned/{sample}/strict90/{hap}_all.bam"


def nearest_gene(pos):
    for g, (a, b) in GENES.items():
        if a - 2000 <= pos <= b + 2000:
            return g
    return min(GENES, key=lambda g: min(abs(pos - GENES[g][0]), abs(pos - GENES[g][1]))) + "~"


def scan():
    chrom, a, b = WIN
    events = []                              # (tissue, pos, svtype, svlen, read, source, mapq)
    reads_at = defaultdict(int)              # tissue -> spanning reads (for VAF/normalisation)
    for sample, hap, tis in DATASETS:
        bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
        for r in bam.fetch(chrom, a, b):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < lp.MAPQ_MIN:
                continue
            if not (a <= r.reference_start < b):
                continue
            reads_at[tis] += 1
            for svtype, pos, svlen in lp.cigar_leads(r):
                if a <= pos < b:
                    events.append((tis, pos, svtype, svlen, r.query_name, "CIGAR", r.mapping_quality))
            for svtype, pos, svlen, mapq, mate in lp.split_leads(r, chrom):
                if mapq >= lp.MAPQ_MIN and a <= pos < b:
                    events.append((tis, pos, svtype, svlen or 0, r.query_name, "SPLIT", mapq))
        bam.close()
    return events, reads_at


def main():
    os.makedirs(OUT, exist_ok=True)
    events, reads_at = scan()
    with open(f"{OUT}/rpp1_events.tsv", "w") as f:
        f.write("tissue\tpos\tsvtype\tsvlen\tsource\tmapq\tnearest_gene\tread\n")
        for tis, pos, svt, svlen, read, src, mapq in sorted(events, key=lambda e: (e[0], e[1])):
            f.write(f"{tis}\t{pos}\t{svt}\t{svlen}\t{src}\t{mapq}\t{nearest_gene(pos)}\t{read}\n")

    print(f"=== RPP1/DM2 cluster {WIN[0]}:{WIN[1]:,}-{WIN[2]:,}  (Col haplotype) ===")
    print(f"spanning reads: pollen={reads_at['pollen']}  leaf={reads_at['leaf']}\n")
    print(f"{'tissue':8}{'type':6}{'all(>=50bp)':>12}{'big(>=1kb)':>11}{'per1k_reads_big':>16}")
    for tis in ("pollen", "leaf"):
        te = [e for e in events if e[0] == tis]
        for svt in ("DEL", "INS", "DUP", "INV", "BND"):
            allc = sum(1 for e in te if e[2] == svt)
            big = sum(1 for e in te if e[2] == svt and abs(e[3] or 0) >= BIG)
            rate = 1000 * big / reads_at[tis] if reads_at[tis] else 0
            if allc:
                print(f"{tis:8}{svt:6}{allc:12}{big:11}{rate:16.2f}")

    # de-novo candidate catalogue: big events, VAF, recurrence in leaf
    print(f"\n=== big (>= {BIG} bp) events — de-novo candidates (pollen, not recurrent in leaf) ===")
    leaf_pos = [e[1] for e in events if e[0] == "leaf" and abs(e[3] or 0) >= BIG]
    big_poll = [e for e in events if e[0] == "pollen" and abs(e[3] or 0) >= BIG]
    print(f"{'pos':>10}{'type':>5}{'svlen':>8}{'src':>7}{'gene':>18}{'leaf_within2kb':>15}")
    for tis, pos, svt, svlen, read, src, mapq in sorted(big_poll, key=lambda e: e[1]):
        nleaf = sum(1 for p in leaf_pos if abs(p - pos) <= 2000)
        flag = "" if nleaf else "  <- DE-NOVO?"
        print(f"{pos:>10}{svt:>5}{svlen:>8}{src:>7}{nearest_gene(pos):>18}{nleaf:>15}{flag}")

    plot(events, reads_at)
    print("\nDONE_RPP1")


def plot(events, reads_at):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, numpy as np
    chrom, a, b = WIN
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    for ax, tis, col in [(axes[0], "pollen", "#d95f0e"), (axes[1], "leaf", "#2c7fb8")]:
        te = [e for e in events if e[0] == tis]
        # all >=50bp as small ticks, big >=1kb as tall lines coloured by type
        tcol = {"DEL": "#c0392b", "INS": "#27ae60", "DUP": "#8e44ad", "INV": "#e67e22", "BND": "#16a085"}
        for _, pos, svt, svlen, read, src, mapq in te:
            big = abs(svlen or 0) >= BIG
            ax.axvline(pos, ymin=0, ymax=(0.9 if big else 0.3), lw=(1.4 if big else 0.4),
                       color=tcol.get(svt, "grey"), alpha=(0.9 if big else 0.25))
        for g, (ga, gb) in GENES.items():
            ax.axvspan(ga, gb, color="gold", alpha=0.35)
            ax.text((ga + gb) / 2, 1.02, g.split("/")[0], ha="center", va="bottom", fontsize=7, rotation=90)
        ax.set_ylim(0, 1.15); ax.set_yticks([])
        ax.set_ylabel(f"{tis}\n(n={reads_at[tis]} reads)")
    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], color=c, lw=2, label=t) for t, c in
           {"DEL": "#c0392b", "INS": "#27ae60", "DUP": "#8e44ad", "INV": "#e67e22", "BND": "#16a085"}.items()]
    axes[0].legend(handles=leg, ncol=5, fontsize=8, loc="upper right")
    axes[0].set_title(f"RPP1/DM2 cluster single-molecule SVs — pollen vs leaf (Col hap)  "
                      f"{chrom}:{a:,}-{b:,}\ntall lines = big (>={BIG}bp) events, short = >=50bp; gold = genes")
    axes[1].set_xlabel(f"{chrom} position (bp)")
    fig.tight_layout()
    fig.savefig(f"{OUT}/rpp1_cluster.png", dpi=130); plt.close(fig)
    print(f"figure -> {OUT}/rpp1_cluster.png")


if __name__ == "__main__":
    main()
