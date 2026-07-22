#!/usr/bin/env python3
"""Exhaustive taxonomy of single-molecule SV / recombination events in a Col×Ler (or MAT×PAT) F1,
and where each class is captured: our non-hybrid SV pipeline (green), CHARLA hybrid/crossover
analysis (blue), a current GAP (red), or a not-a-real-event artefact (grey).
-> docs/mechanism_taxonomy.png"""
import textwrap
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

DOCS = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/docs"
COLR = {"ours": ("#E8F6EF", "#1E8449"), "charla": ("#EAF2FB", "#2471A3"),
        "gap": ("#FDEDEC", "#C0392B"), "art": ("#F2F3F4", "#7F8C8D")}


def box(ax, x, y, w, h, title, body, kind, fs=8.0):
    fc, ec = COLR[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
                                fc=fc, ec=ec, lw=1.5, zorder=2))
    ax.text(x + w / 2, y + h - 0.018, title, ha="center", va="top", fontsize=fs + 0.6,
            fontweight="bold", color=ec, zorder=3)
    ax.text(x + w / 2, y + h - 0.058, "\n".join(textwrap.wrap(body, 42)), ha="center", va="top",
            fontsize=fs - 0.3, color="#222", zorder=3)


def header(ax, x, w, y, text, kind):
    fc, ec = COLR[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, 0.05, boxstyle="round,pad=0.004,rounding_size=0.01",
                                fc=ec, ec="none", zorder=2))
    ax.text(x + w / 2, y + 0.025, text, ha="center", va="center", fontsize=10.5,
            fontweight="bold", color="white", zorder=3)


def main():
    fig, ax = plt.subplots(figsize=(18, 10.5))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.985, "Single-molecule SV / recombination in the Col×Ler F1 — exhaustive mechanism taxonomy",
            ha="center", fontsize=15, fontweight="bold")
    ax.text(0.5, 0.958, "an event joins two genomic loci in one read.  Axis 1 = partner molecule (self/sister · homolog · non-homolog).  "
            "Axis 2 = allelic (same locus) vs non-allelic (offset).  Overlay = compartment (CEN178 satellite · pericentromere · unique arm).",
            ha="center", fontsize=9.5, color="#444", style="italic")

    cols = [0.015, 0.265, 0.515, 0.765]; W = 0.225
    header(ax, cols[0], W, 0.895, "SELF / SISTER  (same haplotype — non-hybrid read)", "ours")
    header(ax, cols[1], W, 0.895, "INTER-HOMOLOG  (Col↔Ler switch — hybrid read)", "charla")
    header(ax, cols[2], W, 0.895, "GAP  (currently missed)", "gap")
    header(ax, cols[3], W, 0.895, "NOT A REAL EVENT  (artefact)", "art")

    # each block = (title, body, kind, HEIGHT); stacked from the top down
    sb = [
        ("1. CEN, non-allelic — UNEQUAL SISTER-CHROMATID EXCHANGE",
         "CEN178 satellite tandem dup/del/expansion (in-register, whole-monomer DEL/INS/DUP). "
         "Donor LOCAL (adjacent monomers) or DISTAL (across the array). intra-chromatid vs inter-sister are "
         "sequence-identical → NOT separable.", "ours", 0.27),
        ("2. ARM, non-allelic — unequal exchange between arm repeats",
         "NAHR between paralogous/inverted arm repeats → INV (inverted repeats), rare DEL/INS. Same haplotype, "
         "same chromosome, offset.", "ours", 0.17),
        ("3. any, ALLELIC — equal exchange",
         "no net gain/loss → NO SV produced (invisible to any caller).", "gap", 0.09),
        ("4. NON-HOMOLOG — ECTOPIC NAHR (same haplotype)",
         "junction to a DIFFERENT chromosome → BND / translocation. Must exclude shared-satellite cross-mapping "
         "(artefacts →).", "ours", 0.15),
    ]
    cb = [
        ("5. same chr, ALLELIC — CROSSOVER / GENE CONVERSION",
         "classic allelic meiotic HR: reciprocal crossover or non-crossover gene conversion. CHARLA hybrid read, "
         "profile Px→Mx at the homologous position.", "charla", 0.26),
        ("6. same chr, NON-ALLELIC — inter-homolog unequal (NAHR)",
         "offset exchange between the two homologs (CEN–CEN or ARM–ARM at non-matching positions). CHARLA hybrid "
         "read with a position offset across the switch.", "charla", 0.19),
        ("7. DIFFERENT chr — ectopic inter-homolog",
         "translocation carrying a haplotype switch (Col-chr1 ↔ Ler-chr3). CHARLA hybrid read, cross-chromosome "
         "profile Px→My.", "charla", 0.15),
    ]
    gb = [
        ("8. HAPLOTYPE-SWITCH inside an INSERTION",
         "the inserted tract is templated from the OTHER homolog (a Ler segment dropped into a Col read). The read "
         "passes strict-90 (still mostly one hap) so the SV is called but its donor haplotype is never checked, and "
         "the sub-segment switch is below CHARLA's block threshold → currently mis-labelled self/sister. "
         "FIX: Col/Ler k-mer check on the inserted bases themselves.", "gap", 0.33),
        ("9. inter-homolog SV on a strict-90-DISCARDED read",
         "a read spanning a true Col↔Ler junction is ~50/50 → fails the strict-90 purity filter and is removed "
         "before SV calling. FIX: also call SVs on the relaxed / un-split reads.", "gap", 0.19),
    ]
    ab = [
        ("10. shared-satellite cross-mapping → false BND",
         "a fragment mis-maps to another centromere (all CEN share CEN178) → inter-CEN BND, not a real junction. "
         "Empirical noise floor; ~11× worse in pollen (shorter reads).", "art", 0.22),
        ("11. homopolymer / CCS quality-decay → false INS",
         "inserted bases are a homopolymer tract or the CCS quality collapses inside the insertion. ~99% of arm / "
         "~80% genome-wide human insertions. Flagged + removed by the insertion QC.", "art", 0.20),
        ("12. satellite mapping ambiguity → false split",
         "an accurate read placed at the wrong monomer copy → spurious split/indel. Residual caveat in deep satellite.",
         "art", 0.15),
    ]

    def place(col, blocks):
        y = 0.87
        for title, body, kind, h in blocks:
            box(ax, col, y - h, W, h, title, body, kind)
            y -= h + 0.013
    place(cols[0], sb); place(cols[1], cb); place(cols[2], gb); place(cols[3], ab)

    # legend
    lx = 0.015
    for i, (k, lab) in enumerate([("ours", "captured by our non-hybrid SV pipeline"),
                                  ("charla", "captured by CHARLA hybrid / crossover analysis"),
                                  ("gap", "current GAP — not yet captured"),
                                  ("art", "not a real event (artefact / invisible)")]):
        fc, ec = COLR[k]
        ax.add_patch(FancyBboxPatch((lx + i * 0.245, 0.005), 0.02, 0.02, boxstyle="round,pad=0.002",
                                    fc=fc, ec=ec, lw=1.5))
        ax.text(lx + i * 0.245 + 0.028, 0.015, lab, va="center", fontsize=8.5, color=ec, fontweight="bold")

    fig.savefig(f"{DOCS}/mechanism_taxonomy.png", dpi=145, bbox_inches="tight")
    print(f"wrote {DOCS}/mechanism_taxonomy.png")


if __name__ == "__main__":
    main()
