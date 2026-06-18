#!/usr/bin/env python3
"""Generate explanatory diagrams for the repo (docs/):
  algorithm_flow.png  — how a read becomes a single-molecule SV call (and what is skipped)
  topology.png        — how sv.classify_splits reads the Part1/Part2 mapping topology
Run with any python that has matplotlib."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

DOCS = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/docs"
os.makedirs(DOCS, exist_ok=True)
BLUE, RED, GREEN, GREY = "#2C6FBB", "#C0392B", "#27AE60", "#777"


def box(ax, x, y, w, h, text, fc="#EFF4FB", ec=BLUE, fs=9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.04",
                                fc=fc, ec=ec, lw=1.4))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, wrap=True)


def arrow(ax, x1, y1, x2, y2, color="#333"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.6, color=color))


def algorithm_flow():
    fig, ax = plt.subplots(figsize=(12, 6.2)); ax.set_xlim(0, 12); ax.set_ylim(0, 6.2); ax.axis("off")
    box(ax, 0.3, 2.6, 2.3, 1.0, "Candidate read\nde≥0.005  |  NM≥50  |  SA tag", "#FEF9E7", "#C8A415", 9)
    # two detector branches
    box(ax, 3.4, 4.0, 3.2, 1.6,
        "leadprov  (Sniffles2)\n• CIGAR  I/D ≥ 50 bp\n• SA split-read alignments", "#EFF4FB", BLUE, 9)
    box(ax, 3.4, 0.7, 3.2, 1.7,
        "split-and-map  (this repo)\n• MD substitution-contrast\n  split point\n• winnowmap re-map 2 frags", "#EAF7EF", GREEN, 9)
    arrow(ax, 2.6, 3.3, 3.4, 4.8); arrow(ax, 2.6, 2.9, 3.4, 1.6)
    # shared classifier
    box(ax, 7.4, 2.6, 2.5, 1.0, "sv.classify_splits\n(Sniffles2 — UNMODIFIED)", "#EFF4FB", BLUE, 9)
    arrow(ax, 6.6, 4.8, 7.6, 3.6); arrow(ax, 6.6, 1.55, 7.6, 2.6)
    box(ax, 10.1, 2.6, 1.7, 1.0, "DEL  DUP\nINV  INS  BND", "#fff", "#333", 9)
    arrow(ax, 9.9, 3.1, 10.1, 3.1)
    box(ax, 10.0, 0.9, 1.9, 1.0, "1 lead =\n1 single-molecule SV", "#FDEDEC", RED, 9)
    arrow(ax, 10.95, 2.6, 10.95, 1.9)
    # skipped
    ax.add_patch(Rectangle((7.3, 4.3), 4.5, 1.2, fc="#f2f2f2", ec=GREY, ls="--", lw=1.2))
    ax.text(9.55, 4.9, "Sniffles stages 2–3:  cluster.py  +  coverage/VAF QC\nSKIPPED — the only steps that need ≥2 reads",
            ha="center", va="center", fontsize=8.5, color=GREY, style="italic")
    ax.text(6, 6.0, "Single-molecule SV calling = Sniffles' per-read classifier, run without its ≥2-read clustering",
            ha="center", fontsize=11, weight="bold")
    fig.savefig(f"{DOCS}/algorithm_flow.png", dpi=140, bbox_inches="tight"); plt.close(fig)


def frag(ax, x, y, w, color, rev=False, label=""):
    ax.add_patch(Rectangle((x, y), w, 0.22, fc=color, ec="none"))
    dx = -0.12 if rev else 0.12
    xa = (x + w - 0.02) if not rev else (x + 0.02)
    ax.annotate("", xy=(xa + dx, y + 0.11), xytext=(xa, y + 0.11),
                arrowprops=dict(arrowstyle="-|>", color="white", lw=1.4))
    if label:
        ax.text(x + w / 2, y + 0.32, label, ha="center", fontsize=7.5)


def topology():
    # rows: DEL, DUP, INV, INS — each shows read Part1/Part2 mapped to a reference line + the rule
    rows = [
        ("DEL", RED, "same strand · big REF gap between Part1→ and Part2→ (query contiguous)"),
        ("DUP", "#E8A20C", "same strand · Part2 maps BEFORE Part1 ends on the reference (overlap)"),
        ("INV", GREEN, "Part2 maps on the OPPOSITE strand to Part1"),
        ("INS", BLUE, "same strand · big QUERY gap, Part1/Part2 adjacent on the reference"),
    ]
    fig, axes = plt.subplots(len(rows), 1, figsize=(12, 7.6))
    for ax, (name, col, rule) in zip(axes, rows):
        ax.set_xlim(0, 12); ax.set_ylim(-0.55, 1.05); ax.axis("off")
        ax.plot([1.1, 7.2], [0.05, 0.05], color="#bbb", lw=2)   # reference (drawing zone x<=7.2)
        ax.text(0.05, 0.45, name, fontsize=13, weight="bold", color=col, va="center")
        if name == "DEL":
            frag(ax, 1.3, 0.45, 1.8, col, label="Part1 →"); frag(ax, 5.2, 0.45, 1.8, col, label="Part2 →")
            ax.plot([3.1, 5.2], [0.05, 0.05], color=col, ls=":", lw=1.6)
        elif name == "DUP":
            frag(ax, 3.6, 0.55, 1.9, col, label="Part2 →"); frag(ax, 2.5, 0.22, 1.9, col, label="Part1 →")
        elif name == "INV":
            frag(ax, 1.6, 0.45, 1.8, col, label="Part1 →"); frag(ax, 4.4, 0.45, 1.8, col, rev=True, label="← Part2")
        elif name == "INS":
            frag(ax, 2.2, 0.45, 1.7, col, label="Part1 →"); frag(ax, 3.95, 0.45, 1.7, col, label="Part2 →")
            ax.annotate("", xy=(3.95, 0.78), xytext=(3.88, 0.78), arrowprops=dict(arrowstyle="<->", color="#333"))
            ax.text(3.05, 0.92, "inserted (not in ref)", ha="center", fontsize=7.5, color="#333")
        ax.text(7.7, 0.35, rule, ha="left", fontsize=9, color="#444", va="center")
    fig.suptitle("How sv.classify_splits infers SV type from Part1/Part2 mapping topology", fontsize=12, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(f"{DOCS}/topology.png", dpi=140, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    algorithm_flow(); topology()
    print(f"wrote {DOCS}/algorithm_flow.png and topology.png")
