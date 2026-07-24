#!/usr/bin/env python3
"""Two simple slide schematics: in-register vs out-of-register satellite events.
CEN178 monomers = head-to-tail ARROWS. A whole-monomer deletion keeps the array in phase
(in-register, unequal-sister-chromatid-HR signature); a partial-monomer deletion breaks a monomer
at the junction and shifts the phase (out-of-register, NHEJ-like).
-> docs/register_schematic.png"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle

DOCS = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/docs"
GREEN, RED, GREY = "#2E7D32", "#C0392B", "#BDBDBD"
W = 0.86           # monomer/arrow body length (of a 1.0 pitch)
HL = 0.30          # arrow head length


def monomer(ax, x, y, color=GREEN, alpha=1.0):
    ax.add_patch(FancyArrow(x + 0.05, y, W, 0, width=0.34, head_width=0.5, head_length=HL,
                            length_includes_head=True, fc=color, ec="none", alpha=alpha, zorder=3))


def half(ax, x, y, side, color=RED):
    """half a monomer: 'L' = tail rectangle (no head), 'R' = head-only stub."""
    if side == "L":
        ax.add_patch(Rectangle((x + 0.05, y - 0.17), 0.45, 0.34, fc=color, ec="none", zorder=3))
    else:
        ax.add_patch(FancyArrow(x + 0.42, y, 0.44, 0, width=0.34, head_width=0.5, head_length=HL,
                                length_includes_head=True, fc=color, ec="none", zorder=3))


def row(ax, x0, y, n, color=GREEN, alpha=1.0):
    for i in range(n):
        monomer(ax, x0 + i, y, color, alpha)


def panel(ax, mode):
    ax.set_xlim(-1.6, 9.6); ax.set_ylim(-1.4, 5.4); ax.axis("off")
    tc = GREEN if mode == "in" else RED
    title = "IN-REGISTER" if mode == "in" else "OUT-OF-REGISTER"
    ax.text(4.0, 5.05, title, ha="center", fontsize=16, fontweight="bold", color=tc)

    # reference array (top): 8 monomers, deleted ones greyed
    ax.text(-1.4, 3.7, "reference\narray", fontsize=9.5, va="center", ha="left", color="#555")
    for i in range(8):
        monomer(ax, i, 3.7, GREY if i in (3, 4) else GREEN, alpha=0.5 if i in (3, 4) else 1.0)
    dx0, dx1 = (3.0, 5.0) if mode == "in" else (3.0, 4.55)
    ax.plot([dx0 + 0.05, dx1 - 0.05], [4.25, 4.25], color=RED, lw=1.4)
    ax.text((dx0 + dx1) / 2, 4.35, "deleted", ha="center", va="bottom", fontsize=8.5, color=RED)

    # middle: the size rule + down arrow
    rule = ("Δ = whole number of monomers   (|Δ| mod 178 ≈ 0)" if mode == "in"
            else "Δ = partial monomer   (|Δ| mod 178 ≠ 0)")
    ax.text(4.0, 2.55, rule, ha="center", fontsize=10, color=tc, fontweight="bold")
    ax.annotate("", xy=(4.0, 1.75), xytext=(4.0, 3.2), arrowprops=dict(arrowstyle="-|>", color="#888", lw=1.5))

    # read (bottom): the two flanks rejoined
    ax.text(-1.4, 1.2, "read\n(deletion)", fontsize=9.5, va="center", ha="left", color="#555")
    if mode == "in":
        row(ax, 0, 1.2, 6, GREEN)                    # 6 whole monomers, all complete, in phase
        ax.axvline(3.0, ymin=0.30, ymax=0.50, color="#333", lw=1.3, ls="--")
        ax.text(3.0, 0.55, "junction at a monomer boundary", ha="center", va="top", fontsize=8.5, color="#333")
        sub = "flanking monomers still tile perfectly → the array stays in phase.\nSignature of unequal sister-chromatid HR."
    else:
        row(ax, 0, 1.2, 3, GREEN)                    # left flank: 3 whole
        half(ax, 3, 1.2, "L", RED); half(ax, 3, 1.2, "R", RED)   # chimeric broken monomer at junction
        for i in range(4, 7):                        # right flank: shifted (out of phase)
            monomer(ax, i + 0.4, 1.2, GREEN, alpha=0.9)
        ax.axvline(3.5, ymin=0.30, ymax=0.50, color=RED, lw=1.3, ls="--")
        ax.text(3.5, 0.55, "junction mid-monomer → broken monomer, phase shifted", ha="center", va="top", fontsize=8.5, color=RED)
        sub = "a monomer is cut at the junction and the downstream array is out of phase.\nNHEJ-like / non-homologous junction."
    ax.text(4.0, -0.7, sub, ha="center", va="center", fontsize=8.6, color="#444")


def main():
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.2))
    panel(axes[0], "in"); panel(axes[1], "out")
    fig.suptitle("CEN178 satellite monomers = head-to-tail arrows (→).  What a deletion does to the register:",
                 fontsize=11, fontweight="bold", y=1.0)
    fig.savefig(f"{DOCS}/register_schematic.png", dpi=160, bbox_inches="tight")
    print(f"wrote {DOCS}/register_schematic.png")


if __name__ == "__main__":
    main()
