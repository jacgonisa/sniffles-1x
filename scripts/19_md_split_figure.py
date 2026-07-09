#!/usr/bin/env python3
"""Illustrate the MD-tag CUSUM split (split-and-map, 03) on a real read -> docs/md_split.png.
Uses the canonical example documented in ALGORITHM.md: the wt_leaf/col Chr1 DUP read whose
clean<->noisy frontier is ~2.9 kb (an earlier global-contrast method cut it at ~1.0 kb).
Run: python scripts/19_md_split_figure.py   (env nextflow_env)."""
import os, pysam, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import bam_path, CEN

DOCS = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/docs"
READ = "m84227_260320_163030_s4/162401934/ccs"
SAMPLE, HAP, CHROM = "wt_leaf", "col", "Chr1"
MIN_FRAG = 1000
WIN = 50  # rolling-window (bp) for the displayed mismatch rate


def profile(read):
    """per-query-position match/mismatch (1=mismatch) from MD via get_aligned_pairs(with_seq)."""
    qpos, mism = [], []
    for q, r, base in read.get_aligned_pairs(with_seq=True):
        if q is None or r is None or base is None:
            continue
        qpos.append(q)
        mism.append(1 if base.islower() else 0)  # with_seq lowercases mismatched ref bases
    return np.array(qpos), np.array(mism, float)


def main():
    lo, hi = CEN[HAP][CHROM]
    bam = pysam.AlignmentFile(bam_path(SAMPLE, HAP), "rb")
    read = next(r for r in bam.fetch(CHROM, lo, hi)
                if r.query_name == READ and not r.is_secondary and not r.is_supplementary)
    qpos, mism = profile(read)
    n = len(qpos)
    mean = mism.mean()
    cum = np.cumsum(mism - mean)
    # candidate cut positions must leave >= MIN_FRAG each side (same gate as best_split)
    ok = (qpos - qpos[0] >= MIN_FRAG) & (qpos[-1] - qpos >= MIN_FRAG)
    i = np.where(ok, np.abs(cum), -1).argmax()
    cut = qpos[i]
    # displayed rolling mismatch rate (bp bins) just for the eye
    kb = qpos / 1000.0
    roll = np.convolve(mism, np.ones(WIN) / WIN, mode="same")
    L, R = mism[:i].mean(), mism[i:].mean()

    fig, ax = plt.subplots(2, 1, figsize=(9, 5.2), sharex=True,
                           gridspec_kw={"height_ratios": [1, 1], "hspace": 0.12})
    ax[0].fill_between(kb, roll * 100, color="#C0392B", alpha=0.25, step="mid")
    ax[0].plot(kb, roll * 100, color="#C0392B", lw=0.8)
    ax[0].axvline(cut / 1000, color="#111", lw=1.6, ls="--")
    ax[0].set_ylabel("mismatch rate\n(%% in %d-bp window)" % WIN)
    ax[0].set_title("MD-tag split-and-map: cut the read at the clean↔noisy frontier (real read, wt_leaf/col Chr1 DUP)")
    yl = ax[0].get_ylim()[1] * 0.7
    lstr = (f"belongs elsewhere\n(out-of-phase)\n~{L*100:.1f}%" if L > R else f"maps cleanly\n~{L*100:.1f}%")
    rstr = (f"maps cleanly\n~{R*100:.1f}%" if L > R else f"belongs elsewhere\n(out-of-phase)\n~{R*100:.1f}%")
    ax[0].annotate(lstr, (kb[i] * 0.5, yl), ha="center", fontsize=9,
                   color="#C0392B" if L > R else "#1E8449")
    ax[0].annotate(rstr, ((kb[-1] + kb[i]) / 2, yl), ha="center", fontsize=9,
                   color="#1E8449" if L > R else "#C0392B")

    ax[1].plot(kb, cum, color="#2980B9", lw=1.4)
    ax[1].axvline(cut / 1000, color="#111", lw=1.6, ls="--")
    ax[1].plot(cut / 1000, cum[i], "o", color="#111", ms=6)
    ax[1].axhline(0, color="#bbb", lw=0.6)
    ax[1].set_ylabel("CUSUM\nΣ (mismatch − mean)")
    ax[1].set_xlabel("query position along the read (kb)")
    ax[1].annotate(f"argmax |CUSUM| = cut @ {cut/1000:.2f} kb",
                   (cut / 1000, cum[i]), (cut / 1000 + 0.4, cum[i]),
                   fontsize=9, color="#111",
                   arrowprops=dict(arrowstyle="->", color="#111"))
    fig.savefig(f"{DOCS}/md_split.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {DOCS}/md_split.png  (n={n} aligned bp, cut @ {cut} bp, L={L:.3f} R={R:.3f})")


if __name__ == "__main__":
    main()
