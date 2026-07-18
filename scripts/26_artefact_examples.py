#!/usr/bin/env python3
"""Step 26 — example panels of the ARTEFACTUAL insertions the QC flags (homopolymer tracts and
CCS-quality collapse). For the current genome profile (SM_GENOME), picks the worst flagged INS from
insertion_qc.tsv and, per example, draws:
  top   — per-base CCS quality along the read around the insertion (insertion shaded; flank vs
          insertion mean-Q lines) -> shows the quality DECAY inside the insertion.
  bottom— the inserted sequence as a coloured base strip with homopolymer runs (>=5 bp) underlined
          -> shows the low-complexity / homopolymer tract.
-> results[_human]/artefact_examples/<organism>_<mode>_<...>.png
Run: python 26_artefact_examples.py   and   SM_GENOME=human python 26_artefact_examples.py"""
import os, csv
import pysam
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from common import OUT, bam_path, GENOME

ORG = "human" if GENOME == "human" else "arabidopsis"
OUTDIR = f"{OUT}/artefact_examples"; os.makedirs(OUTDIR, exist_ok=True)
PAD = 400
BASECOL = {"A": "#2E7D32", "C": "#1565C0", "G": "#F9A825", "T": "#C62828", "N": "#999"}
FLANK = 200


def find_ins(cigar, target):
    q = 0
    for op, ln in cigar:
        if op in (0, 7, 8): q += ln
        elif op == 1:
            if abs(ln - target) <= 3: return q, q + ln
            q += ln
        elif op == 4: q += ln
    return None, None


def hp_runs(s):
    runs = []; i = 0; n = len(s)
    while i < n:
        j = i
        while j + 1 < n and s[j + 1] == s[i]: j += 1
        if j - i + 1 >= 5: runs.append((i, j + 1, s[i]))
        i = j + 1
    return runs


def read_lookup():
    """(sample,hap,chrom,pos) -> read name, from sm_sv_calls.tsv INS rows."""
    d = {}
    for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t"):
        if r["svtype"] == "INS":
            d[(r["sample"], r["hap"], r["chrom"], r["pos"])] = r["read"]
    return d


def pick(path, n_each=2):
    rl = read_lookup()
    rows = [r for r in csv.DictReader(open(path), delimiter="\t") if r["verdict"] == "FLAG"]
    for r in rows:
        r["read"] = rl.get((r["sample"], r["hap"], r["chrom"], r["pos"]), "")
    rows = [r for r in rows if r["read"]]
    decay = sorted((r for r in rows if r["quality_decay"] == "1"), key=lambda r: -float(r["q_contrast"]))[:n_each]
    homo = sorted((r for r in rows if r["low_complexity"] == "1"), key=lambda r: -float(r["hp_frac"]))[:n_each]
    out = []
    seen = set()
    for mode, sel in (("quality_decay", decay), ("homopolymer", homo)):
        for r in sel:
            k = (r["sample"], r["hap"], r["chrom"], r["pos"])
            if k in seen: continue
            seen.add(k); out.append((mode, r))
    return out


def draw(mode, r):
    s, h, chrom, pos, sz = r["sample"], r["hap"], r["chrom"], int(r["pos"]), int(r["ins_bp"])
    bam = pysam.AlignmentFile(bam_path(s, h), "rb")
    rec = None
    for rd in bam.fetch(chrom, max(0, pos - 50000), pos + 50000):
        if rd.query_name != r["read"] or rd.is_unmapped or rd.is_secondary or rd.is_supplementary or not rd.cigartuples:
            continue
        if rd.query_sequence is None: continue
        qs, qe = find_ins(rd.cigartuples, sz)
        if qs is None: continue
        rec = (rd.query_sequence, rd.query_qualities, qs, qe); break
    bam.close()
    if rec is None: return None
    seq, qual, qs, qe = rec
    ins = seq[qs:qe]

    fig, (axq, axb) = plt.subplots(2, 1, figsize=(10, 4.2), gridspec_kw={"height_ratios": [2, 1], "hspace": 0.35})
    # --- CCS quality track ---
    lo, hi = max(0, qs - PAD), min(len(seq), qe + PAD)
    xs = range(lo, hi)
    if qual is not None:
        axq.scatter(list(xs), list(qual[lo:hi]), s=4, c="#263238", alpha=0.5, linewidths=0)
        qi = sum(qual[qs:qe]) / max(qe - qs, 1)
        fl = list(qual[max(0, qs - FLANK):qs]) + list(qual[qe:qe + FLANK]); qf = sum(fl) / max(len(fl), 1)
        axq.axhline(qf, color="#1B5E20", ls="--", lw=1); axq.text(lo, qf + 0.5, f"flank mean Q{qf:.0f}", color="#1B5E20", fontsize=8)
        axq.axhline(qi, color="#C0392B", ls="--", lw=1); axq.text(qe, qi + 0.5, f"insertion mean Q{qi:.0f}", color="#C0392B", fontsize=8)
    axq.axvspan(qs, qe, color="#C0392B", alpha=0.10)
    axq.axvline(qs, color="#C0392B", lw=1, ls=":"); axq.axvline(qe, color="#C0392B", lw=1, ls=":")
    axq.set_ylabel("CCS base quality"); axq.set_xlim(lo, hi); axq.set_ylim(0, max(45, (max(qual[lo:hi]) if qual is not None else 40) + 3))
    axq.set_xlabel("read position (bp)")
    lab = {"quality_decay": "QUALITY DECAY", "homopolymer": "HOMOPOLYMER / low-complexity"}[mode]
    axq.set_title(f"{ORG} · {s} {h} · {chrom}:{pos:,} · {sz} bp insertion — FLAGGED: {lab}\n"
                  f"homopolymer frac {float(r['hp_frac'])*100:.0f}%  ·  entropy {r['entropy']} bits  ·  "
                  f"Q contrast {float(r['q_contrast']):.0f}", fontsize=9)
    # --- inserted-sequence base strip (subsample if long) ---
    disp = ins if len(ins) <= 300 else ins[:300]
    axb.set_xlim(0, len(disp)); axb.set_ylim(0, 1); axb.set_yticks([])
    for i, b in enumerate(disp):
        axb.add_patch(Rectangle((i, 0.25), 1, 0.5, fc=BASECOL.get(b, "#999"), ec="none"))
    for a, bnd, base in hp_runs(disp):
        axb.plot([a, bnd], [0.15, 0.15], color="black", lw=2.5)
    axb.set_xlabel(f"inserted sequence ({'first 300 of ' + str(len(ins)) + ' bp' if len(ins) > 300 else str(len(ins)) + ' bp'})  ·  black bars = homopolymer runs ≥5 bp")
    hnd = [plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=BASECOL[b], ms=8, label=b) for b in "ACGT"]
    axb.legend(handles=hnd, ncol=4, fontsize=7, loc="upper right", framealpha=0.8)
    p = f"{OUTDIR}/{ORG}_{mode}_{s}_{h}_{chrom}_{pos}.png"
    fig.savefig(p, dpi=140, bbox_inches="tight"); plt.close(fig)
    return p


def main():
    picks = pick(f"{OUT}/insertion_qc.tsv")
    print(f"{ORG}: {len(picks)} flagged examples")
    for mode, r in picks:
        p = draw(mode, r)
        if p: print("  wrote", os.path.basename(p))
    print("DONE_ARTEFACT_EXAMPLES")


if __name__ == "__main__":
    main()
