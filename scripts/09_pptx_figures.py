#!/usr/bin/env python3
"""Step 9 — replicate the colleague's SV-analysis slide figures (20260617_SV_analysis.pptx).

Four data plots, from results/sm_sv_calls.tsv (wt_leaf + wt_pollen; artf1 skipped):
  A. genome map per haplotype (Col-HiFi / Ler-HiFi): 5 CEN columns x sample rows,
     each SV at its coordinate, colour=type, size-encoded (dots <100 / 100bp-1kb /
     1-5kb; >=5kb drawn as a horizontal BAR spanning [pos, pos+|svlen|]).
  B. size-binned barplot: facet by type, x=size bin, y=count per MILLION CEN reads,
     fill=sample (col+ler pooled).  -> large SVs enriched in pollen.
  C. log10(width) proportion histograms: facet type x sample, red dashed = 178 bp.
Denominator = # primary reads anchored in CEN per (sample,hap), cached.
-> results/figures/{map_col,map_ler,size_per_million,size_log10}.png + figures_pptx.html
Run with nextflow_env python."""
import os, csv, base64, io, zlib
from collections import defaultdict
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from common import SAMPLES, HAPS, CEN, GROUPS, bam_path, OUT, MONO, refkey

FIGDIR = f"{OUT}/figures"; os.makedirs(FIGDIR, exist_ok=True)
HTML = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/figures_pptx.html"
SAMPLE_ROWS = GROUPS  # wt/cenh3ox × leaf/pollen
SCOL = {"wt_leaf": "#4C9A2A", "cenh3ox_leaf": "#1B5E20", "wt_pollen": "#E8820C", "cenh3ox_pollen": "#B34700"}
# karyogram vertical offset per group (4 lanes around each chromosome bar)
KOFF = {"wt_leaf": 0.42, "cenh3ox_leaf": 0.15, "wt_pollen": -0.15, "cenh3ox_pollen": -0.42}
TYPE4 = ["INS", "DEL", "DUP", "INV"]
COL = {"INS": "#2C6FBB", "DEL": "#D2622B", "DUP": "#E8D44D", "INV": "#3FA45B"}
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
# actual assembly chromosome lengths (.fai) for the karyogram
CHRLEN = {"col": {"Chr1": 32640075, "Chr2": 23012915, "Chr3": 26150667, "Chr4": 22582341, "Chr5": 30170985},
          "ler": {"Chr1": 32485061, "Chr2": 21328600, "Chr3": 27335240, "Chr4": 22700724, "Chr5": 30661135}}
FWD_COL, REV_COL = "#C0392B", "#2C6FBB"   # forward CEN178 array = red, reverse = blue


def load_orient():
    d = {}
    p = f"{OUT}/cen178_orientation.tsv"
    if os.path.exists(p):
        with open(p) as f:
            next(f)
            for ln in f:
                hap, chrom, s, e, st = ln.split()
                d.setdefault((hap, chrom), []).append((int(s), int(e), st))
    return d
BINEDGES = [0, 100, 1000, 10000, 100000, 1000000, float("inf")]
BINLAB = ["<100", "[100-1k)", "[1k-10k)", "[10k-100k)", "[100k-1M)", ">=1M"]


def cen_read_counts():
    """# primary reads anchored in CEN per (sample,hap); cached to tsv."""
    cache = f"{OUT}/cen_read_counts.tsv"
    if os.path.exists(cache):
        out = {}
        with open(cache) as f:
            next(f)
            for ln in f:
                s, h, n = ln.split()
                out[(s, h)] = int(n)
        return out
    out = {}
    for sample, _ in SAMPLES:
        for hap in HAPS:
            bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
            n = 0
            for chrom, (a, b) in CEN[refkey(sample, hap)].items():
                for r in bam.fetch(chrom, a, b):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary:
                        continue
                    if a <= r.reference_start < b:
                        n += 1
            bam.close()
            out[(sample, hap)] = n
    with open(cache, "w") as f:
        f.write("sample\thap\tcen_reads\n")
        for (s, h), n in out.items():
            f.write(f"{s}\t{h}\t{n}\n")
    return out


def load_calls():
    rows = []
    with open(f"{OUT}/sm_sv_calls.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["pos"] = int(r["pos"])
            r["svlen"] = int(r["svlen"]) if r["svlen"] not in ("", "None") else 0
            rows.append(r)
    return rows


def b64(fig):
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=130, bbox_inches="tight"); plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()


def dot_size(a):
    if a < 100: return 6
    if a < 1000: return 16
    return 40   # 1-5kb


def subsample_by_reads(rows, denom):
    """Read-budget match the genome maps: downsample each sample to the same number of
    CEN reads per haplotype (= the smaller sample, i.e. pollen) so visual crowding is fair.
    Deterministic per-read keep (crc32) so all calls of a kept read survive together."""
    N = {h: min(denom[(s, h)] for s in SAMPLE_ROWS) for h in HAPS}
    p = {(s, h): N[h] / denom[(s, h)] for s in SAMPLE_ROWS for h in HAPS}
    out = []
    for r in rows:
        pr = p[(r["sample"], r["hap"])]
        if pr >= 1 or (zlib.crc32(r["read"].encode()) & 0xffffffff) / 2**32 < pr:
            out.append(r)
    return out, N


def fig_map(rows, hap, path, orient, budget=""):
    sub = [r for r in rows if r["hap"] == hap and r["svtype"] in TYPE4]
    nrow, ncol = len(SAMPLE_ROWS), len(CHROMS)
    fig = plt.figure(figsize=(15, 2.5 * nrow + 0.9))
    gs = fig.add_gridspec(nrow + 1, ncol, height_ratios=[0.18] + [1] * nrow, hspace=0.13, wspace=0.18)
    rng = np.random.default_rng(0)
    for j, chrom in enumerate(CHROMS):
        a, b = CEN[hap][chrom]
        # orientation strip (red = forward CEN178 array, blue = reverse)
        axo = fig.add_subplot(gs[0, j])
        for s, e, st in orient.get((hap, chrom), []):
            axo.axvspan(s / 1e6, e / 1e6, color=FWD_COL if st == "+" else REV_COL)
        axo.set_xlim(a / 1e6, b / 1e6); axo.set_ylim(0, 1)
        axo.set_xticks([]); axo.set_yticks([]); axo.set_title(chrom, fontsize=12)
        for i, sample in enumerate(SAMPLE_ROWS):
            ax = fig.add_subplot(gs[i + 1, j])
            for r in (r for r in sub if r["sample"] == sample and r["chrom"] == chrom):
                y = rng.uniform(0, 1)
                a_len = abs(r["svlen"]); c = COL[r["svtype"]]
                if a_len >= 5000 and r["svtype"] != "INS":
                    ax.plot([r["pos"] / 1e6, (r["pos"] + a_len) / 1e6], [y, y], color=c, lw=3, solid_capstyle="butt")
                else:
                    ax.scatter(r["pos"] / 1e6, y, s=dot_size(a_len), color=c, edgecolor="none")
            ax.set_xlim(a / 1e6, b / 1e6); ax.set_ylim(-0.05, 1.05); ax.set_yticks([])
            if i < nrow - 1:
                ax.set_xticklabels([])
            if j == ncol - 1:
                ax.text(1.02, 0.5, sample, transform=ax.transAxes, rotation=270, va="center", fontsize=10)
    fig.text(0.5, 0.04, "Coordinate (Mb)", ha="center", fontsize=12)
    fig.text(0.08, 0.5, f"{'Col' if hap=='col' else 'Ler'}-HiFi", va="center", rotation=90, fontsize=13)
    handles = [Line2D([0], [0], marker='s', color='w', markerfacecolor=COL[t], markersize=10, label=t.lower())
               for t in TYPE4]
    handles += [Line2D([0], [0], marker='o', color='gray', lw=0, markersize=ms / 4, label=lab)
                for ms, lab in [(6, "<100"), (16, "100bp-1kb"), (40, "1-5kb")]]
    handles += [Line2D([0], [0], color='gray', lw=3, label=">=5kb (bar)"),
                Line2D([0], [0], marker='s', color='w', markerfacecolor=FWD_COL, markersize=10, label="CEN178 forward"),
                Line2D([0], [0], marker='s', color='w', markerfacecolor=REV_COL, markersize=10, label="CEN178 reverse")]
    fig.legend(handles=handles, loc="center right", fontsize=9, frameon=False, bbox_to_anchor=(1.0, 0.5))
    fig.suptitle(f"Single-molecule SVs across the centromere — {'Col' if hap=='col' else 'Ler'}-HiFi "
                 f"(read-budget matched: {budget}; top strip = CEN178 orientation)", y=0.99)
    fig.subplots_adjust(left=0.11, right=0.87, bottom=0.11)
    img = b64(fig); open(path, "wb").write(base64.b64decode(img)); return img


def fig_karyogram(rows, hap, path, orient, budget=""):
    """Genome-wide karyogram: full chromosomes, CEN shaded with CEN178 orientation,
    SV events as ticks (leaf above the bar, pollen below), coloured by type."""
    sub = [r for r in rows if r["hap"] == hap and r["svtype"] in TYPE4]
    L = CHRLEN[hap]
    fig, ax = plt.subplots(figsize=(13, 6))
    rng = np.random.default_rng(1)
    for k, chrom in enumerate(CHROMS):
        y = len(CHROMS) - k
        ax.add_patch(Rectangle((0, y - 0.16), L[chrom] / 1e6, 0.32, fc="#f0f0f0", ec="#888", lw=0.7, zorder=1))
        a, b = CEN[hap][chrom]
        # CEN178 orientation inside the centromere band
        for s, e, st in orient.get((hap, chrom), []):
            ax.add_patch(Rectangle((s / 1e6, y - 0.16), (e - s) / 1e6, 0.32,
                                   fc=FWD_COL if st == "+" else REV_COL, ec="none", alpha=0.85, zorder=2))
        ax.text(-0.6, y, chrom, ha="right", va="center", fontsize=11)
        for r in (r for r in sub if r["chrom"] == chrom):
            off = KOFF.get(r["sample"], 0.0)
            jit = rng.uniform(-0.05, 0.05)
            ax.plot([r["pos"] / 1e6], [y + off + jit], marker='|', ms=5, mew=0.8,
                    color=COL[r["svtype"]], zorder=3)
    ax.set_xlim(-2, max(max(v.values()) for v in CHRLEN.values()) / 1e6 + 1)
    ax.set_ylim(0.2, len(CHROMS) + 1.1)
    ax.set_yticks([]); ax.set_xlabel("Coordinate (Mb)")
    ax.text(0.01, 0.99, "lanes top→bottom: WT leaf · CENH3ox leaf · WT pollen · CENH3ox pollen · grey=arm, red/blue=CEN178 fwd/rev",
            transform=ax.transAxes, va="top", fontsize=8, color="#555")
    handles = [Line2D([0], [0], marker='|', color=COL[t], lw=0, markersize=10, mew=2, label=t.lower()) for t in TYPE4]
    ax.legend(handles=handles, loc="lower right", fontsize=9, ncol=4, frameon=False)
    ax.set_title(f"Genome-wide single-molecule SV karyogram — {'Col' if hap=='col' else 'Ler'}-HiFi "
                 f"(read-budget matched: {budget})")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    img = b64(fig); open(path, "wb").write(base64.b64decode(img)); return img


def fig_size_permillion(rows, denom, out="size_per_million", tag=""):
    # pool col+ler; per million CEN reads
    permillion = {s: sum(denom[(s, h)] for h in HAPS) / 1e6 for s in SAMPLE_ROWS}
    cnt = {(s, t): [0] * len(BINLAB) for s in SAMPLE_ROWS for t in TYPE4}
    for r in rows:
        if r["svtype"] not in TYPE4 or r["svlen"] == 0:
            continue
        a = abs(r["svlen"])
        bi = next(k for k in range(len(BINLAB)) if BINEDGES[k] <= a < BINEDGES[k + 1])
        cnt[(r["sample"], r["svtype"])][bi] += 1
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8))
    x = np.arange(len(BINLAB)); w = 0.2
    for ax, t in zip(axes, TYPE4):
        for k, s in enumerate(SAMPLE_ROWS):
            rate = [cnt[(s, t)][bi] / permillion[s] for bi in range(len(BINLAB))]
            ax.bar(x + (k - 1.5) * w, rate, w, color=SCOL[s], label=s)
        ax.set_title(t); ax.set_xticks(x); ax.set_xticklabels(BINLAB, rotation=45, ha="right", fontsize=8)
        ax.set_xlabel("width")
    axes[0].set_ylabel("count per million reads")
    axes[0].legend(fontsize=9, title="sample")
    fig.suptitle(f"SV size spectrum, count per million CEN reads (col+ler pooled){tag}")
    img = b64(fig); open(f"{FIGDIR}/{out}.png", "wb").write(base64.b64decode(img)); return img


def fig_size_log10(rows, out="size_log10", tag=""):
    fig, axes = plt.subplots(len(SAMPLE_ROWS), 4, figsize=(15, 2.6 * len(SAMPLE_ROWS)), squeeze=False)
    line = np.log10(MONO)
    bins = np.linspace(1, 7, 40)
    for j, t in enumerate(TYPE4):
        for i, s in enumerate(SAMPLE_ROWS):
            ax = axes[i][j]
            vals = [np.log10(abs(r["svlen"])) for r in rows
                    if r["sample"] == s and r["svtype"] == t and r["svlen"] and abs(r["svlen"]) >= 10]
            if vals:
                ax.hist(vals, bins=bins, weights=np.ones(len(vals)) / len(vals), color="#888", edgecolor="k", lw=0.3)
            ax.axvline(line, color="#D2622B", ls="--", lw=1.2)
            if i == 0: ax.set_title(t)
            if j == 0: ax.set_ylabel(f"{s}\nProportion", fontsize=9)
            if i == len(SAMPLE_ROWS) - 1: ax.set_xlabel("log10(width)")
            ax.set_xlim(1, 7)
    fig.suptitle(f"SV size distribution (proportion); red dashed = 178 bp CEN178 monomer{tag}")
    img = b64(fig); open(f"{FIGDIR}/{out}.png", "wb").write(base64.b64decode(img)); return img


def load_singletons():
    rows = []
    p = f"{OUT}/singleton_events.tsv"
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            r["svlen"] = int(r["svlen"]) if r["svlen"] not in ("", "None") else 0
            r["pos"] = int(r["pos"])
            rows.append(r)
    return rows


def main():
    rows = load_calls()
    denom = cen_read_counts()
    orient = load_orient()
    print("CEN reads:", {k: v for k, v in denom.items()})
    # read-budget-matched subset for the spatial plots (maps/karyograms) so leaf is not
    # visually denser just from depth; size spectra (per-million) already correct, use full.
    sub, N = subsample_by_reads(rows, denom)
    print(f"read-budget per hap (matched to pollen): {N}")
    print(f"calls after matching: full={len(rows)}  matched={len(sub)}")
    budget = f"{N['col']//1000}k (Col) / {N['ler']//1000}k (Ler) reads per sample"
    mcol = fig_map(sub, "col", f"{FIGDIR}/map_col.png", orient, budget)
    mler = fig_map(sub, "ler", f"{FIGDIR}/map_ler.png", orient, budget)
    kcol = fig_karyogram(sub, "col", f"{FIGDIR}/karyogram_col.png", orient, budget)
    kler = fig_karyogram(sub, "ler", f"{FIGDIR}/karyogram_ler.png", orient, budget)
    # 1x (singleton) genome maps, same read-budget matching
    sing = load_singletons()
    subsing, _ = subsample_by_reads(sing, denom) if sing else ([], None)
    mcol1 = fig_map(subsing, "col", f"{FIGDIR}/map_col_1x.png", orient, budget + " · 1× only") if subsing else ""
    mler1 = fig_map(subsing, "ler", f"{FIGDIR}/map_ler_1x.png", orient, budget + " · 1× only") if subsing else ""
    b7 = fig_size_permillion(rows, denom, "size_per_million", " — all calls")
    b8 = fig_size_log10(rows, "size_log10", " — all calls")
    # 1x (singleton) events only
    singles = load_singletons()
    b7s = fig_size_permillion(singles, denom, "size_per_million_1x", " — 1× (singleton) events only") if singles else ""
    b8s = fig_size_log10(singles, "size_log10_1x", " — 1× (singleton) events only") if singles else ""

    def im(b, cap):
        return f'<figure><img style="max-width:100%" src="data:image/png;base64,{b}"><figcaption>{cap}</figcaption></figure>'
    html = f"""<!doctype html><meta charset=utf-8><title>SV analysis — pptx replication</title>
<style>body{{font-family:-apple-system,Segoe UI,Arial,sans-serif;max-width:1100px;margin:0 auto;padding:24px;color:#1d1d1f}}
h1{{border-bottom:3px solid #C0392B;font-size:23px}}h2{{color:#C0392B;font-size:18px;margin-top:28px}}
figure{{text-align:center;margin:14px 0}}figcaption{{font-size:12.5px;color:#666;font-style:italic}}</style>
<h1>Comprehensive single-molecule SV maps (replicates 20260617_SV_analysis.pptx)</h1>
<p>WT leaf &amp; pollen, both haplotypes. Calls from <code>results/sm_sv_calls.tsv</code> (split-and-map topology, same method
as the slides). artf1 control not yet available. Top strip of each genome map = CEN178 array orientation (red=forward,
blue=reverse), from minimap2 of the 178-bp consensus. Large events (DEL up to ~12.7 Mb, INV up to ~11 Mb, DUP up to ~2.7 Mb).</p>
<h2>1. Genome map — Col-HiFi</h2>{im(mcol, 'Each SV at its CEN coordinate; colour=type; dots sized by class, ≥5 kb as interval-spanning bars. Top strip = CEN178 forward(red)/reverse(blue) arrays. INV events sit at orientation boundaries.')}
<h2>2. Genome map — Ler-HiFi</h2>{im(mler, 'Independent Ler haplotype/assembly, same encoding.')}
<h2>1b. Genome map — Col-HiFi, 1× (singleton) events only</h2>{im(mcol1, 'Same map restricted to support=1 events.') if mcol1 else ''}
<h2>2b. Genome map — Ler-HiFi, 1× (singleton) events only</h2>{im(mler1, 'Ler haplotype, 1× only.') if mler1 else ''}
<h2>3. Genome-wide karyogram — Col-HiFi</h2>{im(kcol, 'Full chromosomes; grey arms, red/blue centromere = CEN178 orientation. SV ticks: leaf above the bar, pollen below. All events are centromeric.')}
<h2>4. Genome-wide karyogram — Ler-HiFi</h2>{im(kler, 'Ler haplotype.')}
<h2>5. Size spectrum — count per million CEN reads (all calls)</h2>{im(b7, 'Per-million-read SV rate by size bin and type (col+ler pooled). Pollen is enriched, especially at 10 kb–Mb scales (DEL/DUP/INV).')}
<h2>6. Size distribution — log10(width) (all calls)</h2>{im(b8, 'Per-facet proportion; red dashed = 178 bp. Pollen carries a heavy large-size tail absent/weak in leaf.')}
<h2>7. Size spectrum — 1× (singleton) events only</h2>{im(b7s, 'Same as §5 but restricted to support=1 events (one read each). The single-molecule size spectrum.')}
<h2>8. Size distribution — 1× (singleton) events only</h2>{im(b8s, 'log10(width) of singleton events; red dashed = 178 bp.')}"""
    open(HTML, "w").write(html)
    print(f"wrote {HTML}"); print("DONE_PPTX_FIGS")


if __name__ == "__main__":
    main()
