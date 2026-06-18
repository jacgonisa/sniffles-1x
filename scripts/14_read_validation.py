#!/usr/bin/env python3
"""Step 14 — per-read validation plots (verify the mapping/split + CEN178 register).

For chosen single-molecule SV reads, draw (style after insertion_origin
plot_del_read_validation.py):
  • REFERENCE band  — genomic window around the event; SV interval marked; CEN178
    monomers (TRASH on the reference) drawn as a register track.
  • READ band       — the read in query coords with each alignment fragment (primary +
    supplementary) coloured; the split / indel junction marked; trapezoids connect each
    read fragment to where it maps on the reference (so you see HOW it is split & mapped).
  • TRASH(read)     — CEN178 monomers on the actual molecule (TRASH on the FULL read),
    coloured by width (178 bp = green), junction line + ±178 bp zone -> is the register
    conserved across the junction?

-> results/read_validation/<sample>_<hap>_<chrom>_<pos>_<svtype>.png
Run with BASE python (trash_py + pysam):
  /home/jg2070/miniforge3/bin/python 14_read_validation.py            # auto-pick examples
  /home/jg2070/miniforge3/bin/python 14_read_validation.py --read <name> --pos <p>
"""
import sys, os, csv, tempfile, argparse, subprocess
from types import SimpleNamespace
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/insertion_origin/trash-py/src")
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts")
import pysam, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import trash_py.pipeline as _tp, trash_py._log as _tl
from common import bam_path, REF, OUT, MONO

OUTDIR = f"{OUT}/read_validation"; os.makedirs(OUTDIR, exist_ok=True)
WIN = "/home/jg2070/miniforge3/envs/nextflow_env/bin/winnowmap"
ST = "/home/jg2070/miniforge3/envs/nextflow_env/bin/samtools"
CEN178 = ("CEN178", "AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAA"
          "GACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAG"
          "AAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG")
FLANK = 2500
FRAGC = ["#1f6fbf", "#e07b1a", "#2e8b57", "#8e44ad"]


def trash(seq):
    if not seq or len(seq) < 250:
        return []
    _tl.configure(quiet=True)
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "s.fasta"); tp = os.path.join(td, "t.fasta")
        open(fa, "w").write(f">r\n{seq}\n"); open(tp, "w").write(f">{CEN178[0]}\n{CEN178[1]}\n")
        try:
            _tp.run_pipeline(SimpleNamespace(fasta=fa, output=td, max_rep_size=250,
                                             min_rep_size=100, templates=tp, processes=1))
        except Exception:
            pass
        rf = os.path.join(td, "s.fasta_repeats.csv"); reps = []
        if os.path.exists(rf):
            for row in csv.DictReader(open(rf)):
                reps.append({"start": int(row["start"]) - 1, "end": int(row["end"]), "width": int(row["width"])})
    return reps


def mono_color(w):
    if w == MONO: return "#2E7D32"
    if abs(w - MONO) <= 2: return "#558B2F"
    if 100 <= w <= 250: return "#F57F17"
    return "#B71C1C"


def get_alignments(bam, chrom, pos, read_id, svlen):
    """All alignments (primary+supp) of the read near the event. Returns list of
    fragments {q_s,q_e,r_s,r_e,strand,primary} in forward-read query coords, plus
    primary read object (for full sequence + CIGAR junction)."""
    frags = []; prim = None
    for r in bam.fetch(chrom, max(0, pos - abs(svlen) - 5000), pos + abs(svlen) + 5000):
        if r.query_name != read_id or r.is_unmapped or r.is_secondary:
            continue
        # query span in forward-read coords
        if r.is_reverse:
            qs = r.infer_read_length() - r.query_alignment_end
        else:
            qs = r.query_alignment_start
        frags.append({"q_s": qs, "q_e": qs + r.query_alignment_length,
                      "r_s": r.reference_start, "r_e": r.reference_end,
                      "strand": "-" if r.is_reverse else "+", "supp": r.is_supplementary})
        if not r.is_supplementary:
            prim = r
    return frags, prim


def best_split(read):
    """MD substitution-contrast split point (same as step 03) in read query coords."""
    qpos, mism = [], []
    for q, r, base in read.get_aligned_pairs(with_seq=True):
        if q is None or r is None or base is None:
            continue
        qpos.append(q); mism.append(1 if base.islower() else 0)
    n = len(qpos)
    if n < 200:
        return None
    pre = [0] * (n + 1)
    for i in range(n):
        pre[i + 1] = pre[i] + mism[i]
    best = None; q0, q1 = qpos[0], qpos[-1]
    for i in range(1, n):
        if qpos[i] - q0 < 1000 or q1 - qpos[i] < 1000:
            continue
        c = abs((pre[n] - pre[i]) / (n - i) - pre[i] / i)
        if best is None or c > best[1]:
            best = (qpos[i], c)
    return best[0] if best and best[1] >= 0.01 else None


def reproduce_split(read, hap):
    """Re-run the split-and-map for this read: split at the contrast point, winnowmap
    re-map both halves, return the 2 fragments (original-read query coords)."""
    qsplit = best_split(read)
    seq = read.query_sequence
    if qsplit is None or seq is None or qsplit < 1000 or len(seq) - qsplit < 1000:
        return None, None
    ref, rep = REF[hap]
    with tempfile.TemporaryDirectory() as td:
        fa = f"{td}/f.fa"; bamf = f"{td}/f.bam"
        open(fa, "w").write(f">A\n{seq[:qsplit]}\n>B\n{seq[qsplit:]}\n")
        cmd = (f"{WIN} -W {rep} -ax map-pb -t 8 {ref} {fa} 2>/dev/null | "
               f"{ST} sort -o {bamf} - && {ST} index {bamf}")
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
        b = pysam.AlignmentFile(bamf, "rb"); frags = []
        for r in b.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            off = 0 if r.query_name == "A" else qsplit
            frags.append({"q_s": off + r.query_alignment_start, "q_e": off + r.query_alignment_end,
                          "r_s": r.reference_start, "r_e": r.reference_end,
                          "strand": "-" if r.is_reverse else "+", "supp": False, "chrom": r.reference_name})
        b.close()
    return (frags, qsplit) if len(frags) == 2 else (None, None)


def cigar_junction(read, pos, svtype, svlen):
    """Read query pos of the indel junction (for CIGAR INS/DEL), else None."""
    if read is None:
        return None
    size = abs(svlen); rpos = read.reference_start; qpos = 0; best = None
    for op, l in read.cigartuples or []:
        if op in (0, 7, 8):
            rpos += l; qpos += l
        elif op == 1:
            if svtype == "INS" and abs(l - size) <= max(15, size * 0.1):
                d = abs(rpos - pos);  best = (d, qpos, l) if best is None or d < best[0] else best
            qpos += l
        elif op == 2:
            if svtype == "DEL" and abs(l - size) <= max(15, size * 0.1):
                d = abs(rpos - pos);  best = (d, qpos, l) if best is None or d < best[0] else best
            rpos += l
        elif op in (4, 5):
            qpos += l
    return best[1] if best else None


def sub_profile(read, readlen, win=151):
    """Per-base smoothed substitution (mismatch) rate of the ORIGINAL linear alignment,
    in read/query coords. Shows where the aligner extended through mismatches."""
    m = np.zeros(readlen); cov = np.zeros(readlen)
    for q, r, base in read.get_aligned_pairs(with_seq=True):
        if q is None or q >= readlen:
            continue
        if r is not None:
            cov[q] = 1
            if base is not None and base.islower():
                m[q] = 1
    k = np.ones(win)
    msum = np.convolve(m, k, "same"); csum = np.convolve(cov, k, "same")
    return np.where(csum > 0, msum / np.maximum(csum, 1.0), np.nan)


def ref_box(ax, f, bx0, bx1, yb, yt, fa, chrom, color, sv_interval=None):
    """Local reference window for one fragment: grey bar, mapped segment highlighted,
    CEN178 monomers, coord label. Returns the box x-range for the connector."""
    flank = 1500
    ws, we = f["r_s"] - flank, f["r_e"] + flank
    span = max(we - ws, 1)
    rxl = lambda r: bx0 + (r - ws) / span * (bx1 - bx0)
    ax.add_patch(mp.Rectangle((bx0, yb), bx1 - bx0, yt - yb, fc="#ECECEC", ec="#999", lw=0.6))
    for rp in trash(fa.fetch(chrom, max(0, ws), we)):
        gx0, gx1 = rxl(ws + rp["start"]), rxl(ws + rp["end"])
        ax.add_patch(mp.Rectangle((gx0, yt - 0.03), gx1 - gx0, 0.03, fc=mono_color(rp["width"]), ec="white", lw=0.2))
    # mapped segment
    ax.add_patch(mp.Rectangle((rxl(f["r_s"]), yb + 0.01), rxl(f["r_e"]) - rxl(f["r_s"]), yt - yb - 0.05,
                              fc=color, alpha=0.85, ec="none"))
    xa0, xa1 = (rxl(f["r_s"]), rxl(f["r_e"])) if f["strand"] == "+" else (rxl(f["r_e"]), rxl(f["r_s"]))
    ax.annotate("", xy=(xa1, (yb + yt) / 2), xytext=(xa0, (yb + yt) / 2),
                arrowprops=dict(arrowstyle="-|>", color="white", lw=1.3))
    if sv_interval:  # deleted ref interval (DEL)
        a, b = sv_interval
        ax.add_patch(mp.Rectangle((rxl(a), yb), rxl(b) - rxl(a), yt - yb, fc="none", ec="#C0392B", lw=1.6, hatch="///"))
    ax.text((bx0 + bx1) / 2, yt + 0.03, f"{chrom}:{f['r_s']/1e6:.3f}–{f['r_e']/1e6:.3f} Mb ({f['strand']})",
            ha="center", fontsize=7.5, color="#333")


def draw(ev, bam, fa, outdir=OUTDIR, prefix=""):
    chrom, pos, svt = ev["chrom"], int(ev["pos"]), ev["svtype"]
    svlen = int(ev["svlen"]); read_id = ev["read"]
    frags, prim = get_alignments(bam, chrom, pos, read_id, svlen)
    if not frags or prim is None:
        print(f"  skip {read_id}: no alignment"); return
    seq = prim.query_sequence; readlen = len(seq)
    # split-and-map events map linearly in the original BAM — reproduce the split & remap
    reproduced = False; qsplit = None
    if "SPLITANDMAP" in ev.get("methods", "") and len(frags) == 1:
        rf, qs = reproduce_split(prim, ev["hap"])
        if rf and all(f["chrom"] == chrom for f in rf):
            frags = rf; reproduced = True; qsplit = qs
    jq = cigar_junction(prim, pos, svt, svlen)
    fsort = sorted(frags, key=lambda f: f["q_s"])
    # the junction in read coords
    if reproduced and qsplit is not None:
        junction_q = qsplit
    elif jq is not None:
        junction_q = jq
    elif len(fsort) >= 2:
        junction_q = fsort[0]["q_e"]
    else:
        junction_q = readlen // 2

    rate = sub_profile(prim, readlen)
    read_reps = trash(seq)
    qx = lambda q: q / readlen

    fig = plt.figure(figsize=(13, 7.8))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1.0, 2.5, 1.25], hspace=0.55)

    # ── Row A: original-mapping substitution profile ───────────────────────────
    axA = fig.add_subplot(gs[0]); axA.set_xlim(0, readlen)
    axA.fill_between(np.arange(readlen), rate * 100, color="#C0392B", alpha=0.55, lw=0)
    axA.axvline(junction_q, color="#222", ls="--", lw=1.3)
    lhs = np.nanmean(rate[:junction_q]) * 100 if junction_q > 0 else 0
    rhs = np.nanmean(rate[junction_q:]) * 100 if junction_q < readlen else 0
    axA.text(0.02, 0.82, f"left of split: {lhs:.2f}% mismatch", transform=axA.transAxes, ha="left", fontsize=7.5, color="#555")
    axA.text(0.98, 0.82, f"right of split: {rhs:.2f}% mismatch", transform=axA.transAxes, ha="right", fontsize=7.5, color="#555")
    axA.set_ylabel("orig. mapping\nsubst. % (151 bp)", fontsize=7.5)
    axA.set_title(f"{ev.get('sample','')} {ev.get('hap','')}  ·  {read_id}  ·  {chrom}:{pos:,}  ·  {svt} {svlen:+,} bp"
                  f"  ·  ORIGINAL linear mapping: mismatches accumulate on one side → the aligner extended instead of "
                  f"splitting (dashed = split point)", fontsize=8)
    axA.spines[["top", "right"]].set_visible(False); axA.set_xticklabels([])

    # ── Row B: read ↔ reference (per-fragment local windows) ───────────────────
    axB = fig.add_subplot(gs[1]); axB.axis("off"); axB.set_xlim(0, 1); axB.set_ylim(0, 1)
    RD_B, RD_T = 0.04, 0.18
    BOX_B, BOX_T = 0.60, 0.78
    axB.add_patch(mp.Rectangle((0, RD_B), 1, RD_T - RD_B, fc="#F5F5F5", ec="#888", lw=0.5))
    axB.text(-0.012, (RD_B + RD_T) / 2, "read\n(query)", ha="right", va="center", fontsize=8)
    axB.text(-0.012, (BOX_B + BOX_T) / 2, "reference\n(per fragment)", ha="right", va="center", fontsize=8)
    n = len(fsort)
    for i, f in enumerate(fsort):
        c = FRAGC[i % len(FRAGC)]
        # read fragment
        axB.add_patch(mp.Rectangle((qx(f["q_s"]), RD_B), qx(f["q_e"] - f["q_s"]), RD_T - RD_B, fc=c, alpha=0.9, ec="none"))
        axB.annotate("", xy=(qx(f["q_e"] if f["strand"] == "+" else f["q_s"]), (RD_B + RD_T) / 2),
                     xytext=(qx(f["q_s"] if f["strand"] == "+" else f["q_e"]), (RD_B + RD_T) / 2),
                     arrowprops=dict(arrowstyle="-|>", color="white", lw=1.2))
        # evenly-spaced reference box for this fragment
        bx0 = 0.02 + i * (0.96 / n); bx1 = bx0 + 0.96 / n - 0.06
        ivl = (pos - abs(svlen), pos) if (svt == "DEL" and not reproduced) else None
        ref_box(axB, f, bx0, bx1, BOX_B, BOX_T, fa, chrom, c, sv_interval=ivl)
        # ribbon read-fragment → its ref box
        axB.add_patch(Polygon([[qx(f["q_s"]), RD_T], [qx(f["q_e"]), RD_T], [bx1, BOX_B], [bx0, BOX_B]],
                              closed=True, fc=c, alpha=0.14, ec="none"))
    # junction line on read
    axB.plot([qx(junction_q)] * 2, [RD_B, RD_T], color="#C0392B", lw=2.5, zorder=6)
    axB.text(qx(junction_q), RD_T + 0.015, f"split / {'Δ' if svt=='DEL' else '+' if svt=='INS' else ''}"
             f"{abs(svlen):,} bp", ha="center", fontsize=7.5, color="#C0392B", fontweight="bold")
    tag = ""
    if reproduced:
        tag = f"  —  split-and-map RE-MAPPED: 2 fragments map {abs(fsort[1]['r_s'] - fsort[0]['r_s'])/1e3:.0f} kb apart"
        if svt == "INV":
            tag += " on opposite strands"
    elif n >= 2:
        tag = "  —  native split read (SA tag from the aligner)"
    axB.text(0.5, 0.93, f"How the read maps & is split{tag}", ha="center", fontsize=8.5, fontweight="bold")

    # ── Row C: TRASH register on read, ±2 kb junction window highlighted ────────
    axt = fig.add_subplot(gs[2]); axt.set_xlim(0, readlen); axt.set_ylim(0, 1)
    axt.axvspan(junction_q - 2000, junction_q + 2000, color="#FFD54F", alpha=0.30, zorder=0)
    axt.axvline(junction_q, color="#C0392B", lw=1.6, ls="--", zorder=4)
    for rp in read_reps:
        axt.add_patch(mp.Rectangle((rp["start"], 0.05), rp["end"] - rp["start"], 0.9,
                                   fc=mono_color(rp["width"]), ec="white", lw=0.4, zorder=2))
    ncen = sum(1 for rp in read_reps if abs(rp["width"] - MONO) <= 2)
    axt.text(0.5, 1.18, f"CEN178 register on the read  ·  ±2 kb around the junction shaded  ·  "
             f"in_register={ev.get('in_register','?')} (rem={ev.get('monomer_rem','?')})  ·  {ncen} monomers",
             transform=axt.transAxes, ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    axt.set_yticks([]); axt.set_ylabel("TRASH\n(read)", fontsize=7.5)
    axt.set_xlabel("position in read (bp)", fontsize=8)
    axt.spines[["top", "right", "left"]].set_visible(False)
    handles = [mp.Patch(fc="#2E7D32", label="178 bp"), mp.Patch(fc="#F57F17", label="100-250 bp"),
               mp.Patch(fc="#B71C1C", label="out of range"), mp.Patch(fc="#FFD54F", alpha=0.5, label="±2 kb of junction")]
    axt.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.65), ncol=4, fontsize=7.5, frameon=False)

    os.makedirs(outdir, exist_ok=True)
    out = f"{outdir}/{prefix}{ev.get('sample','x')}_{ev.get('hap','x')}_{chrom}_{pos}_{svt}.png"
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {out}")
    return out


def pick_examples():
    """A few representative annotated singletons: in-register DEL/INS, out-of-register, split."""
    rows = list(csv.DictReader(open(f"{OUT}/singleton_events_annotated.tsv"), delimiter="\t"))
    for r in rows:
        r["svlen"] = int(r["svlen"]); r["pos"] = int(r["pos"])
    def f(cond, n):
        return [r for r in rows if cond(r)][:n]
    picks = []
    picks += f(lambda r: r["svtype"] == "DEL" and r["in_register"] == "1" and r["confidence"] == "HIGH" and abs(r["svlen"]) >= 356, 2)
    picks += f(lambda r: r["svtype"] == "INS" and r["in_register"] == "1" and r["confidence"] == "HIGH", 1)
    picks += f(lambda r: r["svtype"] in ("DEL", "INS") and r["in_register"] == "0" and (r["in_cen180_array"] in ("1", 1)), 1)
    picks += f(lambda r: "SPLITANDMAP" in r["methods"] and r["svtype"] in ("DUP", "INV", "DEL") and abs(r["svlen"]) >= 5000, 2)
    # dedup
    seen = set(); out = []
    for p in picks:
        k = (p["read"], p["pos"])
        if k not in seen:
            seen.add(k); out.append(p)
    return out


GALLERY = f"{OUT}/validation_gallery"


def gallery():
    """Curated set covering each case -> results/validation_gallery/ + index.html."""
    rows = list(csv.DictReader(open(f"{OUT}/singleton_events_annotated.tsv"), delimiter="\t"))
    for r in rows:
        r["svlen"] = int(r["svlen"]); r["pos"] = int(r["pos"])
    def sel(cond, n):
        return [r for r in rows if cond(r)][:n]
    cats = [
        ("01_inregister_DEL", "In-register deletion (whole CEN178 monomers, HIGH)",
         lambda r: r["svtype"] == "DEL" and r["in_register"] == "1" and r["confidence"] == "HIGH" and abs(r["svlen"]) >= 356, 3),
        ("02_outofregister_DEL", "Out-of-register deletion (in CEN178 array, not a whole monomer)",
         lambda r: r["svtype"] == "DEL" and r["in_register"] == "0" and r["in_cen180_array"] == "1", 2),
        ("03_inregister_INS", "In-register insertion (whole CEN178 monomers)",
         lambda r: r["svtype"] == "INS" and r["in_register"] == "1" and r["confidence"] == "HIGH", 2),
        ("04_splitread_DEL", "Native split-read (aligner made the SA split, Sniffles-style)",
         lambda r: r["methods"] == "SPLITREAD" and r["svtype"] == "DEL", 2),
        ("05_splitandmap_DUP", "split-and-map duplication (we cut & re-map → fragments land apart)",
         lambda r: "SPLITANDMAP" in r["methods"] and r["svtype"] == "DUP" and abs(r["svlen"]) >= 5000, 2),
        ("06_splitandmap_INV", "split-and-map inversion (opposite-strand fragments)",
         lambda r: "SPLITANDMAP" in r["methods"] and r["svtype"] == "INV", 2),
    ]
    os.makedirs(GALLERY, exist_ok=True)
    cache = {}; sections = []
    for prefix, desc, cond, n in cats:
        evs = sel(cond, n)
        print(f"[{prefix}] {len(evs)} events")
        imgs = []
        for ev in evs:
            key = (ev["sample"], ev["hap"])
            if key not in cache:
                cache[key] = (pysam.AlignmentFile(bam_path(*key), "rb"), pysam.FastaFile(REF[ev["hap"]][0]))
            try:
                p = draw(ev, *cache[key], outdir=GALLERY, prefix=prefix + "__")
                if p:
                    imgs.append(os.path.basename(p))
            except Exception as e:
                print(f"   err {ev['read']}: {e}")
        sections.append((desc, imgs))
    html = ["<!doctype html><meta charset=utf-8><title>Read validation gallery</title>",
            "<style>body{font-family:-apple-system,Arial,sans-serif;max-width:1200px;margin:0 auto;padding:20px}"
            "h2{color:#C0392B;border-bottom:1px solid #eee;margin-top:26px}img{max-width:100%;border:1px solid #eee;margin:8px 0}</style>",
            "<h1>Single-molecule SV — read validation gallery</h1>"]
    for desc, imgs in sections:
        html.append(f"<h2>{desc}</h2>")
        for im in imgs:
            html.append(f'<img src="{im}">')
    open(f"{GALLERY}/index.html", "w").write("\n".join(html))
    print(f"wrote {GALLERY}/index.html ({sum(len(i) for _, i in sections)} plots)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--read"); ap.add_argument("--pos", type=int)
    ap.add_argument("--sample"); ap.add_argument("--hap")
    ap.add_argument("--gallery", action="store_true")
    a = ap.parse_args()
    if a.gallery:
        gallery(); print("DONE_VALIDATION"); return
    if a.read:
        rows = [r for r in csv.DictReader(open(f"{OUT}/singleton_events_annotated.tsv"), delimiter="\t")
                if r["read"] == a.read and (a.pos is None or int(r["pos"]) == a.pos)]
        for r in rows:
            r["svlen"] = int(r["svlen"]); r["pos"] = int(r["pos"])
        events = rows
    else:
        events = pick_examples()
    print(f"rendering {len(events)} validation plots")
    cache = {}
    for ev in events:
        key = (ev["sample"], ev["hap"])
        if key not in cache:
            cache[key] = (pysam.AlignmentFile(bam_path(*key), "rb"), pysam.FastaFile(REF[ev["hap"]][0]))
        draw(ev, *cache[key])
    print("DONE_VALIDATION")


if __name__ == "__main__":
    main()
