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
        return None
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
    return frags if len(frags) == 2 else None


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


def draw(ev, bam, fa):
    chrom, pos, svt = ev["chrom"], int(ev["pos"]), ev["svtype"]
    svlen = int(ev["svlen"]); read_id = ev["read"]
    frags, prim = get_alignments(bam, chrom, pos, read_id, svlen)
    if not frags or prim is None:
        print(f"  skip {read_id}: no alignment"); return
    seq = prim.query_sequence; readlen = len(seq)
    # split-and-map events map linearly in the original BAM — reproduce the split & remap
    reproduced = False
    if "SPLITANDMAP" in ev.get("methods", "") and len([f for f in frags if not f["supp"]]) == len(frags) and len(frags) == 1:
        rf = reproduce_split(prim, ev["hap"])
        if rf and all(f["chrom"] == chrom for f in rf):
            frags = rf; reproduced = True
    jq = cigar_junction(prim, pos, svt, svlen)
    rmin = min(f["r_s"] for f in frags); rmax = max(f["r_e"] for f in frags)
    if svt == "DEL" and not reproduced:
        rmin = min(rmin, pos - abs(svlen)); rmax = max(rmax, pos)
    win_s, win_e = rmin - FLANK, rmax + FLANK
    wide = (win_e - win_s) > 60000   # far-apart fragments: skip per-monomer ref track

    # TRASH
    read_reps = trash(seq)
    ref_reps = [] if wide else trash(fa.fetch(chrom, max(0, win_s), win_e))

    qx = lambda q: q / readlen
    rx = lambda r: (r - win_s) / (win_e - win_s)

    fig = plt.figure(figsize=(13, 6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[3.2, 1.0, 1.0], hspace=0.35)
    ax = fig.add_subplot(gs[0]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    REF_T, REF_B = 0.82, 0.70
    RD_T, RD_B = 0.42, 0.30

    # reference band
    ax.add_patch(mp.Rectangle((0, REF_B), 1, REF_T - REF_B, fc="#E8E8E8", ec="#888", lw=0.6))
    ax.text(-0.01, (REF_T + REF_B) / 2, f"{chrom}\nref", ha="right", va="center", fontsize=8)
    # CEN178 on reference (register track sitting just above the ref band)
    for rp in ref_reps:
        gx0 = rx(win_s + rp["start"]); gx1 = rx(win_s + rp["end"])
        ax.add_patch(mp.Rectangle((gx0, REF_T + 0.02), gx1 - gx0, 0.06, fc=mono_color(rp["width"]),
                                  ec="white", lw=0.3))
    # SV interval on reference
    if svt == "DEL":
        dx0, dx1 = rx(pos - abs(svlen)), rx(pos)
        ax.add_patch(mp.Rectangle((dx0, REF_B), dx1 - dx0, REF_T - REF_B, fc="#C0392B", alpha=0.8, ec="#8B0000"))
    else:
        ax.plot([rx(pos)] * 2, [REF_B, REF_T], color="#2980B9", lw=2.5)

    # read band
    ax.add_patch(mp.Rectangle((0, RD_B), 1, RD_T - RD_B, fc="#F5F5F5", ec="#888", lw=0.5))
    ax.text(-0.01, (RD_T + RD_B) / 2, "read\n(query)", ha="right", va="center", fontsize=8)
    for i, f in enumerate(sorted(frags, key=lambda f: f["q_s"])):
        c = FRAGC[i % len(FRAGC)]
        ax.add_patch(mp.Rectangle((qx(f["q_s"]), RD_B), qx(f["q_e"] - f["q_s"]), RD_T - RD_B,
                                  fc=c, ec="none", alpha=0.9))
        ax.annotate("", xy=(qx(f["q_e"] if f["strand"] == "+" else f["q_s"]), (RD_T + RD_B) / 2),
                    xytext=(qx(f["q_s"] if f["strand"] == "+" else f["q_e"]), (RD_T + RD_B) / 2),
                    arrowprops=dict(arrowstyle="-|>", color="white", lw=1.2))
        # connector trapezoid read fragment -> reference span
        ax.add_patch(Polygon([[qx(f["q_s"]), RD_T], [qx(f["q_e"]), RD_T],
                              [rx(f["r_e"]), REF_B], [rx(f["r_s"]), REF_B]],
                             closed=True, fc=c, alpha=0.16, ec="none"))
    # junction marker
    if jq is not None:
        ax.plot([qx(jq)] * 2, [RD_B, RD_T], color="#C0392B", lw=2.5, zorder=6)
        ax.text(qx(jq), RD_T + 0.02, f"{'Δ' if svt=='DEL' else '+'}{abs(svlen):,} bp",
                ha="center", fontsize=8, color="#C0392B", fontweight="bold")

    nfr = len(frags)
    tag = "  ·  split-and-map RE-MAPPED (2 fragments)" if reproduced else ""
    if reproduced and nfr == 2:
        tag += f"  ·  fragments map {abs(frags[1]['r_s'] - frags[0]['r_s'])/1e3:.0f} kb apart"
    ax.set_title(f"{ev.get('sample','')} {ev.get('hap','')}  ·  {read_id}  ·  {chrom}:{pos:,}  ·  {svt} {svlen:+,} bp  "
                 f"·  {nfr} fragment(s){tag}  ·  in_register={ev.get('in_register','?')} "
                 f"(rem={ev.get('monomer_rem','?')})", fontsize=8.5)

    # TRASH(read) detailed register track
    axt = fig.add_subplot(gs[1]); axt.set_xlim(0, readlen); axt.set_ylim(0, 1)
    if jq is not None:
        axt.axvline(jq, color="#C0392B", lw=1.4, ls="--")
        axt.axvspan(jq - MONO, jq + MONO, color="#C0392B", alpha=0.06)
    for rp in read_reps:
        axt.add_patch(mp.Rectangle((rp["start"], 0.05), rp["end"] - rp["start"], 0.9,
                                   fc=mono_color(rp["width"]), ec="white", lw=0.4))
    ncen = sum(1 for rp in read_reps if abs(rp["width"] - MONO) <= 2)
    axt.text(0.99, 0.86, f"{ncen} CEN178 monomers in read", transform=axt.transAxes, ha="right",
             va="top", fontsize=7.5, color="#1B5E20", fontweight="bold")
    axt.set_yticks([]); axt.set_ylabel("TRASH\n(read)", fontsize=7.5)
    axt.set_xlabel("position in read (bp)", fontsize=8)
    axt.spines[["top", "right", "left"]].set_visible(False)

    # legend / register note
    axl = fig.add_subplot(gs[2]); axl.axis("off")
    axl.text(0.0, 0.7, "Register check: if the CEN178 monomers tile continuously up to the red junction and the "
             "deletion/insertion size is a whole number of 178-bp monomers, the array is IN-REGISTER "
             "(unequal sister-chromatid HR). A phase break or non-monomer size = out-of-register.",
             fontsize=8, va="top", wrap=True)
    handles = [mp.Patch(fc="#2E7D32", label="178 bp monomer"), mp.Patch(fc="#F57F17", label="other 100-250 bp"),
               mp.Patch(fc="#B71C1C", label="out of range")]
    axl.legend(handles=handles, loc="lower left", ncol=3, fontsize=8, frameon=False)

    out = f"{OUTDIR}/{ev.get('sample','x')}_{ev.get('hap','x')}_{chrom}_{pos}_{svt}.png"
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {out}")


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--read"); ap.add_argument("--pos", type=int)
    ap.add_argument("--sample"); ap.add_argument("--hap")
    a = ap.parse_args()
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
