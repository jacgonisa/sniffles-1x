#!/usr/bin/env python3
"""Step 20 — insertion ORIGIN tracing: extract the inserted fragment from the read,
map it back to that read's own reference (minimap2 map-hifi), and plot where it came from.

For each of the N largest CIGAR insertions (distinct loci, size >= MINSZ):
  1. pull the read, cut out the inserted bases (the I op) via CIGAR walk
  2. minimap2 the fragment vs REF[refkey]  (primary + up to 8 secondary hits)
  3. classify origin: local tandem (<50 kb, same chrom) / same-CEN / dispersed-CEN178
     (satellite, many hits) / distal / other-chrom, with strand + %identity + distance
  4. per-event panel (chromosome ideogram, acceptor triangle, origin hits, best-hit arc)
     + one genome-wide summary arc plot.
-> results/insertion_origin/*.png  + insertion_origin.tsv
Run with nextflow_env python (pysam + minimap2)."""
import os, csv, subprocess, tempfile
from collections import defaultdict
import numpy as np
import pysam
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from common import OUT, REF, CEN, CHRLEN, bam_path, refkey, GROUPS

MINIMAP2 = "/home/jg2070/minimap2/minimap2"
OUTDIR = f"{OUT}/insertion_origin"; os.makedirs(OUTDIR, exist_ok=True)
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
MINSZ = 1000      # only trace insertions >= 1 kb (origin is meaningful)
N_EVENTS = 10     # how many distinct loci to plot
GLAB = {"wt_leaf": "WT leaf", "cenh3ox_leaf": "CENH3ox leaf",
        "wt_pollen": "WT pollen", "cenh3ox_pollen": "CENH3ox pollen"}


def find_ins_in_query(cigar, target):
    """(q_start,q_end) of the I op whose length ~ target (ref-forward query coords)."""
    q = 0
    for op, ln in cigar:
        if op in (0, 7, 8):
            q += ln
        elif op == 1:
            if abs(ln - target) <= 3:
                return q, q + ln
            q += ln
        elif op == 4:
            q += ln
    return None, None


def pick_events():
    rows = [r for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t")
            if r["svtype"] == "INS" and "CIGAR" in r["methods"]
            and r["svlen"] not in ("", "None") and abs(int(r["svlen"])) >= MINSZ]
    rows.sort(key=lambda r: -abs(int(r["svlen"])))
    picked, seen = [], []
    for r in rows:
        pos = int(r["pos"]); key = (r["sample"], r["hap"], r["chrom"])
        if any(k == key and abs(pos - p) < 5000 for k, p in seen):   # dedupe nearby/recurrent
            continue
        seen.append((key, pos)); picked.append(r)
        if len(picked) >= N_EVENTS:
            break
    return picked


def extract(events):
    """pull inserted sequence for each event; write per-refkey FASTAs."""
    recs = defaultdict(list)   # refkey -> [(seq_id, seq)]
    for i, ev in enumerate(events):
        s, h, chrom, pos, sz = ev["sample"], ev["hap"], ev["chrom"], int(ev["pos"]), abs(int(ev["svlen"]))
        rk = refkey(s, h); ev["rk"] = rk; ev["seq_id"] = f"ins{i:02d}"
        bam = pysam.AlignmentFile(bam_path(s, h), "rb")
        for read in bam.fetch(chrom, max(0, pos - 60000), pos + 60000):
            if read.query_name != ev["read"] or read.is_unmapped or read.cigartuples is None:
                continue
            qs, qe = find_ins_in_query(read.cigartuples, sz)
            if qs is None:
                continue
            recs[rk].append((ev["seq_id"], read.query_sequence[qs:qe]))
            ev["extracted"] = qe - qs
            break
        bam.close()
        ev.setdefault("extracted", 0)
    for rk, rr in recs.items():
        with open(f"{OUTDIR}/frags_{rk}.fa", "w") as f:
            for sid, seq in rr:
                f.write(f">{sid}\n{seq}\n")
    return recs


def map_back(recs):
    """minimap2 each refkey FASTA vs its reference; return seq_id -> [hits]."""
    hits = defaultdict(list)
    for rk, rr in recs.items():
        if not rr:
            continue
        fa = f"{OUTDIR}/frags_{rk}.fa"; paf = f"{OUTDIR}/frags_{rk}.paf"
        subprocess.run([MINIMAP2, "-x", "map-hifi", "-c", "--secondary=yes", "-N", "8",
                        "-t", "8", REF[rk][0], fa], stdout=open(paf, "w"),
                       stderr=subprocess.DEVNULL, check=True)
        for ln in open(paf):
            c = ln.split("\t")
            if len(c) < 12:
                continue
            qn, tname, tstart, tend = c[0], c[5], int(c[7]), int(c[8])
            nmatch, alen, strand = int(c[9]), int(c[10]), c[4]
            hits[qn].append({"chrom": tname, "start": tstart, "end": tend, "strand": strand,
                             "ident": nmatch / max(alen, 1), "alen": alen})
    for qn in hits:
        hits[qn].sort(key=lambda d: -d["alen"] * d["ident"])   # best = longest*identity
    return hits


def classify(ev, hs):
    """origin category from the hits."""
    if not hs:
        return "no_hit", None
    best = hs[0]
    same = [h for h in hs if h["chrom"] == ev["chrom"]]
    cen = CEN.get(ev["rk"], {}).get(ev["chrom"])
    ncen = sum(1 for h in hs if h["chrom"] == ev["chrom"] and cen and cen[0] <= h["start"] < cen[1])
    d = abs(best["start"] - int(ev["pos"])) if best["chrom"] == ev["chrom"] else None
    if len(hs) >= 4 and ncen >= 3:
        return "dispersed_CEN178", best        # satellite: many hits across the array
    if best["chrom"] != ev["chrom"]:
        return "other_chrom", best
    if d is not None and d < 50_000:
        return "local_tandem", best            # nearby duplication (unequal HR / tandem)
    if cen and cen[0] <= best["start"] < cen[1]:
        return "same_CEN", best
    return "distal_same_chrom", best


ARC_COL = {"local_tandem": "#C0392B", "dispersed_CEN178": "#8E44AD", "same_CEN": "#E67E22",
           "distal_same_chrom": "#2980B9", "other_chrom": "#16A085", "no_hit": "#888"}


def panel(ev, hs, cat, best):
    """Zoomed local view: acceptor (insertion site) ▼ + donor fragment(s) ▲ where the
    inserted sequence maps back, with the best-hit arc. x in kb around the event."""
    pos = int(ev["pos"]); sz = abs(int(ev["svlen"]))
    same = [h for h in hs if h["chrom"] == ev["chrom"]]
    xs = [pos] + [h["start"] for h in same] + [h["end"] for h in same]
    lo, hi = (min(xs), max(xs)) if same else (pos - 5000, pos + 5000)
    pad = max((hi - lo) * 0.25, 2500)
    lo -= pad; hi += pad
    fig, ax = plt.subplots(figsize=(9, 3.2))
    ax.axhline(0.5, color="#bbb", lw=1.2, zorder=1)                       # the reference locus
    cen = CEN.get(ev["rk"], {}).get(ev["chrom"])
    in_cen = cen and cen[0] <= pos < cen[1]
    # acceptor (insertion site)
    ax.plot(pos, 0.56, marker="v", ms=15, color="#C0392B", zorder=6)
    ax.text(pos, 0.66, f"insertion site\n{ev['chrom']}:{pos:,}", ha="center", va="bottom",
            fontsize=8.5, color="#C0392B", fontweight="bold")
    # donor fragments (where the inserted sequence maps back)
    for h in same:
        isbest = h is best
        ax.add_patch(Rectangle((h["start"], 0.40), max(h["end"] - h["start"], (hi - lo) * 0.004), 0.05,
                               fc=ARC_COL[cat] if isbest else "#bbb", ec="none", zorder=4))
        if isbest:
            mid = (h["start"] + h["end"]) / 2
            ax.add_patch(FancyArrowPatch((mid, 0.45), (pos, 0.55),
                         connectionstyle="arc3,rad=-0.35", arrowstyle="-|>",
                         mutation_scale=15, lw=2.2, color=ARC_COL[cat], zorder=5))
            ax.text(mid, 0.34, f"origin (donor)\n{h['start']:,} · {h['strand']} · {h['ident']*100:.0f}% id",
                    ha="center", va="top", fontsize=8.5, color=ARC_COL[cat])
    off = ""
    if best and best["chrom"] != ev["chrom"]:
        off = f"  · best hit on {best['chrom']}:{best['start']:,} ({best['ident']*100:.0f}% id)"
    ax.set_xlim(lo, hi); ax.set_ylim(0.15, 0.9); ax.set_yticks([])
    ax.ticklabel_format(axis="x", style="plain")
    ax.set_xticklabels([f"{t/1000:.0f}" for t in ax.get_xticks()])
    ax.set_xlabel(f"{ev['chrom']} position (kb)   ·   {'inside CEN' if in_cen else 'pericentromere'}")
    d = abs(best["start"] - pos) if best and best["chrom"] == ev["chrom"] else None
    dtxt = (f"{d:,} bp away" if d is not None and d < 100000 else
            (f"{d/1000:.0f} kb away" if d is not None else "different chrom"))
    ax.set_title(f"{GLAB[ev['sample']]} · {ev['hap']} · {sz} bp insertion → origin: "
                 f"{cat.replace('_',' ')} ({len(hs)} hit{'s' if len(hs)!=1 else ''}, {dtxt}){off}",
                 fontsize=10)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    p = f"{OUTDIR}/{ev['seq_id']}_{ev['sample']}_{ev['hap']}_{ev['chrom']}_{ev['pos']}.png"
    fig.savefig(p, dpi=140, bbox_inches="tight"); plt.close(fig)
    return p


def summary(events, allhits, cats):
    """genome-wide: acceptor (top) -> origin (bottom) arcs across chromosomes."""
    offs = {}; x = 0
    for c in CHROMS:
        offs[c] = x; x += max(CHRLEN["col"][c], CHRLEN["ler"][c]) / 1e6 + 8
    fig, ax = plt.subplots(figsize=(13, 5))
    for c in CHROMS:
        Ln = CHRLEN["col"][c] / 1e6
        ax.add_patch(Rectangle((offs[c], 2.0), Ln, 0.12, fc="#eee", ec="#999", lw=0.6))
        cen = CEN["col"][c]
        ax.add_patch(Rectangle((offs[c] + cen[0] / 1e6, 2.0), (cen[1] - cen[0]) / 1e6, 0.12, fc="#d8c3e8"))
        ax.text(offs[c] + Ln / 2, 2.25, c, ha="center", fontsize=9)
    for ev in events:
        hs = allhits.get(ev["seq_id"], []); cat = cats.get(ev["seq_id"], "no_hit")
        if not hs or ev["chrom"] not in offs:
            continue
        best = hs[0]; ax = plt.gca()
        ax_x = offs[ev["chrom"]] + int(ev["pos"]) / 1e6
        ax.plot(ax_x, 2.0, marker="v", ms=8, color="#C0392B", zorder=5)
        if best["chrom"] in offs:
            ox = offs[best["chrom"]] + best["start"] / 1e6
            arc = FancyArrowPatch((ax_x, 2.0), (ox, 1.85), connectionstyle="arc3,rad=-0.3",
                                  arrowstyle="-|>", mutation_scale=10, lw=1.6, color=ARC_COL[cat], zorder=4)
            ax.add_patch(arc)
    ax.set_xlim(-3, x); ax.set_ylim(1.4, 2.5); ax.axis("off")
    handles = [plt.Line2D([0], [0], color=ARC_COL[k], lw=2, label=k.replace("_", " ")) for k in ARC_COL if k != "no_hit"]
    ax.legend(handles=handles, loc="lower center", ncol=5, fontsize=9, frameon=False)
    ax.set_title("Insertion origin — inserted fragment mapped back to its reference (▼ = insertion site, arc → origin)")
    p = f"{OUTDIR}/_summary_origin.png"; fig.savefig(p, dpi=140, bbox_inches="tight"); plt.close(fig)
    return p


def main():
    events = pick_events()
    print(f"tracing {len(events)} insertions (>= {MINSZ} bp)")
    recs = extract(events)
    allhits = map_back(recs)
    cats = {}
    with open(f"{OUT}/insertion_origin.tsv", "w") as f:
        f.write("seq_id\tsample\thap\tchrom\tpos\tins_bp\tn_hits\torigin_category\tbest_chrom\tbest_pos\tbest_strand\tbest_ident\tdist_bp\n")
        for ev in events:
            hs = allhits.get(ev["seq_id"], [])
            cat, best = classify(ev, hs); cats[ev["seq_id"]] = cat
            panel(ev, hs, cat, best)
            d = abs(best["start"] - int(ev["pos"])) if best and best["chrom"] == ev["chrom"] else ""
            bc = best["chrom"] if best else ""; bp = best["start"] if best else ""
            bs = best["strand"] if best else ""; bi = f"{best['ident']:.3f}" if best else ""
            f.write(f"{ev['seq_id']}\t{ev['sample']}\t{ev['hap']}\t{ev['chrom']}\t{ev['pos']}\t{abs(int(ev['svlen']))}\t"
                    f"{len(hs)}\t{cat}\t{bc}\t{bp}\t{bs}\t{bi}\t{d}\n")
            print(f"  {ev['seq_id']} {GLAB[ev['sample']]} {ev['hap']} {ev['chrom']}:{ev['pos']} "
                  f"{abs(int(ev['svlen']))}bp -> {cat} ({len(hs)} hits)")
    summary(events, allhits, cats)
    print(f"DONE_INSORIGIN -> {OUTDIR}/")


if __name__ == "__main__":
    main()
