#!/usr/bin/env python3
"""draw_event_trash — the lab's validation-prototype panel (plot_validation_prototype.draw_event)
plus a TRASH CEN178 monomer track drawn as ARROWS (pointing in each monomer's strand), so the
head-to-tail satellite structure and the register across the insertion are visible.

Reuses pv.parse_cigar / build_dotplot / readmer_profile_kmc / load_paf and pv colour/const;
only the figure assembly is re-implemented to add the TRASH row."""
import os, sys, csv, tempfile
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/insertion_origin/scripts")
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/insertion_origin/trash-py/src")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon, FancyArrow
import matplotlib.gridspec as gridspec
from types import SimpleNamespace
import plot_validation_prototype as pv
import trash_py.pipeline as _tp, trash_py._log as _tl

MONO = 178
CEN178 = ("AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAA"
          "GACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAG"
          "AAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG")


def trash_read(seq):
    """TRASH on the full read -> [{start,end,width,strand}] in read query coords."""
    if not seq or len(seq) < 250:
        return []
    _tl.configure(quiet=True)
    reps = []
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "s.fasta"); tp = os.path.join(td, "t.fasta")
        open(fa, "w").write(f">r\n{seq}\n"); open(tp, "w").write(f">CEN178\n{CEN178}\n")
        try:
            _tp.run_pipeline(SimpleNamespace(fasta=fa, output=td, max_rep_size=250,
                                             min_rep_size=100, templates=tp, processes=1))
        except Exception:
            pass
        rf = os.path.join(td, "s.fasta_repeats.csv")
        if os.path.exists(rf):
            for row in csv.DictReader(open(rf)):
                reps.append({"start": int(row["start"]) - 1, "end": int(row["end"]),
                             "width": int(row["width"]), "strand": row.get("strand", "+")})
    return reps


def mono_color(w):
    if abs(w - MONO) <= 2: return "#2E7D32"      # full monomer (green)
    if 100 <= w <= 250:    return "#F9A825"      # partial monomer (amber)
    return "#B71C1C"                              # off-size (red)


def draw_event_trash(ev, seq, quals, rq, np_passes, cigar_str, paf_hits, db_prefix, out_png):
    """Faithful copy of pv.draw_event with a TRASH monomer-ARROW row added under the read."""
    ref_pos = int(ev["ref_pos"]); ins_size = int(ev["ins_size"])
    chrom = ev["chrom"]; label = ev["label"]; call = ev["call"]
    C_L, C_R, C_I, C_SM, C_CL, C_RB, C_CEN = pv.C_L, pv.C_R, pv.C_I, pv.C_SM, pv.C_CL, pv.C_RB, pv.C_CEN

    blocks, insertion = pv.parse_cigar(cigar_str, ref_pos, ins_size)
    total_q = max((b["q_e"] for b in blocks), default=len(seq))
    mblk = [b for b in blocks if b["type"] == "match"]
    if not mblk:
        return
    r_lo = min(b["r_s"] for b in mblk); r_hi = max(b["r_e"] for b in mblk)
    best = paf_hits[0] if paf_hits else None
    pad = max(15_000, int((r_hi - r_lo) * 0.2))
    win_s, win_e = r_lo - pad, r_hi + pad
    if best and best["chrom"] == chrom:
        win_s = min(win_s, best["tstart"] - pad // 2); win_e = max(win_e, best["tend"] + pad // 2)
    win_w = max(win_e - win_s, 1)
    qx = lambda q: q / max(total_q, 1)
    rx = lambda r: max(0., min(1., (r - win_s) / win_w))

    q_xs = np.arange(len(quals)) if quals is not None else np.array([])
    q_ys = np.array(quals, dtype=float) if quals is not None else np.array([])
    r_xs, r_ys = pv.readmer_profile_kmc(seq, label, db_prefix)
    dotmat = pv.build_dotplot(seq)
    monos = trash_read(seq)
    ins_q_s = insertion["q_s"] if insertion else total_q
    ins_q_e = insertion["q_e"] if insertion else total_q

    fig = plt.figure(figsize=(20, 12))
    outer = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[3, 1], wspace=0.08,
                              left=0.07, right=0.97, top=0.91, bottom=0.04)
    left_gs = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=outer[0],
                                               height_ratios=[5, 1.1, 1.2, 1.2], hspace=0.16)
    ax_sr = fig.add_subplot(left_gs[0])
    ax_tr = fig.add_subplot(left_gs[1])   # TRASH monomer arrows (NEW)
    ax_q = fig.add_subplot(left_gs[2])
    ax_rm = fig.add_subplot(left_gs[3])
    ax_dot = fig.add_subplot(outer[1])

    # ---- split-read view (normalized 0..1) ----
    ax_sr.set_xlim(0, 1); ax_sr.set_ylim(0, 1); ax_sr.axis("off")
    TRT, TRB, RDT, RDB, BRT, BRB = 0.88, 0.78, 0.58, 0.48, 0.22, 0.12

    def ref_bar(y_top, y_bot, hi_ranges, color, alpha=0.65):
        ax_sr.add_patch(mpatches.Rectangle((0, y_bot), 1, y_top - y_bot, fc=C_RB, ec="#888", lw=0.5, zorder=1))
        cen = pv.CEN_CENH3OX.get(chrom)
        if cen:
            cx1, cx2 = rx(cen[0]), rx(cen[1])
            if cx2 > cx1 + 0.001:
                ax_sr.add_patch(mpatches.Rectangle((cx1, y_bot), cx2 - cx1, y_top - y_bot, fc=C_CEN, ec="none", zorder=2))
        for rs, re_ in hi_ranges:
            x1, x2 = rx(rs), rx(re_)
            if x2 > x1:
                ax_sr.add_patch(mpatches.Rectangle((x1, y_bot), x2 - x1, y_top - y_bot, fc=color, alpha=alpha, ec="none", zorder=3))

    def trap(x1r, x2r, y_rs, x1f, x2f, y_rf, color, alpha=0.20):
        ax_sr.add_patch(Polygon([[x1r, y_rs], [x2r, y_rs], [x2f, y_rf], [x1f, y_rf]], closed=True, fc=color, alpha=alpha, ec="none", zorder=2))

    ins_ranges = [(best["tstart"], best["tend"])] if best and best["chrom"] == chrom else []
    ref_bar(TRT, TRB, ins_ranges, C_I)
    ax_sr.add_patch(mpatches.Rectangle((0, RDB), 1, RDT - RDB, fc="#F5F5F5", ec="#888", lw=0.5, zorder=1))
    for b in blocks:
        mc = C_L if b["q_e"] <= ins_q_s else C_R
        if b["type"] == "match":
            ax_sr.add_patch(mpatches.Rectangle((qx(b["q_s"]), RDB), qx(b["q_e"]) - qx(b["q_s"]), RDT - RDB, fc=mc, ec="none", alpha=0.92, zorder=4))
        elif b["type"] == "ins_big":
            ax_sr.add_patch(mpatches.Rectangle((qx(b["q_s"]), RDB), qx(b["q_e"]) - qx(b["q_s"]), RDT - RDB, fc=C_I, ec="#8B0000", lw=0.7, hatch="////", alpha=0.85, zorder=5))
        elif b["type"] == "clip":
            ax_sr.add_patch(mpatches.Rectangle((qx(b["q_s"]), RDB), qx(b["q_e"]) - qx(b["q_s"]), RDT - RDB, fc=C_CL, ec="none", alpha=0.6, zorder=4))
    left_r = [(b["r_s"], b["r_e"]) for b in mblk if b["q_e"] <= ins_q_s]
    right_r = [(b["r_s"], b["r_e"]) for b in mblk if b["q_s"] >= ins_q_e]
    ref_bar(BRT, BRB, left_r, C_L); ref_bar(BRT, BRB, right_r, C_R)
    for b in mblk:
        mc = C_L if b["q_e"] <= ins_q_s else C_R
        trap(qx(b["q_s"]), qx(b["q_e"]), RDB, rx(b["r_s"]), rx(b["r_e"]), BRT, mc)
    if insertion and best and best["chrom"] == chrom:
        trap(qx(insertion["q_s"]), qx(insertion["q_e"]), RDT, rx(best["tstart"]), rx(best["tend"]), TRB, C_I)
    kw = dict(transform=ax_sr.transAxes, ha="right", va="center", fontsize=7.5)
    ax_sr.text(-0.01, (TRT + TRB) / 2, f"{chrom}\n(ins maps)", color=C_I, fontweight="bold", **kw)
    ax_sr.text(-0.01, (RDT + RDB) / 2, "Read", color="#222", fontweight="bold", **kw)
    ax_sr.text(-0.01, (BRT + BRB) / 2, f"{chrom}\n(flanks)", color=C_L, fontweight="bold", **kw)
    if best:
        midp = (rx(best["tstart"]) + rx(best["tend"])) / 2
        ax_sr.text(midp, TRB - 0.005, f"{best['chrom']}:{best['tstart']//1000:,}k–{best['tend']//1000:,}k  id={best['identity']:.1f}%", ha="center", va="top", fontsize=6.5, color="#8B0000")
    ncop = round(ins_size / MONO)
    cen_tag = f"  [{ncop}×CEN178]" if abs(ins_size - ncop * MONO) <= 5 else ""
    ax_sr.set_title(f"{label} · {chrom}:{ref_pos:,} · {ins_size:,} bp insertion{cen_tag} · {call} · read {total_q:,} bp",
                    fontsize=9, fontweight="bold", loc="left", pad=4, color="#1F3864")

    # ---- TRASH monomer ARROW track (NEW) ----
    ax_tr.set_xlim(0, total_q); ax_tr.set_ylim(0, 1)
    if insertion:
        ax_tr.axvspan(ins_q_s, ins_q_e, color=C_I, alpha=0.08, zorder=0)
        ax_tr.axvline(ins_q_s, color=C_I, lw=1.0, ls="--", alpha=0.7)
        ax_tr.axvline(ins_q_e, color=C_I, lw=1.0, ls="--", alpha=0.7)
    hh = 0.42
    for m in monos:
        s, e, st = m["start"], m["end"], m["strand"]
        w = max(e - s, 1); col = mono_color(m["width"])
        hl = min(w * 0.35, total_q * 0.004)      # arrow head length
        if st == "-":
            ax_tr.add_patch(FancyArrow(e, 0.5, -(w - hl), 0, width=hh, head_width=hh * 1.5,
                                       head_length=hl, length_includes_head=True, fc=col, ec="none", alpha=0.9, zorder=3))
        else:
            ax_tr.add_patch(FancyArrow(s, 0.5, (w - hl), 0, width=hh, head_width=hh * 1.5,
                                       head_length=hl, length_includes_head=True, fc=col, ec="none", alpha=0.9, zorder=3))
    nfull = sum(1 for m in monos if abs(m["width"] - MONO) <= 2)
    ax_tr.text(0.005, 0.97, f"TRASH CEN178 monomers: {len(monos)} ({nfull} full 178 bp) — arrow = strand",
               transform=ax_tr.transAxes, va="top", fontsize=7, color="#37474F", fontweight="bold")
    ax_tr.set_yticks([]); ax_tr.set_ylabel("TRASH\nCEN178", fontsize=7.5)
    ax_tr.spines[["top", "right", "left"]].set_visible(False)
    ax_tr.tick_params(labelbottom=False)

    # ---- quality ----
    for a in (ax_q, ax_rm):
        a.axvspan(ins_q_s, ins_q_e, color=C_I, alpha=0.06)
        a.axvline(ins_q_s, color=C_I, lw=1.0, ls="--", alpha=0.7)
        if insertion:
            a.axvline(ins_q_e, color=C_I, lw=1.0, ls="--", alpha=0.7)
    if len(q_xs):
        ax_q.scatter(q_xs, q_ys, s=2, color="#263238", alpha=0.45, linewidths=0, rasterized=True)
    ax_q.set_xlim(0, total_q); ax_q.set_ylim(bottom=0)
    ax_q.set_ylabel("CCS quality\n(phred)", fontsize=7.5); ax_q.spines[["top", "right"]].set_visible(False)
    ax_q.tick_params(labelsize=7, labelbottom=False)
    ax_q.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"Q{int(v)}"))
    if rq is not None:
        ax_q.text(0.005, 0.95, f"rq={rq:.4f}" + (f"  np={np_passes}" if np_passes is not None else ""),
                  transform=ax_q.transAxes, fontsize=6.5, va="top", color="#37474F", fontweight="bold")

    # ---- readmer ----
    if len(r_xs):
        ax_rm.scatter(r_xs, r_ys, s=2, color="#1B5E20", alpha=0.45, linewidths=0, rasterized=True)
        ax_rm.set_yscale("log"); ax_rm.set_ylim(bottom=0.8)
    else:
        ax_rm.text(0.5, 0.5, "readmer unavailable", transform=ax_rm.transAxes, ha="center", va="center", fontsize=9, color="#aaa")
    ax_rm.set_xlim(0, total_q)
    ax_rm.set_ylabel(f"Readmer\n(k={pv.K_RM}, log)", fontsize=7.5); ax_rm.set_xlabel("Position in read (bp)", fontsize=8)
    ax_rm.spines[["top", "right"]].set_visible(False); ax_rm.tick_params(labelsize=7)
    xt = np.linspace(0, total_q, 6, dtype=int); ax_rm.set_xticks(xt); ax_rm.set_xticklabels([f"{x//1000}k" for x in xt])

    # ---- dotplot ----
    ax_dot.imshow(dotmat, origin="upper", cmap="hot", aspect="equal", interpolation="nearest")
    if insertion:
        fs = ins_q_s / max(total_q, 1) * pv.DOT_RES; fe = ins_q_e / max(total_q, 1) * pv.DOT_RES
        for v in (fs, fe):
            ax_dot.axhline(v, color="#00E5FF", lw=0.8, ls="--", alpha=0.85); ax_dot.axvline(v, color="#00E5FF", lw=0.8, ls="--", alpha=0.85)
    ax_dot.set_xticks([]); ax_dot.set_yticks([])
    ax_dot.set_title(f"self-dotplot (k={pv.K_DOT})\ncyan = insertion", fontsize=8, fontweight="bold")

    fig.suptitle(f"Read + insertion-origin + TRASH register — {label} ({call}) · insertion {ins_size:,} bp @ {chrom}:{ref_pos:,}\n"
                 f"Top: split-read (▲origin) · TRASH CEN178 arrows · CCS quality · KMC readmer · Right: self-dotplot",
                 fontsize=9.5, fontweight="bold")
    fig.savefig(out_png, dpi=170, bbox_inches="tight"); plt.close(fig)
    return out_png
