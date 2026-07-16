#!/usr/bin/env python3
"""Step 6 — figures + standalone HTML report for the single-molecule SV callset.
Reads results/sm_sv_calls.tsv -> single_molecule_sv/report.html (PNGs inline)."""
import os, csv, base64, io
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import OUT, MONO, GROUPS, genotype

REPORT = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/report.html"
TYPES = ["DEL", "INS", "DUP", "INV", "BND"]
COL = {"DEL": "#C0392B", "INS": "#2980B9", "DUP": "#27AE60", "INV": "#8E44AD", "BND": "#E67E22"}
# the 4 genotype×tissue groups (haplotypes pooled); label + colour
GLAB = {"wt_leaf": "WT leaf", "cenh3ox_leaf": "CENH3ox leaf",
        "wt_pollen": "WT pollen", "cenh3ox_pollen": "CENH3ox pollen"}
GCOL = {"wt_leaf": "#4C9A2A", "cenh3ox_leaf": "#1B5E20",
        "wt_pollen": "#E8820C", "cenh3ox_pollen": "#B34700"}


def load():
    rows = []
    with open(f"{OUT}/sm_sv_calls.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["svlen"] = int(r["svlen"]) if r["svlen"] not in ("", "None") else 0
            rows.append(r)
    return rows


def png(fig):
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=130, bbox_inches="tight"); plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()


def file_b64(path):
    """base64 of a PNG already on disk (the step-09 comprehensive figures), or None."""
    if not os.path.exists(path):
        return None
    return base64.b64encode(open(path, "rb").read()).decode()


def fig_counts(rows):
    """raw calls per group (haps pooled), grouped by SV type."""
    cnt = defaultdict(Counter)
    for r in rows:
        cnt[r["sample"]][r["svtype"]] += 1
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    x = range(len(GROUPS)); w = 0.16
    for i, t in enumerate(TYPES):
        ax.bar([xi + i * w for xi in x], [cnt[g][t] for g in GROUPS], w, label=t, color=COL[t])
    ax.set_xticks([xi + 2 * w for xi in x]); ax.set_xticklabels([GLAB[g] for g in GROUPS])
    ax.set_ylabel("single-molecule SV calls"); ax.legend(ncol=5, fontsize=8)
    ax.set_title("Single-molecule centromere SVs by type — WT vs CENH3ox (raw; haps pooled)")
    return png(fig)


def load_rates():
    rows = []
    p = f"{OUT}/sm_sv_rates.tsv"
    if not os.path.exists(p):
        return rows
    with open(p) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    return rows


def group_rates(rates):
    """pool haplotypes: per group, calls/Mb per type = Σ(per_mb·Mb)/ΣMb."""
    mb = defaultdict(float); cnt = defaultdict(lambda: defaultdict(float))
    for r in rates:
        g = r["sample"]; m = float(r["cen_mapped_mb"]); mb[g] += m
        for t in TYPES:
            cnt[g][t] += float(r.get(f"{t}_per_mb", 0) or 0) * m
    return {g: {t: (cnt[g][t] / mb[g] if mb[g] else 0) for t in TYPES} for g in GROUPS}, mb


def fig_interaction(rates):
    """Everything in one plot: per-Mb SV rate, tissue on x, one line per genotype,
    so the genotype×tissue interaction (CENH3ox lifts leaf most) is visible.
    Thin markers = col/ler haplotypes; thick line = Mb-weighted pooled group rate."""
    gr, mb = group_rates(rates)
    grate = {g: sum(gr[g].values()) for g in GROUPS}
    perhap = {(r["sample"], r["hap"]): float(r["ALL_per_mb"]) for r in rates}
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    X = {"leaf": 0, "pollen": 1}
    for geno, col in (("wt", "#2980B9"), ("cenh3ox", "#C0392B")):
        gs = [g for g in GROUPS if g.startswith(geno)]
        ys = [grate[g] for g in gs]
        xs = [X["leaf" if "leaf" in g else "pollen"] for g in gs]
        ax.plot(xs, ys, "-o", color=col, lw=2.6, ms=9, zorder=3,
                label=("WT" if geno == "wt" else "CENH3ox"))
        for g in gs:
            xx = X["leaf" if "leaf" in g else "pollen"]
            for h in ("col", "ler"):
                if (g, h) in perhap:
                    ax.plot(xx, perhap[(g, h)], "o", color=col, ms=5, alpha=0.35, zorder=2)
            ax.annotate(f"{grate[g]:.2f}", (xx, grate[g]), (xx + 0.03, grate[g] + 0.06),
                        fontsize=9, color=col, fontweight="bold")
    # fold annotations
    for geno, col, yo in (("wt", "#2980B9", -0.28), ("cenh3ox", "#C0392B", 0.18)):
        gl = grate[f"{geno}_leaf"]; gp = grate[f"{geno}_pollen"]
        ax.annotate(f"pollen/leaf = {gp/max(gl,1e-9):.1f}×", (0.5, (gl + gp) / 2 + yo),
                    ha="center", fontsize=8.5, color=col, style="italic")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["leaf", "pollen"], fontsize=11)
    ax.set_ylabel("single-molecule SV calls per Mb (CEN)"); ax.set_xlim(-0.25, 1.35)
    ax.set_ylim(0, max(grate.values()) * 1.2)
    ax.legend(title="genotype", fontsize=10)
    ax.set_title("Genotype × tissue: CENH3ox lifts the leaf (somatic) rate most\n(faint dots = col/ler haplotypes)")
    fig.savefig(f"{OUT}/figures/rate_interaction.png", dpi=140, bbox_inches="tight")  # standalone for slides
    return png(fig)


def fig_rates(rates):
    gr, mb = group_rates(rates)
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    x = range(len(GROUPS)); w = 0.16
    for i, t in enumerate(TYPES):
        ax.bar([xi + i * w for xi in x], [gr[g][t] for g in GROUPS], w, label=t, color=COL[t])
    ax.set_xticks([xi + 2 * w for xi in x]); ax.set_xticklabels([GLAB[g] for g in GROUPS])
    ax.set_ylabel("calls per Mb of CEN-mapped read seq"); ax.legend(ncol=5, fontsize=8)
    ax.set_title("Read-Mb-normalized single-molecule SV rate — WT vs CENH3ox (haps pooled)")
    return png(fig)


def load_qc():
    p = f"{OUT}/read_qc.tsv"
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fig_qc(qc):
    qc = sorted(qc, key=lambda r: (GROUPS.index(r["sample"]) if r["sample"] in GROUPS else 9, r["hap"]))
    labels = [f"{GLAB.get(r['sample'], r['sample'])}\n{r['hap']}" for r in qc]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    x = range(len(qc))
    cols = [GCOL.get(r["sample"], "#888") for r in qc]
    axes[0].bar(x, [float(r["arm_de_pct"]) for r in qc], color=cols)
    axes[0].set_title("Arm de% (≈ seq error)"); axes[0].set_ylabel("% divergence")
    axes[1].bar(x, [float(r["np_med"]) for r in qc], color=cols)
    axes[1].set_title("HiFi passes (np, median)")
    axes[2].bar(x, [float(r["cen_med_kb"]) for r in qc], color=cols)
    axes[2].set_title("CEN read length (median kb)")
    for ax in axes:
        ax.set_xticks(list(x)); ax.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
    fig.suptitle("Read-quality controls per group (WT vs CENH3ox, leaf & pollen)")
    return png(fig)


def fig_sizes(rows):
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8), sharex=True, sharey=True)
    for ax, g in zip(axes, GROUPS):
        for t in ("DEL", "INS", "DUP"):
            sizes = [abs(r["svlen"]) for r in rows if r["sample"] == g and r["svtype"] == t and r["svlen"]]
            if sizes:
                ax.hist(sizes, bins=range(0, 3000, 50), histtype="step", color=COL[t], label=t, linewidth=1.3)
        for m in range(MONO, 3000, MONO):
            ax.axvline(m, color="#bbb", ls=":", lw=0.6)
        ax.set_title(GLAB[g], fontsize=10); ax.set_xlabel("|SV length| (bp)"); ax.set_xlim(0, 3000); ax.legend(fontsize=7)
    axes[0].set_ylabel("count")
    fig.suptitle("Size distribution (dotted = CEN178 178-bp monomer multiples)")
    return png(fig)


def fig_methods(rows):
    m = Counter(r["methods"] for r in rows)
    fig, ax = plt.subplots(figsize=(7, 3.6))
    items = sorted(m.items(), key=lambda kv: -kv[1])
    ax.barh([k for k, _ in items], [v for _, v in items], color="#34495E")
    for i, (_, v) in enumerate(items):
        ax.text(v, i, f" {v}", va="center", fontsize=9)
    ax.set_xlabel("calls"); ax.set_title("Detection method (CIGAR / split-read / split-and-map)")
    ax.invert_yaxis()
    return png(fig)


def fig_routes(path):
    """How many reads each detection route recovered (from source_breakdown.tsv)."""
    if not os.path.exists(path):
        return ""
    rws = list(csv.DictReader(open(path), delimiter="\t"))
    calls = [int(r["calls"]) for r in rws]
    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    y = range(len(rws))
    ax.barh(list(y), calls, color=["#2980B9", "#E67E22", "#27AE60", "#8E44AD"][:len(rws)])
    for i, r in enumerate(rws):
        ax.text(int(r["calls"]), i, f"  {r['calls']}  ({float(r['pct_in_register']):.0f}% in-register)",
                va="center", fontsize=9)
    ax.set_yticks(list(y)); ax.set_yticklabels([r["methods"] for r in rws], fontsize=9)
    ax.set_xlabel("single-molecule SV calls"); ax.invert_yaxis()
    ax.set_xlim(0, max(calls) * 1.45)
    ax.set_title("Reads recovered per detection route (CIGAR vs split-and-map vs split-read)")
    return png(fig)


def main():
    rows = load()
    n = len(rows)
    by_grp = Counter(r["sample"] for r in rows)
    by_type = Counter(r["svtype"] for r in rows)
    stock = sum(int(r["stock_match"]) for r in rows)
    inph = [int(r["in_phase"]) for r in rows if r["in_phase"] not in ("", "None")]
    rates = load_rates()
    qc = load_qc()
    gr, gmb = group_rates(rates)           # per-group per-type calls/Mb, and Mb per group
    grate = {g: sum(gr[g].values()) for g in GROUPS}   # ALL calls/Mb per group
    f1, fr, f2, f3 = fig_counts(rows), fig_rates(rates), fig_sizes(rows), fig_methods(rows)
    finter = fig_interaction(rates)
    fq = fig_qc(qc) if qc else ""

    # dataset-at-a-glance (Arabidopsis / this run) — per group (genotype×tissue) × haplotype
    cenreads = {}
    if os.path.exists(f"{OUT}/cen_read_counts.tsv"):
        for r in csv.DictReader(open(f"{OUT}/cen_read_counts.tsv"), delimiter="\t"):
            cenreads[(r["sample"], r["hap"])] = int(r["cen_reads"])
    cand = {}
    for p in sorted(__import__("glob").glob(f"{OUT}/candidates/*.tsv")):
        b = os.path.basename(p)[:-4]              # e.g. cenh3ox_leaf_col
        s, h = b.rsplit("_", 1); cand[(s, h)] = sum(1 for _ in open(p)) - 1
    calls_sh = Counter((r["sample"], r["hap"]) for r in rows)
    glance = ("<table><tr><th>genotype</th><th>tissue</th><th>hap</th><th>CEN reads</th>"
              "<th>candidate reads (de≥0.005 ∨ NM≥50 ∨ SA)</th><th>SV calls</th></tr>"
              + "".join(
                  f"<tr><td>{genotype(g)}</td><td>{'leaf' if 'leaf' in g else 'pollen'}</td><td>{h}</td>"
                  f"<td>{cenreads.get((g, h), 0):,}</td>"
                  f"<td>{cand.get((g, h), 0):,}</td><td>{calls_sh[(g, h)]:,}</td></tr>"
                  for g in GROUPS for h in ("col", "ler"))
              + f"<tr><td colspan=3><b>total</b></td><td><b>{sum(cenreads.values()):,}</b></td>"
              f"<td><b>{sum(cand.values()):,}</b></td><td><b>{n:,}</b></td></tr></table>")

    qc = sorted(qc, key=lambda r: (GROUPS.index(r["sample"]) if r["sample"] in GROUPS else 9, r["hap"]))
    qc_tbl = ("<table><tr><th>group</th><th>hap</th><th>CEN read len (kb, median)</th>"
              "<th>arm de% (≈ error)</th><th>CEN de%</th><th>np (median passes)</th><th>rq% (median)</th></tr>"
              + "".join(
                  f"<tr><td>{GLAB.get(r['sample'], r['sample'])}</td><td>{r['hap']}</td><td>{float(r['cen_med_kb']):.1f}</td>"
                  f"<td>{float(r['arm_de_pct']):.3f}</td><td>{float(r['cen_de_pct']):.3f}</td>"
                  f"<td>{float(r['np_med']):.0f}</td><td>{float(r['rq_med_pct']):.3f}</td></tr>"
                  for r in qc) + "</table>") if qc else ""

    rate_tbl = ("<table><tr><th>group</th><th>CEN-mapped Mb</th>"
                "<th>ALL / Mb</th><th>DEL / Mb</th><th>INS / Mb</th><th>DUP / Mb</th><th>INV / Mb</th><th>BND / Mb</th></tr>"
                + "".join(
                    f"<tr><td><b>{GLAB[g]}</b></td><td>{gmb[g]:.0f}</td><td><b>{grate[g]:.3f}</b></td>"
                    f"<td>{gr[g]['DEL']:.3f}</td><td>{gr[g]['INS']:.3f}</td><td>{gr[g]['DUP']:.3f}</td>"
                    f"<td>{gr[g]['INV']:.3f}</td><td>{gr[g]['BND']:.3f}</td></tr>" for g in GROUPS)
                + "</table>")

    # per-group QC medians (haps averaged) for the §3 text
    qcg = {}
    for g in GROUPS:
        gr_rows = [r for r in qc if r["sample"] == g]
        if gr_rows:
            qcg[g] = {k: sum(float(r[k]) for r in gr_rows) / len(gr_rows)
                      for k in ("cen_med_kb", "arm_de_pct", "np_med", "rq_med_pct")}

    def qcmp(metric):  # "wt_leaf X vs CENH3ox Y" style helper
        return {g: qcg.get(g, {}).get(metric, 0) for g in GROUPS}

    def im(b, cap):
        return f'<figure><img src="data:image/png;base64,{b}"><figcaption>{cap}</figcaption></figure>'

    def imf(name, cap):  # full-width figure from a step-09 PNG file
        b = file_b64(f"{OUT}/figures/{name}.png")
        return (f'<figure><img style="max-width:100%" src="data:image/png;base64,{b}">'
                f'<figcaption>{cap}</figcaption></figure>') if b else ""

    def imv(name, cap):  # read-validation PNG (step 14)
        b = file_b64(f"{OUT}/read_validation/{name}.png")
        return (f'<figure><img style="max-width:100%" src="data:image/png;base64,{b}">'
                f'<figcaption>{cap}</figcaption></figure>') if b else ""

    # translocations (BND) section (step 15)
    transloc = ""
    tp = f"{OUT}/translocations.tsv"
    if os.path.exists(tp):
        tl = list(csv.DictReader(open(tp), delimiter="\t"))
        catc = Counter(d["category"] for d in tl)
        grpc = Counter(d["sample"] for d in tl)
        ex = "".join(f"<tr><td>{GLAB.get(d['sample'], d['sample'])}</td><td>{d['chrom']}:{d['pos']}</td><td>{d['mate']}</td>"
                     f"<td>{d['category']}</td><td>{d['methods']}</td></tr>" for d in tl[:12])
        bg = {}
        if os.path.exists(f"{OUT}/crossmap_background.tsv"):
            for d in csv.DictReader(open(f"{OUT}/crossmap_background.tsv"), delimiter="\t"):
                bg[d["group"]] = d
        bgtxt = ""
        if bg:
            bgrows = "".join(
                f"<tr><td>{GLAB.get(g, g)}</td><td>{bg[g]['interCEN_BND']}</td><td>{int(bg[g]['cen_reads']):,}</td>"
                f"<td><b>{float(bg[g]['per_million_reads']):.0f}</b></td></tr>" for g in GROUPS if g in bg)
            bgtxt = (f"<div class=box style='background:#FDEDEC;border-left:4px solid #C0392B'>"
                     f"<b>Satellite cross-mapping noise background.</b> Inter-CEN BNDs are fragments demonstrably "
                     f"landing in the <i>wrong</i> centromere (shared CEN178), so their per-million-read rate is an empirical "
                     f"noise floor for single-molecule satellite split calls. It is much higher in pollen than leaf (shorter "
                     f"reads → more ambiguous fragments), so the split-and-map classes (DUP/INV/BND) should be read cautiously; "
                     f"the CIGAR DEL/INS signal (no re-mapping, register-checked) is unaffected. This floor is a property of the "
                     f"satellite, so it lets us compare WT vs CENH3ox on equal footing."
                     f"<table><tr><th>group</th><th>inter-CEN BND</th><th>CEN reads</th><th>/ million reads</th></tr>{bgrows}</table></div>")
        transloc = f"""<h2>11. Translocations (BND) &amp; the cross-mapping noise floor</h2>
<p>A read whose two fragments map to <b>different contigs</b> is classified <b>BND</b> by
<code>sv.classify_splits</code> — the inter-chromosomal / translocation class. <b>{len(tl)} BND calls</b>
({', '.join(f'{GLAB[g]} {grpc[g]}' for g in GROUPS if grpc[g])}); the partner locus is in the <code>mate</code> column.
By mate category: {', '.join(f'{k} {v}' for k, v in catc.items())}. <b>{catc.get('other_CEN',0)} ({100*catc.get('other_CEN',0)//max(len(tl),1)}%)
map to another centromere</b> — all five centromeres share CEN178, so these are overwhelmingly satellite cross-mapping, not
real translocations; {catc.get('unplaced_organellar',0)} hit unplaced/organellar contigs. Mates on other-chromosome arms
({catc.get('other_chrom_arm',0)}) are possible real junctions but unconfirmable from a single read.</p>
{bgtxt}
<table><tr><th>group</th><th>breakpoint</th><th>mate</th><th>category</th><th>method</th></tr>{ex}</table>"""

    # insertion origin (step 20)
    insorig = ""
    iop = f"{OUT}/insertion_origin.tsv"
    if os.path.exists(iop):
        il = list(csv.DictReader(open(iop), delimiter="\t"))
        catc2 = Counter(d["origin_category"] for d in il)
        import glob as _g2
        def iimg(path, cap):
            b = file_b64(path)
            return (f'<figure><img style="max-width:100%" src="data:image/png;base64,{b}">'
                    f'<figcaption>{cap}</figcaption></figure>') if b else ""
        summ = iimg(f"{OUT}/insertion_origin/_summary_origin.png",
                    "Genome-wide: ▼ = insertion site, arc → where the inserted fragment maps back. All trace to their own pericentromere/CEN.")
        # one local-tandem + one dispersed example panel
        expngs = sorted(_g2.glob(f"{OUT}/insertion_origin/ins*.png"))
        expick = expngs[:1]
        for d in il:
            if d["origin_category"] == "dispersed_CEN178":
                m = [p for p in expngs if f"_{d['chrom']}_{d['pos']}." in p]
                expick += m[:1]; break
        panels = "".join(iimg(p, os.path.basename(p)[:-4]) for p in expick)
        # detailed split-read + dotplot + quality/readmer panels (step 21)
        dpngs = sorted(_g2.glob(f"{OUT}/insertion_origin_detailed/*.png"))
        detailed = ""
        if dpngs:
            dp = "".join(iimg(p, os.path.basename(p)[:-4]) for p in dpngs[:2])
            detailed = (f"<h3>14b. Detailed read view — split-read origin · self-similarity dotplot · CCS-quality &amp; k-mer readmer</h3>"
                        f"<p>Per read: the <b>split-read alignment</b> (top bar = where the inserted fragment maps back = origin; "
                        f"middle = the read, blue flanks + red insertion; bottom bar = flanks), the <b>self-similarity dotplot</b> "
                        f"(off-diagonal bands = the internal tandem repeat structure of the inserted copy), and per-base <b>CCS quality</b> "
                        f"+ <b>KMC readmer</b> (dataset k-mer support: high = real/repeated sequence, 1 = sequencing error). "
                        f"Full set: <code>results/insertion_origin_detailed/</code>.</p>{dp}")
        trows = "".join(
            f"<tr><td>{GLAB.get(d['sample'],d['sample'])}</td><td>{d['hap']}</td><td>{d['chrom']}:{int(d['pos']):,}</td>"
            f"<td>{d['ins_bp']}</td><td>{d['origin_category'].replace('_',' ')}</td><td>{d['n_hits']}</td>"
            f"<td>{('%.0f%%'%(float(d['best_ident'])*100)) if d['best_ident'] else '–'}</td>"
            f"<td>{(str(int(d['dist_bp']))+' bp' if d['dist_bp'] not in ('','None') else '–')}</td></tr>" for d in il)
        insorig = f"""<h2>14. Where does the inserted fragment come from? (insertion origin)</h2>
<p>For each large insertion (≥1 kb) we cut the <b>inserted bases out of the read</b> and map that fragment back to the same
reference (minimap2 <code>map-hifi</code>). The location(s) it maps to are the <b>origin</b> of the inserted sequence.
Of the {len(il)} traced, <b>{catc2.get('local_tandem',0)} are local tandem duplications</b> (the fragment is a near-identity
copy of sequence within a few kb — the unequal-sister-chromatid-HR / satellite-expansion signature) and
{catc2.get('dispersed_CEN178',0)} are dispersed CEN178 (the fragment hits many places across the array). None came from a
different chromosome. So these centromeric insertions are <b>locally templated</b>, not captured from elsewhere.</p>
{summ}
{panels}
<table><tr><th>group</th><th>hap</th><th>insertion site</th><th>size (bp)</th><th>origin</th><th>#hits</th><th>best id</th><th>donor distance</th></tr>{trows}</table>
<p class=cap style="font-size:12.5px;color:#666">Full set of per-event panels: <code>results/insertion_origin/</code>.</p>
{detailed}"""

    # CEN vs ARM control (step 16)
    armc = ""
    ap2 = f"{OUT}/arm_control.tsv"
    if os.path.exists(ap2):
        al2 = list(csv.DictReader(open(ap2), delimiter="\t"))
        trow = "".join(
            f"<tr><td>{GLAB.get(d['group'], d['group'])}</td><td>{d['svtype']}</td><td>{d['cen_per_Mreads']}</td>"
            f"<td>{d['arm_per_Mreads']}</td><td><b>{d['enrichment_CEN_over_ARM']}×</b></td></tr>" for d in al2)
        armsm = ""
        sp = f"{OUT}/arm_splitmap_control.tsv"
        if os.path.exists(sp):
            sl = list(csv.DictReader(open(sp), delimiter="\t"))
            srow = "".join(
                f"<tr><td>{GLAB.get(d['group'], d['group'])}</td><td>{d['svtype']}</td><td>{d['cen_per_Mreads']}</td>"
                f"<td>{d['arm_per_Mreads']}</td><td><b>{d['enrichment_CEN_over_ARM']}×</b></td></tr>" for d in sl)
            armsm = f"""<h3>12b. Same control on the split-and-map route</h3>
<table><tr><th>group</th><th>type</th><th>CEN /M reads</th><th>ARM /M reads</th><th>CEN ÷ ARM</th></tr>{srow}</table>
<div class=box style="background:#FDEDEC;border-left:4px solid #C0392B"><b>Reading the split-and-map route:</b>
<b>DEL is robustly real</b> (~55× CEN, and it still occurs in unique arm sequence). <b>BND in the arm is zero</b> —
unique sequence cannot cross-map, confirming split-and-map BND in the CEN is satellite cross-mapping. <b>DUP/INV are
arm-absent (∞), but that is NOT proof they are real</b>: a forced split of unique-sequence fragments structurally cannot
produce a DUP/INV (no overlapping/cross mapping), so the arm gives no usable baseline for these two — they are
satellite-exclusive by construction and remain confounded between true CEN rearrangement and satellite mapping.</div>"""
        armc = f"""<h2>12. Centromere vs chromosome-arm control (is the signal satellite-specific?)</h2>
<p>The same per-read leadprov detection (CIGAR + native split) run on reads anchored in the chromosome <b>arms</b>
(unique sequence) — windows starting <b>5 Mb past the centromere</b> so the CEN↔ARM transition/pericentromere is excluded.
Arms are the no-satellite background. The CEN÷ARM ratio is the centromere-specific enrichment:</p>
<table><tr><th>group</th><th>type</th><th>CEN /M reads</th><th>ARM /M reads</th><th>CEN ÷ ARM</th></tr>{trow}</table>
<div class=box><b>Reading the leadprov route:</b> <b>DEL/INS/DUP are genuinely centromere-specific</b> (several-fold CEN÷ARM in
every group) — real centromere instability, in both WT and CENH3ox. <b>BND is NOT enriched (ratio ≈ 1.0)</b> — same rate in
arm and CEN, a uniform background (justifies treating native-split BND as noise). leadprov INV is not CEN-enriched (arms carry
inverted repeats).</p>{armsm}
<p class=cap style="font-size:12.5px;color:#666"><code>results/arm_control.tsv</code> · <code>results/arm_splitmap_control.tsv</code></p>"""

    # SV source breakdown (step 17)
    srcbreak = ""
    sb = f"{OUT}/source_breakdown.tsv"
    if os.path.exists(sb):
        sbl = list(csv.DictReader(open(sb), delimiter="\t"))
        label = {"CIGAR": "CIGAR (inline indel)", "SPLITREAD": "SPLITREAD (aligner SA split)",
                 "SPLITANDMAP": "SPLITANDMAP (we split + re-map)", "CIGAR+SPLITANDMAP": "CIGAR + SPLITANDMAP (both)"}
        rws = "".join(
            f"<tr><td>{label.get(d['methods'], d['methods'])}</td><td>{d['calls']}</td>"
            f"<td>{d['pct_in_register']}% ({d['in_register']})</td><td>{d['DEL']}</td><td>{d['INS']}</td>"
            f"<td>{d['DUP']}</td><td>{d['INV']}</td><td>{d['BND']}</td></tr>" for d in sbl)
        route_fig = fig_routes(sb)
        # per-group route breakdown (WT vs CENH3ox)
        grp_route_tbl = ""
        gbp = f"{OUT}/source_breakdown_by_group.tsv"
        if os.path.exists(gbp):
            gbl = list(csv.DictReader(open(gbp), delimiter="\t"))
            routes = ["CIGAR", "SPLITREAD", "SPLITANDMAP", "CIGAR+SPLITANDMAP"]
            cell = {(d["group"], d["methods"]): d for d in gbl}
            hdr = "".join(f"<th>{r.replace('CIGAR+SPLITANDMAP','both')}</th>" for r in routes)
            body = "".join(
                "<tr><td><b>" + GLAB[g] + "</b></td>" + "".join(
                    (f"<td>{cell[(g, r)]['calls']} ({cell[(g, r)]['pct_in_register']}%)</td>" if (g, r) in cell else "<td>–</td>")
                    for r in routes) + "</tr>" for g in GROUPS)
            grp_route_tbl = (f"<p><b>Route × group</b> (calls, with in-register %):</p>"
                             f"<table><tr><th>group</th>{hdr}</tr>{body}</table>")
        srcbreak = f"""<h2>13. Where the calls come from (detection route) &amp; in-register by route</h2>
<p>Each call is found by one or more routes: <b>CIGAR</b> (inline I/D in a single alignment), <b>SPLITREAD</b>
(the aligner already split the read via its <code>SA</code> tag — Sniffles-compatible), <b>SPLITANDMAP</b> (we split the
read at the contrast frontier and re-mapped — not visible to stock Sniffles). in-register = whole-CEN178-monomer fraction
(|svlen| mod 178 heuristic; the rigorous TRASH version is §9, singletons).</p>
{im(route_fig, 'Reads recovered per route. CIGAR (inline indel) finds the most and is the most in-register; split-and-map adds ~2000 events invisible to stock Sniffles; split-read is the aligner’s own SA splits.') if route_fig else ''}
<table><tr><th>route</th><th>calls</th><th>in-register</th><th>DEL</th><th>INS</th><th>DUP</th><th>INV</th><th>BND</th></tr>{rws}</table>
{grp_route_tbl}
<div class=box>CIGAR and SPLITANDMAP are both highly in-register. The split point in split-and-map is
found by a <b>CUSUM change-point</b>: walk the read's per-base match/mismatch (0/1), keep a running sum of
(value − mean mismatch rate) — it drifts up through noisy stretches and down through clean ones, so its extreme marks the
clean↔noisy <b>frontier</b>, which is where we cut. This replaced an earlier global-contrast argmax that got pulled to a
front-loaded noisy patch (e.g. cutting at 1.0 kb when the real boundary was 2.95 kb), and it lifted SPLITANDMAP in-register
from ~60% to ~83% — cutting at the true boundary makes the re-mapped fragments align at monomer boundaries. SPLITREAD is
lower (many are large non-monomer junctions / BND). See <code>ALGORITHM.md</code>.</div>"""

    import glob as _glob
    CATLAB = {
        "01_inregister_DEL": "In-register DEL — whole-CEN178-monomer deletion (unequal-sister-chromatid-HR signature)",
        "02_outofregister_DEL": "Out-of-register DEL — junction not on a monomer boundary (NHEJ-like)",
        "03_inregister_INS": "In-register INS — whole-monomer insertion",
        "04_splitread_DEL": "Split-read DEL — the aligner already split the read (its own SA tag)",
        "05_splitandmap_DUP": "Split-and-map DUP — read cut at the CUSUM frontier and re-mapped; the two halves land apart",
        "06_splitandmap_INV": "Split-and-map INV — second fragment re-maps in the opposite orientation",
    }
    gallery = sorted(_glob.glob(f"{OUT}/validation_gallery/*.png"))
    valfigs = sorted(_glob.glob(f"{OUT}/read_validation/*.png"))
    val = ""
    if gallery:
        def gimg(p):
            b = file_b64(p)
            rid = os.path.basename(p)[:-4].split("__", 1)[-1]
            return (f'<figure><img style="max-width:100%" src="data:image/png;base64,{b}">'
                    f'<figcaption>{rid.replace("_", " ")}</figcaption></figure>') if b else ""
        blocks = ""
        cur = None
        for p in gallery:
            cat = os.path.basename(p).split("__", 1)[0]
            if cat != cur:
                blocks += f"<h3>{CATLAB.get(cat, cat)}</h3>"
                cur = cat
            blocks += gimg(p)
        val = f"""<h2>10. Read-level validation — test cases, one figure per read</h2>
<p>Worked examples so the calls can be eyeballed, grouped by category. <b>Each panel is one read.</b> It shows the
<b>reference</b> (SV interval + its CEN178 monomer track), the <b>read</b> in query coordinates with every alignment
fragment and the junction (for split-and-map events the split is <b>reproduced</b> — the read is cut at the CUSUM contrast
point and both halves re-mapped, so you see them land apart), and a <b>TRASH CEN178 track on the read</b> coloured by
monomer width. If the 178-bp monomers tile continuously through the red junction and the size is a whole number of monomers,
the array is in-register. Files: <code>results/validation_gallery/</code>.</p>
{blocks}"""
    elif valfigs:
        items = "".join(imv(os.path.basename(p)[:-4],
                            os.path.basename(p)[:-4].replace("_", " ")) for p in valfigs)
        val = f"""<h2>10. Read-level validation (mapping · split · CEN178 register)</h2>
<p>Per-read sanity plots so the calls can be eyeballed. Each shows the <b>reference</b> (with the SV interval and its CEN178
monomer track), the <b>read</b> in query coordinates with every alignment fragment and the junction (for split-and-map events
the split is <b>reproduced</b> — read cut at the contrast point and both halves re-mapped, so you see them land apart), and a
<b>TRASH track on the read</b> coloured by monomer width. If the 178-bp monomers tile continuously through the red junction and
the size is a whole number of monomers, the array is in-register. Files: <code>results/read_validation/</code>.</p>
{items}"""

    fd = f"{OUT}/figures"
    comp = ""
    if os.path.exists(f"{fd}/size_per_million.png"):
        comp = f"""<h2>6. Comprehensive SV maps (replicates 20260617_SV_analysis.pptx)</h2>
<p>Genome maps per haplotype (≥5 kb SVs drawn as bars spanning the interval), the size spectrum as
<b>count per million CEN reads</b> (col+ler pooled), and log10(width) proportion histograms. Large events
(DEL up to ~12.7 Mb, INV ~11 Mb, DUP ~2.7 Mb) come from the split-read fragments. artf1 control not yet available.
<b>The genome maps &amp; karyograms are read-budget matched</b> (each sample downsampled to the same CEN-read count
per haplotype, ≈25k = pollen's depth) so leaf is not visually denser just from its ~14× higher coverage; even matched,
pollen carries more and longer large-SV bars.</p>
<h3>6a. Genome map — Col-HiFi</h3>{imf('map_col', 'Read-budget matched (leaf downsampled to pollen depth). Each SV at its CEN coordinate; colour=type; dots sized by class, ≥5 kb as interval-spanning bars. Top strip = CEN178 forward(red)/reverse(blue) arrays; INV events sit at orientation boundaries.')}
<h3>6b. Genome map — Ler-HiFi</h3>{imf('map_ler', 'Independent Ler haplotype/assembly, same encoding.')}
<h3>6a′. Genome map — Col-HiFi, 1× (singleton) events only</h3>{imf('map_col_1x', 'Same map restricted to support=1 (single-molecule) events; read-budget matched.')}
<h3>6b′. Genome map — Ler-HiFi, 1× (singleton) events only</h3>{imf('map_ler_1x', 'Ler haplotype, 1× only.')}
<h3>6c. Genome-wide karyogram — Col-HiFi</h3>{imf('karyogram_col', 'Full chromosomes; grey arms, red/blue centromere = CEN178 orientation. SV ticks: leaf above bar, pollen below. All events are centromeric.')}
<h3>6d. Genome-wide karyogram — Ler-HiFi</h3>{imf('karyogram_ler', 'Ler haplotype.')}
<h3>6e. Size spectrum — count per million CEN reads (all calls)</h3>{imf('size_per_million', 'Per-million-read SV rate by size bin and type. Pollen enriched at 10 kb–Mb (DEL/DUP/INV).')}
<h3>6f. Size distribution — log10(width) (all calls)</h3>{imf('size_log10', 'Per-facet proportion; red dashed = 178 bp. Pollen carries a heavy large-size tail absent/weak in leaf.')}
<h3>6g. Size spectrum — 1× (singleton) events only</h3>{imf('size_per_million_1x', 'Same as 6e but restricted to support=1 events (one read each) — the single-molecule size spectrum. Same pollen large-SV enrichment.')}
<h3>6h. Size distribution — 1× (singleton) events only</h3>{imf('size_log10_1x', 'log10(width) of singleton events; red dashed = 178 bp.')}"""

    # recurrence section (step 11)
    recur = ""
    rp = f"{OUT}/recurrent_loci.tsv"
    if os.path.exists(rp):
        rl = list(csv.DictReader(open(rp), delimiter="\t"))
        nfix = sum(1 for d in rl if d["class"] == "FIXED_vs_ref")
        ncand = len(rl) - nfix
        top = "".join(
            f"<tr><td>{d['sample']}</td><td>{d['hap']}</td><td>{d['chrom']}:{d['pos']}</td><td>{d['svtype']}</td>"
            f"<td>{d['support']}</td><td>{d['coverage']}</td><td>{d['vaf']}</td><td>{d['median_size']}</td>"
            f"<td>{'fixed vs ref' if d['class']=='FIXED_vs_ref' else 'hotspot?'}</td></tr>"
            for d in rl[:12])
        recur = f"""<h2>7. Recurrent positions (the vertical lines) — hotspot or artifact?</h2>
<p>A vertical line is many reads carrying the same SV at one coordinate. The discriminator is <b>VAF = supporting
reads / spanning coverage</b>. Of {len(rl)} recurrent loci (≥10 reads), <b>{nfix} are FIXED differences vs the
reference assembly</b> (VAF≥0.30, near-identical whole-CEN178-monomer sizes — 178/356/712/1067/3026 bp): the reads
agree with each other but disagree with the assembly, so these are reference/assembly discrepancies, <b>not somatic
hotspots</b>. They are tallest in deep leaf simply because more reads span them. The remaining {ncand} are low-VAF
recurrent loci — genuine hotspot candidates, though in satellite they cannot be cleanly separated from mapping
ambiguity. Take-home: the prominent vertical lines are mostly fixed assembly differences; filter by VAF before
calling anything a somatic hotspot.</p>
<table><tr><th>sample</th><th>hap</th><th>locus</th><th>type</th><th>support</th><th>cov</th><th>VAF</th><th>median size</th><th>class</th></tr>{top}</table>
<p class=cap style="font-size:12.5px;color:#666">Full list: <code>results/recurrent_loci.tsv</code>.</p>"""

    # support-distribution section (step 12)
    supp = ""
    sd = f"{OUT}/support_distribution.tsv"
    if os.path.exists(f"{OUT}/figures/support_distribution.png") and os.path.exists(sd):
        gcol = "group"  # step 12 emits a per-group column now
        d = defaultdict(lambda: {"n": 0, "one": 0, "mx": 0})
        for r in csv.DictReader(open(sd), delimiter="\t"):
            if r.get("set", "all_reads") != "all_reads":
                continue
            s = int(r["support"]); n = int(r["n_loci"]); t = r[gcol]
            d[t]["n"] += n; d[t]["mx"] = max(d[t]["mx"], s)
            if s == 1: d[t]["one"] += n
        nsing = sum(1 for _ in open(f"{OUT}/singleton_events.tsv")) - 1
        srow = "".join(
            f"<tr><td>{GLAB.get(g, g)}</td><td>{d[g]['n']}</td><td>{d[g]['one']} ({100*d[g]['one']/max(d[g]['n'],1):.0f}%)</td>"
            f"<td>{d[g]['n']-d[g]['one']}</td><td>{d[g]['mx']}×</td></tr>" for g in GROUPS)
        supp = f"""<h2>8. Read-support distribution &amp; the 1× events</h2>
<p>This is <b>not</b> a Sniffles VAF (no allele-frequency model) — it is simply the number of independent reads carrying
the same event (per-read calls clustered into loci within 100 bp). <b>support = 1 means the event is seen in exactly one
read</b> — the exclusively single-molecule (singleton) class, the cleanest somatic candidate. Most loci are singletons in
every group. The high-support tail is <b>coverage-capped</b> (deep leaf reaches &gt;50×, shallow pollen tops out low), so the
right panel read-budget-matches all groups to a common depth for a fair support comparison.
<b>{nsing} singleton (1×) events</b> are written to <code>results/singleton_events.tsv</code>.</p>
<table><tr><th>group</th><th>SV loci</th><th>1× (singleton)</th><th>≥2×</th><th>max support</th></tr>{srow}</table>
{imf('support_distribution', 'Loci binned by supporting-read count (log y); 1× shaded. Left = all reads; right = read-budget matched (all groups downsampled to a common depth/hap).')}"""

    # singleton annotation (step 13) — trustworthiness of the 1x events
    ann = ""
    ap = f"{OUT}/singleton_events_annotated.tsv"
    if os.path.exists(ap):
        al = list(csv.DictReader(open(ap), delimiter="\t"))
        agg = defaultdict(lambda: Counter())
        ireg = defaultdict(int); inarr = defaultdict(int); ndi = defaultdict(int)
        for d in al:
            s = d["sample"]
            agg[s][d["confidence"]] += 1
            if d["svtype"] in ("DEL", "INS"):
                ndi[s] += 1
                ireg[s] += int(d["in_register"]); inarr[s] += int(d["in_cen180_array"] or 0)
        arow = "".join(
            f"<tr><td>{GLAB.get(g, g)}</td><td>{sum(agg[g].values())}</td><td>{agg[g]['HIGH']}</td><td>{agg[g]['MEDIUM']}</td>"
            f"<td>{agg[g]['LOW']}</td><td>{ndi[g]}</td><td>{inarr[g]} ({100*inarr[g]//max(ndi[g],1)}%)</td>"
            f"<td>{ireg[g]} ({100*ireg[g]//max(ndi[g],1)}%)</td></tr>" for g in GROUPS)
        ann = f"""<h2>9. Are the 1× events trustworthy? (per-read TRASH annotation)</h2>
<p>A single read with an SV is either a real somatic molecule or a one-off artifact — support alone cannot tell them apart.
Each 1× event is annotated with orthogonal evidence (<code>results/singleton_events_annotated.tsv</code>): the detector(s)
that found it, the read's divergence (<code>de</code>) and MAPQ, and — the key one — <b>TRASH run on the full read</b>
(the lab's canonical method, <code>analyze_deletions.py</code>). TRASH annotates the CEN178 monomers on the actual molecule;
we then find the monomer immediately <b>left and right of the junction</b>. <b>in_CEN178_array</b> = both flanking monomers
exist (the event is inside a satellite array); <b>in_register</b> = additionally a whole-CEN178-monomer event (|svlen| mod 178
≈ 0) — the unequal-sister-chromatid-HR signature, vs an out-of-phase NHEJ-like junction. Confidence = HIGH/MEDIUM/LOW from the
combination. The key comparison is the <b>in-register (whole-monomer, unequal-HR) fraction in CENH3ox vs WT</b>.</p>
<table><tr><th>group</th><th>1× events</th><th>HIGH</th><th>MEDIUM</th><th>LOW</th><th>DEL/INS</th><th>in CEN178 array</th><th>in-register</th></tr>{arow}</table>
<p class=cap style="font-size:12.5px;color:#666">Sorted, fully annotated list (start with the HIGH rows): <code>results/singleton_events_annotated.tsv</code>.</p>"""

    html = f"""<!doctype html><meta charset=utf-8><title>Single-molecule centromere SVs — WT</title>
<style>body{{font-family:-apple-system,Segoe UI,Arial,sans-serif;max-width:920px;margin:0 auto;padding:26px;color:#1d1d1f;line-height:1.6}}
h1{{border-bottom:3px solid #C0392B;padding-bottom:9px;font-size:24px}}h2{{color:#C0392B;font-size:18px;margin-top:30px}}
figure{{text-align:center;margin:16px 0}}figcaption{{font-size:12.5px;color:#666;font-style:italic}}
table{{border-collapse:collapse;margin:10px 0;font-size:14px}}td,th{{border:1px solid #ddd;padding:4px 10px}}
.box{{background:#EBF5FB;border-left:4px solid #3498DB;padding:10px 15px;margin:12px 0}}</style>
<h1>Single-molecule structural variants in the centromere — WT vs CENH3ox, F1 (Col×Ler)</h1>
<div class=box style="background:#FEF9E7;border-left:4px solid #B7950B"><b>Reference design (important).</b> Each genotype's reads
are mapped to its <b>own centromere baseline</b>, so a call is a somatic deviation from that genotype's assembly — not the
fixed Col-vs-CENH3ox remodelling. The CENH3ox line is <b>CENH3ox-Col × WT-Ler</b>: its Col haplotype carries the remodelled centromere, its Ler haplotype is
<b>wild-type Ler</b>. So WT col → Col-HiFi; <b>CENH3ox col → the CENH3ox line's own assembly (CENH3ox-Col-HiFi)</b>; and
<b>both ler haplotypes → Ler-HiFi, which is the correct baseline for each</b> (the CENH3ox Ler haplotype is genetically WT Ler).
A CENH3ox Ler signal above WT is therefore a <b>trans</b> effect of CENH3 overexpression on the unmodified Ler centromere.
CENH3ox col reads are the winnowmap <code>-ax map-pb --MD</code> alignments to the CENH3ox col-parent assembly
(<code>sv_calling/aligned_matched/</code>), the same reads as the WT-ref col mapping, just against the CENH3ox baseline.</div>
<p>Leaf &amp; pollen, <b>WT and CENH3ox</b>, both haplotypes (col → its genotype's Col assembly, ler → Ler-HiFi; winnowmap). Candidate reads
(de≥0.005 ∨ NM≥50 ∨ SA) were scanned at the <b>single-molecule</b> level by two complementary detectors that
share Sniffles2's own topology classifier (<code>sv.classify_splits</code>, run per-read with no min-support):
(1) <b>CIGAR + split-read leadprov</b> and (2) <b>split-and-map</b> (re-map of substitution-contrast-split fragments).
CENH3 overexpression is expected to destabilize centromeres, so the hypothesis is <b>more single-molecule centromere SVs in
CENH3ox than WT</b>. Because coverage differs, all comparisons are per read-Mb.</p>
<div class=box><b>{n}</b> single-molecule SV calls across 4 groups (haps pooled): {', '.join(f"{GLAB[g]} {by_grp[g]}" for g in GROUPS)} (raw counts; depth differs → see §2).
<b>Per read-Mb</b>: {', '.join(f"{GLAB[g]} {grate[g]:.2f}" for g in GROUPS)} calls/Mb.
By type: {', '.join(f'{t} {by_type[t]}' for t in TYPES if by_type[t])}.
{('CEN178 in-register (whole-monomer indels): %d/%d (%.1f%%).' % (sum(inph), len(inph), 100*sum(inph)/max(len(inph),1))) if inph else ''}</div>
<h2>0. Dataset at a glance</h2>
<p><i>Arabidopsis thaliana</i> F1 hybrid (Col-0 × Ler-0), PacBio HiFi, <b>WT and CENH3ox</b>. Reads are haplotype-split; the
col haplotype is mapped to its <b>own genotype's</b> Col assembly (WT→Col-HiFi, CENH3ox→CENH3ox-Col-HiFi) and the ler haplotype
to Ler-HiFi (see the reference box above). The centromere windows span the CEN178 satellite (178-bp monomer) on the 5
chromosomes. Only reads anchored in those windows are analysed. "Candidate reads" pass the divergence/split gate; each SV call
below comes from one read.</p>
{glance}
<h2>1. Calls by type &amp; group (raw counts)</h2>{im(f1, 'Raw single-molecule SV counts per genotype×tissue (haps pooled). Depth differs between groups (leaf ~500×, pollen ~30×; WT vs CENH3ox also differ), so raw counts are NOT comparable — see §2.')}
<h2>2. Read-Mb-normalized rate — WT vs CENH3ox</h2>
<p>Normalized by Mb of mapped read sequence inside the centromere (sum of aligned bp overlapping the CEN window per primary read),
which removes the depth difference. The comparison of interest is <b>WT vs CENH3ox within each tissue</b>:
<b>CENH3ox leaf {grate['cenh3ox_leaf']:.2f} vs WT leaf {grate['wt_leaf']:.2f} calls/Mb</b>
({grate['cenh3ox_leaf']/max(grate['wt_leaf'],1e-9):.1f}×), and
<b>CENH3ox pollen {grate['cenh3ox_pollen']:.2f} vs WT pollen {grate['wt_pollen']:.2f}</b>
({grate['cenh3ox_pollen']/max(grate['wt_pollen'],1e-9):.1f}×).</p>
{im(fr, 'Calls per Mb of CEN-mapped read sequence, per group. CENH3ox vs WT within each tissue is the key contrast.')}
{rate_tbl}
<p><b>Everything in one plot (the genotype×tissue interaction).</b> Depth-corrected, pollen &gt; leaf in <i>both</i> genotypes
(no true flip), but the <b>tissue bias differs</b>: WT has a strong pollen bias (pollen ≈ 2× leaf), while in CENH3ox the leaf
nearly catches up (pollen ≈ 1.1× leaf). CENH3ox lifts the <b>leaf</b> rate ~4× but the <b>pollen</b> rate only ~2×, i.e. it adds
a large <b>somatic/mitotic</b> instability component (visible in leaf) on top of WT's mostly <b>meiotic</b> (pollen-biased) instability.</p>
{im(finter, 'Per-Mb SV rate, tissue on x, one line per genotype (faint dots = col/ler haplotypes). The WT line rises steeply leaf→pollen; the CENH3ox line is high and nearly flat — the gap between lines (the CENH3ox effect) is largest in leaf.')}
<h2>3. Read-quality controls — is the rate difference an artifact?</h2>
<p>For <b>single-molecule</b> calling each read is an independent sample, so depth does not bias the per-Mb rate. The remaining
ways a group could be inflated are read <i>properties</i>. Per-group medians (haps averaged): CEN read length
{', '.join(f"{GLAB[g]} {qcmp('cen_med_kb')[g]:.0f}kb" for g in GROUPS)}; arm de% (≈ error)
{', '.join(f"{GLAB[g]} {qcmp('arm_de_pct')[g]:.3f}" for g in GROUPS)}; HiFi passes
{', '.join(f"{GLAB[g]} {qcmp('np_med')[g]:.0f}" for g in GROUPS)}. <b>CENH3ox reads are not lower quality than WT</b>
(comparable error rate and passes), so the CENH3ox rate elevation is not a read-quality artifact. The residual caveat
(satellite mapping ambiguity, §11–12) applies equally to all groups.</p>
{im(fq, 'Read-quality controls per group. CENH3ox is comparable to WT in error and passes; length differences are handled by per-Mb normalization.')}
{qc_tbl}
<p class=note style="background:#FDEDEC;border-left:4px solid #C0392B;padding:9px 14px">The one control these cannot rule out is
<b>mapping ambiguity within the satellite array</b> — a perfectly accurate read placed at the wrong monomer copy yields a spurious
split/BND. That depends on array structure, not read quality, and remains the residual caveat.</p>
<h2>4. Size distribution &amp; CEN178 register</h2>{im(f2, 'Indel sizes; dotted lines = 178-bp CEN178 monomer multiples. Peaks on the gridlines = whole-monomer unequal-HR signature.')}
<h2>5. Detection method</h2>{im(f3, 'Which detector found each call. split-and-map adds out-of-phase satellite events missed by CIGAR/split-read.')}
{comp}
{recur}
{supp}
{ann}
{val}
{transloc}
{armc}
{srcbreak}
{insorig}
<h2>Caveat</h2><p>A single ≥50 bp change in deep satellite coverage cannot be fully distinguished from a mapping/sequencing
artifact; the split-and-map re-mapping and the 178-bp register check are the mitigations. Treat single-molecule calls as a
sensitivity ceiling, not a confirmed somatic set.</p>"""
    open(REPORT, "w").write(html)
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
