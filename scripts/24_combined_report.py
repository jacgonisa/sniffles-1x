#!/usr/bin/env python3
"""Combined cross-organism report — all single-molecule SV stats for BOTH Arabidopsis (WT+CENH3ox,
leaf+pollen, Col+Ler) and human (HG002 sperm, MAT+PAT), dissected into CENTROMERE vs ARMS.

- Arabidopsis is centromere-restricted; the CEN-vs-ARM dissection comes from the arm-control steps
  (results/arm_control.tsv = leadprov, results/arm_splitmap_control.tsv = split-and-map): per group ×
  type, calls per million reads in CEN vs ARM + the CEN÷ARM enrichment.
- Human is genome-wide; every call is classified by overlap with the HG002 alpha-satellite CEN bed and
  the ARM bed, then counted per haplotype × type × compartment, with a per-region-Mb density and the
  CEN÷ARM enrichment.
- Insertion-quality QC (homopolymer + CCS-Q) flagged fractions for both.

-> single_molecule_sv/combined_report.html   (standalone; reads results/ and results_human/)."""
import os, csv, base64, io
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv"
A_OUT = f"{ROOT}/results"; H_OUT = f"{ROOT}/results_human"
HBED = "/mnt/ssd-8tb/HUMAN/data/assembly/annotation/cen_arms"
REPORT = f"{ROOT}/combined_report.html"
TYPES = ["DEL", "INS", "DUP", "INV", "BND"]
COL = {"DEL": "#C0392B", "INS": "#2980B9", "DUP": "#27AE60", "INV": "#8E44AD", "BND": "#E67E22"}
AGROUPS = ["wt_leaf", "cenh3ox_leaf", "wt_pollen", "cenh3ox_pollen"]
GLAB = {"wt_leaf": "WT leaf", "cenh3ox_leaf": "CENH3ox leaf", "wt_pollen": "WT pollen",
        "cenh3ox_pollen": "CENH3ox pollen"}
HHAPS = ["MAT", "PAT"]


def load(p):
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def png(fig):
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=130, bbox_inches="tight"); plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()


def im(b, cap):
    return f'<figure><img style="max-width:100%" src="data:image/png;base64,{b}"><figcaption>{cap}</figcaption></figure>'


def load_bed(path):
    iv = defaultdict(list)
    if os.path.exists(path):
        for ln in open(path):
            if ln.strip() and not ln.startswith(("#", "track")):
                c = ln.split("\t"); iv[c[0]].append((int(c[1]), int(c[2])))
    return iv


def bed_mb(path):
    return sum(e - s for ivs in load_bed(path).values() for s, e in ivs) / 1e6


def in_iv(iv, chrom, pos):
    return any(s <= pos < e for s, e in iv.get(chrom, []))


# ---------- Arabidopsis ----------
def arabidopsis():
    rows = load(f"{A_OUT}/sm_sv_calls.tsv")
    cnt = defaultdict(Counter)
    for r in rows:
        cnt[r["sample"]][r["svtype"]] += 1
    ccov = ("<table><tr><th>group</th>" + "".join(f"<th>{t}</th>" for t in TYPES) + "<th>total</th></tr>"
            + "".join("<tr><td><b>" + GLAB[g] + "</b></td>" + "".join(f"<td>{cnt[g][t]}</td>" for t in TYPES)
                      + f"<td><b>{sum(cnt[g].values())}</b></td></tr>" for g in AGROUPS) + "</table>")

    def arm_table(path, label):
        al = load(path)
        if not al:
            return ""
        by = defaultdict(dict)
        for d in al:
            by[d["group"]][d["svtype"]] = d
        rowshtml = ""
        for g in AGROUPS:
            for t in ("DEL", "INS", "DUP", "BND"):
                d = by[g].get(t)
                if not d:
                    continue
                rowshtml += (f"<tr><td>{GLAB[g]}</td><td>{t}</td><td>{d['cen_per_Mreads']}</td>"
                             f"<td>{d['arm_per_Mreads']}</td><td><b>{d['enrichment_CEN_over_ARM']}×</b></td></tr>")
        return (f"<h3>{label}</h3><table><tr><th>group</th><th>type</th><th>CEN /M reads</th>"
                f"<th>ARM /M reads</th><th>CEN÷ARM</th></tr>{rowshtml}</table>")

    n = len(rows)
    return f"""<h2>1. Arabidopsis — WT &amp; CENH3ox, leaf &amp; pollen (centromere-restricted)</h2>
<p><b>{n}</b> single-molecule SV calls in the CEN178 centromeres. Arabidopsis is analysed centromere-restricted, so the
CEN-vs-ARM dissection is a <b>control</b>: the same per-read detectors run on chromosome-arm windows (5 Mb past the CEN,
unique sequence). Rates are calls per million reads in each compartment; CEN÷ARM is the centromere-specific enrichment.</p>
<h3>Calls by type &amp; group (in CEN)</h3>{ccov}
{arm_table(f"{A_OUT}/arm_control.tsv", "CEN vs ARM — CIGAR + native split-read (leadprov)")}
{arm_table(f"{A_OUT}/arm_splitmap_control.tsv", "CEN vs ARM — split-and-remap route")}
<div class=box>DEL/INS/DUP are strongly centromere-enriched (several- to tens-fold); BND is ≈1× (uniform background =
satellite cross-mapping). So the real centromere-instability signal is DEL/INS/DUP, not BND.</div>"""


# ---------- Human ----------
def human():
    rows = load(f"{H_OUT}/sm_sv_calls.tsv")
    cen = {h: load_bed(f"{HBED}/hg002v1.1.{h}.alpha_CEN.bed") for h in HHAPS}
    arm = {h: load_bed(f"{HBED}/hg002v1.1.{h}.all_ARMS.bed") for h in HHAPS}
    cen_mb = {h: bed_mb(f"{HBED}/hg002v1.1.{h}.alpha_CEN.bed") for h in HHAPS}
    arm_mb = {h: bed_mb(f"{HBED}/hg002v1.1.{h}.all_ARMS.bed") for h in HHAPS}
    # classify each call: CEN(alpha) / ARM / other(pericentromere)
    comp = {(h, c, t): Counter() for h in HHAPS for c in ("CEN", "ARM", "other") for t in TYPES}
    tot = defaultdict(Counter)   # (hap, compartment) -> type counts
    for r in rows:
        h, chrom, pos, t = r["hap"], r["chrom"], int(r["pos"]), r["svtype"]
        if r.get("in_cen") == "1" or in_iv(cen[h], chrom, pos):
            k = "CEN"
        elif in_iv(arm[h], chrom, pos):
            k = "ARM"
        else:
            k = "other"
        tot[(h, k)][t] += 1

    def dens(h, k, t, mb):
        return tot[(h, k)][t] / mb if mb else 0

    rowshtml = ""
    for h in HHAPS:
        for t in ("DEL", "INS", "DUP", "BND"):
            cd = dens(h, "CEN", t, cen_mb[h]); ad = dens(h, "ARM", t, arm_mb[h])
            enr = f"{cd/ad:.1f}×" if ad else "∞"
            rowshtml += (f"<tr><td>{h}</td><td>{t}</td><td>{tot[(h,'CEN')][t]}</td><td>{tot[(h,'ARM')][t]}</td>"
                         f"<td>{tot[(h,'other')][t]}</td><td>{cd:.2f}</td><td>{ad:.2f}</td><td><b>{enr}</b></td></tr>")
    nc = sum(sum(tot[(h, 'CEN')].values()) for h in HHAPS)
    na = sum(sum(tot[(h, 'ARM')].values()) for h in HHAPS)
    no = sum(sum(tot[(h, 'other')].values()) for h in HHAPS)
    n = len(rows)
    # figure: CEN vs ARM density by type (MAT+PAT pooled)
    fig, ax = plt.subplots(figsize=(8, 4))
    x = range(len(TYPES)); w = 0.38
    cmb = sum(cen_mb.values()); amb = sum(arm_mb.values())
    cvals = [sum(tot[(h, "CEN")][t] for h in HHAPS) / cmb for t in TYPES]
    avals = [sum(tot[(h, "ARM")][t] for h in HHAPS) / amb for t in TYPES]
    ax.bar([xi - w/2 for xi in x], cvals, w, label=f"alpha-CEN ({cmb:.0f} Mb)", color="#8E44AD")
    ax.bar([xi + w/2 for xi in x], avals, w, label=f"arms ({amb:.0f} Mb)", color="#95a5a6")
    ax.set_yscale("log"); ax.set_xticks(list(x)); ax.set_xticklabels(TYPES)
    ax.set_ylabel("calls per Mb of region (log)"); ax.legend()
    ax.set_title("Human — call density, alpha-satellite CEN vs arms")
    return f"""<h2>2. Human — HG002 sperm, MAT &amp; PAT (genome-wide)</h2>
<p><b>{n}</b> single-molecule SV calls genome-wide; each classified by overlap with the HG002 alpha-satellite CEN bed vs the
ARM bed. <b>{nc} in alpha-CEN, {na} in arms, {no} pericentromere/other.</b> Density = calls per Mb of that region's sequence
(alpha-CEN ≈ {cmb:.0f} Mb, arms ≈ {amb:.0f} Mb), so CEN÷ARM is length-normalized.</p>
{im(png(fig), 'Human call density (per Mb of region, log) — alpha-CEN vs arms, by type. DUP/INS enriched in CEN; BND is where satellite cross-mapping concentrates.')}
<table><tr><th>hap</th><th>type</th><th>CEN n</th><th>ARM n</th><th>other n</th><th>CEN /Mb</th><th>ARM /Mb</th><th>CEN÷ARM</th></tr>{rowshtml}</table>
<div class=box>Only ~{100*nc//max(n,1)}% of genome-wide human calls are in alpha-satellite CEN, but per Mb the CEN is far denser
(DUP/INS especially) — the single-molecule satellite signal stock Sniffles misses. Arms are the unique-sequence background.</div>"""


# ---------- Insertion QC both ----------
def insqc():
    def summ(path, groups, lab):
        ql = load(path)
        if not ql:
            return ""
        agg = defaultdict(lambda: [0, 0, 0, 0])
        keyf = (lambda d: d["sample"]) if groups is AGROUPS else (lambda d: d["hap"])
        for d in ql:
            a = agg[keyf(d)]; a[0] += 1
            a[1] += int(d["low_complexity"]); a[2] += int(d["quality_decay"]); a[3] += (d["verdict"] == "FLAG")
        rr = "".join(f"<tr><td>{GLAB.get(g, g)}</td><td>{agg[g][0]}</td><td>{agg[g][1]}</td><td>{agg[g][2]}</td>"
                     f"<td><b>{agg[g][3]}</b> ({100*agg[g][3]//max(agg[g][0],1)}%)</td></tr>"
                     for g in groups if g in agg)
        return f"<h3>{lab}</h3><table><tr><th>group</th><th>CIGAR INS</th><th>low-cplx</th><th>q-decay</th><th>flagged</th></tr>{rr}</table>"
    return f"""<h2>3. Insertion-sequence quality QC (homopolymer &amp; CCS base-quality) — both organisms</h2>
<p>Every CIGAR insertion's inserted bases + per-base CCS qualities are checked: <b>low-complexity</b> (homopolymer &gt;30% or
entropy &lt;1.2 bits) and <b>quality-decay</b> (insertion ≥5 Q below the ±200 bp flanks). Flagged INS are removed for the
high-confidence set.</p>
{summ(f"{A_OUT}/insertion_qc.tsv", AGROUPS, "Arabidopsis (centromeric insertions, by group)")}
{summ(f"{H_OUT}/insertion_qc.tsv", HHAPS, "Human (genome-wide insertions, by haplotype)")}
{compartment_table()}
{examples_html()}
<div class=box>The <b>CEN-vs-arm contrast is the story.</b> In <b>unique arm sequence almost every</b> single-molecule
insertion is an artefact (<b>99% Arabidopsis, 84% human flagged</b>), whereas in the <b>centromere most are real</b>
(only <b>3% Arabidopsis, 33% human flagged</b>). So a ≥50 bp "insertion" seen in one arm read is essentially always a
homopolymer/STR sequencing error, while centromeric insertions are genuine satellite duplications. This both validates the
CEN signal (esp. the CENH3ox excess) and shows the QC is indispensable outside the centromere.</div>"""


def _fileb64(p):
    return base64.b64encode(open(p, "rb").read()).decode() if os.path.exists(p) else None


MECHLAB = {
    "self_sister_unequal_inreg": ("1", "self/sister unequal exchange, in-register (satellite)", "ours"),
    "self_sister_unequal_offreg": ("1", "self/sister unequal exchange, out-of-register", "ours"),
    "inversion": ("2", "inversion (arm/CEN)", "ours"),
    "ectopic_samehap": ("4", "ectopic (same-hap BND to another chromosome)", "ours"),
    "inter_homolog_crossover": ("5", "inter-homolog crossover / gene conversion", "charla"),
    "inter_homolog_satellite": ("6", "inter-homolog satellite exchange (tandem)", "charla"),
    "ectopic_interhomolog": ("7", "ectopic inter-homolog (non-homologous)", "charla"),
    "artefact_crossmap": ("10", "artefact — satellite cross-mapping BND", "art"),
    "artefact_lowqual": ("11", "artefact — homopolymer / quality INS", "art"),
    "charla_ambiguous": ("12", "artefact — ambiguous / chimeric (CHARLA)", "art"),
}
MCOL = {"ours": "#1E8449", "charla": "#2471A3", "art": "#7F8C8D"}


def _read_mb(path, col):
    d = defaultdict(float)
    for r in load(path):
        d[r["sample"]] += float(r[col])   # sum haplotypes -> per group
    return d


def mechanism_integration_html():
    rows = load(f"{A_OUT}/mechanism_summary.tsv")
    if not rows:
        return ""
    cnt = defaultdict(dict)
    for r in rows:
        cnt[r["mechanism"]][r["group"]] = int(r["count"])
    # denominators: read-Mb ANALYSED by each pipeline. sniffles scanned CEN reads (CEN-mapped-Mb);
    # CHARLA scanned all reads (genome-mapped-Mb). Normalise so leaf/pollen (very different depth) compare.
    cen_mb = _read_mb(f"{A_OUT}/cen_mapped_mb.tsv", "cen_mapped_mb")
    gen_mb = _read_mb(f"{A_OUT}/genome_mapped_mb.tsv", "genome_mapped_mb")  # WT only

    def rate_per_gb(m, g, src):
        n = cnt.get(m, {}).get(g)
        denom = (cen_mb.get(g) if src == "ours" or m == "inter_homolog_satellite" else gen_mb.get(g))
        if n is None or not denom:
            return None
        return n / (denom / 1000.0)     # events per Gb of analysed reads

    order = ["self_sister_unequal_inreg", "self_sister_unequal_offreg", "inversion", "ectopic_samehap",
             "inter_homolog_crossover", "inter_homolog_satellite", "ectopic_interhomolog",
             "artefact_crossmap", "artefact_lowqual", "charla_ambiguous"]
    body = ""
    for m in order:
        if m not in cnt:
            continue
        cls, lab, src = MECHLAB[m]
        denom_lab = "CEN-reads" if (src == "ours" or m == "inter_homolog_satellite") else "genome-reads"
        cells = ""
        for g in AGROUPS:
            r = rate_per_gb(m, g, src); raw = cnt[m].get(g)
            cells += f"<td>{('%.2f' % r) if r is not None else '·'}<br><span style='color:#999;font-size:10px'>{raw if raw is not None else ''}</span></td>"
        body += (f"<tr><td style='color:{MCOL[src]}'><b>{cls}</b></td><td style='color:{MCOL[src]}'>{lab}</td>"
                 f"<td style='font-size:11px'>{denom_lab}</td>{cells}</tr>")
    tbl = ("<table><tr><th>class</th><th>mechanism</th><th>per Gb of…</th>"
           + "".join(f"<th>{GLAB[g]}</th>" for g in AGROUPS) + "</tr>" + body
           + "</table><p class=cap style='font-size:11px;color:#999'>each cell = rate (events per Gb of analysed reads); "
             "small grey number = raw count. · = N/A (CHARLA is WT-only).</p>")

    # figure: two normalised comparisons, WT leaf vs pollen (each in its own denominator)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (m, ttl, col) in zip(axes, [
            ("self_sister_unequal_inreg", "SELF/SISTER satellite unequal exchange\n(sniffles_1x · per Gb CEN-reads)", "#1E8449"),
            ("inter_homolog_crossover", "INTER-HOMOLOG crossover\n(CHARLA · per Gb genome-reads)", "#2471A3")]):
        src = MECHLAB[m][2]
        vals = [rate_per_gb(m, g, src) or 0 for g in ("wt_leaf", "wt_pollen")]
        ax.bar([0, 1], vals, 0.55, color=col)
        for xi, v in zip([0, 1], vals):
            ax.text(xi, v, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
        fold = vals[1] / vals[0] if vals[0] else float("inf")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["WT leaf", "WT pollen"], fontsize=11)
        ax.set_title(ttl, fontsize=9.5); ax.set_ylabel("events per Gb of reads")
        ax.set_ylim(0, max(vals) * 1.25 if max(vals) else 1)
        ax.text(0.5, 0.92, f"pollen/leaf = {fold:.0f}×", transform=ax.transAxes, ha="center",
                fontsize=10, color=col, fontweight="bold")
    fig.suptitle("Depth-normalised: satellite unequal exchange ≈ constitutive; inter-homolog crossover is meiotic (pollen)",
                 fontsize=10.5, y=1.06)
    fig.subplots_adjust(top=0.80, wspace=0.28)
    sat_l = rate_per_gb("self_sister_unequal_inreg", "wt_leaf", "ours"); sat_p = rate_per_gb("self_sister_unequal_inreg", "wt_pollen", "ours")
    cox_l = rate_per_gb("inter_homolog_crossover", "wt_leaf", "charla"); cox_p = rate_per_gb("inter_homolog_crossover", "wt_pollen", "charla")
    return (f"<h2>0b. Integrated mechanism RATES — charla_hifi + sniffles_1x (depth-normalised)</h2>"
            f"<p>Raw counts conflate biology with depth (WT leaf is ~16× deeper in the CEN, ~7× genome-wide, than pollen), so every "
            f"mechanism is normalised by the <b>Mb of reads each pipeline actually analysed</b>: sniffles_1x per <b>CEN-mapped-Mb</b> "
            f"(it scanned CEN reads only), CHARLA per <b>genome-mapped-Mb</b> (it scanned all reads). After normalisation the two "
            f"mechanisms separate cleanly: <b>self/sister satellite unequal exchange is ~constitutive</b> "
            f"(WT leaf {sat_l:.2f} → pollen {sat_p:.2f} per Gb, ~{sat_p/max(sat_l,1e-9):.0f}×), whereas <b>inter-homolog crossover is "
            f"strongly meiotic</b> (WT leaf {cox_l:.2f} → pollen {cox_p:.1f} per Gb, ~{cox_p/max(cox_l,1e-9):.0f}× — the depth-corrected "
            f"pollen enrichment, far larger than the raw ~37× because pollen is much shallower).</p>"
            f"{im(png(fig), 'Depth-normalised rates (per Gb of analysed reads). Satellite unequal exchange (left) is similar in leaf and pollen; inter-homolog crossover (right) is ~270× higher in pollen — the meiotic signal, previously masked by leaf’s greater depth.')}"
            f"{tbl}<p class=cap style='font-size:12px;color:#777'>sniffles_1x satellite events are inherently non-allelic self/sister "
            f"(intra-chromatid vs inter-sister indistinguishable). CHARLA is WT-only (no CENH3ox hybrid-read run).</p>")


def mechanism_html():
    b = _fileb64(f"{ROOT}/docs/mechanism_taxonomy.png")
    if not b:
        return ""
    return (f"<h2>0. Mechanism taxonomy — every kind of event, and where it is captured</h2>"
            f"<p>An event joins two genomic loci in one read. Two axes: <b>partner molecule</b> "
            f"(self/sister · homolog · non-homolog) and <b>allelic vs non-allelic</b> (same locus vs offset), with a "
            f"<b>compartment</b> overlay (CEN178 satellite · pericentromere · unique arm). Colour = where each class is "
            f"captured: <b>our non-hybrid SV pipeline</b> (green: self/sister unequal exchange + ectopic), "
            f"<b>CHARLA hybrid/crossover analysis</b> (blue: inter-homolog), a <b>current gap</b> (red — notably an "
            f"insertion whose donor is templated from the other homolog, and inter-homolog SVs on reads that strict-90 "
            f"discards), or a <b>not-a-real-event artefact</b> (grey).</p>"
            f'<figure><img style="max-width:100%" src="data:image/png;base64,{b}">'
            f"<figcaption>Exhaustive taxonomy. intra-chromatid vs inter-sister are sequence-identical so they are one "
            f"class; inter-homolog lives in CHARLA (hybrid reads); the red gaps are the next things to build.</figcaption></figure>")


def examples_html():
    import glob as _g
    def one(pattern, cap):
        hits = sorted(_g.glob(pattern))
        if not hits:
            return ""
        b = _fileb64(hits[0])
        return im(b, cap) if b else ""
    # one quality-decay + one homopolymer per organism
    a_qd = one(f"{A_OUT}/artefact_examples/arabidopsis_quality_decay_*.png",
               "Arabidopsis — QUALITY DECAY: the inserted bases (shaded) collapse to ~Q7 while the flanks are ~Q40; base strip below.")
    a_hp = one(f"{A_OUT}/artefact_examples/arabidopsis_homopolymer_*.png",
               "Arabidopsis — HOMOPOLYMER: a ~55 bp insertion that is ~98% a single base (poly-G, entropy ≈0.1 bits) — a classic HiFi homopolymer artefact.")
    h_qd = one(f"{H_OUT}/artefact_examples/human_quality_decay_*.png",
               "Human — QUALITY DECAY: CCS quality drops inside the insertion vs the flanks.")
    h_hp = one(f"{H_OUT}/artefact_examples/human_homopolymer_*.png",
               "Human — HOMOPOLYMER / low-complexity inserted tract.")
    if not any((a_qd, a_hp, h_qd, h_hp)):
        return ""
    return (f"<h3>Example flagged insertions (what the artefacts look like)</h3>"
            f"<p>Top panel = per-base <b>CCS quality</b> along the read (insertion shaded red; flank vs insertion mean-Q "
            f"lines). Bottom = the <b>inserted sequence</b> as a coloured base strip, homopolymer runs (≥5 bp) underlined.</p>"
            f"{a_qd}{a_hp}{h_qd}{h_hp}")


def compartment_table():
    p = f"{A_OUT}/insertion_qc_by_compartment.tsv"
    rows = load(p)
    if not rows:
        return ""
    tr = "".join(
        f"<tr><td>{d['organism']}</td><td>{d['compartment']}</td><td>{d['n_ins']}</td>"
        f"<td>{d['low_complexity']}</td><td>{d['quality_decay']}</td>"
        f"<td><b>{d['flagged']}</b> ({d['pct_flagged']}%)</td></tr>" for d in rows)
    return ("<h3>Flagged insertions dissected: CEN vs arms (the key 2×2)</h3>"
            "<table><tr><th>organism</th><th>compartment</th><th>insertions</th><th>low-complexity</th>"
            f"<th>quality-decay</th><th>flagged</th></tr>{tr}</table>"
            "<p class=cap style='font-size:12px;color:#777'>Arabidopsis ARM = CIGAR I≥50 bp scanned in the distal-arm "
            "windows (unique sequence); human CEN/ARM = genome-wide insertions classified by the HG002 alpha-CEN / ARM beds.</p>")


def main():
    html = f"""<!doctype html><meta charset=utf-8><title>Single-molecule SV — Arabidopsis + Human, CEN vs arms</title>
<style>body{{font-family:-apple-system,Segoe UI,Arial,sans-serif;max-width:950px;margin:0 auto;padding:26px;color:#1d1d1f;line-height:1.6}}
h1{{border-bottom:3px solid #16A085;padding-bottom:9px;font-size:23px}}h2{{color:#138D75;font-size:19px;margin-top:30px}}
h3{{font-size:15px;color:#333;margin-top:18px}}figure{{text-align:center;margin:14px 0}}figcaption{{font-size:12.5px;color:#666;font-style:italic}}
table{{border-collapse:collapse;margin:8px 0;font-size:13.5px}}td,th{{border:1px solid #ddd;padding:3px 9px}}
.box{{background:#E8F8F5;border-left:4px solid #16A085;padding:10px 15px;margin:12px 0}}</style>
<h1>Single-molecule structural variants — Arabidopsis &amp; Human, dissected into centromere vs arms</h1>
<p>The cracked Sniffles2 per-read caller (candidate → CIGAR + split-read + split-and-remap → unmodified
<code>sv.classify_splits</code>, no clustering) applied to two datasets. Arabidopsis: WT &amp; CENH3ox, leaf &amp; pollen,
Col &amp; Ler, centromere-restricted. Human: HG002 sperm (BLS0005+BLS0006), MAT &amp; PAT, genome-wide. Full per-dataset
reports: <code>report.html</code> (Arabidopsis) and <code>results_human/report_human.html</code> (human).</p>
{mechanism_html()}
{mechanism_integration_html()}
{arabidopsis()}
{human()}
{insqc()}
<p class=cap style="font-size:12px;color:#777;margin-top:24px">Arabidopsis arm rates = per-million-reads in CEN vs 5-Mb-distal
arm windows (control steps 16/18). Human CEN/arm = per-Mb-of-region density from bed overlap. Different normalizations
(Arabidopsis is CEN-only so has no genome-wide callset), both giving the CEN÷ARM enrichment.</p>"""
    open(REPORT, "w").write(html)
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
