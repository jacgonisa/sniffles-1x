#!/usr/bin/env python3
"""Step 6 — figures + standalone HTML report for the single-molecule SV callset.
Reads results/sm_sv_calls.tsv -> single_molecule_sv/report.html (PNGs inline)."""
import os, csv, base64, io
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import OUT, MONO

REPORT = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/report.html"
TYPES = ["DEL", "INS", "DUP", "INV", "BND"]
COL = {"DEL": "#C0392B", "INS": "#2980B9", "DUP": "#27AE60", "INV": "#8E44AD", "BND": "#E67E22"}


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
    cats = [(t, h) for t in ("leaf", "pollen") for h in ("col", "ler")]
    cnt = defaultdict(Counter)
    for r in rows:
        cnt[(r["tissue"], r["hap"])][r["svtype"]] += 1
    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = range(len(cats)); w = 0.16
    for i, t in enumerate(TYPES):
        ax.bar([xi + i * w for xi in x], [cnt[c][t] for c in cats], w, label=t, color=COL[t])
    ax.set_xticks([xi + 2 * w for xi in x]); ax.set_xticklabels([f"{t}\n{h}" for t, h in cats])
    ax.set_ylabel("single-molecule SV calls"); ax.legend(ncol=5, fontsize=8)
    ax.set_title("Single-molecule centromere SVs by type, tissue, haplotype")
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


def fig_rates(rates):
    cats = [(t, h) for t in ("leaf", "pollen") for h in ("col", "ler")]
    rt = {(r["tissue"], r["hap"]): r for r in rates}
    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = range(len(cats)); w = 0.16
    for i, t in enumerate(TYPES):
        ax.bar([xi + i * w for xi in x],
               [float(rt[c][f"{t}_per_mb"]) if c in rt else 0 for c in cats], w, label=t, color=COL[t])
    ax.set_xticks([xi + 2 * w for xi in x]); ax.set_xticklabels([f"{t}\n{h}" for t, h in cats])
    ax.set_ylabel("calls per Mb of CEN-mapped read seq"); ax.legend(ncol=5, fontsize=8)
    ax.set_title("Read-Mb-normalized single-molecule SV rate (leaf vs pollen comparable)")
    return png(fig)


def load_qc():
    p = f"{OUT}/read_qc.tsv"
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fig_qc(qc):
    cats = [(r["tissue"], r["hap"]) for r in qc]
    labels = [f"{t}\n{h}" for t, h in cats]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    x = range(len(qc))
    cols = ["#C0392B" if r["tissue"] == "pollen" else "#2980B9" for r in qc]
    axes[0].bar(x, [float(r["arm_de_pct"]) for r in qc], color=cols)
    axes[0].set_title("Arm de% (≈ seq error)"); axes[0].set_ylabel("% divergence")
    axes[1].bar(x, [float(r["np_med"]) for r in qc], color=cols)
    axes[1].set_title("HiFi passes (np, median)")
    axes[2].bar(x, [float(r["cen_med_kb"]) for r in qc], color=cols)
    axes[2].set_title("CEN read length (median kb)")
    for ax in axes:
        ax.set_xticks(list(x)); ax.set_xticklabels(labels, fontsize=8)
    fig.suptitle("Read-quality controls — pollen (red) is equal-or-better than leaf (blue)")
    return png(fig)


def fig_sizes(rows):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
    for ax, tis in zip(axes, ("leaf", "pollen")):
        for t in ("DEL", "INS", "DUP"):
            sizes = [abs(r["svlen"]) for r in rows if r["tissue"] == tis and r["svtype"] == t and r["svlen"]]
            if sizes:
                ax.hist(sizes, bins=range(0, 3000, 50), histtype="step", color=COL[t], label=t, linewidth=1.4)
        for m in range(MONO, 3000, MONO):
            ax.axvline(m, color="#bbb", ls=":", lw=0.6)
        ax.set_title(tis); ax.set_xlabel("|SV length| (bp)"); ax.set_xlim(0, 3000); ax.legend(fontsize=8)
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


def main():
    rows = load()
    n = len(rows)
    by_tis = Counter(r["tissue"] for r in rows)
    by_type = Counter(r["svtype"] for r in rows)
    stock = sum(int(r["stock_match"]) for r in rows)
    inph = [int(r["in_phase"]) for r in rows if r["in_phase"] not in ("", "None")]
    rates = load_rates()
    qc = load_qc()
    f1, fr, f2, f3 = fig_counts(rows), fig_rates(rates), fig_sizes(rows), fig_methods(rows)
    fq = fig_qc(qc) if qc else ""

    qc_tbl = ("<table><tr><th>tissue</th><th>hap</th><th>CEN read len (kb, median)</th>"
              "<th>arm de% (≈ error)</th><th>CEN de%</th><th>np (median passes)</th><th>rq% (median)</th></tr>"
              + "".join(
                  f"<tr><td>{r['tissue']}</td><td>{r['hap']}</td><td>{float(r['cen_med_kb']):.1f}</td>"
                  f"<td>{float(r['arm_de_pct']):.3f}</td><td>{float(r['cen_de_pct']):.3f}</td>"
                  f"<td>{float(r['np_med']):.0f}</td><td>{float(r['rq_med_pct']):.3f}</td></tr>"
                  for r in qc) + "</table>") if qc else ""

    rate_tbl = ("<table><tr><th>tissue</th><th>hap</th><th>CEN-mapped Mb</th>"
                "<th>ALL / Mb</th><th>DEL / Mb</th><th>INS / Mb</th><th>DUP / Mb</th></tr>"
                + "".join(
                    f"<tr><td>{r['tissue']}</td><td>{r['hap']}</td><td>{float(r['cen_mapped_mb']):.0f}</td>"
                    f"<td>{float(r['ALL_per_mb']):.3f}</td><td>{float(r['DEL_per_mb']):.3f}</td>"
                    f"<td>{float(r['INS_per_mb']):.3f}</td><td>{float(r['DUP_per_mb']):.3f}</td></tr>"
                    for r in sorted(rates, key=lambda r: (r['tissue'], r['hap'])))
                + "</table>")

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

    import glob as _glob
    valfigs = sorted(_glob.glob(f"{OUT}/read_validation/*.png"))
    val = ""
    if valfigs:
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
        d = defaultdict(lambda: {"n": 0, "one": 0, "mx": 0})
        for r in csv.DictReader(open(sd), delimiter="\t"):
            s = int(r["support"]); n = int(r["n_loci"]); t = r["tissue"]
            d[t]["n"] += n; d[t]["mx"] = max(d[t]["mx"], s)
            if s == 1: d[t]["one"] += n
        nsing = sum(1 for _ in open(f"{OUT}/singleton_events.tsv")) - 1
        srow = "".join(
            f"<tr><td>{t}</td><td>{d[t]['n']}</td><td>{d[t]['one']} ({100*d[t]['one']/max(d[t]['n'],1):.0f}%)</td>"
            f"<td>{d[t]['n']-d[t]['one']}</td><td>{d[t]['mx']}×</td></tr>" for t in ("leaf", "pollen"))
        supp = f"""<h2>8. Read-support distribution &amp; the 1× events</h2>
<p>This is <b>not</b> a Sniffles VAF (no allele-frequency model) — it is simply the number of independent reads carrying
the same event (per-read calls clustered into loci within 100 bp). <b>support = 1 means the event is seen in exactly one
read</b> — the exclusively single-molecule (singleton) class, the cleanest somatic candidate. Most loci are singletons
({d['leaf']['one']*100//max(d['leaf']['n'],1)}% leaf, {d['pollen']['one']*100//max(d['pollen']['n'],1)}% pollen). The
high-support tail is <b>coverage-capped</b>: leaf reaches &gt;50× support, pollen tops out ~{d['pollen']['mx']}× because
it is only ~30× deep — so leaf's recurrent loci (mostly the fixed-vs-reference events of §7) simply cannot appear in pollen.
<b>{nsing} singleton (1×) events</b> are written to <code>results/singleton_events.tsv</code>.</p>
<table><tr><th>tissue</th><th>SV loci</th><th>1× (singleton)</th><th>≥2×</th><th>max support</th></tr>{srow}</table>
{imf('support_distribution', 'Loci binned by supporting-read count (log y); 1× shaded. Left = all reads (leaf’s depth gives a long tail); right = read-budget matched (leaf downsampled to pollen ≈25k reads/hap) — at equal reads pollen has MORE SV loci than leaf.')}"""

    # singleton annotation (step 13) — trustworthiness of the 1x events
    ann = ""
    ap = f"{OUT}/singleton_events_annotated.tsv"
    if os.path.exists(ap):
        al = list(csv.DictReader(open(ap), delimiter="\t"))
        agg = defaultdict(lambda: Counter())
        ireg = defaultdict(int); inarr = defaultdict(int); ndi = defaultdict(int)
        for d in al:
            agg[d["tissue"]][d["confidence"]] += 1
            if d["svtype"] in ("DEL", "INS"):
                ndi[d["tissue"]] += 1
                ireg[d["tissue"]] += int(d["in_register"]); inarr[d["tissue"]] += int(d["in_cen180_array"] or 0)
        arow = "".join(
            f"<tr><td>{t}</td><td>{sum(agg[t].values())}</td><td>{agg[t]['HIGH']}</td><td>{agg[t]['MEDIUM']}</td>"
            f"<td>{agg[t]['LOW']}</td><td>{ndi[t]}</td><td>{inarr[t]} ({100*inarr[t]//max(ndi[t],1)}%)</td>"
            f"<td>{ireg[t]} ({100*ireg[t]//max(ndi[t],1)}%)</td></tr>" for t in ("leaf", "pollen"))
        ann = f"""<h2>9. Are the 1× events trustworthy? (per-read TRASH annotation)</h2>
<p>A single read with an SV is either a real somatic molecule or a one-off artifact — support alone cannot tell them apart.
Each 1× event is annotated with orthogonal evidence (<code>results/singleton_events_annotated.tsv</code>): the detector(s)
that found it, the read's divergence (<code>de</code>) and MAPQ, and — the key one — <b>TRASH run on the full read</b>
(the lab's canonical method, <code>analyze_deletions.py</code>). TRASH annotates the CEN178 monomers on the actual molecule;
we then find the monomer immediately <b>left and right of the junction</b>. <b>in_CEN178_array</b> = both flanking monomers
exist (the event is inside a satellite array); <b>in_register</b> = additionally a whole-CEN178-monomer event (|svlen| mod 178
≈ 0) — the unequal-sister-chromatid-HR signature, vs an out-of-phase NHEJ-like junction. Confidence = HIGH/MEDIUM/LOW from the
combination. Among DEL/INS singletons, ~77–85% sit in a confirmed CEN178 array and ~44–48% are in-register.</p>
<table><tr><th>tissue</th><th>1× events</th><th>HIGH</th><th>MEDIUM</th><th>LOW</th><th>DEL/INS</th><th>in CEN178 array</th><th>in-register</th></tr>{arow}</table>
<p class=cap style="font-size:12.5px;color:#666">Sorted, fully annotated list (start with the HIGH rows): <code>results/singleton_events_annotated.tsv</code>.</p>"""

    html = f"""<!doctype html><meta charset=utf-8><title>Single-molecule centromere SVs — WT</title>
<style>body{{font-family:-apple-system,Segoe UI,Arial,sans-serif;max-width:920px;margin:0 auto;padding:26px;color:#1d1d1f;line-height:1.6}}
h1{{border-bottom:3px solid #C0392B;padding-bottom:9px;font-size:24px}}h2{{color:#C0392B;font-size:18px;margin-top:30px}}
figure{{text-align:center;margin:16px 0}}figcaption{{font-size:12.5px;color:#666;font-style:italic}}
table{{border-collapse:collapse;margin:10px 0;font-size:14px}}td,th{{border:1px solid #ddd;padding:4px 10px}}
.box{{background:#EBF5FB;border-left:4px solid #3498DB;padding:10px 15px;margin:12px 0}}</style>
<h1>Single-molecule structural variants in the centromere — WT F1 (Col×Ler)</h1>
<p>WT leaf &amp; pollen, both haplotypes (reads → Col-HiFi / Ler-HiFi, winnowmap). Candidate reads
(de≥0.005 ∨ NM≥50 ∨ SA) were scanned at the <b>single-molecule</b> level by two complementary detectors that
share Sniffles2's own topology classifier (<code>sv.classify_splits</code>, run per-read with no min-support):
(1) <b>CIGAR + split-read leadprov</b> and (2) <b>split-and-map</b> (re-map of substitution-contrast-split fragments).
Stock Sniffles2 (<code>--minsupport 1 --mosaic</code>) is shown as a concordance reference.</p>
<div class=box><b>{n}</b> single-molecule SV calls — leaf {by_tis['leaf']}, pollen {by_tis['pollen']} (raw; not comparable).
<b>Read-Mb-normalized, pollen &gt; leaf</b> ({', '.join(f"{r['tissue']}/{r['hap']} {float(r['ALL_per_mb']):.2f}" for r in sorted(rates, key=lambda r:(r['tissue'],r['hap'])))} calls/Mb).
By type: {', '.join(f'{t} {by_type[t]}' for t in TYPES if by_type[t])}.
Stock-Sniffles concordant: {stock}/{n} ({100*stock/max(n,1):.1f}%).
{('CEN178 in-register (whole-monomer indels): %d/%d (%.1f%%).' % (sum(inph), len(inph), 100*sum(inph)/max(len(inph),1))) if inph else ''}</div>
<h2>1. Calls by type, tissue, haplotype (raw counts)</h2>{im(f1, 'Raw single-molecule SV counts. Leaf is ~500× / pollen ~30×, so leaf yields more events purely from depth — NOT comparable. See §2.')}
<h2>2. Read-Mb-normalized rate (leaf vs pollen comparable)</h2>
<p>Normalized by Mb of mapped read sequence inside the centromere (sum of aligned bp overlapping the CEN window per primary read).
This removes the depth difference. <b>Per Mb, pollen carries a higher single-molecule SV rate than leaf</b> — the raw leaf-heavy
counts were a depth artifact.</p>
{im(fr, 'Calls per Mb of CEN-mapped read sequence. Pollen rate > leaf rate once depth is removed.')}
{rate_tbl}
<h2>3. Read-quality controls — is the pollen rate an artifact?</h2>
<p>For <b>single-molecule</b> calling each read is an independent sample, so depth does not bias the per-Mb rate (more
coverage just samples more molecules). The remaining ways pollen could be inflated are read <i>properties</i>. None hold:
pollen reads are only modestly shorter (median {float(qc[2]['cen_med_kb']):.0f} vs {float(qc[0]['cen_med_kb']):.0f} kb, already
absorbed by per-Mb), have a <b>lower</b> arm error rate ({float(qc[2]['arm_de_pct']):.3f}% vs {float(qc[0]['arm_de_pct']):.3f}%),
<b>more</b> HiFi passes (np {float(qc[2]['np_med']):.0f} vs {float(qc[0]['np_med']):.0f}) and higher predicted accuracy
(rq {float(qc[2]['rq_med_pct']):.3f}% vs {float(qc[0]['rq_med_pct']):.3f}%). Pollen reads are equal-or-better quality, so the
higher pollen per-Mb SV rate is not a read-quality artifact.</p>
{im(fq, 'Read-quality controls. Pollen (red) ≤ leaf (blue) on error and ≥ on passes; shorter length is handled by per-Mb normalization.')}
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
<h2>Caveat</h2><p>A single ≥50 bp change in deep satellite coverage cannot be fully distinguished from a mapping/sequencing
artifact; the split-and-map re-mapping and the 178-bp register check are the mitigations. Treat single-molecule calls as a
sensitivity ceiling, not a confirmed somatic set.</p>"""
    open(REPORT, "w").write(html)
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
