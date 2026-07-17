#!/usr/bin/env python3
"""Lean human report (SM_GENOME=human) — the cracked single-molecule caller run genome-wide on
HG002 sperm (BLS0005+BLS0006 merged), MAT + PAT. Grouped by haplotype. Sections: counts & size by
type/haplotype, per-mapped-Mb rate, alpha-CEN vs rest overlap, detection route, and the insertion
quality QC (homopolymer + CCS-Q contrast) with the high-confidence filtered counts.
Run: SM_GENOME=human python 06_report_human.py -> results_human/report_human.html"""
import os, csv, base64, io
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import OUT, HAPS

REPORT = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/results_human/report_human.html"
TYPES = ["DEL", "INS", "DUP", "INV", "BND"]
COL = {"DEL": "#C0392B", "INS": "#2980B9", "DUP": "#27AE60", "INV": "#8E44AD", "BND": "#E67E22"}
HAPCOL = {"MAT": "#C0392B", "PAT": "#2980B9"}


def png(fig):
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=130, bbox_inches="tight"); plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()


def im(b, cap):
    return f'<figure><img style="max-width:100%" src="data:image/png;base64,{b}"><figcaption>{cap}</figcaption></figure>' if b else ""


def load(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []


def fig_counts(rows):
    cnt = defaultdict(Counter)
    for r in rows:
        cnt[r["hap"]][r["svtype"]] += 1
    fig, ax = plt.subplots(figsize=(7.5, 4))
    x = range(len(HAPS)); w = 0.16
    for i, t in enumerate(TYPES):
        ax.bar([xi + i * w for xi in x], [cnt[h][t] for h in HAPS], w, label=t, color=COL[t])
    ax.set_xticks([xi + 2 * w for xi in x]); ax.set_xticklabels([f"{h}" for h in HAPS])
    ax.set_ylabel("single-molecule SV calls"); ax.legend(ncol=5, fontsize=8)
    ax.set_title("Human sperm single-molecule SVs by type & haplotype (genome-wide)")
    return png(fig)


def fig_sizes(rows):
    fig, ax = plt.subplots(figsize=(8, 4))
    for t in ("DEL", "INS", "DUP"):
        s = [abs(int(r["svlen"])) for r in rows if r["svtype"] == t and r["svlen"] not in ("", "None") and abs(int(r["svlen"])) <= 5000]
        if s:
            ax.hist(s, bins=range(0, 5000, 100), histtype="step", color=COL[t], label=t, lw=1.4)
    ax.set_xlabel("|SV length| (bp)"); ax.set_ylabel("count"); ax.legend()
    ax.set_title("Size distribution (DEL/INS/DUP, ≤5 kb)")
    return png(fig)


def main():
    rows = load(f"{OUT}/sm_sv_calls.tsv")
    n = len(rows)
    by_hap = Counter(r["hap"] for r in rows)
    by_type = Counter(r["svtype"] for r in rows)
    incen = Counter(r["hap"] for r in rows if r.get("in_cen") == "1")

    # per-mapped-Mb rate (genome-wide) from 07
    rates = load(f"{OUT}/sm_sv_rates.tsv")
    rate_tbl = ""
    if rates:
        rate_tbl = ("<table><tr><th>hap</th><th>mapped Mb</th><th>ALL/Mb</th><th>DEL/Mb</th><th>INS/Mb</th><th>DUP/Mb</th></tr>"
                    + "".join(f"<tr><td>{r['hap']}</td><td>{float(r['cen_mapped_mb']):.0f}</td>"
                              f"<td>{float(r['ALL_per_mb']):.3f}</td><td>{float(r['DEL_per_mb']):.3f}</td>"
                              f"<td>{float(r['INS_per_mb']):.3f}</td><td>{float(r['DUP_per_mb']):.3f}</td></tr>"
                              for r in rates) + "</table>")

    # alpha-CEN vs rest
    cen_tbl = ("<table><tr><th>hap</th><th>total calls</th><th>in alpha-CEN</th><th>rest of genome</th></tr>"
               + "".join(f"<tr><td>{h}</td><td>{by_hap[h]}</td><td>{incen[h]} ({100*incen[h]//max(by_hap[h],1)}%)</td>"
                         f"<td>{by_hap[h]-incen[h]}</td></tr>" for h in HAPS) + "</table>")

    # detection route (17)
    sb = load(f"{OUT}/source_breakdown.tsv")
    src_tbl = ""
    if sb:
        src_tbl = ("<table><tr><th>route</th><th>calls</th><th>DEL</th><th>INS</th><th>DUP</th><th>INV</th><th>BND</th></tr>"
                   + "".join(f"<tr><td>{d['methods']}</td><td>{d['calls']}</td><td>{d['DEL']}</td><td>{d['INS']}</td>"
                             f"<td>{d['DUP']}</td><td>{d['INV']}</td><td>{d['BND']}</td></tr>" for d in sb) + "</table>")

    # insertion QC (23)
    ql = load(f"{OUT}/insertion_qc.tsv")
    insqc = ""
    if ql:
        agg = defaultdict(lambda: [0, 0, 0, 0])
        for d in ql:
            a = agg[d["hap"]]; a[0] += 1
            a[1] += int(d["low_complexity"]); a[2] += int(d["quality_decay"]); a[3] += (d["verdict"] == "FLAG")
        qrow = "".join(f"<tr><td>{h}</td><td>{agg[h][0]}</td><td>{agg[h][1]}</td><td>{agg[h][2]}</td>"
                       f"<td><b>{agg[h][3]}</b> ({100*agg[h][3]//max(agg[h][0],1)}%)</td></tr>" for h in HAPS if h in agg)
        nflag = sum(a[3] for a in agg.values()); nins = sum(a[0] for a in agg.values())
        nhi = sum(1 for _ in open(f"{OUT}/sm_sv_calls_hiconf.tsv")) - 1 if os.path.exists(f"{OUT}/sm_sv_calls_hiconf.tsv") else 0
        insqc = f"""<h2>5. Insertion-sequence quality QC (homopolymer &amp; CCS base-quality)</h2>
<p>For every CIGAR insertion we cut out the inserted bases + their per-base CCS qualities and flag
<b>low-complexity</b> (homopolymer &gt;30% or entropy &lt;1.2 bits) and <b>quality-decay</b> (insertion ≥5 Q below the
±200 bp flanks). <b>{nflag}/{nins} INS flagged</b>; removed for the high-confidence set ({nhi} calls,
<code>results_human/sm_sv_calls_hiconf.tsv</code>).</p>
<table><tr><th>hap</th><th>CIGAR INS</th><th>low-complexity</th><th>quality-decay</th><th>flagged</th></tr>{qrow}</table>"""

    f1, f2 = fig_counts(rows), fig_sizes(rows)
    html = f"""<!doctype html><meta charset=utf-8><title>Human single-molecule SVs — HG002 sperm</title>
<style>body{{font-family:-apple-system,Segoe UI,Arial,sans-serif;max-width:900px;margin:0 auto;padding:26px;color:#1d1d1f;line-height:1.6}}
h1{{border-bottom:3px solid #2980B9;padding-bottom:9px;font-size:23px}}h2{{color:#2471A3;font-size:18px;margin-top:28px}}
figure{{text-align:center;margin:14px 0}}figcaption{{font-size:12.5px;color:#666;font-style:italic}}
table{{border-collapse:collapse;margin:10px 0;font-size:14px}}td,th{{border:1px solid #ddd;padding:4px 10px}}
.box{{background:#EBF5FB;border-left:4px solid #3498DB;padding:10px 15px;margin:12px 0}}</style>
<h1>Single-molecule structural variants — HG002 sperm (genome-wide)</h1>
<p>The cracked Sniffles2 per-read caller (candidate → <b>CIGAR + split-read (SA) leadprov</b> → unmodified
<code>sv.classify_splits</code>, no clustering) run <b>genome-wide</b> on merged BLS0005+BLS0006 sperm HiFi, haplotype-split
(strict-90) MAT + PAT, each mapped to its own HG002 assembly. <i>Note:</i> the split-and-remap detector is deferred for
human — re-mapping fragments requires re-loading the 6 Gb genome per call (impractical); a pre-built minimap2 index would
enable it. So these calls are the CIGAR + native split-read routes only.</p>
<div class=box><b>{n}</b> single-molecule SV calls — MAT {by_hap['MAT']}, PAT {by_hap['PAT']}.
By type: {', '.join(f'{t} {by_type[t]}' for t in TYPES if by_type[t])}.
In alpha-satellite CEN: {sum(incen.values())} ({100*sum(incen.values())//max(n,1)}%).
(Stock Sniffles2 --minsupport 1 on the same data called only ~96 MAT / 29 PAT — see the existing report.)</div>
<h2>1. Calls by type &amp; haplotype</h2>{im(f1, 'Raw single-molecule SV counts per haplotype, genome-wide.')}
<h2>2. Per-mapped-Mb rate</h2>{rate_tbl}
<h2>3. Alpha-satellite CEN vs rest of genome</h2>
<p>Genome-wide calls annotated by overlap with the HG002 alpha-satellite centromere bed. This is where stock Sniffles finds
almost nothing (its clustering discards single-read satellite signal).</p>{cen_tbl}
<h2>4. Size distribution &amp; detection route</h2>{im(f2, 'Indel size distribution (≤5 kb).')}{src_tbl}
{insqc}
<h2>Caveat</h2><p>Single-molecule calls in satellite cannot be separated from mapping/sequencing artefacts by support alone;
the split-and-remap and the §5 insertion QC are mitigations. Register/TRASH annotation (Arabidopsis CEN178-specific) is omitted.</p>"""
    open(REPORT, "w").write(html)
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
