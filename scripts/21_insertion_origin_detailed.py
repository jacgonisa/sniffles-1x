#!/usr/bin/env python3
"""Step 21 — DETAILED insertion-origin panels in the lab's validation-prototype style
(split-read origin view · CCS-quality track · KMC readmer track · self-similarity dotplot).

Reuses insertion_origin/scripts/plot_validation_prototype.py (draw_event, parse_cigar,
build_dotplot, load_paf, readmer) unchanged; we just feed it OUR single-molecule INS calls,
OUR BAMs (per-genotype reference via refkey), minimap2 the inserted fragment back to that
reference for the origin PAF, and reuse the lab's KMC dbs for the readmer track (CENH3ox
datasets; WT shows 'readmer unavailable').

-> results/insertion_origin_detailed/<sample>_<hap>_<label>.png
Run with nextflow_env python (pysam + minimap2 + get_counts_threads)."""
import os, sys, csv, subprocess, tempfile
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/insertion_origin/scripts")
import pysam
import plot_validation_prototype as pv
import insplot
from common import OUT, REF, CEN, bam_path, refkey

MINIMAP2 = "/home/jg2070/minimap2/minimap2"
OUTDIR = f"{OUT}/insertion_origin_detailed"; os.makedirs(OUTDIR, exist_ok=True)
KMC_DIR = "/mnt/ssd-4tb/HIFI_NAMIL/insertion_origin/results/hifi_multi_dataset/kmc_dbs"
MINSZ = 1500          # detailed panels only for the big ones
N = 6
DS = {"cenh3ox_leaf": "CENH3ox_leaf", "cenh3ox_pollen": "CENH3ox_pollen"}  # -> KMC db prefix
HP = {"col": "Col", "ler": "Ler"}
pv.OUT_DIR = OUTDIR    # redirect the lab function's output


def pick_locus(sample, hap, chrom, pos, tol=3000):
    """target one specific INS locus (bypass the size filter) -> [row]."""
    best = None
    for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t"):
        if (r["svtype"] == "INS" and "CIGAR" in r["methods"] and r["sample"] == sample
                and r["hap"] == hap and r["chrom"] == chrom
                and r["svlen"] not in ("", "None") and abs(int(r["pos"]) - pos) <= tol):
            if best is None or abs(int(r["svlen"])) > abs(int(best["svlen"])):
                best = r
    return [best] if best else []


def pick():
    rows = [r for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t")
            if r["svtype"] == "INS" and "CIGAR" in r["methods"]
            and r["svlen"] not in ("", "None") and abs(int(r["svlen"])) >= MINSZ]
    rows.sort(key=lambda r: -abs(int(r["svlen"])))
    out, seen = [], []
    for r in rows:                    # distinct loci, spread across the 4 groups
        pos = int(r["pos"]); k = (r["sample"], r["hap"], r["chrom"])
        if any(kk == k and abs(pos - p) < 5000 for kk, p in seen):
            continue
        seen.append((k, pos)); out.append(r)
        if len({x["sample"] for x in out}) >= 2 and len(out) >= N:
            break
    return out[:N]


def get_read(sample, hap, chrom, pos, read_name, size):
    bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
    try:
        for r in bam.fetch(chrom, max(0, pos - 60000), pos + 60000):
            if r.query_name != read_name or r.is_unmapped or not r.cigartuples:
                continue
            return {"seq": r.query_sequence or "", "quals": list(r.query_qualities or []),
                    "cigar": r.cigarstring or "", "rq": r.get_tag("rq") if r.has_tag("rq") else None,
                    "np": r.get_tag("np") if r.has_tag("np") else None}
    finally:
        bam.close()
    return None


def origin_hits(label, chrom, pos, size, seq, rk):
    """minimap2 the inserted fragment vs REF[rk]; return pv-style hit dicts."""
    with tempfile.TemporaryDirectory() as td:
        fa = f"{td}/x.fa"; paf = f"{td}/x.paf"
        qn = f"{label}|{chrom}|{pos}|{size}bp"
        open(fa, "w").write(f">{qn}\n{seq}\n")
        subprocess.run([MINIMAP2, "-x", "map-hifi", "-c", "--secondary=yes", "-N", "8",
                        "-t", "8", REF[rk][0], fa], stdout=open(paf, "w"),
                       stderr=subprocess.DEVNULL, check=True)
        return [h for h in pv.load_paf(paf) if h["qname"] == qn]


def main():
    if len(sys.argv) >= 5 and sys.argv[1] == "--locus":
        _, _, sample, hap, chrom, pos = sys.argv[:6]
        events = pick_locus(sample, hap, chrom, int(pos))
        print(f"targeted locus {sample} {hap} {chrom}:{pos} -> {len(events)} event")
    else:
        events = pick()
        print(f"detailed panels for {len(events)} insertions (>= {MINSZ} bp)")
    for i, r in enumerate(events):
        s, h, chrom, pos, size = r["sample"], r["hap"], r["chrom"], int(r["pos"]), abs(int(r["svlen"]))
        rk = refkey(s, h); label = f"ins{i:02d}"
        rd = get_read(s, h, chrom, pos, r["read"], size)
        if not rd:
            print(f"  {label}: read not found"); continue
        # inserted sequence for the origin map-back
        b, ins = pv.parse_cigar(rd["cigar"], pos, size)
        if not ins:
            print(f"  {label}: insertion not located in CIGAR"); continue
        ins_seq = rd["seq"][ins["q_s"]:ins["q_e"]]
        hits = origin_hits(label, chrom, pos, size, ins_seq, rk)
        pv.CEN_CENH3OX = CEN[rk]            # centromere band = this event's reference
        db = os.path.join(KMC_DIR, f"{DS.get(s,'NA')}_{HP[h]}")   # WT -> missing -> readmer unavailable
        ev = {"ref_pos": pos, "ins_size": size, "chrom": chrom, "label": label,
              "call": f"{s} {h}", "dataset": s, "haplotype": h}
        print(f"  {label}: {chrom}:{pos} {size}bp, {len(hits)} origin hits")
        out_png = os.path.join(OUTDIR, f"{s}_{h}_{label}.png")
        insplot.draw_event_trash(ev, rd["seq"], rd["quals"] or None, rd["rq"], rd["np"],
                                 rd["cigar"], hits, db, out_png)
    print(f"DONE_INSORIGIN_DETAILED -> {OUTDIR}/")


if __name__ == "__main__":
    main()
