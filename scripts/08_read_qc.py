#!/usr/bin/env python3
"""Step 8 — read-property controls: are pollen single-molecule SVs a read-quality
artifact? Computes, leaf vs pollen, the metrics that could inflate per-read SV calls:
  - read length (mapped CEN, median kb)            -> shorter reads, less opportunity
  - arm de%  (error proxy: divergence in non-repeat arms, ~sequencing error)
  - CEN de%  (divergence in centromere = error + true satellite variation)
  - np (HiFi passes, median) + rq (predicted accuracy) from the SOURCE hifi_reads.bam
    (these PacBio tags are stripped by samtools fastq, so read them pre-mapping)
-> results/read_qc.tsv  (+ printed). All point the same way: pollen reads are
equal-or-better quality, so the higher pollen per-Mb SV rate is not a quality artifact.
Run with nextflow_env python."""
import pysam, statistics as st, sys
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts")
from common import SAMPLES, HAPS, CEN, bam_path, OUT

SRC = {"leaf": "/mnt/ssd-4tb/HIFI_NAMIL/01_f1leaf-wt/data/WT_leaf.hifi_reads.bam",
       "pollen": "/mnt/ssd-4tb/HIFI_NAMIL/03_f1pollen-wt/03_F1pollen_WT-1/data/WT_pollen.hifi_reads.bam"}
ARM = [('Chr1', 2_000_000, 10_000_000), ('Chr2', 8_000_000, 15_000_000)]
CAP = 300000


def de_len(bam, regions, cap=40000):
    de, ln = [], []
    for chrom, a, b in regions:
        for r in bam.fetch(chrom, a, b):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or not (a <= r.reference_start < b):
                continue
            ln.append(r.query_length)
            if r.has_tag("de"):
                de.append(r.get_tag("de"))
            if len(de) >= cap:
                break
        if len(de) >= cap:
            break
    return de, ln


def src_np_rq(path):
    bam = pysam.AlignmentFile(path, "rb", check_sq=False)
    nps, rqs = [], []
    for r in bam:
        if r.has_tag("np"):
            nps.append(r.get_tag("np"))
        if r.has_tag("rq"):
            rqs.append(r.get_tag("rq"))
        if len(nps) >= CAP:
            break
    bam.close()
    return st.median(nps), st.median(rqs) * 100


def main():
    npm = {t: src_np_rq(p) for t, p in SRC.items()}   # tissue -> (np_med, rq_med%)
    rows = []
    for sample, tis in SAMPLES:
        for hap in HAPS:
            bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
            de_arm, _ = de_len(bam, ARM)
            de_cen, ln_cen = de_len(bam, [(c, a, b) for c, (a, b) in CEN[hap].items()])
            bam.close()
            np_med, rq_med = npm[tis]
            rows.append({"sample": sample, "hap": hap, "tissue": tis,
                         "cen_med_kb": st.median(ln_cen) / 1e3,
                         "arm_de_pct": st.mean(de_arm) * 100,
                         "cen_de_pct": st.mean(de_cen) * 100,
                         "np_med": np_med, "rq_med_pct": rq_med})
    cols = ["sample", "hap", "tissue", "cen_med_kb", "arm_de_pct", "cen_de_pct", "np_med", "rq_med_pct"]
    with open(f"{OUT}/read_qc.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    print(f"{'sample':10}{'hap':4}{'CENkb':>7}{'arm_de%':>9}{'cen_de%':>9}{'np':>5}{'rq%':>9}")
    for r in rows:
        print(f"{r['sample']:10}{r['hap']:4}{r['cen_med_kb']:7.1f}{r['arm_de_pct']:9.3f}{r['cen_de_pct']:9.3f}{r['np_med']:5.0f}{r['rq_med_pct']:9.3f}")
    print("DONE_READQC")


if __name__ == "__main__":
    main()
