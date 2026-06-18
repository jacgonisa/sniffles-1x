#!/usr/bin/env python3
"""Step 1 — candidate read selection (colleague's gating).
For each WT BAM, fetch PRIMARY alignments overlapping the centromere and keep a
read if  de >= 0.005  OR  NM >= 50  OR  it has an SA tag.
-> results/candidates/{sample}_{hap}.tsv  (qname chrom pos de nm has_sa)
Run with nextflow_env python (pysam)."""
import os, pysam
from common import SAMPLES, HAPS, CEN, bam_path, OUT

DE_MIN = 0.005
NM_MIN = 50

def main():
    od = f"{OUT}/candidates"; os.makedirs(od, exist_ok=True)
    print(f"{'sample':10} {'hap':4} {'reads_in_cen':>12} {'candidates':>11}")
    for sample, _tis in SAMPLES:
        for hap in HAPS:
            bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
            n_seen = n_cand = 0
            with open(f"{od}/{sample}_{hap}.tsv", "w") as out:
                out.write("qname\tchrom\tpos\tde\tnm\thas_sa\n")
                for chrom, (a, b) in CEN[hap].items():
                    for r in bam.fetch(chrom, a, b):
                        if r.is_unmapped or r.is_secondary or r.is_supplementary:
                            continue
                        if not (a <= r.reference_start < b):   # primary anchored in CEN
                            continue
                        n_seen += 1
                        de = r.get_tag("de") if r.has_tag("de") else 0.0
                        nm = r.get_tag("NM") if r.has_tag("NM") else 0
                        sa = r.has_tag("SA")
                        if de >= DE_MIN or nm >= NM_MIN or sa:
                            n_cand += 1
                            out.write(f"{r.query_name}\t{chrom}\t{r.reference_start}\t{de:.5f}\t{nm}\t{int(sa)}\n")
            print(f"{sample:10} {hap:4} {n_seen:12d} {n_cand:11d}")
            bam.close()
    print("DONE_CANDIDATES")

if __name__ == "__main__":
    main()
