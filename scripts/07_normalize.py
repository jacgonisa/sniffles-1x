#!/usr/bin/env python3
"""Step 7 — read-Mb normalization so leaf (~500x) and pollen (~30x) are comparable.

Denominator = total mapped read sequence inside the centromere (sum over primary
alignments of aligned-bp overlapping the CEN window), in Mb. Rate = calls per Mb.
-> results/cen_mapped_mb.tsv  and  results/sm_sv_rates.tsv  (+ printed table).
Run with nextflow_env python."""
import csv, pysam
from collections import defaultdict
from common import SAMPLES, HAPS, CEN, bam_path, OUT, refkey


def cen_mapped_mb():
    """Mb of mapped read sequence within CEN, per (sample,hap)."""
    mb = {}
    for sample, _tis in SAMPLES:
        for hap in HAPS:
            bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
            bp = 0
            for chrom, (a, b) in CEN[refkey(sample, hap)].items():
                for r in bam.fetch(chrom, a, b):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary:
                        continue
                    # aligned reference bp overlapping the CEN window
                    bp += max(0, min(r.reference_end, b) - max(r.reference_start, a))
            bam.close()
            mb[(sample, hap)] = bp / 1e6
    return mb


def main():
    mb = cen_mapped_mb()
    with open(f"{OUT}/cen_mapped_mb.tsv", "w") as f:
        f.write("sample\thap\tcen_mapped_mb\n")
        for (s, h), v in sorted(mb.items()):
            f.write(f"{s}\t{h}\t{v:.2f}\n")

    counts = defaultdict(lambda: defaultdict(int))   # (sample,hap) -> svtype -> n
    tis = {}
    with open(f"{OUT}/sm_sv_calls.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            k = (r["sample"], r["hap"]); counts[k][r["svtype"]] += 1
            counts[k]["ALL"] += 1; tis[k] = r["tissue"]

    types = ["ALL", "DEL", "INS", "DUP", "INV", "BND"]
    with open(f"{OUT}/sm_sv_rates.tsv", "w") as f:
        f.write("sample\thap\ttissue\tcen_mapped_mb\t" + "\t".join(f"{t}_per_mb" for t in types) + "\n")
        print(f"{'sample':10}{'hap':5}{'tissue':8}{'Mb':>9}{'ALL/Mb':>9}{'DEL/Mb':>9}{'INS/Mb':>9}{'DUP/Mb':>9}")
        for k in sorted(counts):
            s, h = k; m = mb[k]
            rates = [counts[k][t] / m if m else 0 for t in types]
            f.write(f"{s}\t{h}\t{tis[k]}\t{m:.2f}\t" + "\t".join(f"{x:.3f}" for x in rates) + "\n")
            print(f"{s:10}{h:5}{tis[k]:8}{m:9.1f}{rates[0]:9.3f}{rates[1]:9.3f}{rates[2]:9.3f}{rates[3]:9.3f}")
    print("DONE_NORMALIZE")


if __name__ == "__main__":
    main()
