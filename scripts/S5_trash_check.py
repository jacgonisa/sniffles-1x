#!/usr/bin/env python3
"""Synthetic-control step 5 — validate the CIGAR in-register call with the CANONICAL TRASH-py method
(step 13), not just the size test (|svlen| mod 178). For a random sample of real-CEN vs dirty-artefact
CIGAR DEL/INS, run TRASH on the FULL read, find the flanking CEN178 monomers around the junction, and
report in_cen180_array (both flanks are real monomers) and in_register (array AND whole-monomer phase).
Expected: dirty reads are still IN CEN178 arrays (they are real satellite sequence) but are NOT in
register (random-sized glitch indels) — so the 178-bp register is what separates them, confirmed by TRASH.
Run with BASE python (trash_py + pysam): /home/jg2070/miniforge3/bin/python S5_trash_check.py [N]"""
import sys, csv, random, importlib
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts")
mod = importlib.import_module("13_annotate_singletons")
import pysam

ROOT = "/mnt/ssd-4tb/HIFI_NAMIL"
SM = f"{ROOT}/single_molecule_sv"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200


def real_bam(sample, hap):
    if sample.startswith("cenh3ox") and hap == "col":
        return f"{ROOT}/sv_calling/aligned/{sample}/strict90/col_cenh3oxref.bam"
    return f"{ROOT}/sv_calling/aligned/{sample}/strict90/{hap}_all.bam"


def dirty_bam(sample, hap):
    return f"{ROOT}/sv_calling/aligned_synthetic_dirty/{sample}/strict90/{hap}_all.bam"


def sample_calls(path, n):
    rows = [r for r in csv.DictReader(open(path), delimiter="\t")
            if "CIGAR" in r["methods"] and r["svtype"] in ("DEL", "INS")]
    random.Random(7).shuffle(rows)
    return rows[:n]


def annotate(rows, bamfn):
    bams = {}; n = inarr = inreg = insize = 0
    for r in rows:
        key = (r["sample"], r["hap"])
        if key not in bams:
            bams[key] = pysam.AlignmentFile(bamfn(*key), "rb")
        svlen = int(r["svlen"])
        seq, rp, de, mapq = mod.find_read_junction(bams[key], r["chrom"], int(r["pos"]),
                                                   r["svtype"], svlen, r["read"])
        if seq is None or rp is None:
            continue
        reps = mod.run_trash_full(seq, (r["sample"], r["hap"], r["read"]))
        ph = mod.phase(reps, rp, r["svtype"], svlen)
        n += 1
        inarr += ph.get("in_cen180_array", 0)
        inreg += int(bool(ph.get("in_cen180_array", 0) and ph.get("in_phase", 0)))
        insize += int((abs(svlen) % 178 <= 20) or (abs(svlen) % 178 >= 158))  # size test, for contrast
        if n % 100 == 0:
            print(f"    ...{n}")
    return n, inarr, inreg, insize


def main():
    print(f"TRASH-py in-register check (canonical step-13 method), N={N} per group\n")
    print(f"{'group':12}{'n':>5}{'in_CEN178_array':>17}{'in_register(TRASH)':>20}{'in_phase(size)':>16}")
    out = []
    for name, path, bamfn in [("real CEN", f"{SM}/results/sm_sv_calls.tsv", real_bam),
                              ("dirty", f"{SM}/results_synthetic_dirty/sm_sv_calls.tsv", dirty_bam)]:
        rows = sample_calls(path, N)
        n, inarr, inreg, insize = annotate(rows, bamfn)
        if not n:
            print(f"{name:12}{0:>5}  (no reads matched)"); continue
        print(f"{name:12}{n:>5}{f'{100*inarr/n:.0f}%':>17}{f'{100*inreg/n:.0f}%':>20}{f'{100*insize/n:.0f}%':>16}")
        out.append((name, n, round(100 * inarr / n), round(100 * inreg / n), round(100 * insize / n)))
    with open(f"{SM}/results/trash_check.tsv", "w") as f:
        f.write("group\tn\tin_cen178_array_pct\tin_register_trash_pct\tin_phase_size_pct\n")
        for r in out:
            f.write("\t".join(map(str, r)) + "\n")
    print("DONE_TRASH_CHECK")


if __name__ == "__main__":
    main()
