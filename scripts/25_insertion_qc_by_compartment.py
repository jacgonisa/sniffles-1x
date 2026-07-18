#!/usr/bin/env python3
"""Step 25 — insertion-quality QC dissected into CENTROMERE vs ARMS, for both organisms.

Fills the 2x2 (Arabidopsis|human) x (CEN|arm):
  - Arabidopsis CEN: results/insertion_qc.tsv (already CEN-only; step 23).
  - Arabidopsis ARM: NEW — scan the distal-arm windows (CEN_end+5Mb, 3Mb wide, as step 16), find CIGAR
    I ops >= MINSV, extract inserted seq + CCS qualities, apply the SAME homopolymer/quality flags.
  - Human CEN & ARM: results_human/insertion_qc.tsv (genome-wide; step 23) classified by the HG002
    alpha_CEN / all_ARMS beds.
-> results/insertion_qc_by_compartment.tsv  (organism, compartment, n_ins, low_cplx, q_decay, flagged)
Run with nextflow_env python (default SM_GENOME=arabidopsis; also reads the human TSVs)."""
import os, csv, math
from collections import defaultdict, Counter
import pysam
from common import SAMPLES, HAPS, CEN, CHRLEN, bam_path, refkey, OUT

A_OUT = OUT                                    # results/ (arabidopsis)
H_OUT = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/results_human"
HBED = "/mnt/ssd-8tb/HUMAN/data/assembly/annotation/cen_arms"
ARM_BUFFER, ARM_WIN = 5_000_000, 3_000_000
HP_RUN, HP_FRAC_MAX, ENT_MIN, FLANK, Q_MIN, MINSV = 5, 0.30, 1.2, 200, 5.0, 50


def homopoly(s):
    if not s: return 0.0
    covered = 0; i = 0; n = len(s)
    while i < n:
        j = i
        while j + 1 < n and s[j + 1] == s[i]: j += 1
        if j - i + 1 >= HP_RUN: covered += j - i + 1
        i = j + 1
    return covered / n


def entropy(s):
    if not s: return 0.0
    c = Counter(s); n = len(s)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def flags(seq, qual, qs, qe):
    lowc = homopoly(seq) > HP_FRAC_MAX or entropy(seq) < ENT_MIN
    decay = False
    if qual is not None and qe > qs:
        qi = sum(qual[qs:qe]) / (qe - qs)
        fl = list(qual[max(0, qs - FLANK):qs]) + list(qual[qe:qe + FLANK])
        qf = sum(fl) / len(fl) if fl else qi
        decay = (qf - qi) >= Q_MIN
    return lowc, decay


def arm_windows(rk):
    w = {}
    for chrom, (a, b) in CEN[rk].items():
        s = b + ARM_BUFFER; e = min(s + ARM_WIN, CHRLEN[rk][chrom] - 100_000)
        if e - s > 200_000:
            w[chrom] = (s, e)
    return w


def arabidopsis_arm():
    """scan arm windows, QC every CIGAR I op >= MINSV."""
    n = lowc = decay = flag = 0
    for sample, _ in SAMPLES:
        for hap in HAPS:
            rk = refkey(sample, hap)
            bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
            for chrom, (a, b) in arm_windows(rk).items():
                for r in bam.fetch(chrom, a, b):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.cigartuples:
                        continue
                    if not (a <= r.reference_start < b) or r.query_sequence is None:
                        continue
                    seq = r.query_sequence; qual = r.query_qualities; q = 0
                    for op, ln in r.cigartuples:
                        if op in (0, 7, 8): q += ln
                        elif op == 1:
                            if ln >= MINSV:
                                n += 1
                                lc, dc = flags(seq[q:q + ln], qual, q, q + ln)
                                lowc += lc; decay += dc; flag += (lc or dc)
                            q += ln
                        elif op == 4: q += ln
            bam.close()
    return n, lowc, decay, flag


def load_bed(path):
    iv = defaultdict(list)
    if os.path.exists(path):
        for ln in open(path):
            if ln.strip() and not ln.startswith(("#", "track")):
                c = ln.split("\t"); iv[c[0]].append((int(c[1]), int(c[2])))
    return iv


def in_iv(iv, chrom, pos):
    return any(s <= pos < e for s, e in iv.get(chrom, []))


def human_by_compartment():
    cen = {h: load_bed(f"{HBED}/hg002v1.1.{h}.alpha_CEN.bed") for h in ("MAT", "PAT")}
    arm = {h: load_bed(f"{HBED}/hg002v1.1.{h}.all_ARMS.bed") for h in ("MAT", "PAT")}
    agg = defaultdict(lambda: [0, 0, 0, 0])   # compartment -> [n, lowc, decay, flag]
    for d in csv.DictReader(open(f"{H_OUT}/insertion_qc.tsv"), delimiter="\t"):
        h, chrom, pos = d["hap"], d["chrom"], int(d["pos"])
        if in_iv(cen[h], chrom, pos): k = "CEN"
        elif in_iv(arm[h], chrom, pos): k = "ARM"
        else: k = "other"
        a = agg[k]; a[0] += 1
        a[1] += int(d["low_complexity"]); a[2] += int(d["quality_decay"]); a[3] += (d["verdict"] == "FLAG")
    return agg


def arab_cen():
    a = [0, 0, 0, 0]
    for d in csv.DictReader(open(f"{A_OUT}/insertion_qc.tsv"), delimiter="\t"):
        a[0] += 1; a[1] += int(d["low_complexity"]); a[2] += int(d["quality_decay"]); a[3] += (d["verdict"] == "FLAG")
    return a


def main():
    print("scanning Arabidopsis arm windows for insertions ...")
    ac = arab_cen()
    aa = arabidopsis_arm()
    h = human_by_compartment()
    rows = [
        ("arabidopsis", "CEN", *ac),
        ("arabidopsis", "ARM", *aa),
        ("human", "CEN", *h["CEN"]),
        ("human", "ARM", *h["ARM"]),
        ("human", "other(pericen)", *h["other"]),
    ]
    with open(f"{A_OUT}/insertion_qc_by_compartment.tsv", "w") as f:
        f.write("organism\tcompartment\tn_ins\tlow_complexity\tquality_decay\tflagged\tpct_flagged\n")
        for org, comp, n, lc, dc, fl in rows:
            f.write(f"{org}\t{comp}\t{n}\t{lc}\t{dc}\t{fl}\t{100*fl/max(n,1):.1f}\n")
    print(f"\n{'organism':12}{'compartment':16}{'INS':>8}{'low_cplx':>10}{'q_decay':>9}{'flagged':>12}")
    for org, comp, n, lc, dc, fl in rows:
        print(f"{org:12}{comp:16}{n:8d}{lc:10d}{dc:9d}{fl:8d} ({100*fl/max(n,1):.0f}%)")
    print("DONE_INSQC_COMPARTMENT")


if __name__ == "__main__":
    main()
