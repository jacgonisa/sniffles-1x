#!/usr/bin/env python3
"""Step 23 — insertion-sequence quality QC (organism-agnostic).

Some 1x insertions are artefactual: the inserted bases are a low-complexity homopolymer tract,
and/or the CCS base quality collapses inside the insertion relative to its flanks. For every INS
called from a CIGAR I op (so the inserted bases physically exist in one read) we extract the inserted
sequence + its per-base CCS qualities and compute:
  hp_frac    fraction of inserted bases inside homopolymer runs >= HP_RUN bp
  hp_longest longest homopolymer run
  entropy    Shannon entropy over base composition (0-2 bits; low = low complexity)
  q_ins      mean CCS base quality over the inserted bases
  q_flank    mean CCS base quality over +/- FLANK bp of read flanking the insertion
  q_contrast q_flank - q_ins   (>0 = insertion is lower quality than its surroundings)
Flags: low_complexity (hp_frac > HP_FRAC_MAX or entropy < ENT_MIN); quality_decay (q_contrast >= Q_MIN).

-> results/insertion_qc.tsv              per-INS metrics + PASS/FLAG
   results/sm_sv_calls_hiconf.tsv        full callset minus flagged INS (high-confidence set)
Run with nextflow_env python."""
import csv, math
from collections import defaultdict, Counter
import pysam
from common import OUT, bam_path, refkey

HP_RUN = 5            # homopolymer run length that counts as a tract
HP_FRAC_MAX = 0.30    # > this fraction of bases in homopolymer runs -> low complexity
ENT_MIN = 1.2         # Shannon entropy (bits) below this -> low complexity
FLANK = 200           # bp of read flanking the insertion for the quality contrast
Q_MIN = 5.0           # insertion this many Q below flanks -> quality decay


def find_ins(cigar, target):
    q = 0
    for op, ln in cigar:
        if op in (0, 7, 8): q += ln
        elif op == 1:
            if abs(ln - target) <= 3: return q, q + ln
            q += ln
        elif op == 4: q += ln
    return None, None


def homopoly(s):
    if not s: return 0, 0.0
    longest = 0; covered = 0; i = 0; n = len(s)
    while i < n:
        j = i
        while j + 1 < n and s[j + 1] == s[i]: j += 1
        run = j - i + 1
        longest = max(longest, run)
        if run >= HP_RUN: covered += run
        i = j + 1
    return longest, covered / n


def entropy(s):
    if not s: return 0.0
    c = Counter(s); n = len(s)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def main():
    rows = list(csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t"))
    # index INS+CIGAR rows that need a read fetch, grouped by (sample,hap)
    need = defaultdict(list)
    for r in rows:
        if r["svtype"] == "INS" and "CIGAR" in r["methods"] and r["svlen"] not in ("", "None"):
            need[(r["sample"], r["hap"])].append(r)

    qc = {}      # id(row) -> metrics dict
    for (s, h), rr in need.items():
        bam = pysam.AlignmentFile(bam_path(s, h), "rb")
        for r in rr:
            chrom, pos, sz = r["chrom"], int(r["pos"]), abs(int(r["svlen"]))
            found = None
            for rd in bam.fetch(chrom, max(0, pos - 50000), pos + 50000):
                if rd.query_name != r["read"] or rd.is_unmapped or not rd.cigartuples:
                    continue
                if rd.query_sequence is None or rd.is_supplementary or rd.is_secondary:
                    continue
                qs, qe = find_ins(rd.cigartuples, sz)
                if qs is None:
                    continue
                seq = rd.query_sequence[qs:qe]
                qual = rd.query_qualities
                lp, hf = homopoly(seq); ent = entropy(seq)
                if qual is not None and qe > qs:
                    qi = sum(qual[qs:qe]) / (qe - qs)
                    fl = list(qual[max(0, qs - FLANK):qs]) + list(qual[qe:qe + FLANK])
                    qf = sum(fl) / len(fl) if fl else qi
                else:
                    qi = qf = 0.0
                lowc = hf > HP_FRAC_MAX or ent < ENT_MIN
                decay = (qf - qi) >= Q_MIN
                found = {"hp_frac": hf, "hp_longest": lp, "entropy": ent, "q_ins": qi,
                         "q_flank": qf, "q_contrast": qf - qi, "low_complexity": lowc,
                         "quality_decay": decay, "flag": lowc or decay}
                break
            qc[id(r)] = found
        bam.close()

    # per-INS QC table
    cols = ["sample", "hap", "chrom", "pos", "ins_bp", "hp_frac", "hp_longest", "entropy",
            "q_ins", "q_flank", "q_contrast", "low_complexity", "quality_decay", "verdict"]
    n_flag = Counter(); n_lowc = Counter(); n_decay = Counter(); n_ins = Counter()
    with open(f"{OUT}/insertion_qc.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            m = qc.get(id(r))
            if m is None:
                continue
            g = r["sample"]; n_ins[g] += 1
            n_lowc[g] += int(m["low_complexity"]); n_decay[g] += int(m["quality_decay"])
            n_flag[g] += int(m["flag"])
            f.write("\t".join(str(x) for x in [
                r["sample"], r["hap"], r["chrom"], r["pos"], abs(int(r["svlen"])),
                f"{m['hp_frac']:.3f}", m["hp_longest"], f"{m['entropy']:.2f}",
                f"{m['q_ins']:.1f}", f"{m['q_flank']:.1f}", f"{m['q_contrast']:.1f}",
                int(m["low_complexity"]), int(m["quality_decay"]),
                "FLAG" if m["flag"] else "PASS"]) + "\n")

    # high-confidence callset: drop flagged INS, keep everything else
    kept = 0; dropped = 0
    hdr = rows[0].keys() if rows else []
    with open(f"{OUT}/sm_sv_calls_hiconf.tsv", "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            m = qc.get(id(r))
            if m is not None and m["flag"]:
                dropped += 1; continue
            kept += 1
            f.write("\t".join(str(r[c]) for c in hdr) + "\n")

    print("=== insertion QC (CIGAR INS) ===")
    print(f"{'group':16}{'INS':>7}{'low_cplx':>10}{'q_decay':>9}{'FLAGGED':>9}")
    for g in sorted(n_ins):
        print(f"{g:16}{n_ins[g]:7d}{n_lowc[g]:10d}{n_decay[g]:9d}"
              f"{n_flag[g]:9d} ({100*n_flag[g]/max(n_ins[g],1):.0f}%)")
    print(f"high-confidence callset: kept {kept}, dropped {dropped} flagged INS "
          f"-> results/sm_sv_calls_hiconf.tsv")
    print("DONE_INSQC")


if __name__ == "__main__":
    main()
