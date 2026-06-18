#!/usr/bin/env python3
"""Step 13 — annotate the 1x (singleton) events with trustworthiness flags.

A single read carrying an SV (support=1) is either a genuine somatic single molecule or a
one-off mapping/sequencing artifact — support alone cannot tell them apart. We add the
orthogonal evidence:

  methods      which detector(s) found it (CIGAR / SPLITREAD / SPLITANDMAP); multi = corroborated
  read_de      gap-compressed divergence of the read (low = clean read, not a noisy mapping)
  mapq         mapping quality of the read at the locus
  TRASH on the READ window around the breakpoint -> n_monomers, mono_mean_w, cen178 (satellite?)
  in_register  whole-CEN178-monomer event (|svlen|%178 ~ 0) AND TRASH confirms a CEN178 array
               = unequal-sister-chromatid-HR signature (vs out-of-phase NHEJ-like)
  confidence   HIGH / MEDIUM / LOW from the combination above

-> results/singleton_events_annotated.tsv
Run with BASE python (has trash_py + pysam):
  /home/jg2070/miniforge3/bin/python 13_annotate_singletons.py
"""
import sys, os, csv, tempfile
from collections import defaultdict
from types import SimpleNamespace
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/insertion_origin/trash-py/src")
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts")
import pysam
import trash_py.pipeline as _tp, trash_py._log as _tl
from common import bam_path, OUT, MONO, TOL

TPL = ("CEN178", "AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAA"
       "GACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAG"
       "AAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG")
WIN = 1500          # bp each side of the breakpoint extracted from the read for TRASH


def run_trash(seq):
    """Return (n_monomers, mean_width) of CEN178-template repeats in seq."""
    if not seq or len(seq) < 200:
        return 0, 0.0
    _tl.configure(quiet=True)
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "s.fasta"); tp = os.path.join(td, "t.fasta")
        open(fa, "w").write(f">r\n{seq}\n"); open(tp, "w").write(f">{TPL[0]}\n{TPL[1]}\n")
        try:
            _tp.run_pipeline(SimpleNamespace(fasta=fa, output=td, max_rep_size=250,
                                             min_rep_size=100, templates=tp, processes=1))
        except Exception:
            pass
        rf = os.path.join(td, "s.fasta_repeats.csv"); w = []
        if os.path.exists(rf):
            for r in csv.DictReader(open(rf)):
                w.append(int(r["width"]))
    return len(w), (sum(w) / len(w) if w else 0.0)


def read_window(bam, chrom, pos, read_name):
    """Read sequence ±WIN around the breakpoint (read coords), plus de & mapq."""
    for rd in bam.fetch(chrom, max(0, pos - 2000), pos + 2000):
        if rd.query_name != read_name or rd.is_secondary:
            continue
        de = rd.get_tag("de") if rd.has_tag("de") else -1
        rp = None
        for qp, refp in rd.get_aligned_pairs():
            if refp is not None and refp >= pos and qp is not None:
                rp = qp; break
        if rp is None:
            rp = (rd.query_alignment_start + rd.query_alignment_end) // 2
        seq = rd.query_sequence
        return (seq[max(0, rp - WIN):rp + WIN] if seq else None, de, rd.mapping_quality)
    return None, -1, -1


def main():
    sing = list(csv.DictReader(open(f"{OUT}/singleton_events.tsv"), delimiter="\t"))
    # methods per (sample,hap,read,svtype) from the merged calls
    methods = {}
    for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t"):
        methods[(r["sample"], r["hap"], r["read"], r["svtype"])] = r["methods"]

    by_bam = defaultdict(list)
    for e in sing:
        by_bam[(e["sample"], e["hap"])].append(e)

    out = []
    done = 0
    for (sample, hap), events in by_bam.items():
        bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
        for e in events:
            pos = int(e["pos"]); svlen = int(e["svlen"]) if e["svlen"] not in ("", "None") else 0
            seq, de, mapq = read_window(bam, e["chrom"], pos, e["read"])
            nmono, meanw = run_trash(seq) if seq else (0, 0.0)
            cen178 = nmono >= 3 and 165 <= meanw <= 188
            rem = abs(svlen) % MONO if svlen else None
            whole_mono = svlen != 0 and (rem <= TOL or rem >= MONO - TOL)
            in_reg = bool(cen178 and whole_mono and e["svtype"] in ("DEL", "INS", "DUP"))
            meth = methods.get((sample, hap, e["read"], e["svtype"]), "?")
            multi = "+" in meth
            clean = 0 <= de <= 0.01
            # confidence
            if (multi or in_reg) and clean and mapq >= 20:
                conf = "HIGH"
            elif (multi or in_reg or cen178) and mapq >= 10:
                conf = "MEDIUM"
            else:
                conf = "LOW"
            out.append({**e, "methods": meth, "mapq": mapq, "read_de": round(de, 4) if de >= 0 else "",
                        "n_monomers": nmono, "mono_mean_w": round(meanw, 1), "cen178": int(cen178),
                        "monomer_rem": rem if rem is not None else "", "in_register": int(in_reg),
                        "confidence": conf})
            done += 1
            if done % 200 == 0:
                print(f"  ...{done}/{len(sing)}")
        bam.close()

    cols = ["sample", "hap", "tissue", "chrom", "pos", "svtype", "svlen", "methods", "mapq",
            "read_de", "n_monomers", "mono_mean_w", "cen178", "monomer_rem", "in_register",
            "confidence", "read"]
    with open(f"{OUT}/singleton_events_annotated.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for d in sorted(out, key=lambda d: (d["confidence"] != "HIGH", d["confidence"] != "MEDIUM",
                                            d["sample"], d["hap"], d["chrom"], int(d["pos"]))):
            f.write("\t".join(str(d.get(c, "")) for c in cols) + "\n")

    # summary
    from collections import Counter
    print("=== singleton annotation summary ===")
    for tis in ("leaf", "pollen"):
        sub = [d for d in out if d["tissue"] == tis]
        conf = Counter(d["confidence"] for d in sub)
        ireg = sum(d["in_register"] for d in sub)
        c178 = sum(d["cen178"] for d in sub)
        multi = sum(1 for d in sub if "+" in d["methods"])
        print(f"{tis}: n={len(sub)}  HIGH={conf['HIGH']} MED={conf['MEDIUM']} LOW={conf['LOW']}  "
              f"in_register={ireg} ({100*ireg/max(len(sub),1):.0f}%)  cen178={c178}  multi_method={multi}")
    print("DONE_ANNOTATE")


if __name__ == "__main__":
    main()
