#!/usr/bin/env python3
"""Step 13 — annotate the 1x (singleton) events with TRASH-based in-register phase.

Replicates the lab's canonical phase method (insertion_origin/scripts/analyze_deletions.py):
TRASH is run on the FULL read, then for the event junction we find the CEN178 monomer
immediately LEFT and immediately RIGHT of the junction. The event is in a CEN178 array only
if BOTH flanking monomers exist (in_cen180_array); the breakpoint phase position within the
178-bp monomer is bp_pos = (read_pos - last_left_monomer.start) % 178. in_register = whole
monomer (|svlen| % 178 within 20) AND in_cen180_array — the unequal-sister-chromatid-HR
signature. Also records detector(s), read de, MAPQ -> confidence HIGH/MEDIUM/LOW.

-> results/singleton_events_annotated.tsv
Run with BASE python (trash_py + pysam):
  /home/jg2070/miniforge3/bin/python 13_annotate_singletons.py
"""
import sys, os, csv, tempfile
from collections import defaultdict, Counter
from types import SimpleNamespace
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/insertion_origin/trash-py/src")
sys.path.insert(0, "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts")
import pysam
import trash_py.pipeline as _tp, trash_py._log as _tl
from common import bam_path, OUT

MONOMER_BP = 178
PHASE_TOL = 20
CEN178 = ("CEN178", "AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAA"
          "GACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAG"
          "AAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG")

_trash_cache = {}


def run_trash_full(seq, key):
    """Run TRASH on the FULL read sequence; return list of {start,end,width} (0-based)."""
    if key in _trash_cache:
        return _trash_cache[key]
    reps = []
    if seq and len(seq) >= 250:
        _tl.configure(quiet=True)
        with tempfile.TemporaryDirectory() as td:
            fa = os.path.join(td, "s.fasta"); tp = os.path.join(td, "t.fasta")
            open(fa, "w").write(f">r\n{seq}\n"); open(tp, "w").write(f">{CEN178[0]}\n{CEN178[1]}\n")
            try:
                _tp.run_pipeline(SimpleNamespace(fasta=fa, output=td, max_rep_size=250,
                                                 min_rep_size=100, templates=tp, processes=1))
            except Exception:
                pass
            rf = os.path.join(td, "s.fasta_repeats.csv")
            if os.path.exists(rf):
                for row in csv.DictReader(open(rf)):
                    reps.append({"start": int(row["start"]) - 1, "end": int(row["end"]), "width": int(row["width"])})
    _trash_cache[key] = reps
    return reps


def find_read_junction(bam, chrom, ref_pos, svtype, svlen, read_id):
    """Return (read.query_sequence, read_pos_at_junction, de, mapq) or (None,...).
    Locates the I/D CIGAR op of ~|svlen| nearest ref_pos in the read carrying the event."""
    size = abs(svlen)
    lo = max(0, ref_pos - size - 3000)
    for r in bam.fetch(chrom, lo, ref_pos + 3000):
        if r.query_name != read_id or r.is_unmapped or r.is_secondary or not r.cigartuples:
            continue
        de = r.get_tag("de") if r.has_tag("de") else -1
        rpos = r.reference_start; qpos = 0; best = None
        for op, l in r.cigartuples:
            if op in (0, 7, 8):              # M/=/X
                rpos += l; qpos += l
            elif op == 1:                    # I (insertion)
                if svtype == "INS" and abs(l - size) <= max(15, size * 0.1):
                    d = abs(rpos - ref_pos)
                    if best is None or d < best[0]:
                        best = (d, qpos)     # read pos = start of insertion
                qpos += l
            elif op == 2:                    # D (deletion)
                if svtype == "DEL" and abs(l - size) <= max(15, size * 0.1):
                    d = abs(rpos - ref_pos)
                    if best is None or d < best[0]:
                        best = (d, qpos)     # read pos = junction (deletion collapsed)
                rpos += l
            elif op in (4, 5):               # S/H
                qpos += l
        if best is not None:
            return r.query_sequence, best[1], de, r.mapping_quality
        # event present but not a single CIGAR op (split) — still return read for context
        return r.query_sequence, None, de, r.mapping_quality
    return None, None, -1, -1


def phase(reps, read_pos, svtype, svlen):
    """Flanking-monomer phase around the junction (canonical method)."""
    if read_pos is None or not reps:
        return {}
    size = abs(svlen)
    right_anchor = read_pos + size if svtype == "INS" else read_pos
    left = [r for r in reps if r["end"] <= read_pos]
    right = [r for r in reps if r["start"] >= right_anchor]
    if not left or not right:
        return {"in_cen180_array": 0}
    last_L = max(left, key=lambda r: r["end"]); first_R = min(right, key=lambda r: r["start"])
    rem = size % MONOMER_BP
    in_phase = (rem <= PHASE_TOL) or (rem >= MONOMER_BP - PHASE_TOL)
    return {"in_cen180_array": 1, "phase_rem": rem, "in_phase": int(in_phase),
            "bp_pos": (read_pos - last_L["start"]) % MONOMER_BP,
            "n_left": len(left), "n_right": len(right)}


def main():
    sing = list(csv.DictReader(open(f"{OUT}/singleton_events.tsv"), delimiter="\t"))
    methods = {}
    for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t"):
        methods[(r["sample"], r["hap"], r["read"], r["svtype"])] = r["methods"]

    by_bam = defaultdict(list)
    for e in sing:
        by_bam[(e["sample"], e["hap"])].append(e)

    out = []; done = 0
    for (sample, hap), events in by_bam.items():
        bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
        for e in events:
            pos = int(e["pos"]); svlen = int(e["svlen"]) if e["svlen"] not in ("", "None") else 0
            svt = e["svtype"]
            seq, read_pos, de, mapq = find_read_junction(bam, e["chrom"], pos, svt, svlen, e["read"])
            reps = run_trash_full(seq, (sample, hap, e["read"])) if seq else []
            nmono = len(reps)
            widths = [r["width"] for r in reps]
            cen178 = nmono >= 3 and (165 <= (sum(widths) / nmono) <= 188)
            ph = phase(reps, read_pos, svt, svlen) if svt in ("DEL", "INS") else {}
            in_array = ph.get("in_cen180_array", 0)
            in_reg = int(bool(in_array and ph.get("in_phase", 0)))
            meth = methods.get((sample, hap, e["read"], svt), "?")
            multi = "+" in meth
            clean = 0 <= de <= 0.01
            if (multi or in_reg) and clean and mapq >= 20:
                conf = "HIGH"
            elif (multi or in_reg or cen178) and mapq >= 10:
                conf = "MEDIUM"
            else:
                conf = "LOW"
            out.append({**e, "methods": meth, "mapq": mapq, "read_de": round(de, 4) if de >= 0 else "",
                        "n_monomers": nmono, "cen178": int(cen178), "in_cen180_array": in_array,
                        "monomer_rem": ph.get("phase_rem", ""), "bp_pos": ph.get("bp_pos", ""),
                        "in_register": in_reg, "confidence": conf})
            done += 1
            if done % 200 == 0:
                print(f"  ...{done}/{len(sing)}  (TRASH cache {len(_trash_cache)})")
        bam.close()

    cols = ["sample", "hap", "tissue", "chrom", "pos", "svtype", "svlen", "methods", "mapq",
            "read_de", "n_monomers", "cen178", "in_cen180_array", "monomer_rem", "bp_pos",
            "in_register", "confidence", "read"]
    with open(f"{OUT}/singleton_events_annotated.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
        for d in sorted(out, key=lambda d: (order[d["confidence"]], d["sample"], d["hap"], d["chrom"], int(d["pos"]))):
            f.write("\t".join(str(d.get(c, "")) for c in cols) + "\n")

    print("=== singleton annotation summary (TRASH on full read; flanking-monomer phase) ===")
    for samp in sorted({d["sample"] for d in out}):
        sub = [d for d in out if d["sample"] == samp]
        c = Counter(d["confidence"] for d in sub)
        ind = [d for d in sub if d["svtype"] in ("DEL", "INS")]
        ireg = sum(d["in_register"] for d in ind); inarr = sum(d["in_cen180_array"] for d in ind)
        print(f"{samp}: n={len(sub)}  HIGH={c['HIGH']} MED={c['MEDIUM']} LOW={c['LOW']}  | DEL/INS={len(ind)} "
              f"in_CEN178_array={inarr} in_register={ireg} ({100*ireg/max(len(ind),1):.0f}% of DEL/INS)")
    print("DONE_ANNOTATE")


if __name__ == "__main__":
    main()
