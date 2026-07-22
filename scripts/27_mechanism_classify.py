#!/usr/bin/env python3
"""Step 27 — integrate charla_hifi (inter-homolog) + sniffles_1x (self/sister + ectopic) into one
mechanism callset for the Col×Ler F1 (arabidopsis WT leaf/pollen; CENH3ox has no CHARLA).

(A) Classify every sniffles_1x non-hybrid SV call into a taxonomy class:
      self_sister_unequal_inreg   CEN DEL/INS/DUP, in-register (whole-monomer)  [class 1]
      self_sister_unequal_offreg  CEN DEL/INS/DUP, out-of-register (NHEJ-like)   [class 1]
      inversion                   INV                                            [class 2]
      ectopic_samehap             BND to a different-chromosome ARM              [class 4]
      artefact_crossmap           BND to another CEN / unplaced (satellite x-map)[class 10]
      artefact_lowqual            INS flagged homopolymer/quality                [class 11]
(B) Read CHARLA per-sample summary.tsv and map its categories to:
      inter_homolog_crossover     4b-2 + 4b-5 (allelic HR / crossover)           [class 5]
      inter_homolog_satellite     1a-1d (tandem-repeat exchange between homologs)[class 6]
      ectopic_interhomolog        4a  (non_homologous)                           [class 7]
      charla_ambiguous            2,3,4b-1/3/4/6,4c,liftover_failed (artefact)   [class 10/12]
-> results/mechanism_summary.tsv  (source, group, mechanism, class, count)  + printed.
Run with nextflow_env python (default arabidopsis profile)."""
import csv
from collections import defaultdict, Counter
from common import OUT, GROUPS

CHARLA = {   # sample -> CHARLA t50 summary.tsv (WT only; CENH3ox not run through CHARLA)
    "wt_leaf": "/mnt/ssd-4tb/HIFI_NAMIL/01_f1leaf-wt/01-mask_0/output/10-call_recombination_sites/wt_leaf_threshold_50.summary.tsv",
    "wt_pollen": "/mnt/ssd-4tb/HIFI_NAMIL/03_f1pollen-wt/01-mask_0/output/10-call_recombination_sites/wt_pollen_threshold_50.summary.tsv",
}
CHARLA_MAP = {"4b-2": ("inter_homolog_crossover", 5), "4b-5": ("inter_homolog_crossover", 5),
              "1a": ("inter_homolog_satellite", 6), "1b": ("inter_homolog_satellite", 6),
              "1c": ("inter_homolog_satellite", 6), "1d": ("inter_homolog_satellite", 6),
              "4a": ("ectopic_interhomolog", 7)}   # everything else -> charla_ambiguous


def sniffles_mechanism():
    # BND mate category + QC-flagged INS lookups
    bnd = {}
    try:
        for d in csv.DictReader(open(f"{OUT}/translocations.tsv"), delimiter="\t"):
            bnd[(d["sample"], d["hap"], d["chrom"], d["pos"], d["read"])] = d["category"]
    except FileNotFoundError:
        pass
    flagged = set()
    try:
        for d in csv.DictReader(open(f"{OUT}/insertion_qc.tsv"), delimiter="\t"):
            if d["verdict"] == "FLAG":
                flagged.add((d["sample"], d["hap"], d["chrom"], d["pos"]))
    except FileNotFoundError:
        pass

    out = defaultdict(Counter)   # group -> mechanism -> n
    for r in csv.DictReader(open(f"{OUT}/sm_sv_calls.tsv"), delimiter="\t"):
        g, t = r["sample"], r["svtype"]
        key = (r["sample"], r["hap"], r["chrom"], r["pos"])
        if t == "INS" and key in flagged:
            m = ("artefact_lowqual", 11)
        elif t in ("DEL", "INS", "DUP"):
            m = ("self_sister_unequal_inreg", 1) if r["in_phase"] == "1" else ("self_sister_unequal_offreg", 1)
        elif t == "INV":
            m = ("inversion", 2)
        elif t == "BND":
            cat = bnd.get((r["sample"], r["hap"], r["chrom"], r["pos"], r["read"]), "")
            m = ("ectopic_samehap", 4) if cat == "other_chrom_arm" else ("artefact_crossmap", 10)
        else:
            m = ("other", 0)
        out[g][m] += 1
    return out


def charla_mechanism():
    out = defaultdict(Counter)
    for sample, path in CHARLA.items():
        try:
            for d in csv.DictReader(open(path), delimiter="\t"):
                if d["category"] in ("TOTAL", ""):
                    continue
                m = CHARLA_MAP.get(d["category"], ("charla_ambiguous", 12))
                out[sample][m] += int(d["count"])
        except FileNotFoundError:
            print(f"  (CHARLA summary missing: {sample})")
    return out


def main():
    sn = sniffles_mechanism()
    ch = charla_mechanism()
    rows = []
    for g in GROUPS:
        for (m, cls), n in sorted(sn[g].items(), key=lambda x: -x[1]):
            rows.append(("sniffles_1x", g, m, cls, n))
    for g in ("wt_leaf", "wt_pollen"):
        for (m, cls), n in sorted(ch[g].items(), key=lambda x: -x[1]):
            rows.append(("charla_hifi", g, m, cls, n))
    with open(f"{OUT}/mechanism_summary.tsv", "w") as f:
        f.write("source\tgroup\tmechanism\ttaxonomy_class\tcount\n")
        for s, g, m, cls, n in rows:
            f.write(f"{s}\t{g}\t{m}\t{cls}\t{n}\n")

    print("=== integrated mechanism counts (Col×Ler F1; CHARLA = WT only) ===")
    print(f"{'source':13}{'group':16}{'mechanism':30}{'cls':>4}{'count':>8}")
    for s, g, m, cls, n in rows:
        print(f"{s:13}{g:16}{m:30}{cls:>4}{n:>8}")
    print("DONE_MECHANISM")


if __name__ == "__main__":
    main()
