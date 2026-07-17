#!/usr/bin/env python3
"""Shared paths/constants for the single-molecule SV pipeline.

Two genome profiles, selected by env var SM_GENOME (default 'arabidopsis'):
  arabidopsis — WT+CENH3ox, leaf+pollen, Col+Ler; CENTROMERE-restricted (SCAN = CEN).
  human       — HG002 sperm (BLS0005+BLS0006 merged), MAT+PAT; GENOME-WIDE (SCAN = whole contigs),
                CEN kept only as an overlay to annotate which calls fall in alpha-satellite.

Every downstream script imports the same names; the only new abstraction is SCAN[refkey] = the list
of (chrom, start, end) regions to iterate/fetch (was hard-coded to CEN before)."""
import os

ROOT = "/mnt/ssd-4tb/HIFI_NAMIL"
GENOME = os.environ.get("SM_GENOME", "arabidopsis")


def _fai_lengths(fai):
    d = {}
    for ln in open(fai):
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 2:
            d[c[0]] = int(c[1])
    return d


def _bed_intervals(bed):
    d = {}
    for ln in open(bed):
        if not ln.strip() or ln.startswith(("#", "track")):
            continue
        c = ln.split("\t")
        d.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    return d


if GENOME == "human":
    HROOT = "/mnt/ssd-8tb/HUMAN"
    OUT = f"{ROOT}/single_molecule_sv/results_human"
    SAMPLES = [("BLS0005_BLS0006_sperm", "sperm")]
    HAPS = ["MAT", "PAT"]
    GROUPS = ["BLS0005_BLS0006_sperm"]      # one sample; the human report regroups by haplotype
    MONO, TOL = 0, 0                         # CEN178 register is Arabidopsis-only

    def genotype(sample):
        return sample

    def refkey(sample, hap):
        return hap                           # each haplotype -> its own MAT/PAT assembly

    def bam_path(sample, hap):
        return (f"{HROOT}/sv_calling/strict90_sniffles1x_BLS0005_BLS0006/bam/"
                f"{sample}.{hap}.strict90.merged.md.bam")   # calmd'd (MD tags added), see prep_human_md.sh

    _G = f"{HROOT}/data/assembly/genomes"
    _A = f"{HROOT}/data/assembly/annotation/cen_arms"
    REF = {"MAT": (f"{_G}/MAT.fasta", f"{_G}/mat_repetitive_k15.txt"),
           "PAT": (f"{_G}/PAT.fasta", f"{_G}/pat_repetitive_k15.txt")}
    CHRLEN = {"MAT": _fai_lengths(f"{_G}/MAT.fasta.fai"),
              "PAT": _fai_lengths(f"{_G}/PAT.fasta.fai")}
    CEN_IV = {"MAT": _bed_intervals(f"{_A}/hg002v1.1.MAT.alpha_CEN.bed"),
              "PAT": _bed_intervals(f"{_A}/hg002v1.1.PAT.alpha_CEN.bed")}
    CEN = CEN_IV                              # overlay only (genome-wide scan)
    SCAN = {rk: [(c, 0, L) for c, L in CHRLEN[rk].items()] for rk in HAPS}

    def in_phase(svlen):                      # register not defined for human
        return (False, abs(svlen))

else:  # ---- arabidopsis (default) ----
    OUT = f"{ROOT}/single_molecule_sv/results"
    SAMPLES = [("wt_leaf", "leaf"), ("wt_pollen", "pollen"),
               ("cenh3ox_leaf", "leaf"), ("cenh3ox_pollen", "pollen")]
    HAPS = ["col", "ler"]
    GROUPS = ["wt_leaf", "cenh3ox_leaf", "wt_pollen", "cenh3ox_pollen"]
    MONO, TOL = 178, 20

    def genotype(sample):
        return "cenh3ox" if "cenh3ox" in sample else "wt"

    # CENH3ox col reads are mapped to the CENH3ox line's OWN remodelled assembly (its baseline).
    def refkey(sample, hap):
        return "cenh3ox_col" if (sample.startswith("cenh3ox") and hap == "col") else hap

    def bam_path(sample, hap):
        if sample.startswith("cenh3ox") and hap == "col":
            return f"{ROOT}/sv_calling/aligned/{sample}/strict90/col_cenh3oxref.bam"
        return f"{ROOT}/sv_calling/aligned/{sample}/strict90/{hap}_all.bam"

    REF = {
        "col": (f"{ROOT}/01_genomes/Col-HiFi/Col-0.ragtag_scaffolds_with_organellar.fa",
                f"{ROOT}/01_genomes/Col-HiFi/Col-0.ragtag_scaffolds_with_organellar.repetitive_k15.txt"),
        "ler": (f"{ROOT}/01_genomes/Ler-HiFi/Ler-0.ragtag_scaffolds_with_organellar.fa",
                f"{ROOT}/01_genomes/Ler-HiFi/Ler-0.ragtag_scaffolds_with_organellar.repetitive_k15.txt"),
        "cenh3ox_col": (f"{ROOT}/01_genomes/CENH3ox-Col-parent/cenh3ox_col.renamed.fa",
                        f"{ROOT}/01_genomes/CENH3ox-Col-parent/cenh3ox_col_parent.scaffold_with_organellar.repetitive_k15.txt"),
    }
    CEN = {
        "col": {'Chr1': (14841147, 17216861), 'Chr2': (4621558, 6841935), 'Chr3': (13596351, 15826119),
                'Chr4': (5208113, 7982091), 'Chr5': (12402000, 15178500)},
        "ler": {'Chr1': (14538002, 16448388), 'Chr2': (3889714, 5401326), 'Chr3': (13992334, 17050927),
                'Chr4': (5394545, 7765831), 'Chr5': (12296068, 15662532)},
        "cenh3ox_col": {'Chr1': (14836227, 16964795), 'Chr2': (3925257, 6280744), 'Chr3': (13594954, 15842927),
                        'Chr4': (4786803, 7581705), 'Chr5': (12308011, 15067943)},
    }
    CHRLEN = {
        "col": {"Chr1": 32640075, "Chr2": 23012915, "Chr3": 26150667, "Chr4": 22582341, "Chr5": 30170985},
        "ler": {"Chr1": 32485061, "Chr2": 21328600, "Chr3": 27335240, "Chr4": 22700724, "Chr5": 30661135},
        "cenh3ox_col": {"Chr1": 32376493, "Chr2": 22314084, "Chr3": 26166303, "Chr4": 22192233, "Chr5": 30119503},
    }
    # centromere-restricted: iterate the CEN windows themselves
    SCAN = {rk: [(c, a, b) for c, (a, b) in d.items()] for rk, d in CEN.items()}
    CEN_IV = {rk: {c: [iv] for c, iv in d.items()} for rk, d in CEN.items()}

    def in_phase(svlen):
        """Whole-monomer (in-register, unequal-HR) vs out-of-phase (NHEJ-like)."""
        r = abs(svlen) % MONO
        return (r <= TOL or r >= MONO - TOL), r


def in_cen(refkey_, chrom, pos):
    """Does pos fall in a centromere/alpha-satellite interval of this reference? (overlay)"""
    for s, e in CEN_IV.get(refkey_, {}).get(chrom, []):
        if s <= pos < e:
            return True
    return False


os.makedirs(OUT, exist_ok=True)
