#!/usr/bin/env python3
"""Shared paths/constants for the single-molecule centromere SV pipeline.
WT only; leaf + pollen; Col + Ler haplotypes. CEN coords reused from the
existing compartment classifier (sv_compartment_phase.py)."""
import os

ROOT = "/mnt/ssd-4tb/HIFI_NAMIL"
OUT = f"{ROOT}/single_molecule_sv/results"

# (sample, tissue). sample name is genotype_tissue, so `sample` doubles as the
# genotype×tissue group key; genotype(sample) recovers the genotype axis.
SAMPLES = [("wt_leaf", "leaf"), ("wt_pollen", "pollen"),
           ("cenh3ox_leaf", "leaf"), ("cenh3ox_pollen", "pollen")]
HAPS = ["col", "ler"]

# fixed display order for the 4 sample groups (WT then CENH3ox, leaf then pollen)
GROUPS = ["wt_leaf", "cenh3ox_leaf", "wt_pollen", "cenh3ox_pollen"]

def genotype(sample):
    return "cenh3ox" if "cenh3ox" in sample else "wt"

def bam_path(sample, hap):
    return f"{ROOT}/sv_calling/aligned/{sample}/strict90/{hap}_all.bam"

# reference + winnowmap repetitive-kmer file per haplotype (for split-and-map remap)
REF = {
    "col": (f"{ROOT}/01_genomes/Col-HiFi/Col-0.ragtag_scaffolds_with_organellar.fa",
            f"{ROOT}/01_genomes/Col-HiFi/Col-0.ragtag_scaffolds_with_organellar.repetitive_k15.txt"),
    "ler": (f"{ROOT}/01_genomes/Ler-HiFi/Ler-0.ragtag_scaffolds_with_organellar.fa",
            f"{ROOT}/01_genomes/Ler-HiFi/Ler-0.ragtag_scaffolds_with_organellar.repetitive_k15.txt"),
}

# centromere boundaries (same as centromere_sv_pipeline/scripts/sv_compartment_phase.py)
CEN = {
    "col": {'Chr1': (14841147, 17216861), 'Chr2': (4621558, 6841935), 'Chr3': (13596351, 15826119),
            'Chr4': (5208113, 7982091), 'Chr5': (12402000, 15178500)},
    "ler": {'Chr1': (14538002, 16448388), 'Chr2': (3889714, 5401326), 'Chr3': (13992334, 17050927),
            'Chr4': (5394545, 7765831), 'Chr5': (12296068, 15662532)},
}

MONO = 178  # CEN178 monomer
TOL = 20

def in_phase(svlen):
    """Whole-monomer (in-register, unequal-HR) vs out-of-phase (NHEJ-like)."""
    r = abs(svlen) % MONO
    return (r <= TOL or r >= MONO - TOL), r

def in_cen(hap, chrom, pos):
    cc = CEN[hap].get(chrom)
    return cc is not None and cc[0] <= pos < cc[1]

os.makedirs(OUT, exist_ok=True)
