#!/usr/bin/env python3
"""Synthetic-control step 2 — map the simulated HiFi reads back to the FULL assembly with the SAME
winnowmap settings as the real BAMs / step 03 (-ax map-pb --MD). Reads were simulated only from CEN+arm
windows but are mapped to the whole reference, so satellite cross-mapping to other centromeres (a real
artefact) is preserved.
-> sv_calling/aligned_synthetic/<sample>/strict90/<hap>_all.bam
Run: python S2_map_synthetic.py (arabidopsis profile)."""
import os, subprocess, sys
from common import SAMPLES, HAPS, refkey, REF

WIN = "/home/jg2070/miniforge3/envs/nextflow_env/bin/winnowmap"
SAM = "/home/jg2070/miniforge3/envs/nextflow_env/bin/samtools"
SYN = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/synthetic_reads"
SFX = ".dirty" if os.environ.get("SYN_DIRTY") == "1" else ""
ALN = "/mnt/ssd-4tb/HIFI_NAMIL/sv_calling/aligned_synthetic" + ("_dirty" if SFX else "")
THREADS = 16


def one(sample, hap):
    rk = refkey(sample, hap); ref, rep = REF[rk]
    fq = f"{SYN}/{sample}_{hap}{SFX}.hifi.fastq.gz"
    if not os.path.exists(fq):
        print(f"  {sample} {hap}: no fastq, skip"); return
    outdir = f"{ALN}/{sample}/strict90"; os.makedirs(outdir, exist_ok=True)
    bam = f"{outdir}/{hap}_all.bam"
    cmd = (f"{WIN} -W {rep} -ax map-pb --MD -t {THREADS} {ref} {fq} 2>/dev/null "
           f"| {SAM} sort -@4 -o {bam} - && {SAM} index {bam}")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    n = int(subprocess.check_output(f"{SAM} view -c -F0x900 {bam}", shell=True))
    print(f"  {sample} {hap}: mapped {n} primary -> {bam}")


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for sample, _ in SAMPLES:
        for hap in HAPS:
            if only and only != f"{sample}_{hap}":
                continue
            one(sample, hap)
    print("DONE_MAP_SYNTHETIC")


if __name__ == "__main__":
    main()
