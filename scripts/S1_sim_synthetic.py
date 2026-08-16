#!/usr/bin/env python3
"""Synthetic-control step 1 — simulate HiFi reads FROM the assemblies (no variants), mirroring the
real 4-group structure. Reads are simulated only from the analysed regions (CEN windows + arm-control
windows) to keep it light; S2 maps them to the FULL assembly so satellite cross-mapping is preserved.

Simulator: badread (pacbio2021 HiFi model, identity ~99.8% ≈ real de; no junk/chimeras). One streamed
command per (sample,hap) → fastq.gz (no giant intermediates). read length ~ real tissue profile.
-> single_molecule_sv/synthetic_reads/<sample>_<hap>.hifi.fastq.gz
Run inside the readsim conda env (has badread); e.g.:
  conda run -n readsim env SM_GENOME=arabidopsis python S1_sim_synthetic.py [sample_hap]"""
import os, subprocess, sys, glob
from common import SAMPLES, HAPS, refkey, REF, CEN, CHRLEN

BADREAD = "/home/jg2070/miniforge3/envs/readsim/bin/badread"
SAM = "/home/jg2070/miniforge3/envs/nextflow_env/bin/samtools"
SYN = "/mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/synthetic_reads"
ARM_BUFFER, ARM_WIN = 5_000_000, 3_000_000
DEPTH = 15                       # per-Mb rate is depth-invariant; 15x is plenty and fast
LEN = {"leaf": (18000, 3500), "pollen": (14000, 3500)}
# clean (mapping floor only) vs dirty (mapping+sequencing floor, error tail matched to the real reads:
# real de median 0.0007 / mean 0.0015 / p90 0.004  ->  identity mean 99.85 with a real low tail + glitches)
DIRTY = os.environ.get("SYN_DIRTY") == "1"
IDENT = "99.85,100,0.5" if DIRTY else "99.8,100,0.3"
GLITCH = "10000,25,25" if DIRTY else "0,0,0"
SFX = ".dirty" if DIRTY else ""


def regions(rk):
    r = []
    for c, (a, b) in CEN[rk].items():
        r.append(f"{c}:{a}-{b}")
        s = b + ARM_BUFFER; e = min(s + ARM_WIN, CHRLEN[rk][c] - 100_000)
        if e - s > 200_000:
            r.append(f"{c}:{int(s)}-{int(e)}")
    return r


def one(sample, hap, tis):
    rk = refkey(sample, hap); ref = REF[rk][0]
    base = f"{SYN}/{sample}_{hap}"
    regfa = f"{base}.regions.fa"; out = f"{base}{SFX}.hifi.fastq.gz"
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        print(f"  {sample} {hap}: exists, skip"); return
    subprocess.run(f"{SAM} faidx {ref} {' '.join(regions(rk))} > {regfa}",
                   shell=True, check=True, executable="/bin/bash")
    mean, sd = LEN[tis]
    cmd = (f"{BADREAD} simulate --reference {regfa} --quantity {DEPTH}x --length {mean},{sd} "
           f"--identity {IDENT} --error_model pacbio2021 --qscore_model pacbio2021 "
           f"--junk_reads 0 --random_reads 0 --chimeras 0 --glitches {GLITCH} --seed 7 2>/dev/null "
           f"| gzip > {out}")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    n = int(subprocess.check_output(f"zcat {out} | wc -l", shell=True)) // 4
    print(f"  {sample} {hap}: {n} HiFi reads -> {os.path.basename(out)}")
    for f in glob.glob(f"{regfa}*"):
        os.remove(f)


def main():
    os.makedirs(SYN, exist_ok=True)
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for sample, tis in SAMPLES:
        for hap in HAPS:
            if only and only != f"{sample}_{hap}":
                continue
            one(sample, hap, tis)
    print("DONE_SIM")


if __name__ == "__main__":
    main()
