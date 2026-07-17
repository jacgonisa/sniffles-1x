#!/bin/bash
# Add MD tags to the human haplotype BAMs, restricted to the alpha-satellite CEN regions
# (all that step 03 split-and-remap needs). MD is required by best_split's
# get_aligned_pairs(with_seq=True). CEN-only keeps this fast vs calmd on the full 17 GB BAM.
set -euo pipefail
SAM=/home/jg2070/miniforge3/envs/nextflow_env/bin/samtools
BD=/mnt/ssd-8tb/HUMAN/sv_calling/strict90_sniffles1x_BLS0005_BLS0006/bam
G=/mnt/ssd-8tb/HUMAN/data/assembly/genomes
A=/mnt/ssd-8tb/HUMAN/data/assembly/annotation/cen_arms
for hap in MAT PAT; do
  in=$BD/BLS0005_BLS0006_sperm.$hap.strict90.merged.bam
  ref=$G/$hap.fasta
  bed=$A/hg002v1.1.$hap.alpha_CEN.bed
  out=$BD/BLS0005_BLS0006_sperm.$hap.strict90.cen.md.bam
  echo "[$hap] alpha_CEN reads + calmd -> $(basename $out)"
  $SAM view -@4 -b -L "$bed" "$in" | $SAM calmd -@4 -b - "$ref" 2>/dev/null > "$out"
  $SAM index "$out"
  echo "[$hap] done: $($SAM view -c -F0x900 "$out") primary reads"
done
echo PREP_MD_DONE
