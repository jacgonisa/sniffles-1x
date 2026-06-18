#!/bin/bash
# Step 4 — stock Sniffles2 cross-check, forced to single-molecule sensitivity.
# --minsupport 1 --mosaic --minsvlen 50, restricted to the centromere bed.
# Output: results/stock/{sample}_{hap}.vcf   (concordance reference for step 5)
set -euo pipefail
PY=/home/jg2070/miniforge3/envs/nextflow_env/bin
ROOT=/mnt/ssd-4tb/HIFI_NAMIL
OUT=$ROOT/single_molecule_sv/results
mkdir -p "$OUT/stock"

# centromere beds (Col / Ler) — same coords as common.py
cat > "$OUT/stock/cen_col.bed" <<'EOF'
Chr1	14841147	17216861
Chr2	4621558	6841935
Chr3	13596351	15826119
Chr4	5208113	7982091
Chr5	12402000	15178500
EOF
cat > "$OUT/stock/cen_ler.bed" <<'EOF'
Chr1	14538002	16448388
Chr2	3889714	5401326
Chr3	13992334	17050927
Chr4	5394545	7765831
Chr5	12296068	15662532
EOF

declare -A REF=(
  [col]=$ROOT/01_genomes/Col-HiFi/Col-0.ragtag_scaffolds_with_organellar.fa
  [ler]=$ROOT/01_genomes/Ler-HiFi/Ler-0.ragtag_scaffolds_with_organellar.fa
)

for sample in wt_leaf wt_pollen; do
  for hap in col ler; do
    bam=$ROOT/sv_calling/aligned/$sample/strict90/${hap}_all.bam
    echo "[$(date +%T)] sniffles $sample $hap"
    "$PY/sniffles" -i "$bam" -v "$OUT/stock/${sample}_${hap}.vcf" \
      --reference "${REF[$hap]}" \
      --regions "$OUT/stock/cen_${hap}.bed" \
      --minsupport 1 --mosaic --mosaic-af-min 0 --mosaic-include-germline --minsvlen 50 \
      --threads 16 2> "$OUT/stock/${sample}_${hap}.log" || \
      { echo "  sniffles FAILED (see log)"; tail -3 "$OUT/stock/${sample}_${hap}.log"; }
  done
done
echo DONE_STOCK
