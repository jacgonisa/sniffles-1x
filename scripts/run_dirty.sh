#!/bin/bash
set -e
cd /mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts
export SYN_DIRTY=1
PY=/home/jg2070/miniforge3/envs/nextflow_env/bin/python
$PY S1_sim_synthetic.py
$PY S2_map_synthetic.py
export SM_SYNTHETIC=1
for s in 01_candidates 02_leadprov_sm 03_split_and_map 05_merge_classify 07_normalize 16_arm_control 18_arm_splitmap 23_insertion_qc; do
  echo "== $s =="; $PY $s.py
done
echo ALL_DIRTY_DONE
