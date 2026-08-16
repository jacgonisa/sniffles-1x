#!/bin/bash
# Synthetic-control pipeline: same scripts as real, SM_SYNTHETIC=1 -> reads/OUT redirected.
# Requires S1_sim_synthetic.py + S2_map_synthetic.py to have produced aligned_synthetic/ BAMs.
set -e
cd /mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts
export SM_SYNTHETIC=1
PY=/home/jg2070/miniforge3/envs/nextflow_env/bin/python
$PY 01_candidates.py
$PY 02_leadprov_sm.py
$PY 03_split_and_map.py
$PY 05_merge_classify.py
$PY 07_normalize.py
$PY 16_arm_control.py
$PY 18_arm_splitmap.py
$PY 23_insertion_qc.py
echo "ALL_SYNTHETIC_DONE"
