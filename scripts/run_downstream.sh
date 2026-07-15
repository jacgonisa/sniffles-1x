#!/bin/bash
# Downstream pipeline after 03 (split-and-map) finishes, for all 4 samples.
# ponytail: plain sequential run; stops on first failure.
set -e
cd /mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts
PY=/home/jg2070/miniforge3/envs/nextflow_env/bin/python
BASEPY=/home/jg2070/miniforge3/bin/python   # for TRASH steps (13, 14)

$PY 05_merge_classify.py
$PY 07_normalize.py
$PY 08_read_qc.py
$PY 10_cen178_orient.py
$PY 11_recurrence.py
$PY 12_support_distribution.py
$BASEPY 13_annotate_singletons.py
$PY 15_translocations.py
$PY 16_arm_control.py
$PY 17_source_breakdown.py
$PY 18_arm_splitmap.py
$PY 09_pptx_figures.py
$BASEPY 14_read_validation.py
$PY 06_report.py
echo "ALL_DOWNSTREAM_DONE"
