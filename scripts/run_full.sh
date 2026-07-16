#!/bin/bash
# Full re-run after CENH3ox col was re-mapped to its own assembly (CENH3ox-Col-HiFi).
# ponytail: plain sequential; stops on first failure.
set -e
cd /mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts
PY=/home/jg2070/miniforge3/envs/nextflow_env/bin/python
BASEPY=/home/jg2070/miniforge3/bin/python

rm -f ../results/cen_read_counts.tsv          # stale (cenh3ox col counted vs WT-Col); rebuild below
$PY 01_candidates.py
$PY 02_leadprov_sm.py
$PY 03_split_and_map.py
$PY 05_merge_classify.py
$PY -c "import importlib; importlib.import_module('09_pptx_figures').cen_read_counts()"  # rebuild denominators
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
$PY 20_insertion_origin.py
$PY 21_insertion_origin_detailed.py
$PY 06_report.py
echo "ALL_FULL_DONE"
