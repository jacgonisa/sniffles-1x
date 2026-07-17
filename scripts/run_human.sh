#!/bin/bash
# Human post-candidate steps (SM_GENOME=human). Run after 01_candidates finishes.
set -e
cd /mnt/ssd-4tb/HIFI_NAMIL/single_molecule_sv/scripts
export SM_GENOME=human
PY=/home/jg2070/miniforge3/envs/nextflow_env/bin/python
$PY 02_leadprov_sm.py
$PY 03_split_and_map.py        # CEN-only (bam_path_md = calmd cen.md.bam)
$PY 05_merge_classify.py       # no stock VCFs for human -> stock_match=0
$PY 07_normalize.py            # per-mapped-Mb (genome-wide)
$PY 17_source_breakdown.py
$PY 23_insertion_qc.py         # homopolymer + CCS-Q contrast + hiconf filter
$PY 06_report_human.py
echo ALL_HUMAN_DONE
