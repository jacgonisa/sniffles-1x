# Single-molecule centromere SV calling — WT leaf & pollen (Col + Ler)

Detects structural variants supported by a **single read** in the centromeres of
WT F1 (Col×Ler) leaf and pollen, on both haplotypes. Reproduces a colleague's
read-level / mapping-topology approach and additionally runs Sniffles2's own
per-read signal classifier at the single-molecule level.

## Inputs (reused, not regenerated)
- BAMs: `sv_calling/aligned/{wt_leaf,wt_pollen}/strict90/{col,ler}_all.bam`
  (winnowmap `-ax map-pb --MD`, k-mer haplotype-split, strict-90).
- References + winnowmap repetitive-kmer files: `01_genomes/{Col,Ler}-HiFi/`.
- CEN coords / register helper: mirrored in `scripts/common.py` (same values as
  `centromere_sv_pipeline/scripts/sv_compartment_phase.py`).
- Sniffles2 v2.7.5 source: imported from env `nextflow_env`
  (`sniffles.sv.classify_splits` is the shared topology classifier).

## Pipeline (run in order, env `nextflow_env`)
```
python scripts/01_candidates.py      # de>=0.005 | NM>=50 | SA  -> results/candidates/
python scripts/02_leadprov_sm.py     # cracked sniffles per-read -> results/leadprov_sm.tsv
python scripts/03_split_and_map.py   # MD-contrast split + remap -> results/split_and_map.tsv
bash   scripts/04_sniffles_stock.sh  # stock sniffles --minsupport 1 -> results/stock/
python scripts/05_merge_classify.py  # union+dedup+register+concord -> results/sm_sv_calls.tsv
python scripts/07_normalize.py       # read-Mb normalization -> results/sm_sv_rates.tsv
python scripts/08_read_qc.py         # read-quality controls (len/de/np/rq) -> results/read_qc.tsv
python scripts/10_cen178_orient.py   # CEN178 array orientation (minimap2) -> results/cen178_orientation.tsv
python scripts/11_recurrence.py      # recurrent-locus VAF (hotspot vs fixed) -> results/recurrent_loci.tsv
python scripts/12_support_distribution.py  # read-support per locus + 1x events (all-reads & read-budget-matched)
/home/jg2070/miniforge3/bin/python scripts/13_annotate_singletons.py  # TRASH in-register + confidence (BASE python)
python scripts/09_pptx_figures.py    # genome maps (+orientation strip) + karyograms -> results/figures/
python scripts/06_report.py          # -> report.html (run last; embeds 07-13 outputs)
```

## Comprehensive figures (`figures_pptx.html`, step 09)
Replicates `20260617_SV_analysis.pptx` (wt_leaf + wt_pollen, artf1 control not yet available):
genome maps per haplotype (≥5 kb SVs drawn as bars spanning the interval), size-binned
**count per million CEN reads** (col+ler pooled — large SVs enriched in pollen), and
log10(width) proportion histograms with the 178-bp line. Large events (DEL up to ~12.7 Mb,
INV ~11 Mb, DUP ~2.7 Mb) come from the split-read fragments. Denominator cached in
`results/cen_read_counts.tsv`.

`python scripts/02_leadprov_sm.py --selftest` checks the imported classifier
returns DEL/DUP/INV on synthetic split topologies.

## Method notes
- **Candidate filter** (colleague's): primary alignment with de≥0.005, NM≥50, or an SA tag.
- **leadprov (step 2)**: per candidate read, CIGAR I/D ≥50 bp (INLINE) + SA split reads
  fed to `sv.classify_splits` (SPLIT). No clustering / no min-support ⇒ each lead is one
  single-molecule SV. This is Sniffles' own algorithm with the coverage requirement removed.
- **split-and-map (step 3)**: the breakpoint is found from the per-base mismatch profile
  (MD tag) as the query position of maximal left/right substitution-rate contrast (≥0.01);
  the read is split there and both fragments re-mapped with winnowmap. An SV is called when
  both fragments ≥1 kb, MAPQ≥10, ref gap ≥50 bp; topology via the same `classify_splits`.
  Catches out-of-phase satellite events that stay a single linear alignment (no SA).
- **Output** `sm_sv_calls.tsv`: per-read calls with `methods` (CIGAR/SPLITREAD/SPLITMAP),
  CEN178 `in_phase`+`monomer_rem` (whole-monomer = unequal-HR signature), and `stock_match`
  (concordance with stock Sniffles2).

## Caveat
A lone ≥50 bp change in deep satellite coverage cannot be fully separated from a
mapping/sequencing artifact. The split-and-map re-mapping and the 178-bp register check
are mitigations; treat the callset as a sensitivity ceiling, not a confirmed somatic set.
