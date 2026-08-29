# RPP1 / DM2 cluster — preliminary de-novo SV scan in F1 pollen (Col haplotype)

**Question (Ian):** the DM2/RPP1-like TNL cluster on Chr3 is structurally mutable between
accessions — do we see **de novo structural change in the pollen**? Also relevant to the
reviewer's ONT observation that TE/repeat-spanning reads carry "many apparent polymorphisms
with sharp boundaries" (?ONT artefact).

## Locus (Col-HiFi assembly, liftoff annotation)
Cluster window **Chr3:19,150,000–19,460,000** (~310 kb):

| gene | Chr3 coords | strand |
|---|---|---|
| AT3G44400 | 19,179,993–19,185,327 | − |
| AT3G44630 | 19,378,961–19,384,613 | + |
| AT3G44670 | 19,403,010–19,407,575 | + |
| **RPP1 / AT3G44480** | 19,432,477–19,438,182 | + |

## Method
Single-molecule SV detection (same detectors as the main pipeline: CIGAR I/D ≥ 50 bp + SA
split-read → `sv.classify_splits`) on every MAPQ≥10 read in the window. **Pollen = probe for
de-novo events; leaf = deep somatic baseline / paralog-artefact control** (a real de-novo pollen
event should be absent from the ~6× deeper leaf). Col haplotype only (Col-HiFi carries the
annotation and RPP1 is a Col gene). `scripts/rpp1_scan.py` → `results/rpp1/`.

## Findings — the Col cluster looks structurally **stable** in WT pollen

1. **Bulk copy number intact.** Mean depth across the cluster = the adjacent unique flank
   (pollen 58× vs 60×, ratio 0.96; leaf 453× vs 452×, ratio 1.00). No fixed gain/loss of array
   copies.
2. **Uniquely mappable in HiFi.** **0% of reads are MAPQ<10** in the cluster (both tissues) — winnowmap
   places these HiFi reads uniquely even in the TNL array. No paralog-mismap artefact storm, and we
   are not filtering out the interesting reads.
3. **No mutation/divergence hotspot.** Per-read gap-compressed divergence in the cluster ≈ a matched
   control region (pollen 0.00134 vs 0.00145; leaf 0.00161 vs 0.00155). **HiFi does NOT reproduce the
   ONT "extra polymorphisms in repeat-spanning reads"** → supports that being an ONT-specific artefact.
4. **Single-molecule SV rate = genome-wide arm background, not a hotspot.** Cluster DEL rate: pollen
   799/Mread vs arm background ~779; leaf 388 vs ~450. The whole 310 kb window carries only **4 pollen
   / 9 leaf** events (table below); all CIGAR indels are small (<600 bp) at gene edges and also present
   in the deep leaf.
5. **One weak de-novo candidate:** a single pollen read (`m84227…/32375450`) with a **1383 bp inversion
   junction** near AT3G44670 (4653 bp fwd @19,405,776 → 679 bp **rev** @19,411,129, both MAPQ 60), absent
   from leaf. Plus a second single-read 744 bp INV (MAPQ 49) near AT3G44630. **Both are single split
   reads — the artefact-prone class in repeats** (most likely a read traversing an inverted paralog, not
   a de-novo inversion). Need per-read validation before claiming anything.

**Bottom line:** at HiFi single-molecule resolution the Col RPP1/DM2 cluster is structurally quiet in
WT pollen — no bulk CNV, no divergence hotspot, SV rate at arm background, and only single-read
split-based inversion candidates that most likely reflect inverted paralogs. No dramatic de-novo
structural instability jumps out.

## Event catalogue

| tissue | pos | type | svlen | source | mapq | nearest gene |
|---|---|---|---|---|---|---|
| pollen | 19,194,222 | BND | – | SPLIT | 60 | AT3G44400~ |
| pollen | 19,312,843 | INV | 744 | SPLIT | 49 | AT3G44630~ |
| pollen | 19,380,913 | DEL | −54 | CIGAR | 60 | AT3G44630 |
| pollen | 19,410,428 | INV | 1383 | SPLIT | 60 | AT3G44670~ (candidate) |
| leaf | 19,175,989 | DEL | −83 | CIGAR | 60 | AT3G44400~ |
| leaf | 19,282,868 / 922 | INS | 310 / 536 | CIGAR | 60 | AT3G44630~ (one read) |
| leaf | 19,296,856 | DEL | −70 | CIGAR | 60 | AT3G44630~ |
| leaf | 19,426,121 / 19,440,518 / 19,455,092 | BND | – | SPLIT | 60 | RPP1~ |
| leaf | 19,435,198 | DEL | −59 | CIGAR | 60 | RPP1/AT3G44480 |
| leaf | 19,452,667 | INS | 508 | CIGAR | 60 | RPP1~ |

## Statistical test — CIGAR DEL/INS rate (the trustworthy channel)

`scripts/rpp1_cigar_test.py`: CIGAR-only DEL/INS (≥50 bp) per read in the cluster vs **60 matched
310-kb arm windows** tiled across all chromosome arms (pericentromere + telomeres excluded), Col hap.
Poisson-rate statistics.

**(1) Cluster vs surrounding arm regions — NOT significantly different (not a hotspot):**

| tissue | cluster rate /1k | arm median /1k | expected vs observed | Poisson p | percentile |
|---|---|---|---|---|---|
| leaf | 0.78 (6 events / 7725 reads) | 0.52 | 4.4 exp vs 6 obs | **0.56 (ns)** | 73rd |
| pollen | 0.80 (1 event / 1252 reads) | 1.58 | 1.6 exp vs 1 obs | **1.0 (ns)** | 32nd |

The cluster sits inside the arm distribution in both tissues (leaf slightly above median, pollen
below) — no evidence of elevated CIGAR-indel instability at the cluster.

**(2) Leaf vs pollen (Poisson rate-ratio test):**
- **In the cluster: no difference** — leaf 0.78 vs pollen 0.80 /1k, rate ratio 0.97, **p=1.0** (but only
  6 vs 1 events → underpowered at the cluster alone).
- **Genome-wide in arms: pollen is ~2.2× higher — highly significant.** leaf 0.57 vs pollen 1.26 /1k
  (n=247 vs 90 events), rate ratio 0.45, **p=1.5×10⁻⁹**; paired per-window Wilcoxon **p=4.4×10⁻⁶**.
  (Consistent with the pipeline's arm-control step 16, and robust to read length — pollen reads are
  *shorter*, so per-bp the enrichment is larger.)

So: the elevated pollen single-molecule indel rate is a **genome-wide arm phenomenon**, not specific to
RPP1 — and the RPP1 cluster is unremarkable against its own arm background. Figure:
`docs/rpp1_cigar_test.png`.

## Caveats / next steps
- **Col haplotype only.** The **Ler** RPP1 cluster (structurally different, and Ler-HiFi has no liftoff
  annotation yet) needs its own coordinates — map the RPP1 CDS to Ler-HiFi, then rerun. Ler is where
  accession-level structural difference is most likely.
- **WT only** here; CENH3ox pollen not scanned.
- Validate the two INV candidate reads (dotplot the read vs the locus; check whether the reverse
  segment is a paralog elsewhere) before treating as de-novo.
- Split-based INV/BND are unreliable in repeats (established elsewhere in this project); the CIGAR
  indels are the trustworthy channel and they show nothing unusual.
- To push sensitivity for **copy-number** change specifically, a per-read count of RPP1-monomer/TNL
  repeat units (TRASH-style, like the CEN178 register analysis) would detect array expansion/contraction
  a single CIGAR indel can miss.

Figure: `docs/rpp1_cluster.png` (event map, pollen vs leaf).
