# How position, length and breakpoints are calculated

Every call is one read. Coordinates are **0-based** reference positions on the contig the read
(or its main fragment) is anchored to. `svlen` sign follows Sniffles: **negative = deletion**,
**positive = insertion / duplicated / inverted span**. Output columns: `chrom  pos  svtype  svlen
… mate`. The two detectors compute them as follows.

## A. CIGAR (INLINE) indels — `scripts/02_leadprov_sm.py::cigar_leads`

Walk the read's CIGAR, tracking the reference position `pos_ref`. For an `I`/`D` op ≥ 50 bp:

| type | `pos` (reported start) | `svlen` | reference footprint | length source |
|------|------------------------|---------|---------------------|---------------|
| **INS** | `pos_ref` (the ref base at the insertion point) | `+oplen` | none (zero ref span — the inserted bases live in the **read**) | the `I` op length (query) |
| **DEL** | `pos_ref + oplen` (ref base **after** the deletion) | `-oplen` | `[pos + svlen, pos]` = `[pos_ref, pos_ref + oplen]` | the `D` op length (reference) |

This matches Sniffles' `read_iterindels` (DEL `ref_start = pos_ref + oplength`). So for a deletion,
**`pos` is the coordinate immediately downstream of the gap**, and the deleted interval is
`pos − |svlen|` … `pos`. (If you'd prefer `pos` = deletion *start*, it's a one-line change — say so.)

- **end** = `pos` for INS (no ref span); `pos + |svlen|` for DEL (i.e. `pos_ref + oplen`, but recall `pos` already *is* that — the deletion spans `pos−|svlen| … pos`).
- `pos_ref` only advances on `M/=/X/D/N` (ref-consuming ops); `I` and soft-clips advance the read, not the reference.

## B. Split topology (DEL/DUP/INV/INS/BND) — `sv.classify_splits` (Sniffles2, unmodified)

Used by both the native split-read path (`02`, fragments from the aligner's `SA` tag) and the
split-and-map path (`03`, the two fragments we re-map). Fragments are sorted by **query** start;
adjacent pairs (`last`, `curr`) give `pos` and `svlen` directly from Sniffles
([annotated source](sniffles_classify_splits.md)):

| type | condition | `svstart` (=`pos`) | `svlen` |
|------|-----------|--------------------|---------|
| **INS** (fwd) | query gap ≫ ref gap | `curr.ref_start` | `curr.qry_start − last.qry_end` (the query gap) |
| **DEL** (fwd) | ref gap ≫ query gap | `curr.ref_start` | `−(curr.ref_start − last.ref_end)` (the ref gap) |
| **DUP** (fwd) | `curr.ref_start ≤ last.ref_end` (ref overlap) | `curr.ref_start` | `last.ref_end − curr.ref_start` (overlap length) |
| **INV** | fragments on opposite strands | one fragment's `ref_start`/`ref_end` (4 cases) | distance between the two fragments' facing ref ends |
| **BND** | fragments on **different contigs** | breakpoint on the main contig (`a.ref_end`/`ref_start` by strand) | — (no length; a junction, see §C) |

`*.qry_start/qry_end` are the fragment's span in **read** coordinates (forward-read frame; for the
split-and-map fragments the second half is offset by the split point so the gaps are measured in the
original read). `*.ref_start/ref_end` are the fragment's reference span. **end** = `pos + |svlen|`
for DEL/DUP/INV; `pos` for INS.

> Sizes ≥ 50 bp only (`minsvlen_screen = 50`). The `±50 bp` / `≥1 kb` / `MAPQ ≥ 10` gates in
> split-and-map (`03`) are applied to the *fragments* before this step.

## C. Translocations = the BND class — **yes, they are detected**

`classify_splits` returns **BND** whenever a read's two fragments land on **different contigs**.
This is the inter-chromosomal / translocation class. For a BND we store the partner locus in the
**`mate`** column as `mate_contig:mate_ref_start` (from Sniffles' `SVCallBNDInfo`); `svlen` is empty
(a junction has no single length). `pos` is the breakpoint on the read's main (centromere) contig.

**In these centromeres BND is the least trustworthy class** — the 5 centromeres share the same
CEN178 satellite, so a fragment from one CEN can mis-map to another CEN (or to an unplaced/organellar
contig). Of the native split-read BNDs, the mate falls in *another centromere* ~20% of the time and
on *unplaced/organellar* contigs ~14% — both consistent with satellite cross-mapping rather than true
translocations; the rest land on other-chromosome arms (possible real junctions, but unconfirmable
from a single read). Treat single-molecule BND as candidates only, and filter on the `mate` column.
A summary is written to `results/translocations.tsv` by `scripts/15_translocations.py`.
