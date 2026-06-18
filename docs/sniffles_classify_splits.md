# The Sniffles2 code that classifies the topology

This is the **exact, unmodified** function we import and call for every read:
`sv.classify_splits` from **Sniffles2 v2.7.5** (`sniffles/sv.py`, line 620), by Moritz Smolka /
Hermann Romanek, <https://github.com/fritzsedlazeck/Sniffles>, **MIT-licensed**. It is reproduced
here only so the topology logic is visible in this repo — we do not fork or change it; we install
it as a dependency (`sniffles==2.7.5`) and `import` it.

**Where we call it:**
- [`scripts/02_leadprov_sm.py`](../scripts/02_leadprov_sm.py) — `sv.classify_splits(read, leads, CFG, contig)` on a read's primary + `SA` split alignments.
- [`scripts/03_split_and_map.py`](../scripts/03_split_and_map.py) — same call, on the two fragments we produce by splitting and re-mapping.

## How to read it (each branch → one SV type)

The function sorts the read's alignment fragments by **query** position (`qry_start`), then walks
adjacent fragment pairs (`last`, `curr`) and decides the type from their **strand** and their
**reference vs query** coordinate gaps. The decisive branches (the source keeps these as comments):

| branch in the code | condition (plain English) | emits |
|---|---|---|
| `# INS, FWD` / `# INS, REV` | same strand, **query** gap ≫ **reference** gap | `INS` |
| `# DEL, FWD` / `# DEL, REV` | same strand, **reference** gap ≫ **query** gap | `DEL` |
| `# DUP, FWD` / `# DUP, REV` | same strand, `curr` starts **before** `last` ends on the reference (overlap) | `DUP` |
| `# INV  CASE A–D` | the two fragments are on **opposite strands** | `INV` |
| `# BND` | the two fragments are on **different contigs** | `BND` |

`minsvlen_screen` is the size threshold (we set it to 50 bp). The returned leads carry
`svtypes_starts_lens = [(svtype, svstart, svlen), …]`, which our wrappers read out directly — one
entry = one single-molecule SV.

## The source (Sniffles2 2.7.5 · `sniffles/sv.py:620`)

```python
def classify_splits(read, leads, config, main_contig) -> list:
    """
    Determines the SV type of a split read (read with supplementary alignments). Returns (possibly changed) list of leads.
    """
    minsvlen_screen = config.minsvlen_screen
    min_split_len_bnd = config.bnd_min_split_length

    leads.sort(key=lambda ld: ld.qry_start)
    last = leads[0]
    last.svtypes_starts_lens = []
    hints = 0

    if last.qry_start >= config.long_ins_length * 0.5:
        last.svtypes_starts_lens.append(("INS", last.ref_start, None))

    for i in range(1, len(leads)):
        curr = leads[i]
        curr.svtypes_starts_lens = []

        if curr.contig == last.contig:
            rev = (curr.strand == "-")
            fwd = not rev
            if curr.strand == last.strand:
                #
                # INS, DEL, DUP
                #
                if (fwd and (curr.qry_start - last.qry_end) >= minsvlen_screen
                        and (curr.qry_start - last.qry_end) - (curr.ref_start - last.ref_end) >= minsvlen_screen):
                    # INS, FWD
                    svstart = curr.ref_start
                    svlen = (curr.qry_start - last.qry_end)
                    if svlen <= config.dev_seq_cache_maxlen:
                        curr.seq = read.query_sequence[last.qry_end:curr.qry_start]
                    else:
                        curr.seq = None
                    curr.svtypes_starts_lens.append(("INS", svstart, svlen))
                    hints += 1

                elif (rev and (curr.qry_start - last.qry_end) >= minsvlen_screen
                      and (curr.qry_start - last.qry_end) - (last.ref_start - curr.ref_end) >= minsvlen_screen):
                    # INS, REV
                    svstart = last.ref_start
                    svlen = (curr.qry_start - last.qry_end)
                    if svlen <= config.dev_seq_cache_maxlen:
                        curr.seq = read.query_sequence[last.qry_end:curr.qry_start]
                    else:
                        curr.seq = None
                    curr.svtypes_starts_lens.append(("INS", svstart, svlen))
                    hints += 1

                elif (fwd and (curr.ref_start - last.ref_end) >= minsvlen_screen
                      and (curr.ref_start - last.ref_end) - (curr.qry_start - last.qry_end) >= minsvlen_screen):
                    # DEL, FWD
                    svstart = curr.ref_start
                    svlen = (curr.ref_start - last.ref_end)
                    curr.svtypes_starts_lens.append(("DEL", svstart, -svlen))
                    hints += 1

                elif (rev and (last.ref_start - curr.ref_end) >= minsvlen_screen
                      and (last.ref_start - curr.ref_end) - (curr.qry_start - last.qry_end) >= minsvlen_screen):
                    # DEL, REV
                    svstart = last.ref_start
                    svlen = (last.ref_start - curr.ref_end)
                    curr.svtypes_starts_lens.append(("DEL", svstart, -svlen))
                    hints += 1

                elif fwd and curr.ref_start <= last.ref_end:
                    # DUP, FWD
                    svstart = curr.ref_start
                    svlen = (last.ref_end - curr.ref_start)
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("DUP", svstart, svlen))
                        hints += 1

                elif rev and last.ref_start <= curr.ref_end:
                    # DUP, REV
                    svstart = last.ref_start
                    svlen = (curr.ref_end - last.ref_start)
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("DUP", svstart, svlen))
                        hints += 1

            else:
                #
                # INV
                #
                if fwd and curr.ref_start <= last.ref_start:
                    # CASE B
                    svstart = curr.ref_start
                    svlen = last.ref_start - curr.ref_start
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV", svstart, svlen))
                        hints += 1

                elif fwd and curr.ref_start > last.ref_start:
                    # CASE C
                    svstart = last.ref_start
                    svlen = curr.ref_start - last.ref_start
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV", svstart, svlen))
                        hints += 1

                elif rev and curr.ref_end >= last.ref_end:
                    # CASE A
                    svstart = last.ref_end
                    svlen = curr.ref_end - last.ref_end
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV", svstart, svlen))
                        hints += 1

                elif rev and curr.ref_end < last.ref_end:
                    # CASE D
                    svstart = curr.ref_end
                    svlen = last.ref_end - curr.ref_end
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV", svstart, svlen))
                        hints += 1
        else:
            #
            # BND
            #
            if curr.contig == main_contig:
                a, b = curr, last
            else:
                a, b = last, curr

            if a.contig == main_contig and abs(last.qry_end - last.qry_start) >= min_split_len_bnd and abs(curr.qry_end - curr.qry_start) >= min_split_len_bnd:
                is_first = a.qry_start < b.qry_start
                if is_first:
                    if a.strand == "+":
                        svstart = a.ref_end
                    else:
                        svstart = a.ref_start
                else:
                    if a.strand == "+":
                        svstart = a.ref_start
                    else:
                        svstart = a.ref_end
                a.svtypes_starts_lens.append(("BND",
                                              svstart,
                                              SVCallBNDInfo(b.contig,
                                                            b.ref_start,
                                                            is_first,
                                                            a.strand != b.strand)))
                hints += 1
        last = curr

    if not hints and len(leads) > 2:
        # filter out short pseudo BND hints
        left = leads[0]
        leads = [ld for ld in leads if ld.contig == left.contig and ld.strand == left.strand]
        if len(leads) == 2:
            return classify_splits(read, leads, config, main_contig)

    return leads
```

> Verbatim from Sniffles2 2.7.5 (MIT). The only thing this project changes is *what feeds the
> function* (a single read's fragments — including fragments we create by split-and-map) and the
> fact that we **do not** run Sniffles' downstream `cluster.py` / coverage-QC, which is what would
> otherwise require ≥2 reads. See [`../ALGORITHM.md`](../ALGORITHM.md).
