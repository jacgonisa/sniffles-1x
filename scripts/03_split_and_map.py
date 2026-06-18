#!/usr/bin/env python3
"""Step 3 — colleague's split-and-map (the novel detector).

For each candidate read: from the per-base match/mismatch profile (MD tag, via
get_aligned_pairs(with_seq=True)) find the query position where left-vs-right
substitution-rate contrast is maximal. If contrast >= 0.01 split the read there,
re-map both fragments independently with winnowmap, and call an SV when both
fragments are >= 1 kb, MAPQ >= 10, and the inter-fragment ref gap is >= 50 bp.
Topology (DEL/DUP/INV/INS/BND) is decided by Sniffles' sv.classify_splits, run on
the two remapped fragments -> identical logic to step 2, different evidence.

-> results/split_and_map.tsv  (source=SPLITMAP)
Run with nextflow_env python (pysam + winnowmap + samtools).
Catches out-of-phase satellite events that stay a single linear alignment.
"""
import sys, os, subprocess, pysam
from types import SimpleNamespace
sys.path.insert(0, "/home/jg2070/miniforge3/envs/nextflow_env/lib/python3.13/site-packages")
from sniffles import sv
from sniffles.leadprov import Lead
from common import SAMPLES, HAPS, CEN, REF, bam_path, OUT

CFG = SimpleNamespace(minsvlen_screen=50, long_ins_length=2500,
                      bnd_min_split_length=1000, dev_seq_cache_maxlen=0)
WIN = "/home/jg2070/miniforge3/envs/nextflow_env/bin/winnowmap"
ST = "/home/jg2070/miniforge3/envs/nextflow_env/bin/samtools"
THREADS = 16
MIN_FRAG = 1000          # each fragment >= 1 kb
CONTRAST_MIN = 0.01      # min left/right substitution-rate contrast
MAPQ_MIN = 10
MIN_READLEN = 2 * MIN_FRAG + 200


def best_split(read):
    """Return (qsplit, contrast) of the max substitution-rate-contrast point, or None.
    qsplit is a coordinate into read.query_sequence; both sides must have >= MIN_FRAG
    aligned query bases."""
    qpos, mism = [], []
    for q, r, base in read.get_aligned_pairs(with_seq=True):
        if q is None or r is None or base is None:
            continue                       # skip indels / clips
        qpos.append(q)
        mism.append(1 if base.islower() else 0)   # with_seq lowercases mismatches
    n = len(qpos)
    if n < 200:
        return None
    # prefix sum of mismatches; scan interior split points keeping >=MIN_FRAG query bp/side
    pre = [0] * (n + 1)
    for i in range(n):
        pre[i + 1] = pre[i] + mism[i]
    best = None
    q0, q1 = qpos[0], qpos[-1]
    for i in range(1, n):
        if qpos[i] - q0 < MIN_FRAG or q1 - qpos[i] < MIN_FRAG:
            continue
        lrate = pre[i] / i
        rrate = (pre[n] - pre[i]) / (n - i)
        c = abs(rrate - lrate)
        if best is None or c > best[1]:
            best = (qpos[i], c)
    if best is None or best[1] < CONTRAST_MIN:
        return None
    return best


def extract(sample, hap):
    """Write per-read fragment FASTA + metadata; return path, meta dict."""
    cand = set()
    with open(f"{OUT}/candidates/{sample}_{hap}.tsv") as f:
        next(f)
        for ln in f:
            cand.add(ln.split("\t", 1)[0])
    bam = pysam.AlignmentFile(bam_path(sample, hap), "rb")
    fa = f"{OUT}/splitmap/{sample}_{hap}.frags.fa"
    meta = {}                              # idx -> (qname, chrom, qsplit, readlen)
    os.makedirs(f"{OUT}/splitmap", exist_ok=True)
    idx = 0
    with open(fa, "w") as out:
        for chrom, (a, b) in CEN[hap].items():
            for r in bam.fetch(chrom, a, b):
                if r.is_secondary or r.is_supplementary or r.is_unmapped:
                    continue
                if r.query_name not in cand or not (a <= r.reference_start < b):
                    continue
                seq = r.query_sequence
                if seq is None or len(seq) < MIN_READLEN:
                    continue
                bs = best_split(r)
                if bs is None:
                    continue
                qsplit, _c = bs
                if qsplit < MIN_FRAG or len(seq) - qsplit < MIN_FRAG:
                    continue
                meta[idx] = (r.query_name, chrom, qsplit, len(seq))
                out.write(f">{idx}_A\n{seq[:qsplit]}\n>{idx}_B\n{seq[qsplit:]}\n")
                idx += 1
    bam.close()
    print(f"{sample} {hap}: {idx} reads split -> {fa}")
    return fa, meta


def remap(fa, hap):
    ref, rep = REF[hap]
    bam = fa.replace(".fa", ".bam")
    cmd = (f"{WIN} -W {rep} -ax map-pb -t {THREADS} {ref} {fa} 2>/dev/null "
           f"| {ST} sort -@4 -o {bam} - && {ST} index {bam}")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    return bam


def classify(sample, hap, tis, meta, frag_bam, out):
    bam = pysam.AlignmentFile(frag_bam, "rb")
    halves = {}                            # idx -> {'A':aln,'B':aln}
    for r in bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        if r.query_alignment_length < MIN_FRAG or r.mapping_quality < MAPQ_MIN:
            continue
        idx_s, half = r.query_name.rsplit("_", 1)
        halves.setdefault(int(idx_s), {})[half] = r
    bam.close()
    n = 0
    for idx, hv in halves.items():
        if "A" not in hv or "B" not in hv:
            continue
        qname, main_contig, qsplit, _rl = meta[idx]
        A, B = hv["A"], hv["B"]
        leads = [
            Lead(0, qname, A.reference_name, A.reference_start, A.reference_end,
                 A.query_alignment_start, A.query_alignment_end,
                 "-" if A.is_reverse else "+", A.mapping_quality, 0, "SPLITMAP", "?"),
            Lead(0, qname, B.reference_name, B.reference_start, B.reference_end,
                 qsplit + B.query_alignment_start, qsplit + B.query_alignment_end,
                 "-" if B.is_reverse else "+", B.mapping_quality, 0, "SPLITMAP", "?"),
        ]
        for ld in sv.classify_splits(_DUMMY, leads, CFG, main_contig):
            for svtype, svstart, arg in (ld.svtypes_starts_lens or []):
                if svtype == "NOSV":
                    continue
                if svtype == "BND":
                    svlen = ""; mate = f"{arg.mate_contig}:{arg.mate_ref_start}"
                else:
                    svlen = arg; mate = ""
                a, b = CEN[hap][main_contig] if main_contig in CEN[hap] else (0, 0)
                if a <= svstart < b:
                    mq = min(A.mapping_quality, B.mapping_quality)
                    out.write(f"{sample}\t{hap}\t{tis}\t{main_contig}\t{svstart}\t{svtype}\t{svlen}\tSPLITMAP\t{mq}\t{qname}\t{mate}\n")
                    n += 1
    print(f"{sample} {hap}: {n} split-and-map SV calls")


class _Dummy:
    query_sequence = ""
_DUMMY = _Dummy()


def main():
    out = open(f"{OUT}/split_and_map.tsv", "w")
    out.write("sample\thap\ttissue\tchrom\tpos\tsvtype\tsvlen\tsource\tmapq\tread\tmate\n")
    for sample, tis in SAMPLES:
        for hap in HAPS:
            fa, meta = extract(sample, hap)
            if not meta:
                continue
            frag_bam = remap(fa, hap)
            classify(sample, hap, tis, meta, frag_bam, out)
    out.close(); print("DONE_SPLITMAP")


if __name__ == "__main__":
    main()
