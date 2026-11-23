#!/usr/bin/env python3
#-----------------------------------------------------------------------------
# peek.py
#-----------------------------------------------------------------------------
# Author:  Reza Mozafari
# Date:    2026-11-22
# Repo:    https://github.com/rmozafari/peptide-ngs-counterselection
#
# Summary: Sanity-check a single FASTQ before committing to a full run: flank
#          recovery, read orientation, insert-length distribution, stop-codon
#          rate and positive-control detection.
#-----------------------------------------------------------------------------
"""
peek.py - sanity-check a single real FASTQ before committing to a full run.

Reads the first N records of one file and reports whether the library actually
looks the way Library.docx says it does: are the flanks present, which
orientation, what insert lengths, is the positive control there.

    python3 peek.py SOME_R1_001.fastq.gz
    python3 peek.py SOME_R1_001.fastq.gz 500000
"""
import gzip, sys, os
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from avb3_ngs import (LEFT_FLANK, RIGHT_FLANK, CONTROL_DNA, CONTROL_PEPTIDE,
                      reverse_complement, translate, INSERT_LENGTHS)

path = sys.argv[1]
limit = int(sys.argv[2]) if len(sys.argv) > 2 else 200000

ANCHORS = [8, 14, 24]
hits = {a: 0 for a in ANCHORS}
orient = Counter()
lengths = Counter()
peptides = Counter()
read_lengths = Counter()
control = 0
n = 0

print(f"Reading up to {limit:,} reads from {os.path.basename(path)}\n")

with gzip.open(path, "rt", errors="replace") as fh:
    for i, line in enumerate(fh):
        if i % 4 != 1:
            continue
        n += 1
        if n > limit:
            break
        seq = line.strip().upper()
        read_lengths[len(seq)] += 1

        for a in ANCHORS:
            la, ra = LEFT_FLANK[-a:], RIGHT_FLANK[:a]
            for s, o in ((seq, "fwd"), (reverse_complement(seq), "rev")):
                p = s.find(la)
                if p == -1:
                    continue
                q = s.find(ra, p + a)
                if q == -1:
                    continue
                L = q - (p + a)
                if L in INSERT_LENGTHS:
                    hits[a] += 1
                    if a == 14:
                        orient[o] += 1
                        lengths[L] += 1
                        ins = s[p + a:q]
                        pep = translate(ins)
                        peptides[pep] += 1
                        if ins == CONTROL_DNA:
                            control += 1
                    break

n = min(n, limit)
print(f"Read length: {', '.join(f'{k} nt ({v/n*100:.0f}%)' for k, v in read_lengths.most_common(3))}")
print(f"Amplicon is {len(LEFT_FLANK) + 30 + len(RIGHT_FLANK)} nt at most, so a read "
      f"of {max(read_lengths)} nt {'spans it' if max(read_lengths) >= 99 else 'DOES NOT span it'}\n")

print("Flank recovery by anchor length:")
for a in ANCHORS:
    print(f"  {a:>2} nt/side: {hits[a]:>8,} / {n:,}  ({hits[a]/n*100:5.1f}%)")

if hits[14] == 0:
    print("\n*** NO READS MATCHED ***")
    print("The constant flanks are not where Library.docx says. Before anything else,")
    print("look at a few raw reads and find the real flanking sequence:\n")
    with gzip.open(path, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                print("   ", line.strip()[:120])
            if i > 20:
                break
    sys.exit(1)

print(f"\nOrientation: {dict(orient)}  "
      f"({orient['rev']/(orient['fwd']+orient['rev'])*100:.1f}% reverse-complement)")

print("\nInsert length distribution:")
tot = sum(lengths.values())
for L in sorted(INSERT_LENGTHS):
    c = lengths.get(L, 0)
    bar = "#" * int(c / tot * 40) if tot else ""
    print(f"  {L:>2} nt = {INSERT_LENGTHS[L]:>2} aa: {c:>8,} ({c/tot*100:5.1f}%) {bar}")

stops = sum(v for k, v in peptides.items() if "*" in k)
# NNK restricts the third base to G or T, so TAA and TGA are impossible and TAG
# is the only stop: exactly 1 of 32 codons. The chance a peptide contains at
# least one therefore grows with length, 1 - (31/32)^L, and is 27% for a 10-mer.
# Compare against the observed length mix, not against a flat "few percent".
exp_stop = sum((lengths.get(nt, 0) / tot) * (1 - (31 / 32) ** aa)
               for nt, aa in INSERT_LENGTHS.items()) if tot else 0
print(f"\nStop codons: {stops:,} of {tot:,} ({stops/tot*100:.1f}%)")
print(f"  expected for this length mix: {exp_stop*100:.1f}%   "
      f"[NNK: TAG is 1 of 32 codons, so 1-(31/32)^L per peptide]")
if exp_stop:
    ratio = (stops / tot) / exp_stop
    verdict = ("as expected" if 0.8 <= ratio <= 1.2 else
               "LOWER than expected - selection may have depleted truncated clones"
               if ratio < 0.8 else
               "HIGHER than expected - check for frameshifts or synthesis errors")
    print(f"  observed/expected = {ratio:.2f}  ({verdict})")

print(f"\nPositive control {CONTROL_PEPTIDE} / M{CONTROL_PEPTIDE}")
print(f"  exact DNA {CONTROL_DNA}: {control:,} reads")
pep_hits = peptides.get(CONTROL_PEPTIDE, 0)
print(f"  peptide level (any encoding):      {pep_hits:,} reads")
if pep_hits == 0:
    wrong = peptides.get("RGDWAGC", 0)
    print("  NOT FOUND in this sample. Either it is rarer than this subsample,")
    print("  or the control construct does not carry the same flanks.")
    print(f"  (the incorrect 'RGDWAGC' from the original brief: {wrong:,} reads)")

print(f"\nTop 10 peptides in this subsample:")
for pep, c in peptides.most_common(10):
    print(f"  {('M'+pep):<14} {len(pep):>2} aa  {c:>7,}  ({c/tot*100:.3f}%)")
