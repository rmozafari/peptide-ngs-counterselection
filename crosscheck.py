#!/usr/bin/env python3
#-----------------------------------------------------------------------------
# crosscheck.py
#-----------------------------------------------------------------------------
# Author:  Reza Mozafari
# Date:    2026-11-22
# Repo:    https://github.com/rmozafari/peptide-ngs-counterselection
#
# Summary: Independent re-implementation of the extraction and counting, written
#          against the library spec rather than against avb3_ngs.py, for
#          comparing counts between two implementations that share no code.
#-----------------------------------------------------------------------------
"""
crosscheck.py - deliberately naive, independent re-implementation.

Written against the library spec rather than against avb3_ngs.py: R1 only, no
pairing, full-length flank matching, regex instead of index arithmetic, its own
codon table built by hand. If two implementations that share no code agree on
the counts, the extraction logic is probably right.
"""
import gzip, re, sys, json
from collections import Counter

LEFT = "TAATACGACTCACTATAGGGACAATTACTATTTACAATTACAATG"
RIGHT = "GGTAGTGGTACGTCTGGCTCCAGC"
PAT = re.compile(LEFT + r"([ACGT]{18}|[ACGT]{21}|[ACGT]{24}|[ACGT]{30})" + RIGHT)

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {}
i = 0
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            CODE[b1 + b2 + b3] = AA[i]; i += 1

def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

def tr(d):
    return "".join(CODE.get(d[j:j+3], "X") for j in range(0, len(d), 3))

counts = Counter()
for path in sys.argv[1:]:
    with gzip.open(path, "rt") as fh:
        for n, line in enumerate(fh):
            if n % 4 != 1:
                continue
            s = line.strip()
            m = PAT.search(s) or PAT.search(rc(s))
            if m:
                p = tr(m.group(1))
                if "*" not in p and "X" not in p:
                    counts[p] += 1

print(json.dumps({"total": sum(counts.values()),
                  "unique": len(counts),
                  "top": counts.most_common(12)}))
