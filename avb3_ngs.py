#!/usr/bin/env python3
#-----------------------------------------------------------------------------
# avb3_ngs.py
#-----------------------------------------------------------------------------
# Author:  Reza Mozafari
# Date:    2026-11-22
# Repo:    https://github.com/rmozafari/peptide-ngs-counterselection
#
# Summary: NGS analysis for mRNA display peptide selection: extraction, translation,
#          counting and enrichment across three selection rounds, for pooled
#          libraries of different randomized-region lengths.
#-----------------------------------------------------------------------------
"""
avb3_ngs.py - NGS analysis for mRNA-display peptide selection (integrin avb3)

Adapted from the architecture of rmozafari/illumina-peptide-ngs for the new
library design: variable 6 / 8 / 10 x NNK inserts between fixed flanks, plus a
7-aa positive control, across rounds R5 / R8 / R10.

Key differences from the original pipeline:
  * variable insert length (18 / 21 / 24 / 30 nt) instead of a fixed 30 nt
  * long constant flanks, no terminal cysteine
  * paired-end aware: R1 and R2 are read in lockstep and must agree
  * single streaming pass, no multi-GB intermediate text dumps
  * two comparisons (R8/R5 and R10/R8) instead of one

Dependencies: Python 3.8+ stdlib. openpyxl and matplotlib optional
(for the Excel workbook and the plots).
"""

import argparse
import csv
import gzip
import json
import math
import os
import sys
from collections import Counter, defaultdict
from itertools import product

# --------------------------------------------------------------------------
# Library definition
# --------------------------------------------------------------------------

LEFT_FLANK = "TAATACGACTCACTATAGGGACAATTACTATTTACAATTACAATG"   # ends with the ATG
RIGHT_FLANK = "GGTAGTGGTACGTCTGGCTCCAGC"                        # GSGTSGSS linker

# Accepted insert lengths in nt -> peptide length in aa
INSERT_LENGTHS = {18: 6, 21: 7, 24: 8, 30: 10}
RANDOM_LENGTHS = (6, 8, 10)          # the NNK libraries
CONTROL_LENGTH = 7                   # the positive control is the odd one out

# Positive control, taken from the construct in Library.docx:
#   ATG GGC AGG GGC GAT TGG GCC GGC  =  M  G  R  G  D  W  A  G
# The insert after the initiating ATG is 21 nt = GRGDWAG.
CONTROL_DNA = "GGCAGGGGCGATTGGGCCGGC"
CONTROL_PEPTIDE = "GRGDWAG"

ROUNDS = ("R5", "R8", "R10")

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

_codons = ["".join(c) for c in product("TCAG", repeat=3)]
_aas = ("FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG")
GENETIC_CODE = dict(zip(_codons, _aas))

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(dna):
    return "".join(GENETIC_CODE.get(dna[i:i + 3], "X") for i in range(0, len(dna), 3))


# --------------------------------------------------------------------------
# Read parsing
# --------------------------------------------------------------------------

def iter_fastq(path):
    """Yield (sequence, quality) from a plain or gzipped FASTQ file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="ascii", errors="replace") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\n").upper()
            plus = fh.readline()
            qual = fh.readline().rstrip("\n")
            if not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"Malformed FASTQ record in {path} near: {header.strip()!r}")
            yield seq, qual


class Extractor:
    """Pull the variable insert out of a single read, in either orientation."""

    def __init__(self, anchor_len=14, min_q=0, enforce_nnk=False):
        if anchor_len > len(RIGHT_FLANK):
            anchor_len = len(RIGHT_FLANK)
        self.left_anchor = LEFT_FLANK[-anchor_len:]
        self.right_anchor = RIGHT_FLANK[:anchor_len]
        self.min_q = min_q
        self.enforce_nnk = enforce_nnk

    def _try_one(self, seq, qual):
        """Return (insert_dna, status). status is a short QC string."""
        pos = seq.find(self.left_anchor)
        if pos == -1:
            return None, "no_left_anchor"
        start = pos + len(self.left_anchor)
        rpos = seq.find(self.right_anchor, start)
        if rpos == -1:
            return None, "no_right_anchor"
        length = rpos - start
        if length not in INSERT_LENGTHS:
            return None, "bad_insert_length"
        insert = seq[start:rpos]
        if "N" in insert:
            return None, "ambiguous_base"
        if self.min_q:
            window = qual[start:rpos]
            if min(ord(c) - 33 for c in window) < self.min_q:
                return None, "low_quality"
        if self.enforce_nnk and length != 21:
            # NNK: third base of every codon must be G or T
            if any(insert[i + 2] not in "GT" for i in range(0, length, 3)):
                return None, "not_nnk"
        return insert, "ok"

    def extract(self, seq, qual):
        """Try the read forward, then reverse-complemented."""
        insert, status = self._try_one(seq, qual)
        if insert is not None:
            return insert, "fwd", status
        rc_seq = reverse_complement(seq)
        rc_qual = qual[::-1]
        insert, status_rc = self._try_one(rc_seq, rc_qual)
        if insert is not None:
            return insert, "rev", status_rc
        # report the more informative of the two failures
        return None, None, status if status != "no_left_anchor" else status_rc


# --------------------------------------------------------------------------
# Per-round counting
# --------------------------------------------------------------------------

def count_round(r1_paths, r2_paths, extractor, strict_pairs=False,
                keep_dna=True, progress_every=1000000, label=""):
    """
    Stream one round. R1 and R2 are read in lockstep: Illumina writes mates in
    the same order, so record i of R1 pairs with record i of R2.

    A fragment is counted ONCE. If both mates yield an insert they must agree,
    which suppresses single-base sequencing errors. If only one mate yields an
    insert the fragment is still counted unless --strict-pairs is set.
    """
    peptide_counts = Counter()
    dna_counts = Counter()
    pep_to_dna = defaultdict(Counter) if keep_dna else None

    qc = Counter()
    length_hist = Counter()

    r1_iter = _chain_files(r1_paths)
    r2_iter = _chain_files(r2_paths) if r2_paths else None

    n = 0
    while True:
        rec1 = next(r1_iter, None)
        rec2 = next(r2_iter, None) if r2_iter else None
        if rec1 is None and rec2 is None:
            break
        if rec1 is None or (r2_iter and rec2 is None):
            qc["mate_file_length_mismatch"] += 1
            break

        n += 1
        if progress_every and n % progress_every == 0:
            print(f"    {label}: {n:,} fragments, {qc['counted']:,} usable",
                  file=sys.stderr, flush=True)

        ins1, orient1, st1 = extractor.extract(*rec1)
        if rec2 is not None:
            ins2, orient2, st2 = extractor.extract(*rec2)
        else:
            ins2, orient2, st2 = None, None, "no_mate"

        insert = None
        if ins1 is not None and ins2 is not None:
            if ins1 == ins2:
                insert = ins1
                qc["pair_concordant"] += 1
            else:
                qc["pair_discordant"] += 1
                qc["rejected"] += 1
                continue
        elif ins1 is not None or ins2 is not None:
            if strict_pairs and rec2 is not None:
                qc["rejected_single_mate"] += 1
                qc["rejected"] += 1
                continue
            insert = ins1 if ins1 is not None else ins2
            qc["single_mate"] += 1
        else:
            qc[f"fail_{st1}"] += 1
            qc["rejected"] += 1
            continue

        if orient1 == "rev" or (ins1 is None and orient2 == "rev"):
            qc["orientation_reverse"] += 1
        else:
            qc["orientation_forward"] += 1

        peptide = translate(insert)
        if "*" in peptide:
            qc["stop_codon"] += 1
            qc["rejected"] += 1
            continue
        if "X" in peptide:
            qc["untranslatable"] += 1
            qc["rejected"] += 1
            continue

        qc["counted"] += 1
        length_hist[len(peptide)] += 1
        peptide_counts[peptide] += 1
        if keep_dna:
            dna_counts[insert] += 1
            pep_to_dna[peptide][insert] += 1

    qc["total_fragments"] = n
    return {
        "peptides": peptide_counts,
        "dna": dna_counts,
        "pep_to_dna": pep_to_dna,
        "qc": qc,
        "length_hist": length_hist,
    }


def _chain_files(paths):
    for p in paths:
        for rec in iter_fastq(p):
            yield rec


# --------------------------------------------------------------------------
# Enrichment maths
# --------------------------------------------------------------------------

def fold_change(count_num, total_num, count_den, total_den):
    """
    Frequency ratio with an add-one pseudocount on both counts. The pseudocount
    keeps the ratio finite for peptides absent from one round and damps the
    wild ratios you otherwise get from single-read peptides.
    """
    f_num = (count_num + 1) / (total_num + 1)
    f_den = (count_den + 1) / (total_den + 1)
    return f_num / f_den


def build_table(round_data, min_total=1):
    """
    Merge the per-round counters into one row per peptide.

    Two normalizations are produced, because the three NNK libraries are pooled
    in one tube and are NOT equally represented. Measured in R5 the mix is
    roughly 95% 10-mers, 3% 6-mers, 2% 8-mers.

      freq_global  = count / all valid reads in the round
      freq_inlib   = count / reads of the SAME peptide length in that round

    If the pooling ratio drifts between rounds - entirely plausible, since
    longer peptides often win these selections - then freq_global moves for
    every peptide in a minority library even when nothing happened to that
    peptide. Worked example: an 8-mer holding a steady 5% of the 8-mer pool,
    while that pool shrinks from 1.8% to 0.7% of the tube, shows a 2.6-fold
    "depletion" in freq_global and no change at all in freq_inlib.

    freq_inlib answers "did this peptide beat its own competitors", which is
    the question the selection actually posed, so classification uses it.
    The positive control is a spike-in rather than a member of a competing
    library, so for it freq_inlib is defined as freq_global.
    """
    totals = {r: sum(round_data[r]["peptides"].values()) for r in ROUNDS}

    # per-round, per-length denominators
    len_totals = {r: Counter() for r in ROUNDS}
    for r in ROUNDS:
        for pep, c in round_data[r]["peptides"].items():
            len_totals[r][len(pep)] += c

    peptides = set()
    for r in ROUNDS:
        peptides.update(round_data[r]["peptides"])

    rows = []
    for pep in peptides:
        c = {r: round_data[r]["peptides"].get(pep, 0) for r in ROUNDS}
        if sum(c.values()) < min_total:
            continue

        dna_counter = Counter()
        for r in ROUNDS:
            p2d = round_data[r]["pep_to_dna"]
            if p2d and pep in p2d:
                dna_counter.update(p2d[pep])
        top_dna = dna_counter.most_common(1)[0][0] if dna_counter else ""

        # How many INDEPENDENT library members encoded this peptide?
        #
        # Each physical library member is one DNA molecule. A peptide recovered
        # from many distinct encodings was enriched as many independent events,
        # which PCR jackpotting cannot manufacture - jackpotting amplifies one
        # molecule, not a peptide sequence across different DNA.
        #
        # The raw variant count is NOT trustworthy alone: a synonymous
        # sequencing error creates a new "encoding" that is not a new founder,
        # and at high depth error alone fabricates dozens. So also report how
        # many encodings clear a read floor, and what share of reads sits in
        # the single commonest one. A peptide whose reads are concentrated in
        # one encoding is one clone however many variants are listed.
        dna_total = sum(dna_counter.values())
        n_ge5 = sum(1 for v in dna_counter.values() if v >= 5)
        n_ge20 = sum(1 for v in dna_counter.values() if v >= 20)
        top_share = (dna_counter.most_common(1)[0][1] / dna_total) if dna_total else 0.0

        f = {r: (c[r] / totals[r] * 100 if totals[r] else 0.0) for r in ROUNDS}
        fold_85 = fold_change(c["R8"], totals["R8"], c["R5"], totals["R5"])
        fold_108 = fold_change(c["R10"], totals["R10"], c["R8"], totals["R8"])

        length = len(pep)

        # in-library denominators; the spike-in control uses the global one
        is_control = (pep == CONTROL_PEPTIDE)
        lt = {r: (totals[r] if is_control else len_totals[r].get(length, 0))
              for r in ROUNDS}
        fi = {r: (c[r] / lt[r] * 100 if lt[r] else 0.0) for r in ROUNDS}
        fold_85_il = fold_change(c["R8"], lt["R8"], c["R5"], lt["R5"])
        fold_108_il = fold_change(c["R10"], lt["R10"], c["R8"], lt["R8"])

        rows.append({
            "freq_inlib_pct_R5": fi["R5"], "freq_inlib_pct_R8": fi["R8"],
            "freq_inlib_pct_R10": fi["R10"],
            "fold_inlib_R8_vs_R5": fold_85_il,
            "log2fc_inlib_R8_vs_R5": math.log2(fold_85_il),
            "fold_inlib_R10_vs_R8": fold_108_il,
            "log2fc_inlib_R10_vs_R8": math.log2(fold_108_il),
            "peptide": pep,
            "peptide_with_M": "M" + pep,
            "length_aa": length,
            "library": "positive_control" if pep == CONTROL_PEPTIDE
                       else ("control_length_7" if length == CONTROL_LENGTH else f"{length}AA"),
            "dna_sequence": top_dna,
            "n_dna_variants": len(dna_counter),
            "n_dna_variants_ge5": n_ge5,
            "n_dna_variants_ge20": n_ge20,
            "top_dna_read_share": top_share,
            "count_R5": c["R5"], "count_R8": c["R8"], "count_R10": c["R10"],
            "freq_pct_R5": f["R5"], "freq_pct_R8": f["R8"], "freq_pct_R10": f["R10"],
            "fold_R8_vs_R5": fold_85,
            "log2fc_R8_vs_R5": math.log2(fold_85),
            "fold_R10_vs_R8": fold_108,
            "log2fc_R10_vs_R8": math.log2(fold_108),
        })
    return rows, totals, len_totals


def classify(rows, min_count_r8=200, min_count_r5=5, up=2.0, down=0.5):
    """
    Label the two populations Yousef asked for.

    The read-count floor matters more than the fold thresholds: with one sample
    per round there is no error model, so a peptide seen a handful of times
    carries no usable fold change whatever ratio it produces.
    """
    for row in rows:
        if (row["count_R8"] < min_count_r8 or row["count_R5"] < min_count_r5
                or row["library"] == "positive_control"):
            row["class"] = ""
            continue
        enriched = row["fold_inlib_R8_vs_R5"] >= up
        if not enriched:
            row["class"] = ""
        elif row["fold_inlib_R10_vs_R8"] <= down:
            row["class"] = "enriched_then_depleted"      # lost on avb5 counterselection
        else:
            row["class"] = "enriched_maintained"         # avb3-specific candidate
    return rows


# --------------------------------------------------------------------------
# Conserved positions / motifs
# --------------------------------------------------------------------------

def conserved_positions(rows, length, top_n=500, key="count_R8"):
    """Per-position amino-acid frequency and Shannon information for the top peptides."""
    subset = [r for r in rows if r["length_aa"] == length and r["library"] == f"{length}AA"]
    subset.sort(key=lambda r: r[key], reverse=True)
    subset = subset[:top_n]
    if not subset:
        return []

    out = []
    for pos in range(length):
        counter = Counter(r["peptide"][pos] for r in subset)
        total = sum(counter.values())
        entropy = -sum((v / total) * math.log2(v / total) for v in counter.values() if v)
        consensus, consensus_n = counter.most_common(1)[0]
        row = {
            "position": pos + 1,
            "consensus": consensus,
            "consensus_pct": consensus_n / total * 100,
            "bits": math.log2(20) - entropy,       # information content vs uniform
            "n_peptides": total,
        }
        for aa in AMINO_ACIDS:
            row[aa] = counter.get(aa, 0) / total * 100
        out.append(row)
    return out


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

MAIN_COLUMNS = [
    "peptide", "peptide_with_M", "length_aa", "library", "dna_sequence",
    "n_dna_variants", "n_dna_variants_ge5", "n_dna_variants_ge20",
    "top_dna_read_share",
    "count_R5", "count_R8", "count_R10",
    "freq_pct_R5", "freq_pct_R8", "freq_pct_R10",
    "freq_inlib_pct_R5", "freq_inlib_pct_R8", "freq_inlib_pct_R10",
    "fold_inlib_R8_vs_R5", "log2fc_inlib_R8_vs_R5",
    "fold_inlib_R10_vs_R8", "log2fc_inlib_R10_vs_R8",
    "fold_R8_vs_R5", "log2fc_R8_vs_R5",
    "fold_R10_vs_R8", "log2fc_R10_vs_R8",
    "class",
]


def write_csv(path, rows, columns):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def control_report(round_data, totals):
    """
    Track the positive control at both peptide and DNA level, separately from
    the randomized libraries because its length differs.
    """
    report = []
    for r in ROUNDS:
        pep_count = round_data[r]["peptides"].get(CONTROL_PEPTIDE, 0)
        dna_count = round_data[r]["dna"].get(CONTROL_DNA, 0)
        total = totals[r]
        all_7mer = sum(v for k, v in round_data[r]["peptides"].items() if len(k) == 7)
        report.append({
            "round": r,
            "peptide": CONTROL_PEPTIDE,
            "peptide_with_M": "M" + CONTROL_PEPTIDE,
            "expected_dna": CONTROL_DNA,
            "count_peptide_level": pep_count,
            "count_exact_dna": dna_count,
            "count_synonymous_dna": pep_count - dna_count,
            "freq_pct": pep_count / total * 100 if total else 0.0,
            "total_valid_reads": total,
            "all_7aa_reads": all_7mer,
        })
    return report


def write_excel(path, rows, control, qc_rows, motifs, totals, excel_cap):
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, PatternFill, Alignment
        from openpyxl.utils import get_column_letter
    except ImportError:
        print("  openpyxl not installed - skipping Excel workbook "
              "(CSV outputs are complete). Install with: pip3 install openpyxl")
        return False

    wb = Workbook()
    head_font = Font(bold=True, color="FFFFFF")
    head_fill = PatternFill("solid", fgColor="2F5597")

    def add_sheet(title, data, columns, note=None):
        ws = wb.create_sheet(title[:31])
        start = 1
        if note:
            ws.cell(row=1, column=1, value=note).font = Font(italic=True, color="808080")
            start = 3
        for j, col in enumerate(columns, 1):
            cell = ws.cell(row=start, column=j, value=col)
            cell.font = head_font
            cell.fill = head_fill
            cell.alignment = Alignment(horizontal="center")
        for i, row in enumerate(data, start + 1):
            for j, col in enumerate(columns, 1):
                ws.cell(row=i, column=j, value=row.get(col))
        ws.freeze_panes = ws.cell(row=start + 1, column=1)
        for j, col in enumerate(columns, 1):
            width = max(len(str(col)) + 2, 12)
            if col in ("peptide", "peptide_with_M", "dna_sequence", "expected_dna"):
                width = 34
            ws.column_dimensions[get_column_letter(j)].width = width
        return ws

    # Summary sheet
    ws = wb.active
    ws.title = "Summary"
    ws["A1"] = "mRNA display NGS - integrin avb3 selection"
    ws["A1"].font = Font(bold=True, size=14)
    summary = [
        ("Rounds", "R5 (reference), R8 (avb3 selection), R10 (after avb5 counterselection)"),
        ("Libraries", "6, 8, 10 x NNK randomized; 7-aa positive control analyzed separately"),
        ("Left flank", LEFT_FLANK),
        ("Right flank", RIGHT_FLANK),
        ("Positive control peptide", CONTROL_PEPTIDE + "  (M-prefixed: M" + CONTROL_PEPTIDE + ")"),
        ("Positive control DNA", CONTROL_DNA),
        ("", ""),
        ("Valid peptide reads R5", totals["R5"]),
        ("Valid peptide reads R8", totals["R8"]),
        ("Valid peptide reads R10", totals["R10"]),
        ("Unique peptides", len(rows)),
        ("", ""),
        ("Frequency", "percent of valid peptide reads within that round"),
        ("Fold change", "frequency ratio with add-one pseudocount on both counts"),
        ("Caveat", "single sample per round: no replicates, so no error model and no p-values"),
        ("Caveat", "peptides with low R8 counts carry no usable fold change"),
        ("Caveat", "M is an in silico annotation, not evidence of Met incorporation"),
    ]
    for i, (k, v) in enumerate(summary, 3):
        ws.cell(row=i, column=1, value=k).font = Font(bold=True)
        ws.cell(row=i, column=2, value=v)
    ws.column_dimensions["A"].width = 26
    ws.column_dimensions["B"].width = 95

    add_sheet("Positive_control", control, list(control[0].keys()) if control else [],
              note="Positive control tracked separately - its length (7 aa) differs from the randomized libraries.")

    ranked = sorted(rows, key=lambda r: r["count_R8"], reverse=True)
    truncated = len(ranked) > excel_cap
    add_sheet("All_peptides", ranked[:excel_cap], MAIN_COLUMNS,
              note=(f"Top {excel_cap:,} of {len(ranked):,} peptides by R8 count. "
                    "Full table in all_peptides.csv." if truncated else None))

    for length in RANDOM_LENGTHS:
        subset = [r for r in ranked if r["library"] == f"{length}AA"]
        add_sheet(f"{length}AA", subset[:excel_cap], MAIN_COLUMNS)

    for name, label in (("enriched_then_depleted", "Enriched_then_depleted"),
                        ("enriched_maintained", "Enriched_maintained")):
        subset = [r for r in ranked if r.get("class") == name]
        subset.sort(key=lambda r: r["fold_inlib_R8_vs_R5"], reverse=True)
        note = ("Up from R5 to R8, then lost after avb5 counterselection - consistent with cross-reactivity."
                if name == "enriched_then_depleted"
                else "Up from R5 to R8 and retained in R10 - candidate avb3-specific binders.")
        add_sheet(label, subset[:excel_cap], MAIN_COLUMNS, note=note)

    for length, data in motifs.items():
        if data:
            add_sheet(f"Conserved_{length}AA", data, list(data[0].keys()),
                      note=f"Per-position amino acid composition (%) across the top enriched {length}-aa peptides. "
                           f"'bits' is information content; higher means more conserved.")

    add_sheet("QC", qc_rows, list(qc_rows[0].keys()) if qc_rows else [])

    wb.save(path)
    return True


def make_plots(outdir, rows, round_data, totals, control):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not installed - skipping plots. Install with: pip3 install matplotlib")
        return []

    written = []

    # 1. Read accounting and insert-length distribution
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].bar(ROUNDS, [totals[r] for r in ROUNDS], color="#2F5597")
    axes[0].set_ylabel("valid peptide reads")
    axes[0].set_title("Usable reads per round")
    for i, r in enumerate(ROUNDS):
        axes[0].text(i, totals[r], f"{totals[r]:,}", ha="center", va="bottom", fontsize=8)

    width = 0.25
    for i, r in enumerate(ROUNDS):
        hist = round_data[r]["length_hist"]
        tot = sum(hist.values()) or 1
        lengths = sorted(INSERT_LENGTHS.values())
        xs = [x + (i - 1) * width for x in range(len(lengths))]
        axes[1].bar(xs, [hist.get(L, 0) / tot * 100 for L in lengths], width, label=r)
    axes[1].set_xticks(range(len(sorted(INSERT_LENGTHS.values()))))
    axes[1].set_xticklabels([f"{L} aa" for L in sorted(INSERT_LENGTHS.values())])
    axes[1].set_ylabel("% of valid reads")
    axes[1].set_title("Library composition by insert length")
    axes[1].legend()
    fig.tight_layout()
    p = os.path.join(outdir, "qc_reads_and_lengths.png")
    fig.savefig(p, dpi=150); plt.close(fig); written.append(p)

    # 2. The money plot: R8/R5 vs R10/R8
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    plotted = [r for r in rows if r["count_R8"] >= 5 and r["library"].endswith("AA")]
    groups = {
        "other": ("#C8C8C8", 6, 0.35),
        "enriched_then_depleted": ("#C00000", 22, 0.85),
        "enriched_maintained": ("#2E7D32", 22, 0.85),
    }
    for name, (color, size, alpha) in groups.items():
        sel = [r for r in plotted if (r.get("class") or "other") == name]
        if not sel:
            continue
        ax.scatter([r["log2fc_inlib_R8_vs_R5"] for r in sel],
                   [r["log2fc_inlib_R10_vs_R8"] for r in sel],
                   s=size, c=color, alpha=alpha, linewidths=0,
                   label=f"{name} (n={len(sel)})")
    ctrl = [r for r in rows if r["library"] == "positive_control"]
    if ctrl:
        ax.scatter([ctrl[0]["log2fc_inlib_R8_vs_R5"]], [ctrl[0]["log2fc_inlib_R10_vs_R8"]],
                   s=140, marker="*", c="#FF8C00", edgecolors="black",
                   linewidths=0.6, label="positive control", zorder=5)
    ax.axhline(0, color="black", lw=0.8); ax.axvline(0, color="black", lw=0.8)
    ax.axvline(1, color="gray", lw=0.6, ls="--"); ax.axhline(-1, color="gray", lw=0.6, ls="--")
    ax.set_xlabel("log2 fold change  R8 / R5   (avb3 selection, within-library)")
    ax.set_ylabel("log2 fold change  R10 / R8   (after avb5 counterselection, within-library)")
    ax.set_title("Enrichment then counterselection\n(peptides with R8 count >= 5)")
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    p = os.path.join(outdir, "enrichment_scatter.png")
    fig.savefig(p, dpi=150); plt.close(fig); written.append(p)

    # 3. Frequency trajectories of the strongest candidates
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, cls, title in (
            (axes[0], "enriched_maintained", "Retained in R10 (avb3-specific candidates)"),
            (axes[1], "enriched_then_depleted", "Lost in R10 (likely cross-reactive)")):
        sel = sorted([r for r in rows if r.get("class") == cls],
                     key=lambda r: r["freq_pct_R8"], reverse=True)[:15]
        for r in sel:
            ax.plot(ROUNDS, [r["freq_pct_R5"], r["freq_pct_R8"], r["freq_pct_R10"]],
                    marker="o", lw=1.2, ms=4, alpha=0.8, label=r["peptide_with_M"])
        if ctrl:
            c = ctrl[0]
            ax.plot(ROUNDS, [c["freq_pct_R5"], c["freq_pct_R8"], c["freq_pct_R10"]],
                    marker="*", ms=13, lw=2, color="#FF8C00", label="control " + c["peptide_with_M"])
        ax.set_yscale("log")
        ax.set_ylabel("frequency (% of valid reads, log scale)")
        ax.set_title(title)
        if sel or ctrl:
            ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    p = os.path.join(outdir, "frequency_trajectories.png")
    fig.savefig(p, dpi=150); plt.close(fig); written.append(p)

    # 4. Rank-abundance: how converged is each round?
    fig, ax = plt.subplots(figsize=(7, 5))
    for r in ROUNDS:
        counts = sorted(round_data[r]["peptides"].values(), reverse=True)[:10000]
        tot = totals[r] or 1
        ax.plot(range(1, len(counts) + 1), [c / tot * 100 for c in counts], lw=1.4, label=r)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("peptide rank"); ax.set_ylabel("frequency (%)")
    ax.set_title("Rank-abundance: convergence across rounds")
    ax.legend()
    fig.tight_layout()
    p = os.path.join(outdir, "rank_abundance.png")
    fig.savefig(p, dpi=150); plt.close(fig); written.append(p)

    # 5. Positive control
    if control:
        fig, ax = plt.subplots(figsize=(6, 4.5))
        vals = [c["freq_pct"] for c in control]
        ax.bar(ROUNDS, vals, color="#FF8C00")
        for i, c in enumerate(control):
            ax.text(i, c["freq_pct"], f"{c['count_peptide_level']:,} reads",
                    ha="center", va="bottom", fontsize=8)
        ax.set_ylabel("frequency (%)")
        ax.set_title(f"Positive control M{CONTROL_PEPTIDE} across rounds")
        fig.tight_layout()
        p = os.path.join(outdir, "positive_control.png")
        fig.savefig(p, dpi=150); plt.close(fig); written.append(p)

    return written


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="NGS analysis for mRNA-display selection against integrin avb3 (R5/R8/R10).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    for r in ROUNDS:
        ap.add_argument(f"--{r.lower()}", nargs="+", required=True,
                        metavar="FASTQ", help=f"{r} read-1 FASTQ(.gz) file(s)")
        ap.add_argument(f"--{r.lower()}-r2", nargs="*", default=[],
                        metavar="FASTQ", help=f"{r} read-2 FASTQ(.gz) file(s) (optional)")
    ap.add_argument("--output", "-o", default="ngs_results", help="output directory")
    ap.add_argument("--anchor-len", type=int, default=14,
                    help="nt of constant flank used as the anchor on each side")
    ap.add_argument("--min-q", type=int, default=0,
                    help="minimum Phred score across the insert (0 = off)")
    ap.add_argument("--enforce-nnk", action="store_true",
                    help="require G or T at the third base of every randomized codon")
    ap.add_argument("--strict-pairs", action="store_true",
                    help="require both mates to yield the same insert")
    ap.add_argument("--min-count-r8", type=int, default=200,
                    help="minimum R8 read count before a peptide can be called a candidate")
    ap.add_argument("--min-count-r5", type=int, default=5,
                    help="minimum R5 read count; below this the R8/R5 ratio is not a "
                         "measurement, because a peptide seen once or twice is a sampling "
                         "coin-flip rather than an observed frequency")
    ap.add_argument("--up", type=float, default=2.0, help="fold-change threshold for enrichment R8/R5")
    ap.add_argument("--down", type=float, default=0.5, help="fold-change threshold for depletion R10/R8")
    ap.add_argument("--min-total", type=int, default=2,
                    help="drop peptides with fewer than this many reads summed over all rounds")
    ap.add_argument("--motif-top-n", type=int, default=500,
                    help="peptides per length used for the conserved-position analysis")
    ap.add_argument("--excel-cap", type=int, default=50000,
                    help="maximum rows written to any Excel sheet")
    ap.add_argument("--skip-dna", action="store_true",
                    help="do not keep DNA-level counts (lower memory, loses synonymous-codon info)")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--progress-every", type=int, default=1000000)
    args = ap.parse_args()

    os.makedirs(args.output, exist_ok=True)
    extractor = Extractor(anchor_len=args.anchor_len, min_q=args.min_q,
                          enforce_nnk=args.enforce_nnk)

    print("Library configuration")
    print(f"  left anchor   {extractor.left_anchor}  (last {args.anchor_len} nt of the constant 5' flank)")
    print(f"  right anchor  {extractor.right_anchor}  (first {args.anchor_len} nt of the constant 3' flank)")
    print(f"  insert sizes  {', '.join(f'{k} nt = {v} aa' for k, v in sorted(INSERT_LENGTHS.items()))}")
    print(f"  control       {CONTROL_PEPTIDE} / M{CONTROL_PEPTIDE}   DNA {CONTROL_DNA}")
    print()

    round_data = {}
    for r in ROUNDS:
        r1 = getattr(args, r.lower())
        r2 = getattr(args, r.lower() + "_r2")
        print(f"Counting {r}")
        for p in r1:
            print(f"  R1  {p}")
        for p in r2:
            print(f"  R2  {p}")
        round_data[r] = count_round(r1, r2, extractor,
                                    strict_pairs=args.strict_pairs,
                                    keep_dna=not args.skip_dna,
                                    progress_every=args.progress_every,
                                    label=r)
        qc = round_data[r]["qc"]
        total = qc["total_fragments"] or 1
        print(f"  {qc['total_fragments']:,} fragments -> {qc['counted']:,} usable "
              f"({qc['counted'] / total * 100:.1f}%)")
        if qc.get("pair_discordant"):
            print(f"  {qc['pair_discordant']:,} discarded: mates disagreed")
        print()

    print("Building merged table")
    rows, totals, len_totals = build_table(round_data, min_total=args.min_total)
    rows = classify(rows, min_count_r8=args.min_count_r8,
                    min_count_r5=args.min_count_r5, up=args.up, down=args.down)
    rows.sort(key=lambda x: x["count_R8"], reverse=True)
    print(f"  {len(rows):,} unique peptides")

    control = control_report(round_data, totals)

    print("\n  Library composition (% of valid reads):")
    print("    length " + "".join(f"{r:>12}" for r in ROUNDS))
    for L in sorted(INSERT_LENGTHS.values()):
        cells = "".join(f"{len_totals[r].get(L,0)/totals[r]*100:>11.2f}%" if totals[r]
                        else f"{'-':>12}" for r in ROUNDS)
        print(f"    {L:>3} aa {cells}")
    print("    (in-library frequencies correct for drift in these proportions)\n")

    qc_rows = []
    for r in ROUNDS:
        row = {"round": r}
        row.update(round_data[r]["qc"])
        for L in sorted(INSERT_LENGTHS.values()):
            row[f"reads_{L}aa"] = round_data[r]["length_hist"].get(L, 0)
        row["unique_peptides"] = len(round_data[r]["peptides"])
        row["unique_dna_variants"] = len(round_data[r]["dna"])
        qc_rows.append(row)
    all_qc_keys = []
    for row in qc_rows:
        for k in row:
            if k not in all_qc_keys:
                all_qc_keys.append(k)
    qc_rows = [{k: row.get(k, 0) for k in all_qc_keys} for row in qc_rows]

    motifs = {L: conserved_positions(rows, L, top_n=args.motif_top_n) for L in RANDOM_LENGTHS}

    print("Writing output")
    write_csv(os.path.join(args.output, "all_peptides.csv"), rows, MAIN_COLUMNS)
    for L in RANDOM_LENGTHS:
        subset = [r for r in rows if r["library"] == f"{L}AA"]
        write_csv(os.path.join(args.output, f"peptides_{L}aa.csv"), subset, MAIN_COLUMNS)
    write_csv(os.path.join(args.output, "positive_control.csv"), control,
              list(control[0].keys()) if control else [])
    for name in ("enriched_then_depleted", "enriched_maintained"):
        subset = sorted([r for r in rows if r.get("class") == name],
                        key=lambda r: r["fold_inlib_R8_vs_R5"], reverse=True)
        write_csv(os.path.join(args.output, f"{name}.csv"), subset, MAIN_COLUMNS)
        print(f"  {name}: {len(subset):,} peptides")
    for L, data in motifs.items():
        if data:
            write_csv(os.path.join(args.output, f"conserved_positions_{L}aa.csv"),
                      data, list(data[0].keys()))

    with open(os.path.join(args.output, "QC_summary.json"), "w") as fh:
        json.dump({
            "rounds": {r: dict(round_data[r]["qc"]) for r in ROUNDS},
            "valid_peptide_reads": totals,
            "unique_peptides": len(rows),
            "length_distribution": {r: dict(round_data[r]["length_hist"]) for r in ROUNDS},
            "settings": vars(args),
            "library": {
                "left_flank": LEFT_FLANK, "right_flank": RIGHT_FLANK,
                "left_anchor": extractor.left_anchor, "right_anchor": extractor.right_anchor,
                "control_peptide": CONTROL_PEPTIDE, "control_dna": CONTROL_DNA,
            },
        }, fh, indent=2)

    xlsx = os.path.join(args.output, "avb3_NGS_results.xlsx")
    if write_excel(xlsx, rows, control, qc_rows, motifs, totals, args.excel_cap):
        print(f"  {xlsx}")

    if not args.no_plots:
        for p in make_plots(args.output, rows, round_data, totals, control):
            print(f"  {p}")

    print("\nPositive control")
    for c in control:
        print(f"  {c['round']}: {c['count_peptide_level']:,} reads "
              f"({c['freq_pct']:.4f}%), exact DNA {c['count_exact_dna']:,}")
    if all(c["count_peptide_level"] == 0 for c in control):
        print("  NOT DETECTED - check that the control construct carries the same flanks.")

    print(f"\nDone. Results in {args.output}/")


if __name__ == "__main__":
    main()
