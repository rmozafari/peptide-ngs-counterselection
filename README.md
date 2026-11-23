# peptide-ngs-counterselection

Multi-round enrichment and counterselection analysis for mRNA display peptide NGS:
variable-length NNK libraries, paired-end error suppression, and separate tracking
of a positive control whose length differs from the randomized libraries.

Companion to [`illumina-peptide-ngs`](https://github.com/rmozafari/illumina-peptide-ngs),
which handles fixed-length libraries and a single enrichment comparison. This repo
shares no code with it; the design is a descendant, not a fork.

## What it does

Takes three rounds of Illumina FASTQ from a pooled mRNA display selection and
returns ranked candidate peptides with counts, frequencies and fold changes.

| round | role |
|---|---|
| first | reference round for enrichment |
| second | after selection on the target |
| third | after counterselection against a related target |

Two populations are separated: peptides that rise between the first two rounds
and then fall in the third (lost on counterselection), and peptides that rise
and hold (candidate target-specific binders).

Author: Reza Mozafari. Licensed MIT.


Adapted from the architecture of `rmozafari/illumina-peptide-ngs` for the new
library design. Validated against simulated data with known ground truth before
being pointed at real reads.

## Files

| File | Purpose |
|---|---|
| `avb3_ngs.py` | the analysis — FASTQ in, Excel/CSV/plots out |
| `simulate_fastq.py` | generates synthetic R5/R8/R10 with planted peptides and known behaviour |
| `validate.py` | scores analysis output against the simulator's ground truth |
| `crosscheck.py` | deliberately naive independent re-implementation, for count comparison |

## Requirements

Python 3.8+. The analysis runs on the standard library alone; `openpyxl` adds the
Excel workbook and `matplotlib` adds the plots.

```bash
pip3 install openpyxl matplotlib
```

## Running it on the real data

Put the six Dropbox files in one folder. Note the naming trap: in
`..._R10_S59_L004_R1_001.fastq.gz` the **first** `R10` is the selection round and
`_R1_001` is the sequencing read. Round and read are both "R".

```bash
cd ~/Desktop/avb3_ngs

python3 avb3_ngs.py \
  --r5  2533TCLT4_4_2160378671_R5_S57_L004_R1_001.fastq.gz \
  --r5-r2  2533TCLT4_4_2160378671_R5_S57_L004_R2_001.fastq.gz \
  --r8  2533TCLT4_4_2160378671_R8_S58_L004_R1_001.fastq.gz \
  --r8-r2  2533TCLT4_4_2160378671_R8_S58_L004_R2_001.fastq.gz \
  --r10 2533TCLT4_4_2160378671_R10_S59_L004_R1_001.fastq.gz \
  --r10-r2 2533TCLT4_4_2160378671_R10_S59_L004_R2_001.fastq.gz \
  --output results
```

Expect roughly 30–60 minutes for ~4.2 GB of gzipped input on a laptop, single core,
and 3–6 GB of RAM depending on library diversity. If memory is tight add
`--skip-dna`, which drops the DNA-variant table and roughly halves peak usage.

### Useful options

| Flag | Effect |
|---|---|
| `--anchor-len N` | nt of constant flank matched on each side (default 14) |
| `--strict-pairs` | discard a fragment unless both mates give the same insert |
| `--min-q N` | minimum Phred across the insert |
| `--enforce-nnk` | require G/T at the third base of every randomized codon |
| `--min-count-r8 N` | read floor before a peptide can be called a candidate (default 10) |
| `--up` / `--down` | fold-change thresholds (defaults 2.0 and 0.5) |
| `--skip-dna` | lower memory, loses synonymous-codon information |

Run once with defaults, then once with `--min-q 20 --enforce-nnk --strict-pairs`,
and compare. Numbers that move a lot under stricter settings were not solid.

## Design decisions worth knowing

**Short anchors, not full flanks.** The script matches the last 14 nt of the 5'
flank and the first 14 nt of the 3' flank rather than all 69 nt of constant
sequence. Every constant base you require to be error-free is a base that can
throw the read away. At a 0.1% per-base error rate, demanding all 69 nt costs
about 6.7% of reads; 28 nt costs about 2.8%. Measured on simulated data at four
anchor lengths, observed yield tracked `0.999^n` to within 0.1% — so the loss is
entirely explained by flank errors, and the shorter anchor is strictly better
without any loss of specificity (a 14-mer still occurs by chance once in 2.7×10^8 nt).

**Paired-end used for error suppression, not for depth.** The amplicon is 99 nt,
so a 150-cycle R1 already spans the whole insert — merging is unnecessary. Instead
R1 and R2 are read in lockstep and the fragment is counted once, requiring the two
mates to agree. A sequencing error in the insert makes them disagree and the
fragment is dropped. This costs a few percent of reads and removes a corresponding
number of spurious one-off peptides, which matters because the long tail of
single-read peptides is exactly where false enrichment lives.

**Insert length is discovered, not assumed.** The script finds the 3' anchor and
measures the gap, accepting 18/21/24/30 nt. Off-size inserts are counted in QC
rather than silently dropped, so indels and synthesis errors are visible.

**Pseudocounts on fold changes.** `fold = ((c_num+1)/(T_num+1)) / ((c_den+1)/(T_den+1))`.
Keeps ratios finite for peptides absent from a round and damps the implausible
ratios you get from single-read peptides.

**Frequencies are percent of valid peptide reads in that round.** Peptides dropped
by `--min-total` remain in the denominator, so frequencies in the table sum to
slightly under 100%. That is intentional — the denominator is the library, not the
subset you kept.

## Outputs

`avb3_NGS_results.xlsx` with sheets: Summary, Positive_control, All_peptides,
6AA, 8AA, 10AA, Enriched_then_depleted, Enriched_maintained, Conserved_6/8/10AA, QC.
The same content as CSV (full, uncapped), plus `QC_summary.json` and five PNGs.

Columns: peptide, peptide_with_M, length_aa, library, dna_sequence (most common
encoding), n_dna_variants, counts and frequencies for R5/R8/R10, fold and log2
fold for R8/R5 and R10/R8, class.

## Interpretation limits

These are not script bugs; they are properties of the experiment.

- **One sample per round.** No replicates means no error model, no p-values, and
  no principled threshold. Any fold-change cutoff, including the default, is a
  convention. Treat the output as a ranked list, not a hit list.
- **Depth vs diversity.** A 10×NNK library has ~10^15 theoretical members; a few
  million reads samples essentially none of it. Most peptides will have counts in
  the low single digits, where an "8-fold enrichment" is one read becoming eight.
  The `--min-count-r8` floor matters far more than the fold threshold.
- **Loss between R8 and R10 is consistent with avb5 binding, not proof of it.**
  PCR bias, round-to-round bottlenecks, and differences in sequencing depth all
  produce the same signature.
- **`peptide_with_M` is an annotation, not an observation.** The initiating ATG is
  in the constant flank and is not sequenced as part of the variable region.
  Under reprogrammed initiation the actual N-terminal residue may not be methionine.

## Validation performed

`simulate_fastq.py` writes paired FASTQ with ten planted peptides at defined
frequencies — three that rise then fall, four that rise and stay, two flat, plus
the positive control — over a nested background (R10 ⊂ R8 ⊂ R5), with 50% of
fragments reverse-complemented, 0.1% per-base error, adapter readthrough, stop
codons and 2% junk reads.

```bash
python3 simulate_fastq.py sim_fastq 200000
python3 avb3_ngs.py --r5 sim_fastq/SIM_R5_S1_L004_R1_001.fastq.gz --r5-r2 sim_fastq/SIM_R5_S1_L004_R2_001.fastq.gz \
                    --r8 sim_fastq/SIM_R8_S1_L004_R1_001.fastq.gz --r8-r2 sim_fastq/SIM_R8_S1_L004_R2_001.fastq.gz \
                    --r10 sim_fastq/SIM_R10_S1_L004_R1_001.fastq.gz --r10-r2 sim_fastq/SIM_R10_S1_L004_R2_001.fastq.gz \
                    --output sim_results --min-total 1
python3 validate.py sim_fastq/ground_truth.json sim_results
```

All ten planted peptides are recovered and classified correctly; translation,
orientation handling, stop-codon exclusion, length assignment, frequency
normalization and the M-prefix column all check out. `crosscheck.py` is an
independent re-implementation sharing no code; its counts agree once the
difference in anchor length is accounted for.

## Author

Reza Mozafari — https://github.com/rmozafari

Released under the MIT License. If you use this in published work, a citation
of the repository is appreciated.
