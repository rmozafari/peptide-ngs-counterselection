#!/usr/bin/env python3
#-----------------------------------------------------------------------------
# simulate_fastq.py
#-----------------------------------------------------------------------------
# Author:  Reza Mozafari
# Date:    2026-11-22
# Repo:    https://github.com/rmozafari/peptide-ngs-counterselection
#
# Summary: Generate synthetic paired FASTQ with known ground truth, for validating
#          the analysis before it is pointed at real reads.
#-----------------------------------------------------------------------------
"""
simulate_fastq.py - synthetic paired-end FASTQ with known ground truth.

Builds R5/R8/R10 the way the real experiment should look: a diverse NNK
background, a handful of peptides planted at defined frequencies with defined
behaviour across rounds, the positive control at a known level, plus the
nuisances that break naive parsers (reverse-complemented reads, sequencing
errors, adapter readthrough, stop codons, off-target junk).

Writes the planted truth to ground_truth.json so the analysis can be scored
against it rather than eyeballed.
"""
import gzip, json, os, random, sys
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from avb3_ngs import (LEFT_FLANK, RIGHT_FLANK, CONTROL_DNA, CONTROL_PEPTIDE,
                      GENETIC_CODE, reverse_complement, translate)

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTG"
READ_LEN = 150
ERROR_RATE = 0.001          # per-base substitution
RC_FRACTION = 0.5           # half the fragments sequenced in reverse orientation
JUNK_FRACTION = 0.02        # reads with no recognizable construct

# NNK codons: third base restricted to G or T
NNK_CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "GT"]
NNK_SENSE = [c for c in NNK_CODONS if GENETIC_CODE[c] != "*"]


def random_insert(n_codons, allow_stop=False):
    pool = NNK_CODONS if allow_stop else NNK_SENSE
    return "".join(random.choice(pool) for _ in range(n_codons))


def mutate(seq, rate):
    out = []
    for b in seq:
        if random.random() < rate:
            out.append(random.choice([x for x in "ACGT" if x != b]))
        else:
            out.append(b)
    return "".join(out)


def build_amplicon(insert):
    return LEFT_FLANK + insert + RIGHT_FLANK


def write_pair(r1_fh, r2_fh, idx, amplicon):
    frag = amplicon + ADAPTER
    if random.random() < RC_FRACTION:
        frag_r1 = reverse_complement(amplicon) + ADAPTER
        frag_r2 = amplicon + ADAPTER
    else:
        frag_r1 = frag
        frag_r2 = reverse_complement(amplicon) + ADAPTER

    for fh, seq in ((r1_fh, frag_r1), (r2_fh, frag_r2)):
        seq = (seq + "G" * READ_LEN)[:READ_LEN]
        seq = mutate(seq, ERROR_RATE)
        qual = "I" * READ_LEN
        fh.write(f"@SIM:{idx}\n{seq}\n+\n{qual}\n")


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "sim_fastq"
    reads_per_round = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    random.seed(20261121)
    os.makedirs(out_dir, exist_ok=True)

    # ---- planted peptides with defined behaviour -------------------------
    # frequencies are fractions of the round's fragments
    planted = []

    def plant(name, n_codons, f5, f8, f10, insert=None, peptide=None):
        ins = insert or random_insert(n_codons)
        planted.append({
            "name": name,
            "insert_dna": ins,
            "peptide": peptide or translate(ins),
            "length_aa": len(ins) // 3,
            "freq": {"R5": f5, "R8": f8, "R10": f10},
        })

    # positive control: present throughout, modestly enriched
    plant("positive_control", 7, 0.0020, 0.0060, 0.0055,
          insert=CONTROL_DNA, peptide=CONTROL_PEPTIDE)

    # up R5->R8, lost R8->R10  (cross-reactive: binds avb5 too)
    plant("cross_reactive_10mer", 10, 0.00020, 0.0120, 0.00090)
    plant("cross_reactive_8mer",  8,  0.00030, 0.0080, 0.00120)
    plant("cross_reactive_6mer",  6,  0.00015, 0.0045, 0.00060)

    # up R5->R8, retained in R10  (avb3-specific candidates)
    plant("specific_10mer_A", 10, 0.00020, 0.0150, 0.0220)
    plant("specific_10mer_B", 10, 0.00010, 0.0060, 0.0075)
    plant("specific_8mer",     8, 0.00025, 0.0090, 0.0110)
    plant("specific_6mer",     6, 0.00020, 0.0040, 0.0048)

    # flat: no selection, should fall in neither class
    plant("neutral_10mer", 10, 0.0010, 0.0011, 0.0009)
    plant("neutral_8mer",   8, 0.0010, 0.0010, 0.0012)

    # make sure every planted insert is unique
    seen = set()
    for p in planted:
        assert p["insert_dna"] not in seen, "duplicate planted insert"
        seen.add(p["insert_dna"])

    truth = {"reads_per_round": reads_per_round, "planted": planted,
             "error_rate": ERROR_RATE, "rc_fraction": RC_FRACTION,
             "junk_fraction": JUNK_FRACTION}

    # background library composition (fractions of the non-planted remainder)
    background_mix = {6: 0.30, 8: 0.35, 10: 0.35}
    # convergence: later rounds are less diverse, so draw background from a
    # shrinking pool to mimic a bottlenecked population
    pool_size = {"R5": 60000, "R8": 12000, "R10": 6000}

    # Pools must be NESTED: selection acts on the population that already
    # exists, so R8 is a subset of R5 and R10 a subset of R8. Drawing three
    # independent random pools would make every R8 peptide absent from R10
    # and label the entire background "enriched then depleted".
    pool_r5 = []
    for n_codons, share in background_mix.items():
        k = max(1, int(pool_size["R5"] * share))
        pool_r5 += [random_insert(n_codons) for _ in range(k)]
    random.shuffle(pool_r5)
    pool_r8 = pool_r5[:pool_size["R8"]]
    pool_r10 = pool_r8[:pool_size["R10"]]
    pools = {"R5": pool_r5, "R8": pool_r8, "R10": pool_r10}

    for rnd in ("R5", "R8", "R10"):
        r1_path = os.path.join(out_dir, f"SIM_{rnd}_S1_L004_R1_001.fastq.gz")
        r2_path = os.path.join(out_dir, f"SIM_{rnd}_S1_L004_R2_001.fastq.gz")
        counts = {p["name"]: 0 for p in planted}
        n_junk = 0
        n_stop = 0

        with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
            for i in range(reads_per_round):
                roll = random.random()

                # junk reads: no construct at all
                if roll < JUNK_FRACTION:
                    seq = "".join(random.choice("ACGT") for _ in range(READ_LEN))
                    q = "I" * READ_LEN
                    f1.write(f"@SIM:{i}\n{seq}\n+\n{q}\n")
                    f2.write(f"@SIM:{i}\n{seq}\n+\n{q}\n")
                    n_junk += 1
                    continue

                # planted peptides
                cum = JUNK_FRACTION
                chosen = None
                for p in planted:
                    cum += p["freq"][rnd]
                    if roll < cum:
                        chosen = p
                        break

                if chosen is not None:
                    write_pair(f1, f2, i, build_amplicon(chosen["insert_dna"]))
                    counts[chosen["name"]] += 1
                    continue

                # background, with occasional stop codons from real NNK
                if random.random() < 0.04:
                    n_codons = random.choices([6, 8, 10], weights=[3, 3.5, 3.5])[0]
                    insert = random_insert(n_codons, allow_stop=True)
                    if "*" in translate(insert):
                        n_stop += 1
                else:
                    insert = random.choice(pools[rnd])
                write_pair(f1, f2, i, build_amplicon(insert))

        truth.setdefault("observed_planted_fragments", {})[rnd] = counts
        truth.setdefault("junk_fragments", {})[rnd] = n_junk
        print(f"{rnd}: {reads_per_round:,} fragments -> {r1_path}, {r2_path}")

    with open(os.path.join(out_dir, "ground_truth.json"), "w") as fh:
        json.dump(truth, fh, indent=2)
    print(f"\nGround truth -> {os.path.join(out_dir, 'ground_truth.json')}")


if __name__ == "__main__":
    main()
