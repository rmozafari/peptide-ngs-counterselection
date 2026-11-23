#!/usr/bin/env python3
#-----------------------------------------------------------------------------
# validate.py
#-----------------------------------------------------------------------------
# Author:  Reza Mozafari
# Date:    2026-11-22
# Repo:    https://github.com/rmozafari/peptide-ngs-counterselection
#
# Summary: Score analysis output against the simulator's ground truth.
#-----------------------------------------------------------------------------
"""validate.py - score the pipeline output against the simulator's ground truth."""
import csv, json, os, sys

truth_path = sys.argv[1] if len(sys.argv) > 1 else "sim_fastq/ground_truth.json"
results_dir = sys.argv[2] if len(sys.argv) > 2 else "sim_results"

truth = json.load(open(truth_path))
rows = {r["peptide"]: r for r in csv.DictReader(open(os.path.join(results_dir, "all_peptides.csv")))}
qc = json.load(open(os.path.join(results_dir, "QC_summary.json")))

EXPECT = {
    "positive_control":    "none",
    "cross_reactive_10mer": "enriched_then_depleted",
    "cross_reactive_8mer":  "enriched_then_depleted",
    "cross_reactive_6mer":  "enriched_then_depleted",
    "specific_10mer_A":     "enriched_maintained",
    "specific_10mer_B":     "enriched_maintained",
    "specific_8mer":        "enriched_maintained",
    "specific_6mer":        "enriched_maintained",
    "neutral_10mer":        "none",
    "neutral_8mer":         "none",
}

print(f"{'planted':<22} {'len':>3} {'found':>5} {'R5':>7} {'R8':>7} {'R10':>7} "
      f"{'fc85':>7} {'fc108':>7}  {'class':<24} {'ok'}")
print("-" * 116)

failures = []
for p in truth["planted"]:
    pep = p["peptide"]
    row = rows.get(pep)
    expect = EXPECT[p["name"]]
    if row is None:
        print(f"{p['name']:<22} {p['length_aa']:>3} {'NO':>5}")
        failures.append(f"{p['name']}: not found in output")
        continue

    obs = {r: int(row[f"count_{r}"]) for r in ("R5", "R8", "R10")}
    exp_frag = {r: truth["observed_planted_fragments"][r][p["name"]] for r in ("R5", "R8", "R10")}
    got = row["class"] or "none"
    ok = (got == expect)

    # recovery rate: counted reads vs fragments actually written
    rec = [obs[r] / exp_frag[r] if exp_frag[r] else 1.0 for r in ("R5", "R8", "R10")]

    print(f"{p['name']:<22} {p['length_aa']:>3} {'yes':>5} "
          f"{obs['R5']:>7,} {obs['R8']:>7,} {obs['R10']:>7,} "
          f"{float(row['fold_R8_vs_R5']):>7.2f} {float(row['fold_R10_vs_R8']):>7.2f}  "
          f"{got:<24} {'PASS' if ok else 'FAIL'}")

    if not ok:
        failures.append(f"{p['name']}: expected {expect}, got {got}")
    # Reads carrying a sequencing error inside the insert are correctly
    # discarded (mates disagree), so recovery is ~95%, not 100%. Tolerance has
    # to widen for small counts or binomial noise alone trips it.
    import math as _m
    for rname, frac, n in zip(("R5", "R8", "R10"), rec, [exp_frag[r] for r in ("R5","R8","R10")]):
        if n < 5:
            continue
        sd = _m.sqrt(max(n, 1) * 0.05 * 0.95) / max(n, 1)
        if frac < 0.95 - 4 * sd - 0.02:
            failures.append(f"{p['name']} {rname}: read recovery {frac:.1%} "
                            f"below 4 SD of expected 95% (n={n})")
    if row["peptide_with_M"] != "M" + pep:
        failures.append(f"{p['name']}: M-prefix column wrong")
    if p["name"] == "positive_control":
        if row["library"] != "positive_control":
            failures.append("control not flagged as positive_control")
        if row["dna_sequence"] != p["insert_dna"]:
            failures.append(f"control DNA mismatch: {row['dna_sequence']} != {p['insert_dna']}")

print()

# translation spot check
import importlib.util
spec = importlib.util.spec_from_file_location("avb3", os.path.join(os.path.dirname(os.path.abspath(__file__)), "avb3_ngs.py"))
avb3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(avb3)
for p in truth["planted"]:
    if avb3.translate(p["insert_dna"]) != p["peptide"]:
        failures.append(f"translation mismatch for {p['name']}")
print(f"Translation check: {len(truth['planted'])} planted inserts re-translated"
      f" -> {'all match' if not any('translation' in f for f in failures) else 'MISMATCH'}")

# control DNA -> peptide, independently
assert avb3.translate("ATG" + avb3.CONTROL_DNA) == "M" + avb3.CONTROL_PEPTIDE
print(f"Control construct ATG+{avb3.CONTROL_DNA} translates to "
      f"M{avb3.CONTROL_PEPTIDE}  (matches Library.docx codon annotation)")

# no stop codons leaked through
stops = [p for p in rows if "*" in p]
if stops:
    failures.append(f"{len(stops)} peptides with stop codons leaked into output")
print(f"Stop-codon leakage: {len(stops)} peptides")

# length sanity
bad_len = [p for p, r in rows.items() if len(p) != int(r["length_aa"])]
if bad_len:
    failures.append(f"{len(bad_len)} rows where length_aa disagrees with the peptide")
lens = sorted({int(r['length_aa']) for r in rows.values()})
print(f"Peptide lengths present: {lens}")
if lens != [6, 7, 8, 10]:
    failures.append(f"unexpected length set {lens}")

# frequencies sum to 100 per round
min_total = qc["settings"]["min_total"]
for r in ("R5", "R8", "R10"):
    s = sum(float(row[f"freq_pct_{r}"]) for row in rows.values())
    total = qc["valid_peptide_reads"][r]
    counted = sum(int(row[f"count_{r}"]) for row in rows.values())
    print(f"Frequency sum {r}: {s:.4f}%  "
          f"({counted:,}/{total:,} reads in retained rows; remainder is peptides "
          f"below --min-total {min_total}, which stay in the denominator)")
    # the sum must equal the fraction of reads in retained rows, exactly
    if abs(s - counted / total * 100) > 1e-6:
        failures.append(f"{r}: frequency column inconsistent with counts/total")
    if min_total <= 1 and abs(s - 100) > 1e-6:
        failures.append(f"{r}: with no min-total filter frequencies must sum to 100, got {s:.4f}")

# orientation handling
for r in ("R5", "R8", "R10"):
    q = qc["rounds"][r]
    fwd, rev = q.get("orientation_forward", 0), q.get("orientation_reverse", 0)
    frac = rev / (fwd + rev) if (fwd + rev) else 0
    print(f"Orientation {r}: {frac:.1%} reverse-complement (simulated 50%)")
    if not 0.4 < frac < 0.6:
        failures.append(f"{r}: reverse-orientation fraction {frac:.1%} off target")

print()
if failures:
    print(f"{len(failures)} FAILURE(S):")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print("ALL CHECKS PASSED")
