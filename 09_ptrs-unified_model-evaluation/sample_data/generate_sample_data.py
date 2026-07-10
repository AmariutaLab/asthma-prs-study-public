#!/usr/bin/env python3
"""
generate_sample_data.py — synthetic fixture for the 09 evaluation demo.

Produces the *structure* the module's notebooks consume (per-feature PTRS CSVs
with `Sample_ID, Keep_Vector, asthma` + a PRS table), but with FULLY SYNTHETIC
values — no individual-level study data. A weak signal is planted so the
cross-modal integrated model beats PRS-alone, illustrating the headline result.

Layout written under sample_data/:
    ptrs_results-concat_17CT_gacrs_train_test/cd4_naive_results.csv     (train: GACRS)
    ptrs_results-concat_39_gacrs_train_test/Esophagus_Mucosa_results.csv
    ptrs_results-concat_17CT_camp/cd4_naive_results.csv                 (test: CAMP, balanced)
    ptrs_results-concat_39_camp/Esophagus_Mucosa_results.csv
    prs_demo.csv                                                        (Sample_ID, PRS)

Deterministic (seed=42) so run_demo.py reproduces the committed expected_outputs.
"""
import os
import numpy as np
import pandas as pd

RNG = np.random.RandomState(42)
HERE = os.path.dirname(os.path.abspath(__file__))

# effect sizes: cd4_naive (CT) + Esophagus_Mucosa (tissue) are orthogonal signals;
# PRS carries weaker signal on its own — so integration > PRS alone.
B_CD4, B_ESO, B_PRS = 1.1, 1.0, 0.6

def make_cohort(prefix, n_case, n_ctrl):
    n = n_case + n_ctrl
    cd4 = RNG.normal(0, 1, n)
    eso = RNG.normal(0, 1, n)
    prs = RNG.normal(0, 1, n)
    logit = B_CD4*cd4 + B_ESO*eso + B_PRS*prs + RNG.normal(0, 1.4, n)
    if n_case is not None:            # force a balanced case/control count
        order = np.argsort(-logit)
        y = np.zeros(n, int); y[order[:n_case]] = 1
    ids = [f"{prefix}_{i:04d}" for i in range(n)]
    return pd.DataFrame({"Sample_ID": ids, "cd4": cd4, "eso": eso, "prs": prs, "asthma": y})

def write_feature(df, dirname, feature, col):
    d = os.path.join(HERE, dirname); os.makedirs(d, exist_ok=True)
    out = pd.DataFrame({"Sample_ID": df.Sample_ID, "Keep_Vector": df[col], "asthma": df.asthma})
    out.to_csv(os.path.join(d, f"{feature}_results.csv"), index=False)

# --- GACRS train+test cohort (imbalanced, like the real data) ---
gacrs = make_cohort("GACRS", n_case=140, n_ctrl=260)
write_feature(gacrs, "ptrs_results-concat_17CT_gacrs_train_test", "cd4_naive", "cd4")
write_feature(gacrs, "ptrs_results-concat_39_gacrs_train_test",  "Esophagus_Mucosa", "eso")

# --- CAMP-only test cohort (balanced 64 vs 64) ---
camp = make_cohort("CAMP", n_case=64, n_ctrl=64)
write_feature(camp, "ptrs_results-concat_17CT_camp", "cd4_naive", "cd4")
write_feature(camp, "ptrs_results-concat_39_camp",  "Esophagus_Mucosa", "eso")

# --- PRS table (both cohorts) ---
prs = pd.concat([gacrs[["Sample_ID", "prs"]], camp[["Sample_ID", "prs"]]])
prs.columns = ["Sample_ID", "PRS"]
prs.to_csv(os.path.join(HERE, "prs_demo.csv"), index=False)

print(f"GACRS train+test: {len(gacrs)} ({gacrs.asthma.sum()} cases) | "
      f"CAMP test: {len(camp)} ({camp.asthma.sum()} cases)")
print("wrote synthetic per-feature PTRS CSVs + prs_demo.csv under sample_data/")
