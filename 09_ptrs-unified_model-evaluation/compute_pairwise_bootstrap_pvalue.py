"""
Pairwise bootstrap AUC / delta-AUC / empirical one-sided P-value across all
Figure 4B/4C-cited models on CAMP-only Balanced.

For each model, reconstruct the 100 per-iteration CAMP-Balanced AUCs from the
model's persisted per-sample predictions using the same RNG the notebooks use:
  seed = iteration_index in 0..99
  np.random.RandomState(seed).choice(cases, 64, replace=False) + all_controls

Then compute, for every ordered pair (A, B):
  Delta_AUC_i           = AUC_A_i - AUC_B_i  (paired, same iteration)
  mean_Delta_AUC        = mean over 100 iterations
  empirical_one_sided_P = fraction of iterations with Delta_AUC <= 0

Writes a long-format CSV to OUTPUT_ROOT/pairwise_comparisons/pairwise_bootstrap_auc.csv.

Usage: python compute_pairwise_bootstrap_pvalue.py
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import product
from sklearn.metrics import roc_auc_score

# ============================================================
# Config
# ============================================================
OUTPUT_ROOT = Path("/Users/nancyh/Desktop/asthma-prs-study-fresh/09_ptrs-unified_model-evaluation/data/predictions")
BAL_N       = 64
N_BOOT      = 100
SEED_BASE   = 0  # per-iteration seed = SEED_BASE + i

ARTIFACT = OUTPUT_ROOT / "pairwise_comparisons"
ARTIFACT.mkdir(parents=True, exist_ok=True)
OUT_CSV       = ARTIFACT / "pairwise_bootstrap_auc.csv"
OUT_MATRIX    = ARTIFACT / "per_iteration_auc_matrix.csv"
OUT_MEANS     = ARTIFACT / "model_mean_auc.csv"

# ============================================================
# Model registry — (label, per-sample-preds source, filter)
# Each entry loads CAMP-only per-sample predictions for one model.
# ============================================================
def _load_generic_long(path, filters, score_col="score", y_col="y_true", id_col="sample_id"):
    df = pd.read_csv(path)
    for col, val in filters.items():
        df = df[df[col] == val]
    return df[[id_col, y_col, score_col]].rename(columns={id_col: "sample_id", y_col: "y_true", score_col: "score"})


def _load_per_feature_csv(path):
    """`predictions/best_consistent__<Feature>__<Model>__CAMP_only.csv` files
    from meta_model_{tissue,ct}: columns Sample_ID, y_true, y_pred."""
    df = pd.read_csv(path)
    return df.rename(columns={"Sample_ID": "sample_id", "y_pred": "score"})[["sample_id", "y_true", "score"]]


MODELS = []

# --- PRS baselines (from prscs_evaluation) ---
prs_pred = OUTPUT_ROOT / "prscs_evaluation" / "prs_predictions.csv"
for method in ("PRS-CS", "PRS-CSx"):
    MODELS.append({
        "label":  f"{method} alone (PRS+PCs)",
        "source": prs_pred,
        "loader": lambda p, m=method: _load_generic_long(
            p, {"method": m, "config": "PRS + Ancestry PCs", "eval_set": "CAMP-only"}),
    })

# --- Best FOCUS PTRS per modality (from meta_model best-valid picks) ---
# Read consistent_features.csv to get the best classifier per feature.
for mv, feat in (("tissue", "Esophagus_Mucosa"), ("ct", "cd4_naive")):
    cf = pd.read_csv(OUTPUT_ROOT / f"meta_model_{mv}" / "consistent_features.csv")
    row = cf[cf.Feature == feat].iloc[0]
    model_str = row.Model  # e.g. "Gradient Boosting"
    safe = model_str.replace(" ", "_").replace("(", "").replace(")", "").replace(".", "_")
    per_file = OUTPUT_ROOT / f"meta_model_{mv}" / "predictions" / f"best_consistent__{feat}__{safe}__CAMP_only.csv"
    MODELS.append({
        "label":  f"FOCUS {mv} — {feat} ({model_str})",
        "source": per_file,
        "loader": _load_per_feature_csv,
    })

# --- Unified PTRS alone (top-by-AUC classifier per modality, from meta_model) ---
uni_tissue = OUTPUT_ROOT / "meta_model_tissue" / "predictions" / "unified_ptrs_long.csv"
uni_ct     = OUTPUT_ROOT / "meta_model_ct"     / "predictions" / "unified_ptrs_long.csv"
for mv, path, top_method in (("tissue", uni_tissue, "RF (GridSearch)"),
                              ("ct",     uni_ct,     "Gradient Boosting")):
    MODELS.append({
        "label":  f"Unified PTRS {mv} — {top_method}",
        "source": path,
        "loader": lambda p, m=top_method: _load_generic_long(
            p, {"method": m, "cohort": "CAMP_only"}),
    })

# --- Uni-modal integrated PRS × PTRS (top-by-AUC classifier per (modality, PRS)) ---
uni_int_tissue = OUTPUT_ROOT / "integrated_ptrs_prs_tissue" / "integrated_predictions.csv"
uni_int_ct     = OUTPUT_ROOT / "integrated_ptrs_prs_ct"     / "integrated_predictions.csv"
uni_int_specs = [
    ("tissue", uni_int_tissue, "PRS-CS",  "Random Forest (tuned)"),
    ("tissue", uni_int_tissue, "PRS-CSx", "Gradient Boosting (tuned)"),
    ("ct",     uni_int_ct,     "PRS-CS",  "Random Forest (tuned)"),
    ("ct",     uni_int_ct,     "PRS-CSx", "Gradient Boosting (tuned)"),
]
for mv, path, prs, method in uni_int_specs:
    MODELS.append({
        "label":  f"Unified {mv} + {prs} — {method}",
        "source": path,
        "loader": lambda p, m=method, pr=prs: _load_generic_long(
            p, {"method": m, "prs_type": pr, "eval_set": "CAMP-only"}),
    })

# --- Cross-modal (top RF-tuned for each PRS variant) ---
xm_direct = OUTPUT_ROOT / "integrated_ptrs_prs_combined" / "direct_predictions.csv"
for prs in ("PRS-CS", "PRS-CSx"):
    MODELS.append({
        "label":  f"Cross-modal + {prs} — Random Forest (tuned)",
        "source": xm_direct,
        "loader": lambda p, pr=prs: _load_generic_long(
            p, {"method": "Random Forest (tuned)", "prs_type": pr, "eval_set": "CAMP-only"}),
    })

# ============================================================
# Bootstrap: compute per-iteration AUC for each model
# ============================================================
def bootstrap_auc_series(y_true, scores, n_boot=N_BOOT, bal_n=BAL_N):
    y = np.asarray(y_true).astype(int)
    s = np.asarray(scores)
    cases = np.where(y == 1)[0]
    ctrls = np.where(y == 0)[0]
    aucs = np.zeros(n_boot)
    for i in range(n_boot):
        rng = np.random.RandomState(SEED_BASE + i)
        sub = np.concatenate([rng.choice(cases, bal_n, replace=False), ctrls])
        aucs[i] = roc_auc_score(y[sub], s[sub])
    return aucs


# ============================================================
# Load + reconstruct per-iteration AUC series for every model
# ============================================================
per_model_auc = {}
per_model_mean = []
print(f"{'#':>2} {'model':60s} {'N':>5} {'cases':>6} {'AUC_mean':>9} {'AUC_SD':>9}")
print("-" * 96)
for idx, m in enumerate(MODELS):
    df = m["loader"](m["source"])
    df = df.dropna(subset=["score", "y_true"]).drop_duplicates("sample_id")
    if df.empty:
        print(f"{idx:>2} {m['label']:60s} SKIPPED (no rows)")
        continue
    aucs = bootstrap_auc_series(df["y_true"].values, df["score"].values)
    per_model_auc[m["label"]] = aucs
    per_model_mean.append({
        "model":     m["label"],
        "N":         len(df),
        "n_cases":   int((df["y_true"] == 1).sum()),
        "AUC_mean":  float(aucs.mean()),
        "AUC_SD":    float(aucs.std()),
    })
    print(f"{idx:>2} {m['label']:60s} {len(df):>5} {int((df['y_true']==1).sum()):>6} "
          f"{aucs.mean():>9.4f} {aucs.std():>9.4f}")

means_df = pd.DataFrame(per_model_mean).sort_values("AUC_mean", ascending=False).reset_index(drop=True)
means_df.to_csv(OUT_MEANS, index=False)
print(f"\nSaved -> {OUT_MEANS}  ({len(means_df)} models)")

# Also save the per-iteration AUC matrix (models × 100 iterations)
matrix_df = pd.DataFrame(per_model_auc)
matrix_df.index.name = "iteration"
matrix_df.to_csv(OUT_MATRIX)
print(f"Saved -> {OUT_MATRIX}  ({matrix_df.shape[0]} iterations × {matrix_df.shape[1]} models)")

# ============================================================
# Pairwise ΔAUC + empirical one-sided P
#   Delta_AUC_i = AUC_A_i - AUC_B_i (A is 'y-axis' model — the one hypothesized to be larger)
#   P_one_sided = fraction of iterations with Delta_AUC <= 0
# ============================================================
pairs = []
labels = list(per_model_auc.keys())
for A, B in product(labels, labels):
    if A == B:
        continue
    dA = per_model_auc[A] - per_model_auc[B]
    n_le0 = int((dA <= 0).sum())
    pairs.append({
        "model_A":                    A,
        "model_B":                    B,
        "mean_AUC_A":                 float(per_model_auc[A].mean()),
        "mean_AUC_B":                 float(per_model_auc[B].mean()),
        "mean_delta_AUC_AminusB":     float(dA.mean()),
        "n_iter_A_gt_B":              int((dA > 0).sum()),
        "n_iter_A_lt_B":              int((dA < 0).sum()),
        "n_iter_A_eq_B":              int((dA == 0).sum()),
        "empirical_one_sided_P":      n_le0 / N_BOOT,
        "significant_at_0_05":        (n_le0 / N_BOOT) < 0.05,
        "significant_at_0_01":        (n_le0 / N_BOOT) < 0.01,
    })

pairs_df = (pd.DataFrame(pairs)
              .sort_values(["mean_AUC_A", "mean_delta_AUC_AminusB"], ascending=[False, False])
              .reset_index(drop=True))
pairs_df.to_csv(OUT_CSV, index=False)
print(f"Saved -> {OUT_CSV}  ({len(pairs_df)} ordered pairs)")

# ============================================================
# Sanity: reproduce the paper's cited P = 0.13 for XM PRS-CSx RF vs Esophagus_Mucosa GB
# ============================================================
XM = next(l for l in labels if "Cross-modal + PRS-CSx" in l)
ESO = next((l for l in labels if l.startswith("FOCUS tissue")), None)
if ESO:
    row = pairs_df[(pairs_df.model_A == XM) & (pairs_df.model_B == ESO)].iloc[0]
    print(f"\n=== Sanity check ===")
    print(f"  A = {XM}")
    print(f"  B = {ESO}")
    print(f"  mean ΔAUC (A - B):      {row['mean_delta_AUC_AminusB']:+.4f}   (paper: +0.036)")
    print(f"  empirical one-sided P:  {row['empirical_one_sided_P']:.3f}      (paper: 0.13)")
