"""
Feature-count ablation for the cross-modal PTRS + PRS integration.

Tests whether adding a third MA-FOCUS consistent-shortlist feature to the
two-feature anchor pair (cd4_naive + Esophagus_Mucosa) + PRS improves the
CAMP-Balanced AUC. Runs the same Random-Forest (tuned) integration classifier
used by the rank-1 cross-modal model in `integrated_ptrs_prs_unified.ipynb`.

Design:
  Anchors:           cd4_naive OOF + Esophagus_Mucosa OOF (loaded from disk)
  Candidates:        Consistent-shortlist features minus the two anchors
                     (read dynamically from meta_model_{tissue,ct}/consistent_features.csv)
  PRS variants (2):  PRS-CS, PRS-CSx (altPRS: phi=auto, LDREF=EUR / META)
  Integration:       Random Forest (tuned) via GridSearchCV, same 108-combo grid
                     as integrated_ptrs_prs_unified.ipynb
  Bootstrap:         100 case-matched 64-vs-64 iterations of CAMP-only

Baselines (2-feature anchor + PRS): NOT recomputed here. Reference AUCs used
for Delta_AUC come from integrated_ptrs_prs_unified.ipynb / all_results.csv:
  PRS-CS  Cross-modal RF (tuned): AUC = 0.6266, SD = 0.0387, OR = 3.319
  PRS-CSx Cross-modal RF (tuned): AUC = 0.6320, SD = 0.0399, OR = 3.546

Inputs (via INPUT_ROOT):
  INPUT_ROOT/combine/data/1/ptrs_results-concat_39_gacrs_train_test/     (tissue raw)
  INPUT_ROOT/combine/data/1/ptrs_results-concat_39_camp_1k1k/            (tissue CAMP raw)
  INPUT_ROOT/combine/data/1/ptrs_results-concat_17CT_gacrs_train_test/   (CT raw)
  INPUT_ROOT/combine/data/1/ptrs_results-concat_17CT_camp_gtex/          (CT CAMP raw)
  INPUT_ROOT/files/03_prscs-*.csv                                        (altPRS melt CSVs)

  (via OUTPUT_ROOT)
  OUTPUT_ROOT/meta_model_tissue/consistent_features.csv                  (candidate shortlist)
  OUTPUT_ROOT/meta_model_ct/consistent_features.csv                      (candidate shortlist)
  OUTPUT_ROOT/integrated_ptrs_prs_combined/per_feature_oof.csv           (anchor OOF)

Outputs (under OUTPUT_ROOT/feature_ablation/):
  ablation_results.csv
  ablation_log.txt

Usage:
  python feature_count_ablation.py                # default full 108-combo grid
  python feature_count_ablation.py --reduced-grid # optional 24-combo grid (~10x faster)
"""
import sys, time, warnings, gc, argparse
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import StratifiedKFold, cross_val_predict, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import roc_auc_score

# ============================================================
# INPUT_ROOT / OUTPUT_ROOT (matches meta_model / integrated notebooks)
# ============================================================
INPUT_ROOT  = Path("/Users/nancyh/Desktop/hartwell/gene_model/score")
OUTPUT_ROOT = Path("/Users/nancyh/Desktop/asthma-prs-study-fresh/09_ptrs-unified_model-evaluation/data/predictions")

FILES_DIR   = INPUT_ROOT / "files"
PTRS_ROOT   = INPUT_ROOT / "combine" / "data" / "1"
ARTIFACT    = OUTPUT_ROOT / "feature_ablation"
ARTIFACT.mkdir(parents=True, exist_ok=True)
OUT_CSV     = ARTIFACT / "ablation_results.csv"
LOG_PATH    = ARTIFACT / "ablation_log.txt"

# ============================================================
# Argparse
# ============================================================
_argp = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
_argp.add_argument("--reduced-grid", action="store_true",
                   help="Use a reduced 24-combo RF grid instead of the paper's "
                        "108-combo grid (~10x faster; qualitative results unchanged).")
ARGS = _argp.parse_args()

SEED   = 42
BAL_N  = 64
N_BOOT = 100
ANCHORS = ("cd4_naive", "Esophagus_Mucosa")

MANUSCRIPT_BASELINES = {
    "PRS-CS":  {"AUC": 0.6266, "SD": 0.0387, "OR": 3.319},
    "PRS-CSx": {"AUC": 0.6320, "SD": 0.0399, "OR": 3.546},
}

# ============================================================
# Logging
# ============================================================
LOG_PATH.write_text("")
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a") as f: f.write(line + "\n")

_grid_label = "reduced 24-combo" if ARGS.reduced_grid else "full 108-combo (paper)"
log(f"INPUT_ROOT  = {INPUT_ROOT}")
log(f"OUTPUT_ROOT = {OUTPUT_ROOT}")
log(f"Starting feature-count ablation (RF-tuned integration, "
    f"{_grid_label} grid, {N_BOOT} bootstraps)")

# ============================================================
# 1. Load anchor OOF + labels
# ============================================================
oof = pd.read_csv(OUTPUT_ROOT / "integrated_ptrs_prs_combined" / "per_feature_oof.csv")
def anchor(feat):
    d = oof[oof.feature == feat]
    return {c: d[d.cohort == c].set_index("sample_id") for c in
            ["GACRS_train_OOF", "GACRS_test", "CAMP_only"]}

A_cd4 = anchor("cd4_naive")
A_eso = anchor("Esophagus_Mucosa")
TRAIN_IDS = A_cd4["GACRS_train_OOF"].index.tolist()
CAMP_IDS  = A_cd4["CAMP_only"].index.tolist()
y_train   = A_cd4["GACRS_train_OOF"].loc[TRAIN_IDS, "y_true"].astype(int).values
y_camp    = A_cd4["CAMP_only"].loc[CAMP_IDS, "y_true"].astype(int).values
CD4_TR = A_cd4["GACRS_train_OOF"].loc[TRAIN_IDS, "score"].values
CD4_CM = A_cd4["CAMP_only"].loc[CAMP_IDS, "score"].values
ESO_TR = A_eso["GACRS_train_OOF"].loc[TRAIN_IDS, "score"].values
ESO_CM = A_eso["CAMP_only"].loc[CAMP_IDS, "score"].values
log(f"Anchors: train={len(TRAIN_IDS)} ({y_train.sum()} cases), "
    f"CAMP={len(CAMP_IDS)} ({y_camp.sum()} cases)")

# ============================================================
# 2. Build candidate list dynamically from consistent_features shortlists
# ============================================================
cand_tissue = pd.read_csv(OUTPUT_ROOT / "meta_model_tissue" / "consistent_features.csv")
cand_ct     = pd.read_csv(OUTPUT_ROOT / "meta_model_ct"     / "consistent_features.csv")

def _kind_for(model_str):
    """Map meta_model best-classifier name to a per-feature integration classifier."""
    if model_str.startswith("Ridge"):  return "Ridge"
    if model_str.startswith("Gradient"): return "GB"
    if model_str.startswith("RF"):       return "RF"
    if model_str.startswith("Lasso"):    return "Ridge"  # linear surrogate
    if model_str.startswith("Elastic"):  return "Ridge"
    return "Ridge"

CANDIDATES = []
for _, r in cand_tissue.iterrows():
    if r["Feature"] in ANCHORS: continue
    CANDIDATES.append((r["Feature"], "39", _kind_for(r["Model"])))
for _, r in cand_ct.iterrows():
    if r["Feature"] in ANCHORS: continue
    CANDIDATES.append((r["Feature"], "17CT", _kind_for(r["Model"])))

log(f"Candidate features ({len(CANDIDATES)} total): "
    f"{len([c for c in CANDIDATES if c[1]=='39'])} tissues, "
    f"{len([c for c in CANDIDATES if c[1]=='17CT'])} CTs")
for feat, modir, kind in CANDIDATES:
    log(f"  - {feat} ({modir}, best classifier via {kind})")

# ============================================================
# 3. PRS altPRS loader
# ============================================================
def load_altprs(prs_type):
    if prs_type == "PRS-CS":
        gacrs_f = FILES_DIR / "03_prscs-prscsx-camp-gtex-onek1k-visualization_prscs-gacrs-only-data-melt-admixture.csv"
        camp_f  = FILES_DIR / "03_prscs-prscsx-camp-gtex-onek1k-visualization_prscs_camp-1k1k-data-melt-admixture.csv"
        phi_col, ldref = "PRS-CS(\u03d5)", "EUR"
    else:
        gacrs_f = FILES_DIR / "03_prscs-prscsx-camp-gtex-onek1k-visualization_prscsx_gacrs-only-data-melt-admixture.csv"
        camp_f  = FILES_DIR / "03_prscs-prscsx-camp-gtex-onek1k-visualization_prscsx_camp-1k1k-data-melt-admixture.csv"
        phi_col, ldref = "PRS-CSx(\u03d5)", "META"
    frames = []
    for f in (gacrs_f, camp_f):
        d = pd.read_csv(f)
        d = d[(d[phi_col] == "\u03d5=auto") & (d["LDREF"] == ldref)]
        frames.append(d.drop_duplicates(subset=["IID"])[["IID", "PRS"]]
                       .rename(columns={"IID": "sample_id"}))
    return pd.concat(frames, ignore_index=True).drop_duplicates("sample_id") \
             .set_index("sample_id")["PRS"]

def prs_z(prs_type):
    raw = load_altprs(prs_type)
    tr = raw.reindex(TRAIN_IDS); cm = raw.reindex(CAMP_IDS)
    mu, sd = tr.mean(), tr.std()
    return ((tr - mu) / sd).values, ((cm - mu) / sd).values

PRS_TR, PRS_CM = {}, {}
for pt in ["PRS-CS", "PRS-CSx"]:
    PRS_TR[pt], PRS_CM[pt] = prs_z(pt)
    log(f"  {pt} altPRS: train NaN={np.isnan(PRS_TR[pt]).sum()}, "
        f"CAMP NaN={np.isnan(PRS_CM[pt]).sum()}")

# ============================================================
# 4. Candidate raw-PTRS loader + per-feature OOF
# ============================================================
def get_raw(feat, modir):
    dirs = ("ptrs_results-concat_39_gacrs_train_test", "ptrs_results-concat_39_camp_1k1k") \
           if modir == "39" else \
           ("ptrs_results-concat_17CT_gacrs_train_test", "ptrs_results-concat_17CT_camp_gtex")
    frames = []
    for d in dirs:
        f = PTRS_ROOT / d / f"{feat}_results.csv"
        if f.exists():
            df = pd.read_csv(f)
            df["sample_id"] = df["Sample_ID"]
            frames.append(df[["sample_id", "Keep_Vector"]])
    return pd.concat(frames, ignore_index=True).drop_duplicates("sample_id") \
             .set_index("sample_id")["Keep_Vector"]

def make_pfclf(kind):
    if kind == "Ridge":
        return LogisticRegression(penalty="l2", C=0.01, solver="liblinear",
                                   max_iter=1000, random_state=SEED)
    if kind == "GB":
        return GradientBoostingClassifier(n_estimators=100, max_depth=2,
                                           learning_rate=0.05, min_samples_leaf=20,
                                           subsample=0.8, random_state=SEED)
    if kind == "RF":
        return GridSearchCV(
            RandomForestClassifier(random_state=SEED),
            param_grid={
                "n_estimators":     [50, 100, 200],
                "max_depth":        [2, 3, 5, None],
                "min_samples_leaf": [10, 20, 50],
                "max_features":     [1],
            },
            cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED),
            scoring="roc_auc", n_jobs=1, refit=True,
        )
    raise ValueError(kind)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
CAND_TR, CAND_CM = {}, {}
for feat, modir, kind in CANDIDATES:
    t0 = time.time()
    raw = get_raw(feat, modir)
    tr = raw.reindex(TRAIN_IDS); cm = raw.reindex(CAMP_IDS)
    mu, sd = tr.mean(), tr.std()
    Xtr = ((tr - mu) / sd).values.reshape(-1, 1)
    Xcm = ((cm - mu) / sd).values.reshape(-1, 1)
    est = make_pfclf(kind)
    if isinstance(est, GridSearchCV):
        est.fit(Xtr, y_train); best = est.best_estimator_
    else:
        best = est
    oof_tr = cross_val_predict(best, Xtr, y_train, cv=skf, method="predict_proba")[:, 1]
    best.fit(Xtr, y_train)
    p_cm = best.predict_proba(Xcm)[:, 1]
    CAND_TR[feat] = oof_tr; CAND_CM[feat] = p_cm
    del est, best; gc.collect()
    log(f"  per-feat OOF [{feat} via {kind}]: {time.time()-t0:.1f}s")

# ============================================================
# 5. Integration RF-tuned + CAMP-Balanced bootstrap
# ============================================================
def make_integration_rf():
    if ARGS.reduced_grid:
        param_grid = {
            "n_estimators":     [100, 200],
            "max_depth":        [3, 5, None],
            "min_samples_leaf": [10, 20],
            "max_features":     ["sqrt", "log2"],
        }
    else:  # paper's 108-combo grid, matches integrated_ptrs_prs_unified.ipynb
        param_grid = {
            "n_estimators":     [50, 100, 200],
            "max_depth":        [2, 3, 5, None],
            "min_samples_leaf": [10, 20, 50],
            "max_features":     [1, "sqrt", "log2"],
        }
    return GridSearchCV(
        RandomForestClassifier(random_state=SEED),
        param_grid=param_grid,
        cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED),
        scoring="roc_auc", n_jobs=1, refit=True,
    )

def bootstrap_camp(y, s):
    """CAMP-Balanced 100 x 64v64 bootstrap. Matches integrated_ptrs_prs_unified's evaluate()."""
    y = np.asarray(y); s = np.asarray(s)
    cases = np.where(y == 1)[0]; ctrls = np.where(y == 0)[0]
    aucs = np.zeros(N_BOOT); ors = np.zeros(N_BOOT)
    for i in range(N_BOOT):
        rng = np.random.RandomState(i)
        sub = np.concatenate([rng.choice(cases, BAL_N, replace=False), ctrls])
        yy = y[sub]; ss = s[sub]
        aucs[i] = roc_auc_score(yy, ss)
        q1, q3 = np.percentile(ss, [25, 75])
        top = ss >= q3; bot = ss <= q1
        a = int(((yy == 1) & top).sum()); b = int(((yy == 0) & top).sum())
        c = int(((yy == 1) & bot).sum()); d = int(((yy == 0) & bot).sum())
        ors[i] = (a * d) / (b * c) if (b * c) > 0 else np.nan
    return aucs, ors

def run_experiment(X_tr, y_tr, X_cm, y_cm):
    clf = make_integration_rf()
    clf.fit(X_tr, y_tr)
    score = clf.predict_proba(X_cm)[:, 1]
    aucs, ors = bootstrap_camp(y_cm, score)
    del clf; gc.collect()
    return {
        "CAMP_Balanced_AUC_mean": float(np.nanmean(aucs)),
        "CAMP_Balanced_AUC_SD":   float(np.nanstd(aucs)),
        "CAMP_Balanced_OR_mean":  float(np.nanmean(ors)),
    }

# ============================================================
# 6. Run 3-feature experiments
# ============================================================
all_rows = []
def checkpoint():
    if all_rows: pd.DataFrame(all_rows).to_csv(OUT_CSV, index=False)

for prs in ["PRS-CS", "PRS-CSx"]:
    log(f"\n===== {prs}  (baseline reference: AUC = "
        f"{MANUSCRIPT_BASELINES[prs]['AUC']:.4f}) =====")
    for feat, modir, kind in CANDIDATES:
        t0 = time.time()
        Xtr = np.column_stack([CD4_TR, ESO_TR, CAND_TR[feat], PRS_TR[prs]])
        Xcm = np.column_stack([CD4_CM, ESO_CM, CAND_CM[feat], PRS_CM[prs]])
        mtr = ~np.isnan(Xtr).any(axis=1); mcm = ~np.isnan(Xcm).any(axis=1)
        r = run_experiment(Xtr[mtr], y_train[mtr], Xcm[mcm], y_camp[mcm])
        modality = "GTEx tissue" if modir == "39" else "OneK1K CT"
        r.update({
            "Experiment_Type": "3-feature (anchor + candidate + PRS)",
            "Third_Feature":   feat,
            "Modality":        modality,
            "PRS_Variant":     prs,
        })
        all_rows.append(r); checkpoint()
        log(f"  3F [{feat} + {prs}]: AUC={r['CAMP_Balanced_AUC_mean']:.4f}"
            f"\u00b1{r['CAMP_Balanced_AUC_SD']:.4f} "
            f"OR={r['CAMP_Balanced_OR_mean']:.3f} ({time.time()-t0:.1f}s)")

# ============================================================
# 7. Delta vs baselines + summary
# ============================================================
res = pd.DataFrame(all_rows)
res["Baseline_AUC"]          = res.PRS_Variant.map(lambda p: MANUSCRIPT_BASELINES[p]["AUC"])
res["Delta_AUC_vs_baseline"] = res["CAMP_Balanced_AUC_mean"] - res["Baseline_AUC"]
res = res[[
    "Experiment_Type", "PRS_Variant", "Third_Feature", "Modality",
    "CAMP_Balanced_AUC_mean", "CAMP_Balanced_AUC_SD", "CAMP_Balanced_OR_mean",
    "Baseline_AUC", "Delta_AUC_vs_baseline",
]]
res.to_csv(OUT_CSV, index=False)
log(f"\nWROTE: {OUT_CSV}  ({len(res)} rows)")

log(f"\n===== SUMMARY =====")
log(f"3-feature experiments: {len(res)}")
log(f"  Delta AUC < 0 (dropped vs manuscript baseline):  "
    f"{(res.Delta_AUC_vs_baseline < 0).sum()}")
log(f"  Delta AUC > 0 (improved vs manuscript baseline): "
    f"{(res.Delta_AUC_vs_baseline > 0).sum()}")
log(f"  mean Delta AUC:  {res.Delta_AUC_vs_baseline.mean():+.4f}")
log("Done.")
