#!/usr/bin/env python3
"""
run_demo.py — self-contained smoke test for the 09 integrated PRS + PTRS model.

Reproduces the module's headline computation on the committed SYNTHETIC fixture
(sample_data/) — no individual-level study data. It mirrors
`integrated_ptrs_prs_unified.ipynb`:

  1. per-feature OOF transform  — calibrate each feature's PTRS (Keep_Vector)
     through its own classifier (out-of-fold on train, fit-then-predict on test)
  2. Direct integration          — [cd4_naive_OOF, Esophagus_Mucosa_OOF, PRS]
     -> final classifier
  3. CAMP-only Balanced bootstrap — 64 vs 64 x 100 resamples -> AUC mean +/- std
     + odds ratio, compared against PRS-alone (the headline delta-AUC)

Train cohort = synthetic "GACRS train+test"; test cohort = synthetic "CAMP-only
Balanced". Deterministic (fixed seeds) so it reproduces expected_outputs/.

Requires: scikit-learn, pandas, numpy  (the `test` conda env has these).
Run:  cd 09_ptrs-unified_model-evaluation && python run_demo.py
"""
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
SD   = os.path.join(HERE, "sample_data")
OUT  = os.path.join(HERE, "demo_run"); os.makedirs(OUT, exist_ok=True)

def load(feature, ct_dir, ts_dir):
    ct = pd.read_csv(os.path.join(SD, ct_dir, f"{feature}_results.csv"))
    return ct.rename(columns={"Keep_Vector": feature})[["Sample_ID", feature, "asthma"]]

def build(split):  # 'gacrs_train_test' or 'camp'
    cd4 = load("cd4_naive",       f"ptrs_results-concat_17CT_{split}", None)
    eso = load("Esophagus_Mucosa", f"ptrs_results-concat_39_{split}",  None)
    df  = cd4.merge(eso[["Sample_ID", "Esophagus_Mucosa"]], on="Sample_ID")
    prs = pd.read_csv(os.path.join(SD, "prs_demo.csv"))
    return df.merge(prs, on="Sample_ID")

train = build("gacrs_train_test")
test  = build("camp")
FEATS = ["cd4_naive", "Esophagus_Mucosa"]
yq = train.asthma.values

def per_feature_oof(clf_factory):
    """OOF prob on train, fit-then-predict prob on test, per feature."""
    cv = StratifiedKFold(5, shuffle=True, random_state=0)
    tr, te = {}, {}
    for f in FEATS:
        Xtr = train[[f]].values; Xte = test[[f]].values
        tr[f] = cross_val_predict(clf_factory(), Xtr, yq, cv=cv, method="predict_proba")[:, 1]
        te[f] = clf_factory().fit(Xtr, yq).predict_proba(Xte)[:, 1]
    return tr, te

rf = lambda: RandomForestClassifier(n_estimators=200, max_depth=3, random_state=0)
lr = lambda: LogisticRegression(max_iter=1000)
tr_oof, te_oof = per_feature_oof(rf)

def score_columns(which):
    """Return dict: model_name -> per-sample predicted score on CAMP test."""
    Xtr_int = np.column_stack([tr_oof[f] for f in FEATS] + [train.PRS.values])
    Xte_int = np.column_stack([te_oof[f] for f in FEATS] + [test.PRS.values])
    out = {}
    out["PRS alone"] = test.PRS.values
    for name, mk in [("Cross-modal + PRS (RF)", rf), ("Cross-modal + PRS (LR)", lr)]:
        out[name] = mk().fit(Xtr_int, yq).predict_proba(Xte_int)[:, 1]
    # PTRS-only (no PRS) reference
    Xtr_p = np.column_stack([tr_oof[f] for f in FEATS])
    Xte_p = np.column_stack([te_oof[f] for f in FEATS])
    out["Cross-modal PTRS (no PRS, RF)"] = rf().fit(Xtr_p, yq).predict_proba(Xte_p)[:, 1]
    return out

scores = score_columns("camp")
yt = test.asthma.values
cases = np.where(yt == 1)[0]; ctrls = np.where(yt == 0)[0]

def balanced_bootstrap_auc(pred, n_boot=100, n=64):
    aucs = []
    for i in range(n_boot):
        r = np.random.RandomState(i)
        idx = np.concatenate([r.choice(cases, n, replace=True), r.choice(ctrls, n, replace=True)])
        aucs.append(roc_auc_score(yt[idx], pred[idx]))
    return np.mean(aucs), np.std(aucs)

def odds_ratio(pred):
    z = (pred - pred.mean()) / (pred.std() + 1e-12)
    b = LogisticRegression(max_iter=1000).fit(z.reshape(-1, 1), yt).coef_[0, 0]
    return float(np.exp(b))

rows = []
for name, pred in scores.items():
    m, s = balanced_bootstrap_auc(pred)
    rows.append({"Model": name, "CAMP_bal_AUC": round(m, 4),
                 "AUC_std": round(s, 4), "OR_per_SD": round(odds_ratio(pred), 3)})
res = pd.DataFrame(rows).sort_values("CAMP_bal_AUC", ascending=False).reset_index(drop=True)
base = res.loc[res.Model == "PRS alone", "CAMP_bal_AUC"].iloc[0]
res["dAUC_vs_PRS"] = (res.CAMP_bal_AUC - base).round(4)

res.to_csv(os.path.join(OUT, "summary_records.csv"), index=False)
print("\n=== CAMP-only Balanced (64v64 x 100 bootstrap) — synthetic demo ===")
print(res.to_string(index=False))
rank1 = res.iloc[0]
print(f"\nRank-1: {rank1.Model}  AUC={rank1.CAMP_bal_AUC}+/-{rank1.AUC_std}  "
      f"OR={rank1.OR_per_SD}  dAUC vs PRS={rank1.dAUC_vs_PRS}")
print(f"\nwrote {os.path.join(OUT, 'summary_records.csv')}")

# regression check against committed reference
exp = os.path.join(SD, "expected_outputs", "summary_records.csv")
if os.path.exists(exp):
    a = res.set_index("Model").CAMP_bal_AUC
    b = pd.read_csv(exp).set_index("Model").CAMP_bal_AUC
    d = (a - b).abs().max()
    print(f"\nmax |dAUC| vs committed reference: {d:.2e}  -> {'OK' if d < 1e-9 else 'DIFF'}")
