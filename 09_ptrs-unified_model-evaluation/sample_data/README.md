# `sample_data/` — synthetic fixture for the 09 evaluation demo

These files let anyone run [`../run_demo.py`](../run_demo.py) end-to-end (the
cross-modal integrated PRS + PTRS model) without any individual-level study data.

> [!IMPORTANT]
> **Everything here is FULLY SYNTHETIC.** The real per-sample PTRS/PRS
> predictions in `../data/` are individual-level and are **not** used by the
> demo. This fixture only reproduces the *structure* the notebooks expect
> (`Sample_ID, Keep_Vector, asthma` per-feature CSVs + a PRS table), with values
> drawn from a random model that plants a weak signal so the integrated model
> beats PRS-alone — illustrating the headline result, not reproducing it.

## Layout (matches the real `data/` schema)

| Path | Cohort | Feature |
|------|--------|---------|
| `ptrs_results-concat_17CT_gacrs_train_test/cd4_naive_results.csv` | GACRS train+test | cd4_naive (cell type) |
| `ptrs_results-concat_39_gacrs_train_test/Esophagus_Mucosa_results.csv` | GACRS train+test | Esophagus_Mucosa (tissue) |
| `ptrs_results-concat_17CT_camp/cd4_naive_results.csv` | CAMP-only Balanced (64v64) | cd4_naive |
| `ptrs_results-concat_39_camp/Esophagus_Mucosa_results.csv` | CAMP-only Balanced (64v64) | Esophagus_Mucosa |
| `prs_demo.csv` | both | per-sample `PRS` |

Each per-feature CSV has columns `Sample_ID, Keep_Vector, asthma` — the three
columns the evaluation actually consumes (the real files carry extra unused
columns like `Sign_Vector`, `V1`, `IID_base`).

## Generating

```bash
cd 09_ptrs-unified_model-evaluation/sample_data
python generate_sample_data.py      # deterministic, seed=42
```

`generate_sample_data.py` draws `cd4_naive`, `Esophagus_Mucosa`, and `PRS` from
independent normals and sets the case/control label from
`logit = 1.1·cd4 + 1.0·eso + 0.6·PRS + noise` — so the two PTRS features carry
orthogonal signal and PRS adds a weaker independent component (mirroring *why*
cross-modal integration helps). The GACRS cohort is imbalanced (140/400 cases);
the CAMP test set is balanced 64 vs 64.

## Expected outputs

`expected_outputs/summary_records.csv` — the committed reference produced by
`run_demo.py`: CAMP-Balanced bootstrap AUC (± std), odds ratio, and ΔAUC vs
PRS-alone for PRS-alone, cross-modal PTRS (no PRS), and the two integrated
models. Deterministic — a re-run matches it exactly (max |ΔAUC| = 0).

The AUCs here are *higher* than the manuscript's ~0.63 because the planted signal
is deliberately clean; what the demo reproduces is the **ordering and mechanism**
(integrated > PTRS-only > PRS-alone), not the paper's exact numbers.

## Full pipeline vs demo

`run_demo.py` is a compact stand-in for `integrated_ptrs_prs_unified.ipynb`
(per-feature OOF → Direct integration → CAMP-Balanced bootstrap). The full
notebooks additionally run the 7-classifier zoo, both altPRS variants
(PRS-CS/PRS-CSx), the TWAS-P+T thresholds, and the summary-table / figure /
pairwise-P steps against the real `data/` — those require the individual-level
inputs and are not part of this public demo.
