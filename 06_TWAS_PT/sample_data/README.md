# `sample_data/` — demo fixture for TWAS P+T LD-clumping

These files let anyone run [`../run_demo.sh`](../run_demo.sh) end-to-end (the
tissue-track gene-selection method) without the full GTEx predicted-expression
release. **No individual-level study data is included.**

Scope: **Whole_Blood · 12 chr1 genes** (the same genes as the 05 TWAS demo).

| File | Real or fake? | Provenance |
|------|---------------|------------|
| `designmat_Whole_Blood_demo.RData` (object `df1`) | **REAL public data** | GTEx v8 predicted-expression matrix (`designmat_Whole_Blood_v8_320EUR_double`), subset to the 12 demo gene columns (489 samples × 12 genes). |
| `Marginal_alphas_NEW_TWAS_Whole_Blood.demo.txt.gz` | **REAL (derived)** | TWAS Z-scores for the 12 genes, taken from the [`05_twas-analysis/`](../../05_twas-analysis/) demo FUSION run. One column, row-aligned with the transcripts file. |
| `TranscriptsInWhole_BloodModel.demo.txt` | derived | the 12 gene ids (versioned), row-aligned with the `df1` columns. |
| `TranscriptsInWhole_BloodModel_keep.demo.txt` | derived | heritability keep-list (all 12 genes for the demo). |
| `expected_outputs/*.txt` | reference run | committed gene lists at the four p-thresholds + their per-chromosome splits. |

Gene annotation is **not duplicated here** — the runner reuses the copy already
committed at `../07_focus-finemapping/data/gene_annotation.txt.gz`.

## Regenerating the fixture

`df1` is subset from the full GTEx designmat to the 12 demo gene columns; the
TWAS Z-scores are copied from the 05 demo output; the transcripts / keep lists
are the matching gene ids. Full inputs come from the GTEx v8 / TCSC release
(see the module README).

## Cell-type (OneK1K) track

Not included — the OneK1K cell-type predicted-expression CSVs are not available
to redistribute here. See the module README, *Cell-type (OneK1K) track — pending*.
