# `sample_data/` — demo fixture for ma-FOCUS fine-mapping

These files let anyone run [`../run_demo.sh`](../run_demo.sh) end-to-end
(`focus import` + `focus finemap`, tissue track) without the full 12 GB GTEx
weight panel or the controlled-access study cohorts. **No individual-level study
data is included.**

Scope: **1 tissue (Whole_Blood) · chromosome 1 · 12 genes** — the same genes as
the [`05_twas-analysis/`](../../05_twas-analysis/) demo.

| File | Real or fake? | Provenance |
|------|---------------|------------|
| `weights/META_Whole_Blood/*.wgt.RDat` (12) | **REAL public data** | GTEx v8 EUR FUSION expression-weight models, Whole_Blood, 12 chr1 genes (same as the 05 demo). |
| `twas_whole_blood.demo.pos` | derived | FOCUS `.pos` (`PANEL WGT ID CHR P0 P1 N`); `WGT` is **relative to this directory** (`weights/…`) because `focus import` resolves weight paths against the `.pos` file's location. |
| `sumstat/chr1_formatted_meta_analysis.sumstats.mod2` | **REAL (derived)** | asthma meta-GWAS in FOCUS munged format (`SNP Z A2 A1 CHR BP N`), built from the 05 demo sumstats + the LD-reference base positions. Real Z-scores. |
| `expected_outputs/focus_result_demo.focus.tsv` | reference run | committed FOCUS output. The 90% credible set is **SDF4** (pip≈0.85) + **TNFRSF18** — the two highest-\|Z\| genes in the window. |

The LD reference (`1000G.EUR.1`) is **reused** from
`../../05_twas-analysis/sample_data/LDREF/` (not duplicated here).

## Result interpretation

Fine-mapping is per tissue, so PIPs sum to ~1 within the region and the credible
set collects the genes carrying the TWAS signal. Here SDF4 and TNFRSF18 (TWAS
Z ≈ −5.2 and −4.5 in the 05 demo) dominate; the other 10 window genes stay at
low PIP. This is a genuine — if tiny — reproduction of the method, not a
pre-baked table.

## Regenerating the fixture

Weights are copied from the 05 fixture; the `.mod2` is the 05 demo sumstats with
`CHR`/`BP` added from the LD-reference `.bim` and a representative `N`; the
`.pos` lists the 12 genes with `WGT` paths relative to `sample_data/`. Full
inputs: GTEx FUSION weights + the asthma sumstats from
[`01_meta-analysis/`](../../01_meta-analysis/).

## Cell-type (OneK1K) track

Not included — the OneK1K cell-type FUSION weight panels are not available to
redistribute here. The committed `../data/focus_credset_gene_tissue_17CT.tsv`
still lets the downstream cell-type gene-list + figure steps run. See the module
README, *Cell-type (OneK1K) track — pending*.
