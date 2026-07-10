# `sample_data/` — demo fixture for the FUSION TWAS step

These files let anyone run [`../run_demo.sh`](../run_demo.sh) end-to-end (FUSION
`assoc_test` on the **tissue track**) without the full 12 GB GTEx weight panel or
the controlled-access study cohorts. **No individual-level study data is included.**

Scope of the demo: **1 tissue (Whole_Blood) · chromosome 1 · 12 genes** in the
first ~2.7 Mb of chr1. This is a smoke test of the pipeline mechanics, not a
result.

| File | Real or fake? | Provenance |
|------|---------------|------------|
| `weights/META_Whole_Blood/*.wgt.RDat` (12) | **REAL public data** | GTEx v8 European-ancestry FUSION expression-weight models (`v8_320EUR` meta panel), Whole_Blood, for 12 chr1 genes (NOC2L, PLEKHN1, HES4, ISG15, …). 11.5 KB each. |
| `LDREF/1000G.EUR.1.{bed,bim,fam}` | **REAL public data** | 1000 Genomes Project Phase 3 EUR LD reference (FUSION `LDREF` release), chr1 restricted to 0–3.5 Mb (978 SNPs). |
| `demo_sumstats_chr1.txt` | **REAL (derived)** | Asthma meta-GWAS summary statistics (`SNP Z A2 A1`) produced by [`01_meta-analysis/`](../../01_meta-analysis/), sliced to the SNPs in the demo LD reference. Real association statistics; just row-subset. |
| `twas_whole_blood.demo.pos` | derived | FUSION `.pos` table (`PANEL WGT ID CHR P0 P1 N`) over the 12 demo genes; `WGT` paths point at `weights/`. |
| `expected_outputs/twas_whole_blood_chr1.dat` | reference run | Committed FUSION output for regression comparison. The `TWAS.Z` values match the full-panel production run to <1e-2 (the reduced LD window does not change results for these genes). |

## Regenerating the fixture

Carved from the study's local FUSION assets:

```bash
# from the full whole-blood .pos, take the first 12 chr1 genes that have a
# local .wgt.RDat, copy those weights, and write a demo .pos pointing at them
# (WGT = META_Whole_Blood/<gene>.wgt.RDat).

# LD reference: one chromosome, restricted to the demo window
plink --bfile 1000G.EUR.1 --chr 1 --from-bp 1 --to-bp 3500000 \
      --make-bed --out LDREF/1000G.EUR.1

# sumstats: rows of formatted_meta_analysis.txt whose SNP is in the demo .bim
```

Public sources for the full inputs: GTEx FUSION weights and the `LDREF` panel
from the [FUSION website](http://gusevlab.org/projects/fusion/); asthma sumstats
from [`01_meta-analysis/`](../../01_meta-analysis/).

## Cell-type (OneK1K) track

Not included — the OneK1K cell-type eQTL weight panels are not available to
redistribute here. See the top-level module README, *Cell-type (OneK1K) track —
pending*.
