# `sample_data/` — demo fixture for PRS-CS / PRS-CSx construction

These files let anyone run [`../run_demo.sh`](../run_demo.sh) end-to-end (PRS-CS
+ PRS-CSx + `plink --score`) without the official ~4.5 GB LD reference panel or
the controlled-access study cohorts. **No individual-level study data is
included.**

Scope: **chromosome 22 · 250 HapMap3 SNPs · 2 LD blocks**. A smoke test of the
pipeline mechanics, not a real risk score.

| File | Real or fake? | Provenance |
|------|---------------|------------|
| `reference/ldblk_1kg_eur/ldblk_1kg_chr22.hdf5` | **REAL public data** | LD reference built from 1000 Genomes Project Phase 3 **EUR** chr22 genotypes (503 samples): PLINK `--r` correlation matrices for 250 HapMap3 SNPs (16–25 Mb), split into 2 blocks, in the PRS-CS `blk_*/{ldblk,snplist}` HDF5 layout. |
| `reference/ldblk_1kg_eur/snpinfo_1kg_hm3` | **REAL (derived)** | PRS-CS SNP info (`CHR SNP BP A1 A2 MAF`) for those SNPs, from the 1000G `.bim`/`.frq`. |
| `reference/snpinfo_mult_1kg_hm3` | **REAL (derived)** | PRS-CSx multi-population SNP info (`… FRQ_{AFR,AMR,EAS,EUR,SAS} FLP_…`). This is an **EUR-only demo**, so every `FRQ` column carries the EUR MAF and `FLP=1`; only the EUR panel is provided. |
| `demo_sumstats_chr22.txt` | **FAKE — demonstration only** | GWAS summary statistics (`SNP A1 A2 BETA SE`) with **simulated** effect sizes (random draws; 8 planted "causal" SNPs), *not* real asthma associations. The format matches the [`01_meta-analysis/`](../../01_meta-analysis/) output; in production this file is that real meta-GWAS. |
| `target/target_chr22.{bed,bim,fam}` | **REAL public data** | 1000 Genomes EUR chr22 subset (100 samples, same 250 SNPs) — a public stand-in for the controlled-access validation cohort (GACRS/CAMP). |
| `expected_outputs/` | reference run | committed PRS-CS + PRS-CSx posterior effect files and the PRS `.profile`. MCMC is seed-fixed (`--seed 9500`), so a re-run reproduces these **exactly** (max \|Δbeta\| = 0). |

## Why the LD panel is safe to ship

The official PRS-CS panels (`ldblk_1kg_eur`, ~4.5 GB) can't be redistributed
here, so this fixture **rebuilds a tiny one from public 1000G chr22 genotypes**
in the exact format PRS-CS/PRS-CSx expect (`snpinfo_*` + per-chr `ldblk_*.hdf5`
with `blk_N/{ldblk,snplist}`). It covers only the 250 demo SNPs — enough for the
tools to run, far too small for real prediction.

## Regenerating the fixture

```
plink --bfile <1000G_EUR_chr22> --extract <hm3 chr22> --maf 0.05 --thin-count 250 \
      --make-bed --out ref250
plink --bfile ref250 --freq --out ref250
plink --bfile ref250 --r square --out ref250
# then: pack ref250.{bim,frq,ld} into the ldblk HDF5 + snpinfo files, simulate
# the fake sumstats, and subset 100 samples as the target (see build.py).
```

Full inputs for a real run: the official [PRS-CS LD panels](https://github.com/getian107/PRScs#getting-started)
and the asthma meta-GWAS from [`01_meta-analysis/`](../../01_meta-analysis/).
