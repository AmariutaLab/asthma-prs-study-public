# `sample_data/` — demo inputs for the ancestry pipeline

These files let anyone run [`../run_demo.sh`](../run_demo.sh) end-to-end without
access to the controlled-access study cohorts. **No individual-level study data
is included here.**

| File | Real or fake? | Provenance |
|------|---------------|------------|
| `demo_1kg.chr22.{bed,bim,fam}` | **REAL public data** | 1000 Genomes Project Phase 3, European-ancestry (EUR) samples, chromosome 22 only, downsampled to **100 samples / 6,000 SNPs** (biallelic ACGT SNPs, `--thin-count`, seed 7). This stands in for the "study cohort merged with the 1000 Genomes reference" fileset that Step 01 expects. |
| `w_hm3_chr22_subset.snplist` | **REAL public data** | HapMap3 SNP list (`w_hm3.snplist`, Alkes Group, [Zenodo 7773502](https://zenodo.org/records/7773502)), subset to the chr22 rsIDs present in the demo genotypes. Header (`SNP A1 A2`) preserved. |
| `expected_outputs/` | reference run | Committed outputs of `run_demo.sh`: the `.filt-hm3.{bed,bim,fam}` fileset, ADMIXTURE `.{1,2,3}.Q` ancestry proportions, and `cv_error.txt`. CV error favours *K*=1 (single-ancestry demo), as expected. |

> [!NOTE]
> The PLINK steps are bitwise-reproducible (the `.filt-hm3` fileset is identical
> across PLINK 1.9 b7.2–b7.7). ADMIXTURE at *K*=1 and *K*=2 is stable, but at the
> over-specified *K*=3 on this single-ancestry demo the likelihood surface is flat,
> so `cv_error`/`.3.Q` can vary slightly run-to-run — that variation is ADMIXTURE
> convergence, not a PLINK-version effect.

## Why single-ancestry?

The public 1000 Genomes subset used here is EUR-only, so ADMIXTURE
cross-validation error will favour *K*=1. This is intentional: the demo is a
**smoke test of the pipeline mechanics**, not a biological ancestry result. Point
the same commands at a genuinely multi-ancestry fileset (e.g. a study cohort
merged with the full 1000 Genomes panel across all chromosomes) to recover real
population structure.

## Regenerating the demo data

The demo genotypes were carved from a full 1000 Genomes chr22 PLINK fileset:

```bash
# SRC = 1000 Genomes EUR chr22 PLINK prefix (.bed/.bim/.fam)
plink --bfile "$SRC" \
      --keep first_100_samples.txt \
      --snps-only just-acgt --biallelic-only \
      --thin-count 6000 --seed 7 \
      --make-bed --out demo_1kg.chr22

# w_hm3_chr22_subset.snplist = header + rows of w_hm3.snplist whose rsID is on chr22
```

A public 1000 Genomes fileset can be obtained from the
[PLINK 2.0 resources page](https://www.cog-genomics.org/plink/2.0/resources) or
the [1000 Genomes FTP](https://www.internationalgenome.org/data).
