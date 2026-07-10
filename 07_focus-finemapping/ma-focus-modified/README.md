## ma-focus (modified for asthma-prs-study)

This is a vendored, lightly modified copy of [`mancusolab/ma-focus`](https://github.com/mancusolab/ma-focus) at upstream commit `f3a38fd` (pyfocus v0.803, 2024-04-26), used by [`07_focus-finemapping/`](../).

### Why we modified it

A subset of our FUSION weight files (the meta-analyzed bulk-tissue panels — the `META_*` directories such as `META_Whole_Blood/`, `META_Adipose_Subcutaneous/`, etc.) store only a single META model column in `wgt.matrix` and do **not** include a `cv.performance` object. Stock `pyfocus.models.convert.import_fusion()` assumes every `.wgt.RDat` carries the FUSION multi-model layout (`top1` / `blup` / `lasso` / `enet`) plus a `cv.performance` table for model selection, so `focus import` raises on the META panels.

### What changed

A single file: **`pyfocus/models/convert.py`**, function `import_fusion`. The modification adds a branch inside the per-gene loop:

1. Detect whether the loaded `.wgt.RDat` has a `cv.performance` object **and** a multi-column `wgt.matrix`.
2. **If yes** — original CV-based model selection (pick the non-`top1` method with the best R², fall back to `top1`). Behavior unchanged for standard FUSION panels.
3. **If no** (single-column META layout) — take the only column directly, label `method="meta"`, and leave `cv.R2` / `cv.R2.pval` as `None`.

Also added `robj.r("rm(list = ls())")` at the top of each per-gene iteration so leftover R objects from the previous gene cannot leak into the next when consecutive weight files have different schemas (e.g., a META file immediately after a FUSION file).

The remaining diff vs upstream is cosmetic (dict comprehensions, dropped redundant comments, regex escape fix). The full diff is:

```bash
diff -u <(git -C <fresh-clone> show f3a38fd:pyfocus/models/convert.py) pyfocus/models/convert.py
```

### How we use it

Install in a Python 3.8 conda env (same as the upstream instructions):

```bash
conda create -n ma-focus python=3.8
conda activate ma-focus
cd 07_focus-finemapping/ma-focus-modified
pip install .
pip install mygene 'rpy2==3.5.12'
```

After `pip install .`, the `focus` CLI behaves like stock ma-focus; the `import` subcommand now silently accepts META-only weight files.

### What's NOT in this directory

Only the ma-focus source tree (the Python package, CLI launcher, setup files, LICENSE). The analysis-side artifacts that lived alongside it in our working tree on Expanse — weight directories (`B_intermediate/`, `META_*/`, `v8_allEUR_*/`, …), per-gene SQLite DBs (`db_17CTnew/`, `db_39new/`, …), FOCUS result files, slurm logs, core dumps — are not version-controlled. Re-run the [`07_focus-finemapping/`](../) pipeline to regenerate them.

### Upstream

Track upstream changes against:

- Repo: https://github.com/mancusolab/ma-focus
- Fork base: `f3a38fd` (Merge PR #24 from `update_zl`, 2024-04-26)
- License: MIT (preserved — see `LICENSE`)

---

Original upstream README follows.

---

[![Github](https://img.shields.io/github/stars/mancusolab/sushie?style=social)](https://github.com/mancusolab/sushie)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Project generated with PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)

FOCUS & MA-FOCUS
=====
FOCUS (Fine-mapping Of CaUsal gene Sets) is software to fine-map transcriptome-wide association study statistics at genomic risk regions. The software takes as input summary GWAS data along with eQTL weights and outputs a credible set of _genes_ to explain observed genomic risk.

MA-FOCUS (Multi-Ancestry Fine-mapping Of CaUsal gene Sets) is an extension of FOCUS that leverages summary GWAS data with eQTL weights from multiple ancestries to increase the precision of credible sets for causal genes.

```diff
- We detest usage of our software or scientific outcome to promote racial discrimination.
```

FOCUS is described in:

> [Probabilistic fine-mapping of transcriptome-wide association studies](https://www.nature.com/articles/s41588-019-0367-1). Nicholas Mancuso, Malika K. Freund, Ruth Johnson, Huwenbo Shi, Gleb Kichaev, Alexander Gusev, and Bogdan Pasaniuc. ***Nature Genetics*** 51, 675-682 (2019).

MA-FOCUS is described in:

> [Multi-ancestry fine-mapping improves precision to identify causal genes in transcriptome-wide association studies](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00306-8). Zeyun Lu<sup>\*</sup>, Shyamalika Gopalan<sup>\*</sup>, Dong Yuan, David V. Conti, Bogdan Pasaniuc, Alexander Gusev, Nicholas Mancuso. ***American Journal of Human Genetics*** VOLUME 109, ISSUE 8, P1388-1404, AUGUST 04, 2022.

<sup>\*</sup> indicates equal contribution

**We greatly appreciate users cite both papers if they use our software.**

  [**Installation**](#installation)
  | [**Example**](#get-started-with-example)
  | [**Version History**](#version-history)
  | [**Support**](#support)
  | [**Other Software**](#other-software)

## Installation

Before installing ma-focus, we recommend users to create a new conda environment:

``` bash
conda create -n ma-focus python=3.8
```

Then activate the environment:

``` bash
conda activate ma-focus
```

Users can download the latest repository and then use `pip`:

``` bash
git clone https://github.com/mancusolab/ma-focus.git
cd ma-focus
pip install .
```


## Get Started with Example

### Data Preparation
We will use the white blood count (WBC) GWAS measured in European and African ancestry individuals from [Chen et al. 2020](10.1016/j.cell.2020.06.045) as an example. The eQTL weights are computed from LCLs of American European and African individuals in GENOA study ([Shang et al. 2020](10.1016/j.ajhg.2020.03.002)). The LD reference is from the 1000 Genomes Project.

Users can download all the testing data from [here](https://www.mancusolab.com/ma-focus).

### FOCUS Example

Here is an example of how to perform fine-mapping :

    focus finemap EA_WBC_new_munged_chr22.tsv.gz ./1000GP3/EUR/1000G.EUR.QC.allelesAligned.22 ea_fusion_genoa.db --chr 22 --prior-prob "gencode38" --locations 38:EUR --out ./focus_result

This command will scan `EA_WBC_new_munged_chr22.tsv.gz` for risk regions `38:EUR` generated by LDetect on GRCh38 for European ancestry and then perform TWAS+fine-mapping using LD estimated from plink-formatted `1000G.EUR.QC.allelesAligned.22` and eQTL weights from `ea_fusion_genoa.db` GRCh38.

### MA-FOCUS Example
Here is an example of how to perform multi-ancestry fine-mapping from European and African ancestry (just use `:` to concatenate):

    focus finemap EA_WBC_new_munged_chr22.tsv.gz:AA_WBC_new_munged_chr22.tsv.gz ./1000GP3/EUR/1000G.EUR.QC.allelesAligned.22:./1000GP3/AFR/1000G.AFR.QC.allelesAligned.22 ea_fusion_genoa.db:aa_fusion_genoa.db --chr 22 --prior-prob "gencode38" --locations 38:EUR-AFR --out ./ma-focus_result

Note: `--locations` is a new required parameter compared to versions prior to v0.8. This parameter specifies the genomic regions to be fine-mapped. We recommend to use independent regions for single ancestry or multiple ancestries.

Please see the [wiki](https://github.com/mancusolab/ma-focus/wiki) for more details on how to use focus, ma-focus and links to database files.

## Version History

| Version | Description                                                                                                                                                                                                                                             |
|---------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 0.3 | Initial release. More to come soon.                                                                                                                                                                                                                     |
 | 0.4 | Added FUSION import support.                                                                                                                                                                                                                            |
 | 0.5 | Plotting sorts genes based on tx start. Various bug fixes that limited the number of queried SNPs and plotting when using newer matplotlib.                                                                                                             |
| 0.6 | Fixed bug where only one of the two alleles was reversed complemented breaking alignment. For now these instances are dropped. Added option `--use-ens-id` for FUSION import to indicate the main model label is an Ensembl ID rather than HGNC symbol. |
 | 0.6.5 | Fixed bug in newer versions of matplotlib not accepting string for colormaps. Fixed legend bug in plot. Fixed bug that mismatched string and category when supplying custom locations.                                                                  |
 |0.6.10 | Fixed bug where weight database allele mismatch with GWAS broke inference.                                                                                                                                                                              |
 | 0.8 | Added MA-FOCUS support. Added GWAS imputation using imp-G. Added additional choice for prior probability for causal genes (number of genes in the risk regions).                                                                                        |
| 0.801 | Added gencode_map_v38 and multiple LD block files in GRCh38. Fixed prior_prob bugs.                                                                                                                                                                     |
| 0.802 | Fix the bug that .gitignore includes *.tsv so that gencode files couldn't be pushed to github.                                                                                                                                                          |
| 0.803 | Fix the bugs in pandas/numpy outdated version, intercept.                                                                                                                                                                                               |
| 0.9 | Convert to pyscaffold style, fix bugs and typos, and update documentation and read me.                                                                                                                                                                  |

## Software and support

If users have any questions or comments, please contact Zeyun Lu (<zeyunlu@usc.edu>) and Nicholas Mancuso (<nmancuso@usc.edu>).

Feel free to use other software developed by [Mancuso
Lab](https://www.mancusolab.com/):

-   [SuShiE](https://github.com/mancusolab/sushie): a scalable variational approach to perform SNP fine-mapping for molecular data across multiple ancestries.
-   [SuSiE-PCA](https://github.com/mancusolab/susiepca): a scalable Bayesian variable selection technique for sparse principal component analysis
-   [twas_sim](https://github.com/mancusolab/twas_sim): a Python software to simulate [TWAS](https://www.nature.com/articles/ng.3506) statistics.
-   [FactorGo](https://github.com/mancusolab/factorgo): a scalable variational factor analysis model that learns pleiotropic factors from GWAS summary statistics.
-   [HAMSTA](https://github.com/tszfungc/hamsta): a Python software to  estimate heritability explained by local ancestry data from  admixture mapping summary statistics.
-   [Traceax](https://github.com/tszfungc/traceax): a Python library to perform stochastic trace estimation for linear operators.

------------------------------------------------------------------------

This project has been set up using PyScaffold 4.1.1. For details and
usage information on PyScaffold see <https://pyscaffold.org/>.
