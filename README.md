# A multi-ancestry framework integrating polygenic and transcriptomic risk for childhood-onset asthma prediction

**Contributors**:

* Nan Huang (n5huang@ucsd.edu)
* Michelle Franc Ragsac (mragsac@ucsd.edu)
* Tiffany Amariuta (tamariutabartell@ucsd.edu)

This repository contains the code and associated notes for the pediatric asthma polygenic risk score (PRS) project.

## General Analysis Steps

| Step | Directory | Description |
|------|-----------|-------------|
| 1 | [`01_meta-analysis/`](01_meta-analysis/) | IVW fixed-effects meta-analysis of TAGC + GBMI European-ancestry asthma GWAS summary statistics |
| 2 | [`02_data-harmonization/`](02_data-harmonization/) | Data harmonization |
| 3 | [`03_ancestry-analysis/`](03_ancestry-analysis/) | ADMIXTURE ancestry analysis |
| 4 | [`04_prscs-prscsx-construction/`](04_prscs-prscsx-construction/) | PRS-CS / PRS-CSx construction |
| 5 | [`05_twas-analysis/`](05_twas-analysis/) | FUSION TWAS marginal-alpha computation per (chr, tissue) + sorted_gene / `.pos` artifacts consumed by FOCUS and PTRS construction |
| 7 | [`07_focus-finemapping/`](07_focus-finemapping/) | ma-FOCUS TWAS fine-mapping → per-(chr, tissue) credible-set gene lists for step 8 |
| 8 | [`08_ptrs-construction/`](08_ptrs-construction/) | PTRS construction (per-(chr, tissue) scoring + cross-chromosome concat) |
| 9 | [`09_ptrs-unified_model-evaluation/`](09_ptrs-unified_model-evaluation/) | Unified PTRS exploration + PTRS/PRS integration notebooks |
