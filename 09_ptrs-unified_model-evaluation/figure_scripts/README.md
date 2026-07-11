# Figure 4 — build scripts

Scripts that render the manuscript's **Figure 4** (integrated PRS + PTRS model
evaluation). Rendered outputs live in [`../figures/`](../figures/); the
per-sample prediction inputs are in [`../data/predictions/`](../data/) (this
module).

## Build pipeline (dependency chain)

Three scripts, run in order:

| Script | Reads | Writes |
|---|---|---|
| **`build_figure4_panels.py`** | per-sample prediction CSVs under `data/predictions/` (10 models; anchors resolved dynamically via `resolve_focus_anchor` / `resolve_twas_pt_anchor` from `best_valid_per_feature.csv`) | `figure4_panel{B_violins,C_pairwise,D_or}.{png,pdf}` + `figure4_revised_data.csv` (10-row per-model AUC/OR bootstrap summary) |
| **`build_figure4_combined.py`** | `figure4_panel{B,C,D}.png` | `figure4_combined.{png,pdf}` (B+C+D only) |
| **`build_figure4_pdf.py`** | `panelA.png` (curated schematic) + `figure4_panel{B,C,D}.png` | `Figure4.pdf` + `Figure4_preview.png` (final publication figure, A+B+C+D) |

```
data/predictions/*  ──build_figure4_panels.py──►  panel B/C/D + figure4_revised_data.csv
                                                        │
panelA.png ──────────────────────────┐                 ├──build_figure4_combined.py──► figure4_combined.{png,pdf}
                                      └──build_figure4_pdf.py──►  Figure4.pdf (+ preview)
```

### Configuring paths

All three resolve their working directory from the **`FIG4_ROOT`** environment
variable (default = the original lab path); override it to point at your own
copy:

```bash
export FIG4_ROOT=/path/to/your/figures-root      # contains data/predictions/ + figures/
python build_figure4_panels.py && python build_figure4_combined.py && python build_figure4_pdf.py
```

### Inputs NOT in this repository

- **`data/predictions/*`** — per-sample model predictions are **individual-level**
  data; they are kept only in the private repository, so `build_figure4_panels.py`
  cannot be run from the public repo. `figure4_revised_data.csv` is the aggregate
  (non-individual) summary it produces.
- **`panelA.png`** — the curated Panel A schematic (drawn externally); needed by
  `build_figure4_pdf.py`.

Requires: `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`, `Pillow`.
