"""
Count FOCUS gene-context pairs retained under a PIP >= 0.1 threshold WITHOUT the
credible-set filter.

Companion to post_process.py / post_process_17CT.py, which build the
`focus_credset_gene_tissue_{39,17CT}.tsv` files by:
    walk *.focus.tsv -> drop NULL.MODEL -> split ens_gene_id into (symbol, tissues)
    -> filter in_cred_set_pop1 == 1 -> drop_duplicates((tissues, symbol))

Here we do the SAME aggregation/cleaning but:
    * DO NOT filter on in_cred_set_pop1
    * keep gene-context pairs with pips_pop1 >= 0.1
so we can see how many pairs / unique genes survive a pure PIP >= 0.1 cut.

For reference the script also reproduces the credible-set count and the
credset AND PIP>=0.1 count (the latter is what the manuscript's "retained"
tissue number reflects).

Run:
    python count_pip_no_credset.py
"""
import os
from pathlib import Path
import pandas as pd

FOCUS_DIR = Path(os.environ.get("FOCUS_DIR", "/expanse/lustre/projects/ddp412/n5huang/FOCUS"))
PIP_THRESHOLD = 0.1

MODALITIES = [
    ("TISSUE (GTEx, 39)",       FOCUS_DIR / "results_39new"),
    ("CELLTYPE (OneK1K, 17)",   FOCUS_DIR / "results_17CTnew"),
]


def load_focus_dir(results_dir):
    """Aggregate every *.focus.tsv under results_dir, mirroring post_process.py's
    cleaning EXCEPT the credible-set filter. Returns the full (uncut) DataFrame."""
    frames = []
    for root, _dirs, files in os.walk(results_dir):
        for fn in files:
            if not fn.endswith(".focus.tsv"):
                continue
            try:
                df = pd.read_csv(os.path.join(root, fn), sep="\t")
            except pd.errors.EmptyDataError:
                # empty FOCUS output (tracked in FOCUS/empty_focus_files*.txt) -> skip
                continue
            # drop NULL.MODEL rows (same as post_process.py)
            df = df[df["ens_gene_id"].notna()]
            df = df[df["ens_gene_id"] != "NULL.MODEL"]
            # ens_gene_id == "ENSG..._tissue" -> symbol = ENSG id, tissues = context
            df[["symbol", "tissues"]] = df["ens_gene_id"].str.split("_", n=1, expand=True)
            frames.append(df)
    if not frames:
        raise RuntimeError(f"No .focus.tsv files under {results_dir}")
    out = pd.concat(frames, ignore_index=True)
    out["pips_pop1"] = pd.to_numeric(out["pips_pop1"], errors="coerce")
    return out


def summarize(label, results_dir):
    full = load_focus_dir(results_dir)
    n_files_rows = len(full)

    # --- reproduce the credible-set file (sanity check) -------------------
    cred = full[full["in_cred_set_pop1"] == 1].drop_duplicates(
        subset=["tissues", "symbol"], keep="first")
    # credset AND PIP>=0.1 (what post_process.py prints as "pip filtered")
    cred_pip = cred[cred["pips_pop1"] >= PIP_THRESHOLD]

    # --- PIP >= 0.1 WITHOUT the credible-set filter -----------------------
    # keep pairs with any block at PIP>=threshold, then dedup to unique pairs
    pip_only = full[full["pips_pop1"] >= PIP_THRESHOLD].drop_duplicates(
        subset=["tissues", "symbol"], keep="first")

    print(f"=== {label} : {results_dir.name} ===")
    print(f"  raw rows (post-clean, no filter)         : {n_files_rows:,}")
    print(f"  [sanity] credible-set pairs (dedup)      : {len(cred):,}   "
          f"unique genes: {cred['symbol'].nunique():,}")
    print(f"  [ref]    credset AND PIP>={PIP_THRESHOLD} pairs      : {len(cred_pip):,}   "
          f"unique genes: {cred_pip['symbol'].nunique():,}")
    print(f"  >>> PIP>={PIP_THRESHOLD}, NO credset filter, pairs   : {len(pip_only):,}   "
          f"unique genes: {pip_only['symbol'].nunique():,}")
    print()
    return pip_only


if __name__ == "__main__":
    for label, d in MODALITIES:
        summarize(label, d)
