import os
import pandas as pd

# Base directory containing tissue subfolders
base_directory = "${PROJECT_DIR}/FOCUS/results_17CTnew/"

all_processed_dfs = []

for root, dirs, files in os.walk(base_directory):
    for file in files:
        if file.endswith(".focus.tsv"):
            file_path = os.path.join(root, file)
            print(f"Processing file: {file_path}")

            try:
                df = pd.read_csv(file_path, sep="\t")

                # ---------------------------
                # Remove NULL.MODEL rows
                # ---------------------------
                df = df[df["ens_gene_id"].notna()]
                df = df[df["ens_gene_id"] != "NULL.MODEL"]

                # ---------------------------
                # Split gene + tissue
                # ENSGxxxx_tissue_name
                # ---------------------------
                df[["symbol", "tissues"]] = (
                    df["ens_gene_id"]
                    .str.split("_", n=1, expand=True)
                )

                # ---------------------------
                # Filter: credible set only
                # ---------------------------
                df = df[df["in_cred_set_pop1"] == 1]

                all_processed_dfs.append(df)

            except Exception as e:
                print(f"Error processing {file_path}: {e}")

# ---------------------------
# Combine everything
# ---------------------------
if not all_processed_dfs:
    raise RuntimeError("No .focus.tsv files processed.")

final_combined_df = pd.concat(all_processed_dfs, ignore_index=True)

# ---------------------------
# Deduplicate by (tissue, gene)
# ---------------------------
final_unique = (
    final_combined_df
    .drop_duplicates(subset=["tissues", "symbol"], keep="first")
)

print("Combined shape (cred set only):", final_unique.shape)


# ---------------------------
# Filter: PIP threshold
# ---------------------------
final_pip = final_unique[final_unique["pips_pop1"] >= 0.1]

print("Combined shape (pip filtered)):", final_pip.shape)


# ---------------------------
# Save final output
# ---------------------------
output_path = "${PROJECT_DIR}/FOCUS/focus_credset_gene_tissue_17CT.tsv"
final_unique.to_csv(output_path, sep="\t", index=False)

print(f"Final output written to:\n{output_path}")
