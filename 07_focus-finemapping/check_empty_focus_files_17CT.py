import os
import pandas as pd

base_directory = "${PROJECT_DIR}/FOCUS/results_17CTnew/"

empty_files = []
non_empty_files = []
error_files = []

for root, dirs, files in os.walk(base_directory):
    for file in files:
        if file.endswith(".focus.tsv"):
            file_path = os.path.join(root, file)

            try:
                # Check file size first (fast)
                if os.path.getsize(file_path) == 0:
                    empty_files.append(file_path)
                    continue

                # Try reading header only
                df = pd.read_csv(file_path, sep="\t", nrows=1)

                if df.shape[1] == 0:
                    empty_files.append(file_path)
                else:
                    non_empty_files.append(file_path)

            except pd.errors.EmptyDataError:
                # Classic "No columns to parse from file"
                empty_files.append(file_path)

            except Exception as e:
                error_files.append((file_path, str(e)))

# ---------------------------
# Summary
# ---------------------------
print("====== FOCUS file integrity check (17CT) ======")
print(f"Total .focus.tsv files checked: {len(empty_files) + len(non_empty_files) + len(error_files)}")
print(f"Empty / unreadable files:      {len(empty_files)}")
print(f"Readable (non-empty) files:    {len(non_empty_files)}")
print(f"Other errors:                  {len(error_files)}")

# ---------------------------
# Save lists
# ---------------------------
with open("empty_focus_files_17CT.txt", "w") as f:
    for x in empty_files:
        f.write(x + "\n")

with open("error_focus_files_17CT.txt", "w") as f:
    for path, err in error_files:
        f.write(f"{path}\t{err}\n")

print("\nLists written to:")
print("  - empty_focus_files_17CT.txt")
print("  - error_focus_files_17CT.txt")
