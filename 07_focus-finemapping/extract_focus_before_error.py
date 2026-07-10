import glob

ERROR_KEY = "Unable to allocate"
FINEMAP_KEY = "focus finemap"

OUTPUT_FILE = "focus_finemap_before_errors.txt"

with open(OUTPUT_FILE, "w") as out:
    for slurm_file in sorted(glob.glob("slurm-*.out")):
        with open(slurm_file) as f:
            lines = f.readlines()

        error_indices = [
            i for i, line in enumerate(lines)
            if ERROR_KEY in line
        ]

        if not error_indices:
            continue

        for err_idx in error_indices:
            # walk backwards to find most recent "focus finemap"
            fm_start = None
            for i in range(err_idx, -1, -1):
                if FINEMAP_KEY in lines[i]:
                    fm_start = i
                    break

            if fm_start is None:
                continue

            # capture full finemap command block
            fm_block = []
            for i in range(fm_start, len(lines)):
                if lines[i].strip() == "":
                    break
                fm_block.append(lines[i])

            # write output
            out.write("=" * 80 + "\n")
            out.write(f"SLURM FILE: {slurm_file}\n")
            out.write(f"ERROR LINE: {lines[err_idx]}")
            out.write("\n--- LAST focus finemap COMMAND ---\n")
            out.writelines(fm_block)
            out.write("\n--- LOG CONTEXT (after command, before error) ---\n")
            out.writelines(lines[fm_start + len(fm_block):err_idx])
            out.write("\n")

print(f"Saved results to {OUTPUT_FILE}")

