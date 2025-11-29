
#!/usr/bin/env python3
import pandas as pd
import os
import sys

# === File paths ===
base_dir = os.getcwd()
sra_file = os.path.join(base_dir, "SraRunTable.csv")
geo_file = os.path.join(base_dir, "GSE181294_conditions.csv")

# === Check files ===
if not os.path.exists(sra_file):
    sys.exit("❌ SraRunTable.csv not found. Please download it first.")
if not os.path.exists(geo_file):
    sys.exit("❌ GSE181294_conditions.csv not found. Please download or generate it from GEO.")
# === Read both tables ===
sra = pd.read_csv(sra_file)
geo = pd.read_csv(geo_file)

# Clean quotes if present
geo["Sample"] = geo["Sample"].astype(str).str.replace('"', '').str.strip()
sra["Sample"] = sra["Sample"].astype(str).str.replace('"', '').str.strip()

# === Merge on GSM Sample ID ===
merged = pd.merge(sra, geo, on="Sample", how="left")

# === Clean and summarize ===
print(" Merge complete.")
print(f"Total samples: {len(merged)}")
print(" Condition distribution:")
print(merged["Condition"].value_counts(dropna=False))
print("Grade distribution:")
print(merged["Grade"].value_counts(dropna=False))

# === Save merged table ===
out_file = os.path.join(base_dir, "merged_SRA_GEO.csv")
merged.to_csv(out_file, index=False)
print(f"\n Saved merged metadata to: {out_file}")

# === Extract SRR IDs for a specific group ===
target_group = merged[
    (merged["Condition"] == "tumor") &
    (merged["Grade"] == "high")
]["Run"].dropna().unique().tolist()

if not target_group:
    print("\n No SRR IDs found for Condition='tumor' and Grade='high'.")
else:
    print(f"\n Found {len(target_group)} SRR IDs for high-grade tumor samples:")
    for srr in target_group:
        print(" ", srr)

    # Save to file
    with open(os.path.join(base_dir, "selected_SRR_ids.txt"), "w") as f:
        for srr in target_group:
            f.write(srr + "\n")
    print("\n Saved SRR list to selected_SRR_ids.txt")
