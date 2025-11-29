#!/usr/bin/env python3
import os
import pandas as pd
import subprocess
import sys

# === Define all paths (organized once at top) ===
base_dir = os.getcwd()
fastq_dir = os.path.join(base_dir, "fastq_output")
results_dir = os.path.join(base_dir, "results")
salmon_index = os.path.join(base_dir, "salmon_index")
gtf_uncompressed = os.path.join(base_dir, "gencode.v43.annotation.gtf")
gtf_gz = os.path.join(base_dir, "gencode.v43.annotation.gtf.gz")
events_dir = os.path.join(base_dir, "events_ioe")
suppa_out_dir = os.path.join(base_dir, "suppa_results")
tpm_out = os.path.join(base_dir, "TPM_SUPPA2_final.tsv")

# === Create required directories ===
for d in [fastq_dir, results_dir, events_dir, suppa_out_dir]:
    os.makedirs(d, exist_ok=True)

# === SRA sample IDs: read all from SraRunTable.txt ===
# sep=None + engine="python" lets pandas auto-detect the delimiter (tab, comma, etc.)
sra_table = pd.read_csv("SraRunTable.txt", sep=None, engine="python")

print("Columns in SraRunTable.txt:")
print(list(sra_table.columns))

# Try to find a column whose name (stripped, case-insensitive) matches 'Run'
run_cols = [c for c in sra_table.columns if c.strip().lower() == "run"]

if not run_cols:
    sys.exit("ERROR: Could not find a 'Run' column in SraRunTable.txt. Please inspect the printed column names.")

run_col = run_cols[0]  # first matching column

ids = sra_table[run_col].astype(str).tolist()

print("Found", len(ids), "SRA IDs:")
print(ids)

# === Build Salmon index ===
subprocess.run("salmon index --gencode -t gencode.v43.transcripts.fa.gz -i salmon_index -p 8 -k 31", shell=True)

# === Step 1–4: prefetch → fasterq-dump → salmon quant ===
for sra in ids:
    print(f"\n===== Processing {sra} =====")
    os.system(f"prefetch -v {sra} -O {fastq_dir}/ ")
    os.system(f"fasterq-dump {sra} --split-3 -O {fastq_dir}/")
    os.system(f"salmon quant -i {salmon_index} -l A -r {fastq_dir}/{sra}.fastq -p 8 --validateMappings -o {results_dir}/{sra}_quant")
    os.system(f"rm -fr {fastq_dir}/*")

# === Collect all quant.sf files ===
quant_files = sorted([
    os.path.join(results_dir, d, "quant.sf")
    for d in os.listdir(results_dir)
    if os.path.isdir(os.path.join(results_dir, d)) and d.endswith("_quant")
])

if not quant_files:
    sys.exit("No quant.sf files found in results/. Run Salmon quant first.")

print(f"Found {len(quant_files)} quant.sf files")

# === Build TPM matrix ===
dfs = []
for f in quant_files:
    sample = os.path.basename(os.path.dirname(f)).replace("_quant", "")
    df = pd.read_csv(f, sep="\t", usecols=["Name", "TPM"])
    df["Name"] = df["Name"].str.split("|").str[0]
    dfs.append(df.rename(columns={"TPM": sample}))

merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on="Name", how="outer")

# Sort transcripts for consistency
merged = merged.sort_values(by="Name")

# === Correct SUPPA2 format ===
merged = merged.set_index("Name")  #most important part for SUPPA2

# Write header (sample names only) and data (no "Name" column)
header_line = "\t".join(merged.columns)
with open(tpm_out, "w") as f:
    f.write(header_line + "\n")
    merged.to_csv(f, sep="\t", header=False)

print("TPM_SUPPA2_final.tsv created successfully")

# === Check GTF file ===
if os.path.exists(gtf_uncompressed):
    gtf_file = gtf_uncompressed
elif os.path.exists(gtf_gz):
    gtf_file = gtf_gz
else:
    sys.exit("GTF file not found (gencode.v43.annotation.gtf or .gz).")

# === SUPPA2: generate IOE events ===
generate_cmd = [
    "suppa.py", "generateEvents",
    "-i", gtf_file,
    "-o", os.path.join(events_dir, "events"),
    "-f", "ioe",
    "-e", "SE", "SS", "MX", "RI", "FL"
]
print("Generating IOE event files...")
subprocess.run(generate_cmd, check=True)
# === SUPPA2: calculate PSI values ===
print("Calculating PSI per event...")

ri_ioe = os.path.join(events_dir, "events_RI_strict.ioe")

if os.path.exists(ri_ioe):
    subprocess.run([
        "suppa.py", "psiPerEvent",
        "-i", ri_ioe,
        "-e", tpm_out,
        "-o", os.path.join(suppa_out_dir, "RI.psi")
    ], check=True)

print("SUPPA2 PSI calculation complete.")
