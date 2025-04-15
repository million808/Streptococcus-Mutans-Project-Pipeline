# Author: Max Balter
# Script: download.py
# Purpose: Download top structural matches from AlphaFold for uncharacterized proteins based on Foldseek results and save summary.

import os  # For handling file paths and directories
import pandas as pd  # For reading and processing .m8 Foldseek result files
import requests  # For downloading AlphaFold PDB models

# === CONFIGURABLE PARAMETERS ===

# Only the three result files from the AlphaFold databases
foldseek_result_files = [
    "prostt5_uncharacterized_results/results_afdb50.m8",
    "prostt5_uncharacterized_results/results_afdb_proteome.m8",
    "prostt5_uncharacterized_results/results_afdb_swissprot.m8"
]

# Folder to save downloaded models
output_folder = "top_hit_structures"
os.makedirs(output_folder, exist_ok=True)

# Summary output file path
summary_output_path = "prostt5_uncharacterized_results/downloaded_af_models_summary.csv"

# Keep only the top hit (highest bitscore) per query
top_n = 1

# Foldseek .m8 column layout
column_names = [
    "query", "target", "q_len", "t_len", "evalue", "bitscore",
    "q_start", "q_end", "t_start", "t_end", "aln_len", "ident", "q_cov", "t_cov"
]

# === STEP 1: Load and merge all Foldseek result files ===

all_hits = pd.DataFrame()

for file in foldseek_result_files:
    print(f"📂 Loading: {file}")
    df = pd.read_csv(file, sep='\t', header=None, names=column_names)
    all_hits = pd.concat([all_hits, df], ignore_index=True)

# === STEP 2: Keep only the top N hits per query based on bitscore ===

all_hits_sorted = all_hits.sort_values(by="bitscore", ascending=False)
top_hits = all_hits_sorted.groupby("query").head(top_n)

# === STEP 3: Filter to only AlphaFold targets ===

af_top_hits = top_hits[top_hits["target"].str.startswith("AF-")].copy()

# === STEP 4: Download AlphaFold structures and track success ===

downloaded_records = []

def download_alphafold_model(model_id, query_id, bitscore):
    output_path = os.path.join(output_folder, f"{model_id}.pdb")

    # Skip if already exists
    if os.path.exists(output_path):
        print(f"✅ Already exists: {model_id}")
        downloaded_records.append((query_id, model_id, bitscore))
        return

    url = f"https://alphafold.ebi.ac.uk/files/{model_id}.pdb"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            with open(output_path, "w") as f:
                f.write(response.text)
            print(f"⬇️ Downloaded: {model_id}")
            downloaded_records.append((query_id, model_id, bitscore))
        else:
            print(f"❌ Failed to download {model_id} (HTTP {response.status_code})")
    except Exception as e:
        print(f"❌ Error downloading {model_id}: {e}")

# Loop through each row of the filtered hits
for _, row in af_top_hits.iterrows():
    query_id = row["query"]
    model_id = row["target"]
    bitscore = row["bitscore"]
    download_alphafold_model(model_id, query_id, bitscore)

# === STEP 5: Save summary CSV ===

summary_df = pd.DataFrame(downloaded_records, columns=[
    "Query UniProt ID", "AlphaFold Target Model", "Foldseek Bitscore"
])

summary_df.to_csv(summary_output_path, index=False)
print(f"📄 Summary written to: {summary_output_path}")
