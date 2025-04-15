# Author: Max Balter
# Script: download.py
# Purpose: Download the most structurally similar and statistically significant AlphaFold matches
#          for significantly regulated uncharacterized proteins using Foldseek results.

import os           # I use this to manage folders and file paths.
import pandas as pd # I use pandas to load, sort, and filter the Foldseek hit tables.
import requests     # I use requests to download PDB files from AlphaFold’s public server.

# === STEP 1: Define input files and output folders ===

# These are my Foldseek search results against the three major AlphaFold databases.
# Each file contains uncharacterized proteins and their structural hits.
foldseek_result_files = [
    "prostt5_uncharacterized_results/results_afdb50.m8",
    "prostt5_uncharacterized_results/results_afdb_proteome.m8",
    "prostt5_uncharacterized_results/results_afdb_swissprot.m8"
]

# This is where I want to store the downloaded PDB structure files from AlphaFold.
output_folder = "top_hit_structures"
os.makedirs(output_folder, exist_ok=True)  # I make sure the folder exists or create it.

# This is the final summary file that will list each downloaded structure and its metadata.
summary_output_path = "prostt5_uncharacterized_results/downloaded_af_models_summary.csv"

# === STEP 2: Define the column layout of the .m8 Foldseek output files ===

# These column names are based on the actual structure of my .m8 files.
# Column 11 is E-value, and column 12 is bitscore.
column_names = [
    "query",      # 0 – UniProt ID of the uncharacterized protein
    "target",     # 1 – AlphaFold model matched
    "score",      # 2 – Some alignment score or % identity
    "q_len",      # 3 – Query length
    "t_len",      # 4 – Target length
    "gap_open",   # 5 – Gap opens
    "q_start",    # 6 – Start of alignment in query
    "q_end",      # 7 – End of alignment in query
    "t_start",    # 8 – Start in target
    "t_end",      # 9 – End in target
    "evalue",     # 10 – Match confidence
    "bitscore"    # 11 – Structural similarity score
]

# === STEP 3: Load and combine all result files ===

# I initialize an empty dataframe to store hits from all three databases.
all_hits = pd.DataFrame()

# I load each .m8 file, assign the correct column names, and tag the source database.
for file in foldseek_result_files:
    print(f"📂 Loading: {file}")
    df = pd.read_csv(file, sep='\t', header=None, names=column_names)
    df["source_db"] = os.path.basename(file)  # I add a column to trace where the hit came from
    all_hits = pd.concat([all_hits, df], ignore_index=True)

# === STEP 4: Apply filtering criteria ===

# I first filter for only confident hits — those with an E-value less than 1e-2.
evalue_cutoff = 1e-2
filtered_hits = all_hits[all_hits["evalue"] < evalue_cutoff]

# Then I sort the confident hits by bitscore in descending order,
# so the most structurally similar hits appear first.
filtered_hits_sorted = filtered_hits.sort_values(by="bitscore", ascending=False)

# I keep only the best structural match per query protein.
top_hits = filtered_hits_sorted.drop_duplicates(subset="query", keep="first")

# === STEP 5: Keep only AlphaFold targets (ignore other database hits) ===

# I only want to download AlphaFold predictions, so I filter targets starting with "AF-"
af_top_hits = top_hits[top_hits["target"].str.startswith("AF-")].copy()

# === STEP 6: Download AlphaFold models and track what I download ===

# I will log all successful downloads here for the final summary file.
downloaded_records = []

# This function downloads a PDB structure for a given AlphaFold target model.
def download_alphafold_model(model_id, query_id, evalue, bitscore, source_db):
    output_path = os.path.join(output_folder, f"{model_id}.pdb")

    # I skip the download if the file already exists locally.
    if os.path.exists(output_path):
        print(f"✅ Already exists: {model_id}")
        downloaded_records.append((query_id, model_id, evalue, bitscore, source_db))
        return

    # I generate the direct URL to fetch the PDB model from EBI's AlphaFold server.
    url = f"https://alphafold.ebi.ac.uk/files/{model_id}.pdb"

    try:
        response = requests.get(url)
        if response.status_code == 200:
            with open(output_path, "w") as f:
                f.write(response.text)
            print(f"⬇️ Downloaded: {model_id}")
            downloaded_records.append((query_id, model_id, evalue, bitscore, source_db))
        else:
            print(f"❌ Failed to download {model_id} (HTTP {response.status_code})")
    except Exception as e:
        print(f"❌ Error downloading {model_id}: {e}")

# I loop through all AlphaFold hits and download each structure.
print(f"🔍 Found {len(af_top_hits)} AlphaFold hits to download...")
for _, row in af_top_hits.iterrows():
    download_alphafold_model(
        model_id=row["target"],
        query_id=row["query"],
        evalue=row["evalue"],
        bitscore=row["bitscore"],
        source_db=row["source_db"]
    )

# === STEP 7: Write the download summary to a CSV ===

# I now convert my download records into a DataFrame and export to CSV.
summary_df = pd.DataFrame(downloaded_records, columns=[
    "Query UniProt ID",
    "AlphaFold Target Model",
    "E-value",
    "Bitscore",
    "Source Database"
])

summary_df.to_csv(summary_output_path, index=False)
print(f"📄 Summary written to: {summary_output_path}")
