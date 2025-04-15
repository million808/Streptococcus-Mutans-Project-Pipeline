# summarize.py

import os
import pandas as pd

# Define paths
structure_folder = "top_hit_structures"
results_folder = "prostt5_uncharacterized_results"
summary_output = os.path.join(results_folder, "top_hits_summary_final.csv")

# Load categorized hits
with open(os.path.join(results_folder, "hits_upregulated_uncharacterized.txt")) as f:
    upregulated = set(line.split('\t')[0] for line in f if line.strip())

with open(os.path.join(results_folder, "hits_downregulated_uncharacterized.txt")) as f:
    downregulated = set(line.split('\t')[0] for line in f if line.strip())

with open(os.path.join(results_folder, "hits_significant_uncharacterized.txt")) as f:
    significant = set(line.split('\t')[0] for line in f if line.strip())

# Detect actual downloaded structures
downloaded_targets = set()
for file in os.listdir(structure_folder):
    if file.endswith(".pdb"):
        parts = file.split("_", 1)
        if len(parts) == 2:
            target_id = parts[1].replace(".pdb", "")
            downloaded_targets.add(target_id)

print(f"🔎 Detected {len(downloaded_targets)} downloaded AlphaFold structures.")

# Read and parse each .m8 file to find matching hits
hits = []
for file in os.listdir(results_folder):
    if file.endswith(".m8"):
        path = os.path.join(results_folder, file)
        with open(path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 12:
                    query, target, bitscore = parts[0], parts[1], float(parts[11])
                    if target in downloaded_targets:
                        category = (
                            "Upregulated" if query in upregulated else
                            "Downregulated" if query in downregulated else
                            "Significant" if query in significant else
                            "Uncategorized"
                        )
                        hits.append({
                            "query_id": query,
                            "target_id": target,
                            "bitscore": bitscore,
                            "category": category
                        })

# Create and save DataFrame
df_hits = pd.DataFrame(hits)

if df_hits.empty:
    print("⚠️ No matches found between Foldseek results and downloaded structures.")
else:
    df_hits_sorted = df_hits.sort_values(by="bitscore", ascending=False)
    df_hits_sorted.to_csv(summary_output, index=False)
    print("✅ Final summary written to:", summary_output)
