import os
import requests

# === CONFIGURABLE PARAMETERS ===
m8_file = "prostt5_results/results_afdb_proteome.m8"  # <-- I can change this
output_folder = "top_hit_structures"
top_n = 10  # <-- How many top hits to process

# === Ensure output directory exists ===
os.makedirs(output_folder, exist_ok=True)

# === Read and sort the .m8 file by bit score (column 12) ===
print(f"Reading and sorting {m8_file}...")
with open(m8_file) as f:
    lines = [line.strip().split('\t') for line in f]

# Sort by bit score (column 12) in descending order
sorted_lines = sorted(lines, key=lambda x: float(x[11]), reverse=True)
top_hits = sorted_lines[:top_n]

print(f"Found {len(top_hits)} top hits.")

# === Download each target structure ===
for i, hit in enumerate(top_hits):
    query_id, target_id = hit[0], hit[1]
    structure_file = os.path.join(output_folder, f"{i+1:02d}_{target_id}.pdb")

    # Try to interpret the target_id
    if target_id.upper().startswith("AF-"):
        # Likely AlphaFold ID
        uniprot_id = target_id.split("-")[1]
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    elif len(target_id) == 4 and target_id.isalnum():
        # Likely PDB ID
        url = f"https://files.rcsb.org/download/{target_id.upper()}.pdb"
    else:
        print(f"[{i+1}] Skipping unknown format: {target_id}")
        continue

    # Download structure
    print(f"[{i+1}] Downloading {target_id} from: {url}")
    response = requests.get(url)
    if response.status_code == 200:
        with open(structure_file, 'w') as out:
            out.write(response.text)
        print(f"✅ Saved to {structure_file}")
    else:
        print(f"❌ Failed to download {target_id} (status {response.status_code})")

print("🎉 Done downloading top hit structures!")
