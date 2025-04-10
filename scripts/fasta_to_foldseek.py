# Author: Max Balter
# Script: 04_submit_fasta_to_foldseek_api.py (ProstT5 edition, updated to scan all databases)

# ======================= OVERVIEW ==========================
# I originally designed this script to query Foldseek using PDB structures.
# However, many of my differentially expressed proteins don't have available structures.
# To resolve this, I'm switching to using Foldseek's CLI tools along with ProstT5 modeling.
# ProstT5 allows me to create a structural database from FASTA sequences.
# I will run Foldseek locally to search my custom database against all local Foldseek databases I downloaded.
# Now I want to split the results into upregulated, downregulated, and all significant proteins,
# based on their original accession groupings. I no longer filter only for "uncharacterized proteins"
# because I want to functionally annotate and interpret **all** differential expression hits.

import os  # I use this module to execute system commands and build file paths across my script.

# === STEP 1: Define my input and output structure ===
fasta_file = "ua5v7_DE_proteins.fasta"  # This is the FASTA file that I previously generated with DE protein sequences.
query_db = "ua5v7_structural_db"        # This will be the output database name when I run Foldseek's createdb with ProstT5.
results_dir = "prostt5_results"         # I will save all results to this directory so I can keep everything organized.
prost_model = "weights"                 # This is the folder where I downloaded the ProstT5 model.
database_dir = "/mnt/classes/biol_594_694/group6/databases"  # This is the path to all my usable Foldseek databases.

# === STEP 2: Load accession lists for category filtering ===
def load_accessions(path):
    # I open the file safely and load each line as a set element.
    # I use sets because they allow for fast lookup and remove duplicate entries.
    if os.path.exists(path):
        with open(path, 'r') as f:
            return set(line.strip() for line in f)  # I strip each line of whitespace and collect all in a set.
    return set()  # If the file doesn't exist, I return an empty set to avoid crashing the script.

# I now load my upregulated, downregulated, and all significant UniProt IDs from previous analyses.
upregulated = load_accessions("upregulated_accessions.txt")
downregulated = load_accessions("downregulated_accessions.txt")
all_significant = load_accessions("all_significant_accessions.txt")

# === STEP 3: Define temporary working directory for Foldseek ===
tmp_dir = os.path.join(results_dir, "tmp")  # I create a folder path for temporary Foldseek files.
os.makedirs(tmp_dir, exist_ok=True)  # I make sure this folder exists (creates it if missing).

# === STEP 4: Create ProstT5 database from FASTA ===
print("[Step 1] Creating ProstT5 structural database from FASTA...")
cmd_createdb = f"foldseek createdb {fasta_file} {query_db} --prostt5-model {prost_model}"  # I build the Foldseek command using ProstT5.
os.system(cmd_createdb)  # I run the command in the system shell.

# === STEP 5: Run Foldseek search against all local structural databases ===
print("[Step 2] Running Foldseek easy-search against all local databases...")

# I will now prepare four lists to hold hits categorized by significance class.
all_hits = []  # This will hold every result from all databases.
up_hits = []   # This will hold only hits where the query came from the upregulated list.
down_hits = [] # This will hold only hits where the query came from the downregulated list.
sig_hits = []  # This will hold only hits that were in the significant group.

# I go through each file in the database directory and filter out the valid database names by checking for the .dbtype extension.
for db_file in os.listdir(database_dir):
    if db_file.endswith(".dbtype") and not db_file.startswith(query_db):  # I ignore my own createdb file.
        db_name = db_file.replace(".dbtype", "")  # I remove the extension to get the usable Foldseek DB name.
        db_path = os.path.join(database_dir, db_name)  # I create the full path to the database.
        db_results_m8 = os.path.join(results_dir, f"results_{db_name}.m8")  # I define a name for the result file from this DB.
        print(f"\n▶ Searching against {db_name}...")  # I print which database I'm working on.

        # I now build and execute the Foldseek easy-search command.
        cmd_search = f"foldseek easy-search {query_db} {db_path} {db_results_m8} {tmp_dir}"
        os.system(cmd_search)  # I run the search using the shell.

        # === STEP 6: Parse and classify results from this search ===
        if os.path.exists(db_results_m8):  # I confirm the result file exists before parsing.
            with open(db_results_m8, 'r') as f:
                for line in f:
                    all_hits.append(line.strip())  # I collect the raw line in the overall hits list.
                    query_id = line.split('\t')[0].split('|')[1] if '|' in line.split('\t')[0] else line.split('\t')[0]  # I extract the UniProt ID.
                    if query_id in upregulated:
                        up_hits.append(line.strip())  # I check and store hits in the upregulated group.
                    if query_id in downregulated:
                        down_hits.append(line.strip())  # I do the same for downregulated hits.
                    if query_id in all_significant:
                        sig_hits.append(line.strip())  # I also collect all significant hits for separate analysis.

# === STEP 7: Save categorized results ===
def save_hits(hit_list, filename):
    # I open the desired output file and write all lines from the hit list to it.
    with open(os.path.join(results_dir, filename), 'w') as f:
        for hit in hit_list:
            f.write(hit + '\n')  # I write one hit per line.

# I now save each result set to its own file for downstream analysis.
save_hits(all_hits, "hits_all.txt")
save_hits(up_hits, "hits_upregulated.txt")
save_hits(down_hits, "hits_downregulated.txt")
save_hits(sig_hits, "hits_significant.txt")

# I finish by printing a simple summary showing how many hits were saved.
print(f"✅ Done. Saved {len(all_hits)} total hits.")
print(f"  - {len(up_hits)} upregulated")
print(f"  - {len(down_hits)} downregulated")
print(f"  - {len(sig_hits)} all significant")
