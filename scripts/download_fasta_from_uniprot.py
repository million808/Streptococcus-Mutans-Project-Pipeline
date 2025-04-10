# Author: Max Balter
# Script: download_fasta_from_uniprot.py

# ======================= OVERVIEW ==========================
# This script is part of my data organization workflow.
# I use this to download FASTA sequences for all differentially expressed proteins
# identified in the UA 5v7 strain of Streptococcus mutans.
# The goal is to convert a list of UniProt accession IDs into a usable FASTA file
# that I can feed into AlphaFold, FoldSeek, or other structure/function prediction tools.

# What this script does:
# - Reads a .txt file containing UniProt accessions (one per line)
# - Queries the UniProt API for each sequence
# - Appends the results to a single output file: ua5v7_DE_proteins.fasta

import requests  # I use this to send GET requests to the UniProt API.
import time      # I use time.sleep to avoid overloading the API.

# === USER SETTINGS ===
input_file = "all_significant_accessions.txt"  # I can change this to upregulated or downregulated as needed.
output_fasta = "ua5v7_DE_proteins.fasta"       # This is the combined FASTA file I'll generate.

# === STEP 1: Load all UniProt accession IDs ===
# I open the input file and read in each accession ID.
with open(input_file, 'r') as f:
    accessions = [line.strip() for line in f if line.strip()]  # I strip blank lines just in case.

print(f"Loaded {len(accessions)} accessions from {input_file}")

# === STEP 2: Fetch FASTA from UniProt API ===
# I open my output file in write mode so I can add one FASTA at a time.
with open(output_fasta, 'w') as outfile:
    for i, acc in enumerate(accessions):
        # I construct the UniProt API URL to get FASTA format
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        response = requests.get(url)

        if response.status_code == 200:
            outfile.write(response.text)  # I write the full FASTA entry to the combined file.
            print(f"[{i+1}/{len(accessions)}] Downloaded: {acc}")
        else:
            print(f"[{i+1}/{len(accessions)}] Failed: {acc} (status {response.status_code})")

        time.sleep(0.5)  # I pause briefly to avoid hitting UniProt's rate limits.

print(f"Finished writing combined FASTA to {output_fasta}")

