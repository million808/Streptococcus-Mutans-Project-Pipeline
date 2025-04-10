# README: Foldseek Differential Expression Structural Search Pipeline

### Author: Max Balter

---

## Overview
This project presents a fully automated structural bioinformatics workflow I developed to uncover the structural characteristics of differentially expressed (DE) proteins. In most differential expression studies, researchers stop at identifying which genes or proteins are significantly up- or downregulated. However, understanding **how** these proteins may function differently—especially in terms of their **3D structure**—requires a deeper layer of analysis.

To address that, I built a pipeline that starts with a spreadsheet of DE proteins and ends with structurally annotated, categorized, and downloadable 3D models. This is made possible by integrating:
- **ProstT5**, which predicts structures from FASTA sequences using transformer-based models.
- **Foldseek**, which performs rapid and sensitive structural similarity searches across enormous databases such as AlphaFold DB and ESMAtlas.

With this pipeline, I can examine not just sequence changes but also potential **functional impacts** based on structure, enabling further tasks such as active site mapping, druggability prediction, or structure-guided annotation using tools like ActSeek or CLEAN.

Overall, I used ProstT5 to predict the 3D structures of my differentially expressed proteins, starting from the FASTA sequences I retrieved from UniProt. Then, I used Foldseek to search those predicted structures against large databases like the AlphaFold Protein Structure Database. Foldseek returned a ranked list of structurally similar proteins, which were written into .m8 output files. These hits represent proteins with the closest 3D conformation to my queries, based on structural alignment scoring. Finally, my download.py script parsed these .m8 files, extracted the top hits by score, and downloaded their corresponding structure files directly from the AlphaFold repository. So ultimately, I was retrieving real AlphaFold-predicted structures that are statistically the most similar in shape to the proteins I found to be differentially expressed.
---

## Objective
The goal is to:
- Preprocess expression data and group DE proteins.
- Fetch full protein sequences via UniProt.
- Predict 3D structure with ProstT5.
- Search against Foldseek databases for structural matches.
- Group hits by expression category.
- Download structure files (PDB) for the top hits.

---

## Scripts and Their Roles

### **1. `preprocessing_FoldSeek_input.py`**
**Purpose:**
- Preprocess a spreadsheet of DE results (CSV/Excel).
- Extract and group UniProt accessions.

**What it does:**
- Parses the input spreadsheet.
- Separates accessions into:
  - `upregulated_accessions.txt`
  - `downregulated_accessions.txt`
  - `all_significant_accessions.txt`
- These groupings are used later to tag hits by DE category.

---

### **2. `download_fasta_from_uniprot.py`**
**Purpose:**
- Query UniProt’s API to retrieve the full protein sequences for each DE protein.

**What it does:**
- Reads accession IDs from `all_significant_accessions.txt`.
- Queries UniProt and retrieves sequences.
- Saves all results in a combined FASTA file: `ua5v7_DE_proteins.fasta`.

---

### **3. `fasta_to_foldseek.py`**
**Purpose:**
- Predict protein structures using ProstT5.
- Search predicted structures using Foldseek.
- Categorize and group results.

**What it does step-by-step:**
1. **Creates a structural DB** from FASTA using ProstT5: `ua5v7_structural_db/`
2. **Runs Foldseek easy-search** against all indexed databases in `databases/`
3. **Generates `.m8` result files** for each DB:
   - `results_afdb_proteome.m8`, `results_pdb.m8`, etc.
4. **Categorizes hits** by matching query accessions to:
   - `upregulated_accessions.txt`
   - `downregulated_accessions.txt`
   - `all_significant_accessions.txt`
5. **Writes grouped results** to:
   - `hits_upregulated.txt`
   - `hits_downregulated.txt`
   - `hits_significant.txt`
   - `hits_all.txt`

---

### **4. `download.py`**  
**Location:** `scripts_WIP/prostt5_results/`

**Purpose:**
- Download AlphaFold structure files for top Foldseek hits.

**What it does:**
- Opens a `.m8` result file (like `results_afdb_proteome.m8`).
- Sorts hits by bit score.
- Selects the top 10 hits.
- Downloads AlphaFold `.pdb` files by target UniProt ID.
- Saves structures in `top_hit_structures/`

---

## Required Files Before Running
- `ua5v7_DE_proteins.fasta`
- `weights/` folder containing ProstT5 GGUF model.
- `all_significant_accessions.txt`, `upregulated_accessions.txt`, etc.
- Foldseek databases (indexed `.dbtype`) in `databases/`

---

## Output Summary

| Output | Description |
|--------|-------------|
| `ua5v7_structural_db/` | ProstT5 structural prediction database |
| `prostt5_results/*.m8` | Foldseek search result files |
| `hits_*.txt` | Grouped DE protein hits by category |
| `top_hit_structures/*.pdb` | Downloaded AlphaFold PDB files from top matches |

---

## Interpretation of `.m8` Files
Each line in `.m8` files contains:
```
query_id   target_id   fident   alnlen   mismatch   gapopen   qstart   qend   tstart   tend   evalue   bitscore
```
This lets me:
- Identify which **DE protein** matched which structure.
- Assess **alignment quality** (higher bitscore = better).
- Sort and prioritize hits for functional follow-up.

---

## Running the Pipeline Step-by-Step
```bash
# Step 1: Group DE accessions
python preprocessing_FoldSeek_input.py

# Step 2: Download FASTA sequences from UniProt
python download_fasta_from_uniprot.py

# Step 3: Run ProstT5 and Foldseek
python fasta_to_foldseek.py

# Step 4: Download AlphaFold structures for top hits
cd scripts_WIP/prostt5_results/
python download.py
```

---

## Notes
- Run scripts in order.
- You can inspect `.m8` files manually or with pandas to explore matches.
- `top_hit_structures/*.pdb` are ready for visualization (PyMOL, ChimeraX) or active site matching (ActSeek).


---

## Acknowledgements
- **Foldseek**: van Kempen, M., Kim, S.S., Tumescheit, C., Mirdita, M., et al. (2023). Fast and accurate protein structure search with Foldseek. *Nature Biotechnology*. doi:10.1038/s41587-023-01773-0
- **ChatGPT**: Used for occasional debugging guidance and readability suggestions. All code logic, development, and scientific interpretation were conducted independently.
