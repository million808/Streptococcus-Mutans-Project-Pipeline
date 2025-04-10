# Author: Max Balter
# Debugging done by ChatGPT

# ======================= OVERVIEW ==========================
# In this script, I am preparing input files for downstream Foldseek analysis
# by extracting protein accessions from a differential expression dataset
# specific to the UA5v7 strain of Streptococcus mutans.

# What does "differentially expressed" mean?
# A protein is considered differentially expressed if its abundance significantly
# changes between two experimental conditions. In our case, we're comparing UA5v7
# grown at pH 7 versus pH 5. We use a p-value threshold (commonly < 0.05) to determine
# if this change is statistically significant.

# What is "upregulated"?
# A protein is upregulated if its expression level increases under the condition of interest.
# Here, this means a log2 fold change greater than 1 when comparing pH 7 to pH 5.

# What is "downregulated"?
# A protein is downregulated if its expression level decreases under the condition of interest.
# Here, this means a log2 fold change less than -1 when comparing pH 7 to pH 5.

# The script extracts three groups of proteins:
#   1. All significant proteins (p < 0.05 regardless of direction)
#   2. Significantly upregulated proteins (log2FC > 1, p < 0.05)
#   3. Significantly downregulated proteins (log2FC < -1, p < 0.05)
# It outputs three text files, each containing UniProt accession IDs for use with Foldseek.

import pandas as pd  # I import pandas to handle the spreadsheet and manipulate tabular data easily.

# === Step 1: Load the Excel file ===
# I define the path to the Excel file and the name of the sheet that contains the UA5v7-specific data.
excel_file = "/mnt/classes/biol_594_694/group6/datasets/diffExpData.xlsx"  # Full path to the uploaded Excel file.
sheet_name = "UA 5v7 (Main)"            # This is the worksheet focused specifically on UA5v7 strain comparisons.

# I use pandas to load the relevant sheet into a DataFrame so I can filter and analyze the data.
df = pd.read_excel(excel_file, sheet_name=sheet_name)

# === Step 2: Clean and filter the data ===
# Some rows may have missing data or headers — I need to filter them out before continuing.
# I keep only rows that contain:
# - A valid UniProt accession in the "Accession" column
# - A valid p-value in the "P value " column (note the space)
# - A valid log2 fold change in the "Ratio(log2)" column
df = df[df['Accession'].notna() & df['P value '].notna() & df['Ratio(log2)'].notna()]

# Just to be safe, I ensure the numeric columns are correctly typed as floats
# in case Excel formatted them as strings. This also catches any non-numeric values.
df['P value '] = pd.to_numeric(df['P value '], errors='coerce')
df['Ratio(log2)'] = pd.to_numeric(df['Ratio(log2)'], errors='coerce')

# === Step 3: Apply filtering criteria ===
# Now I define three groups based on differential expression logic.

# Group 1: All significantly differentially expressed proteins (p < 0.05)
# This includes both upregulated and downregulated proteins.
significant_df = df[df['P value '] < 0.05]
all_significant = significant_df['Accession'].unique()  # I extract just the unique accession IDs

# Group 2: Upregulated proteins — log2FC > 1 and significant
# These proteins increased in abundance at pH 7.
upregulated_df = df[(df['P value '] < 0.05) & (df['Ratio(log2)'] > 1)]
upregulated = upregulated_df['Accession'].unique()

# Group 3: Downregulated proteins — log2FC < -1 and significant
# These proteins decreased in abundance at pH 7.
downregulated_df = df[(df['P value '] < 0.05) & (df['Ratio(log2)'] < -1)]
downregulated = downregulated_df['Accession'].unique()

# === Step 4: Write each list to a .txt file ===
# I define a helper function to save each list of accessions to a file.
def save_list(accessions, filename):
    with open(filename, 'w') as f:
        for acc in accessions:
            f.write(str(acc).strip() + '\n')  # I write one accession per line
    print(f"Saved {len(accessions)} accessions to {filename}")

# I save all three groups to their respective output files
save_list(all_significant, "all_significant_accessions.txt")
save_list(upregulated, "upregulated_accessions.txt")
save_list(downregulated, "downregulated_accessions.txt")

# Final message to confirm everything worked
print("Completed Organization")
