import argparse
import pandas as pd
from ete3 import NCBITaxa
import warnings

# Suppress merged taxid warnings
warnings.filterwarnings("ignore", message="taxid .* was translated into .*")

# Set up CLI arguments
parser = argparse.ArgumentParser(description="Add or update taxonomic ranks using taxonomy IDs")
parser.add_argument("input_file", help="Input CSV file with a 'taxonomy_id' column")
parser.add_argument("--output", default="output_with_ranks.csv", help="Output CSV file")
args = parser.parse_args()

# Load input file
df = pd.read_csv(args.input_file)

# Ensure 'taxonomy_id' column exists
if 'taxonomy_id' not in df.columns:
    raise ValueError("Input file must have a 'taxonomy_id' column.")

# Convert taxonomy_id to numeric, keep NaNs for missing values
df['taxonomy_id'] = pd.to_numeric(df['taxonomy_id'], errors='coerce')

# Initialize taxonomy
ncbi = NCBITaxa()

# Prepare lists
ranks = []
unresolved = []

# Process each taxonomy_id
for idx, row in df.iterrows():
    rank = None
    tax_id = row['taxonomy_id']
    name = row['scientific_name']

    # Try taxonomy ID first
    if pd.notna(tax_id):
        try:
            lineage = ncbi.get_rank(ncbi.get_lineage(int(tax_id)))
            rank = lineage.get(int(tax_id), None)
        except Exception:
            pass

    # Fallback to scientific name
    if rank is None and name and name.lower() != 'nan':
        try:
            taxid_lookup = ncbi.get_name_translator([name])
            if name in taxid_lookup:
                fallback_id = taxid_lookup[name][0]
                lineage = ncbi.get_rank(ncbi.get_lineage(fallback_id))
                rank = lineage.get(fallback_id, None)
        except Exception:
            pass

    if rank is None:
        unresolved.append((tax_id, name))

    ranks.append(rank)


# Assign the rank column (always matches original row count)
df['rank'] = ranks

# Save output
df.to_csv(args.output, index=False)

# Optional: report missing tax IDs
if unresolved:
    print("\n⚠️ WARNING: These taxonomy IDs could not be resolved or have no rank:")
    for tid,name in sorted(set(unresolved)):
        print(f" - taxonomy_id: {tid}, scientific_name: {name}")
