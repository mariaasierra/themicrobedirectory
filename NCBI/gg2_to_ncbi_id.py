import csv
from ete3 import NCBITaxa

ncbi = NCBITaxa()

# Input files
taxonomy_file = "2024.09.taxonomy.id.tsv"
expected_ids_file = "exp_ids.txt"  # One ID per line from tree$tip.label
output_file = "gg2_ncbi_mapped_v2.csv"
missing_ids_file = "missing_from_taxonomy.csv"

# Load expected tip labels
with open(expected_ids_file) as f:
    expected_ids = set(line.strip() for line in f if line.strip())

# Load taxonomy data
taxonomy_data = {}
with open(taxonomy_file, "r") as tsv:
    reader = csv.DictReader(tsv, delimiter="\t")
    for row in reader:
        taxonomy_data[row["Feature ID"].strip()] = row["Taxon"].strip()

# Set up cache and output
name_to_taxid = {}

def clean_taxonomy(taxon_string):
    return [part.split("__")[-1].strip() for part in taxon_string.split(";") if part]

def find_taxid(taxa_list):
    for name in reversed(taxa_list):
        if name in name_to_taxid:
            return name, name_to_taxid[name]
        try:
            result = ncbi.get_name_translator([name])
            if result and name in result:
                taxid = result[name][0]
                name_to_taxid[name] = taxid
                return name, taxid
        except Exception:
            continue
    return None, None

# Write output
with open(output_file, "w", newline="") as out, open(missing_ids_file, "w", newline="") as missing:
    writer = csv.writer(out)
    missing_writer = csv.writer(missing)

    writer.writerow(["Greengenes2_ID", "Taxon", "Matched_Taxon", "NCBI_TaxID", "Status"])
    missing_writer.writerow(["Greengenes2_ID", "Status"])

    for gg_id in expected_ids:
        taxon = taxonomy_data.get(gg_id)

        if taxon is None:
            missing_writer.writerow([gg_id, "Not found in taxonomy file"])
            continue

        if not taxon or taxon.lower() in {"", "na", "n/a", "none"}:
            writer.writerow([gg_id, taxon, "", "", "No taxonomy string"])
            continue

        taxa = clean_taxonomy(taxon)
        matched_name, taxid = find_taxid(taxa)

        if taxid:
            writer.writerow([gg_id, taxon, matched_name, taxid, "Matched"])
        else:
            writer.writerow([gg_id, taxon, "", "", "No match found"])

