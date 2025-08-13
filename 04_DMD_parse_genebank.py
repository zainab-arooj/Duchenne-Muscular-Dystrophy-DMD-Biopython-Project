from Bio import SeqIO

# Parse GenBank file
genbank_file = "DMD_sequence.gb"
record = SeqIO.read(genbank_file, "genbank")

# Initialize lists for features
gene_features = [f for f in record.features if f.type == "gene"]
cds_features = [f for f in record.features if f.type == "CDS"]
exon_features = [f for f in record.features if f.type == "exon"]
# Handle 'note' as a list for promoter filtering
promoter_features = [f for f in record.features if any(note.startswith("Dp427m") for note in f.qualifiers.get("note", []))]

# Print results
# Gene features
print("DMD Gene Features:")
for feature in gene_features:
    gene_name = feature.qualifiers.get("gene", ["N/A"])[0]
    synonyms = feature.qualifiers.get("gene_synonym", ["N/A"])
    db_xrefs = feature.qualifiers.get("db_xref", ["N/A"])
    print(f"  Gene: {gene_name}, Synonyms: {synonyms}, DB_Xrefs: {db_xrefs}")

# CDS features (Dp427m isoform)
print("\nCDS Features (Dp427m isoform):")
for feature in cds_features:
    if "Dp427m" in feature.qualifiers.get("product", [""])[0]:
        product = feature.qualifiers.get("product", ["N/A"])[0]
        translation = feature.qualifiers.get("translation", ["N/A"])[0][:50]
        print(f"  Product: {product}")
        print(f"  Translation (first 50 aa): {translation}")

# Exon count
print(f"\nTotal Exons: {len(exon_features)}")

# Promoter features (Dp427m)
print("\nPromoter Features (Dp427m):")
if promoter_features:
    for feature in promoter_features:
        note = ", ".join(feature.qualifiers.get("note", ["N/A"]))
        location = feature.location
        print(f"  Note: {note}, Location: {location}")
else:
    print("  No Dp427m promoter annotations found.")

# Sequence length
print(f"\nSequence Length: {len(record.seq)} bases")