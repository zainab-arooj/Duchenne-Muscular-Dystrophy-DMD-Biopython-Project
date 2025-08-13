from Bio import SeqIO

# Parse GenBank file
genbank_file = "DMD_sequence.gb"
genbank_record = SeqIO.read(genbank_file, "genbank")

# Extract gene feature (assuming first 'gene' feature)
gene_features = [f for f in genbank_record.features if f.type == "gene"]
if gene_features:
    gene = gene_features[0]
    print("DMD Gene Synonyms:", gene.qualifiers.get('gene_synonym', 'N/A'))
    print("Database Cross-References:", gene.qualifiers.get('db_xref', 'N/A'))

# Extract CDS feature and protein sequence
cds_features = [f for f in genbank_record.features if f.type == "CDS"]
if cds_features:
    cds = cds_features[0]
    print("Dystrophin Protein Sequence (first 50 aa):", cds.qualifiers.get('translation', ['N/A'])[0][:50])
    print("Product:", cds.qualifiers.get('product', 'N/A'))
    print("Note:", cds.qualifiers.get('note', 'N/A'))

# Count total exons (hard to count manually on NCBI)
exons = [f for f in genbank_record.features if f.type == "exon"]
print(f"Total Exons in DMD: {len(exons)}")