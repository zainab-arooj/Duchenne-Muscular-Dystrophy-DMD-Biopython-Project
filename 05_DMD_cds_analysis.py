from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt

# Load GenBank file
genbank_file = "DMD_sequence.gb"
record = SeqIO.read(genbank_file, "genbank")
dmd_sequence = record.seq

# Extract CDS for Dp427m
cds_feature = None
cds_sequence = None
for feature in record.features:
    if feature.type == "CDS" and "Dp427m" in feature.qualifiers.get("product", [""])[0]:
        cds_feature = feature
        cds_sequence = feature.extract(dmd_sequence)
        break

# Process CDS
if cds_sequence:
    print(f"CDS Sequence (first 100 bases): {cds_sequence[:100]}")
    print(f"CDS Length: {len(cds_sequence)} bases")
    dmd_rna = cds_sequence.transcribe()
    print(f"DMD RNA Sequence (first 100 bases): {dmd_rna[:100]}")
    dmd_protein = dmd_rna.translate(to_stop=True)
    print(f"DMD Protein Sequence (first 50 aa): {dmd_protein[:50]}")
    print(f"Protein Length: {len(dmd_protein)} amino acids")
else:
    print("No Dp427m CDS found in GenBank file")
    exit()

# Map CDS to exons
exon_features = [f for f in record.features if f.type == "exon"]
cds_exon_lengths = []
for exon in exon_features:
    # Initialize length as 0 (non-coding exon)
    exon_length = 0
    # Check if exon overlaps with CDS location
    if cds_feature and exon.location.start < cds_feature.location.end and exon.location.end > cds_feature.location.start:
        # Calculate overlap with CDS
        start = max(exon.location.start, cds_feature.location.start)
        end = min(exon.location.end, cds_feature.location.end)
        if end > start:
            exon_length = end - start
    cds_exon_lengths.append(exon_length)

# Filter exons contributing to CDS
cds_exon_lengths = [length for length in cds_exon_lengths if length > 0]
exon_numbers = list(range(1, len(cds_exon_lengths) + 1))

# Print exon-wise CDS lengths
print("\nExon-wise CDS Lengths (non-zero only):")
for i, length in enumerate(cds_exon_lengths, 1):
    if length > 0:
        print(f"Exon {i}: {length} bases")

# Create bar plot
plt.figure(figsize=(10, 5))
plt.bar(exon_numbers, cds_exon_lengths, color="#1f77b4")
plt.xlabel("Exon Number")
plt.ylabel("CDS Length (bases)")
plt.title("Dp427m CDS Length per Exon (NG_012232.1)")
plt.show()