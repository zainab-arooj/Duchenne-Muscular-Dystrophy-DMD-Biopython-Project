from Bio import SeqIO
from Bio.Seq import Seq

# Load GenBank file for CDS
genbank_file = "DMD_sequence.gb"
record = SeqIO.read(genbank_file, "genbank")
dmd_sequence = record.seq
print(f"DMD Gene Sequence (first 100 bases): {dmd_sequence[:100]}")

# Extract CDS sequence
cds_sequence = None
for feature in record.features:
    if feature.type == "CDS":
        cds_sequence = feature.extract(dmd_sequence)
        break
if cds_sequence:
    print(f"CDS Sequence (first 100 bases): {cds_sequence[:100]}")
    dmd_rna = cds_sequence.transcribe()
    print(f"DMD RNA Sequence (first 100 bases): {dmd_rna[:100]}")
    dmd_protein = dmd_rna.translate(to_stop=True)
    print(f"DMD Protein Sequence (partial): {dmd_protein[:50]}")
else:
    print("No CDS found in GenBank file")

# Sequence properties
print(f"Length of DMD sequence: {len(dmd_sequence)} bases")
print(f"GC content: {(dmd_sequence.count('G') + dmd_sequence.count('C')) / len(dmd_sequence) * 100:.2f}%")