from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt

# Load GenBank file
genbank_file = "DMD_sequence.gb"
try:
    record = SeqIO.read(genbank_file, "genbank")
except FileNotFoundError as e:
    print(f"Error: {e}")
    exit()

# Extract Dp427m protein sequence
dp427m_cds = None
for feature in record.features:
    if feature.type == "CDS" and "Dp427m" in feature.qualifiers.get("product", [""])[0]:
        dp427m_cds = feature.extract(record.seq)
        break

if dp427m_cds is None:
    print("Error: Dp427m CDS not found")
    exit()

dp427m_protein = dp427m_cds.transcribe().translate(to_stop=True)

# Analyze protein properties
protein_analyzer = ProteinAnalysis(str(dp427m_protein))
aa_composition = protein_analyzer.count_amino_acids()
molecular_weight = protein_analyzer.molecular_weight()
isoelectric_point = protein_analyzer.isoelectric_point()
instability_index = protein_analyzer.instability_index()

# Print results
print(f"Protein Length: {len(dp427m_protein)} amino acids")
print(f"Molecular Weight: {molecular_weight:.2f} Da")
print(f"Isoelectric Point: {isoelectric_point:.2f}")
print(f"Instability Index: {instability_index:.2f}")
print(f"Amino Acid Composition: {aa_composition}")

# Plot amino acid composition
aa_names = list(aa_composition.keys())
aa_counts = list(aa_composition.values())
plt.figure(figsize=(10, 5))
plt.bar(aa_names, aa_counts, color="#1f77b4")
plt.xlabel("Amino Acid")
plt.ylabel("Count")
plt.title("Amino Acid Composition of Dp427m Protein (NG_012232.1)")
plt.show()