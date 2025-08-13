from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt

# Load GenBank file
genbank_file = "DMD_sequence.gb"
record = SeqIO.read(genbank_file, "genbank")

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

# Simple secondary structure prediction (rule-based, e.g., based on amino acid propensities)
helix_aas = {'A', 'E', 'L', 'M', 'Q', 'R'}
sheet_aas = {'F', 'I', 'T', 'V', 'W', 'Y'}
structure_counts = {'helix': 0, 'sheet': 0, 'coil': 0}
for aa in dp427m_protein:
    if aa in helix_aas:
        structure_counts['helix'] += 1
    elif aa in sheet_aas:
        structure_counts['sheet'] += 1
    else:
        structure_counts['coil'] += 1

# Print results
print(f"Protein Length: {len(dp427m_protein)} amino acids")
print(f"Secondary Structure Counts: {structure_counts}")

# Plot secondary structure distribution
plt.figure(figsize=(8, 5))
plt.bar(structure_counts.keys(), structure_counts.values(), color="#1f77b4")
plt.xlabel("Secondary Structure")
plt.ylabel("Residue Count")
plt.title("Secondary Structure Distribution of Dp427m Protein (NG_012232.1)")
plt.show()