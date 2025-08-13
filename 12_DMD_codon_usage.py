from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from collections import Counter

# Load GenBank file
genbank_file = "DMD_sequence.gb"
try:
    record = SeqIO.read(genbank_file, "genbank")
except FileNotFoundError as e:
    print(f"Error: {e}")
    exit()

# Extract Dp427m CDS
dp427m_cds = None
for feature in record.features:
    if feature.type == "CDS" and "Dp427m" in feature.qualifiers.get("product", [""])[0]:
        dp427m_cds = feature.extract(record.seq)
        break

if dp427m_cds is None:
    print("Error: Dp427m CDS not found")
    exit()

# Ensure CDS length is divisible by 3
if len(dp427m_cds) % 3 != 0:
    print(f"Warning: CDS length {len(dp427m_cds)} is not divisible by 3")
    dp427m_cds = dp427m_cds[:-(len(dp427m_cds) % 3)]

# Count codons
codons = [str(dp427m_cds[i:i+3]) for i in range(0, len(dp427m_cds), 3)]
codon_counts = Counter(codons)
codon_frequencies = {codon: count / len(codons) * 100 for codon, count in codon_counts.items()}

# Print results
print(f"Dp427m CDS Length: {len(dp427m_cds)} bases")
print(f"Number of Codons: {len(codons)}")
print(f"Codon Frequencies (%): {codon_frequencies}")

# Plot codon frequencies
codon_names = sorted(codon_frequencies.keys())
codon_values = [codon_frequencies[codon] for codon in codon_names]
plt.figure(figsize=(15, 5))
plt.bar(codon_names, codon_values, color="#1f77b4")
plt.xlabel("Codon")
plt.ylabel("Frequency (%)")
plt.title("Codon Usage in Dp427m CDS (NG_012232.1)")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()