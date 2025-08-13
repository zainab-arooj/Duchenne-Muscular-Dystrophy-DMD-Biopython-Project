from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt

# Load GenBank file
genbank_file = "DMD_sequence.gb"
record = SeqIO.read(genbank_file, "genbank")

# Extract Dp427m CDS
dp427m_cds = None
for feature in record.features:
    if feature.type == "CDS" and "Dp427m" in feature.qualifiers.get("product", [""])[0]:
        dp427m_cds = feature.extract(record.seq)
        break

if dp427m_cds is None:
    print("Error: Dp427m CDS not found")
    exit()

# Calculate 6-mer frequencies
k = 6
kmers = [str(dp427m_cds[i:i+k]) for i in range(len(dp427m_cds) - k + 1)]
kmer_counts = Counter(kmers)
kmer_frequencies = {kmer: count / len(kmers) * 100 for kmer, count in kmer_counts.items()}

# Get top 10 6-mers
top_kmers = sorted(kmer_frequencies.items(), key=lambda x: x[1], reverse=True)[:10]
top_kmer_names = [kmer for kmer, freq in top_kmers]
top_kmer_freqs = [freq for kmer, freq in top_kmers]

# Print results
print(f"Dp427m CDS Length: {len(dp427m_cds)} bases")
print(f"Number of 6-mers: {len(kmers)}")
print(f"Top 10 6-mers: {top_kmers}")

# Plot top 10 6-mers
plt.figure(figsize=(10, 5))
plt.bar(top_kmer_names, top_kmer_freqs, color="#1f77b4")
plt.xlabel("6-mer")
plt.ylabel("Frequency (%)")
plt.title("Top 10 6-mer Frequencies in Dp427m CDS (NG_012232.1)")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()