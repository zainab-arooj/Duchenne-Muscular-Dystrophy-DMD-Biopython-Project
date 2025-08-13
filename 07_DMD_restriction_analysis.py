from Bio import SeqIO
from Bio.Restriction import EcoRI, BamHI, Analysis
import matplotlib.pyplot as plt

# Load FASTA file
fasta_file = "DMD_sequence.fasta"
record = SeqIO.read(fasta_file, "fasta")
dmd_sequence = record.seq

# Load GenBank file for exon locations
genbank_file = "DMD_sequence.gb"
gb_record = SeqIO.read(genbank_file, "genbank")

# Restriction enzyme analysis
rb = Analysis([EcoRI, BamHI], dmd_sequence)
eco_sites = rb.with_sites()[EcoRI]
bam_sites = rb.with_sites()[BamHI]

# Map cut sites to exons
exon_features = [f for f in gb_record.features if f.type == "exon"]
eco_exon_counts = [0] * len(exon_features)
bam_exon_counts = [0] * len(exon_features)
for i, exon in enumerate(exon_features):
    start, end = exon.location.start, exon.location.end
    eco_exon_counts[i] = sum(start <= site <= end for site in eco_sites)
    bam_exon_counts[i] = sum(start <= site <= end for site in bam_sites)

# Print results
print(f"EcoRI Cut Sites: {len(eco_sites)}")
print(f"BamHI Cut Sites: {len(bam_sites)}")
print("\nEcoRI Cut Sites per Exon (non-zero only):")
for i, count in enumerate(eco_exon_counts, 1):
    if count > 0:
        print(f"Exon {i}: {count} sites")
print("\nBamHI Cut Sites per Exon (non-zero only):")
for i, count in enumerate(bam_exon_counts, 1):
    if count > 0:
        print(f"Exon {i}: {count} sites")

# Plot cut site distribution
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
exon_numbers = list(range(1, len(exon_features) + 1))
ax1.bar(exon_numbers, eco_exon_counts, color="#1f77b4", label="EcoRI")
ax1.set_ylabel("Cut Sites")
ax1.set_title("EcoRI Cut Sites per Exon (NG_012232.1)")
ax1.legend()
ax2.bar(exon_numbers, bam_exon_counts, color="#ff7f0e", label="BamHI")
ax2.set_xlabel("Exon Number")
ax2.set_ylabel("Cut Sites")
ax2.set_title("BamHI Cut Sites per Exon (NG_012232.1)")
ax2.legend()
plt.tight_layout()
plt.show()