from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt

# Load FASTA file
fasta_file = "DMD_sequence.fasta"
record = SeqIO.read(fasta_file, "fasta")
dmd_sequence = record.seq

# Sliding window GC content
window_size = 1000
step_size = 500
gc_contents = []
positions = []
for i in range(0, len(dmd_sequence) - window_size + 1, step_size):
    window = dmd_sequence[i:i + window_size]
    gc_contents.append(gc_fraction(window) * 100)
    positions.append(i + window_size // 2)

# Motif analysis
tata_motif = "TATAAA"
cpg_motif = "CG"
tata_count = str(dmd_sequence).count(tata_motif)
cpg_count = str(dmd_sequence).count(cpg_motif)
tata_positions = [i for i in range(len(dmd_sequence) - len(tata_motif) + 1) if str(dmd_sequence[i:i + len(tata_motif)]) == tata_motif]
cpg_positions = [i for i in range(len(dmd_sequence) - len(cpg_motif) + 1) if str(dmd_sequence[i:i + len(cpg_motif)]) == cpg_motif]

# Calculate motif density (motifs per 1000 bp)
tata_density = [0] * len(gc_contents)
cpg_density = [0] * len(gc_contents)
for pos in tata_positions:
    window_idx = min(pos // step_size, len(gc_contents) - 1)
    tata_density[window_idx] += 1
for pos in cpg_positions:
    window_idx = min(pos // step_size, len(gc_contents) - 1)
    cpg_density[window_idx] += 1

# Print results
print(f"Sequence Length: {len(dmd_sequence)} bases")
print(f"Overall GC Content: {gc_fraction(dmd_sequence) * 100:.2f}%")
print(f"TATA Box (TATAAA) Count: {tata_count}")
print(f"First 5 TATA Box Positions: {tata_positions[:5]}")
print(f"CpG Dinucleotide Count: {cpg_count}")
print(f"GC Content Windows (first 5): {gc_contents[:5]}")

# Plot GC content and motif density
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax1.plot(positions, gc_contents, color="#1f77b4")
ax1.set_ylabel("GC Content (%)")
ax1.set_title("GC Content Variation in NG_012232.1 (1000 bp windows)")
ax2.bar(positions, tata_density, width=step_size / 2, color="#ff7f0e", alpha=0.5, label="TATA Box")
ax2.bar([p + step_size / 2 for p in positions], cpg_density, width=step_size / 2, color="#2ca02c", alpha=0.5, label="CpG")
ax2.set_xlabel("Position (bp)")
ax2.set_ylabel("Motif Count per Window")
ax2.set_title("Motif Density in NG_012232.1")
ax2.legend()
plt.tight_layout()
plt.show()