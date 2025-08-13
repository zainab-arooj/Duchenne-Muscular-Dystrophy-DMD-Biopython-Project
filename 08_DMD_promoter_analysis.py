from Bio import SeqIO
import matplotlib.pyplot as plt

# Load GenBank file
genbank_file = "DMD_sequence.gb"
record = SeqIO.read(genbank_file, "genbank")
dmd_sequence = record.seq

# Find Dp427m promoter region (1000 bp upstream of exon 1)
promoter_seq = None
tss = None
for feature in record.features:
    if feature.type == "exon" and feature.qualifiers.get("number", [""])[0] == "1":
        tss = feature.location.start
        promoter_seq = dmd_sequence[max(0, tss - 1000):tss]
        break

if promoter_seq is None:
    print("Exon 1 not found in GenBank file")
    exit()

# Motif analysis
tata_motif = "TATAAA"
caat_motif = "CCAAT"
tata_positions = [i - len(promoter_seq) for i in range(len(promoter_seq) - len(tata_motif) + 1) if str(promoter_seq[i:i + len(tata_motif)]) == tata_motif]
caat_positions = [i - len(promoter_seq) for i in range(len(promoter_seq) - len(caat_motif) + 1) if str(promoter_seq[i:i + len(caat_motif)]) == caat_motif]

# Print results
print(f"Promoter Region Length: {len(promoter_seq)} bases")
print(f"TATA Box (TATAAA) Count: {len(tata_positions)}")
print(f"TATA Box Positions (relative to TSS): {tata_positions}")
print(f"CAAT Box (CCAAT) Count: {len(caat_positions)}")
print(f"CAAT Box Positions (relative to TSS): {caat_positions}")

# Plot motif positions
plt.figure(figsize=(10, 5))
plt.scatter(tata_positions, [1] * len(tata_positions), color="#ff7f0e", label="TATA Box", marker="^")
plt.scatter(caat_positions, [0.5] * len(caat_positions), color="#2ca02c", label="CAAT Box", marker="o")
plt.axvline(0, color="black", linestyle="--", label="TSS")
plt.xlabel("Position Relative to TSS (bp)")
plt.ylabel("Motif Type")
plt.title("Promoter Motifs in Dp427m Upstream Region (NG_012232.1)")
plt.yticks([0.5, 1], ["CAAT", "TATA"])
plt.legend()
plt.show()