from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt

# Load GenBank file
genbank_file = "DMD_sequence.gb"
try:
    record = SeqIO.read(genbank_file, "genbank")
except FileNotFoundError as e:
    print(f"Error: {e}")
    exit()

# Extract Dp427m promoter region (1000 bp upstream of exon 1)
promoter_seq = None
tss = None
for feature in record.features:
    if feature.type == "exon" and feature.qualifiers.get("number", [""])[0] == "1":
        tss = feature.location.start
        promoter_seq = record.seq[max(0, tss - 1000):tss]
        break

if promoter_seq is None:
    print("Error: Exon 1 not found")
    exit()

# Calculate CpG ratios and GC content in sliding windows
window_size = 200
step_size = 50
positions = []
cpg_ratios = []
gc_contents = []
for i in range(0, len(promoter_seq) - window_size + 1, step_size):
    window = promoter_seq[i:i + window_size]
    c_count = window.count("C")
    g_count = window.count("G")
    cpg_count = str(window).count("CG")
    expected_cpg = (c_count * g_count) / window_size if c_count * g_count > 0 else 0
    cpg_ratio = cpg_count / expected_cpg if expected_cpg > 0 else 0
    cpg_ratios.append(cpg_ratio)
    gc_contents.append(gc_fraction(window) * 100)
    positions.append(i + window_size // 2 - 1000)  # Relative to TSS

# Identify CpG islands
cpg_islands = [(pos, ratio, gc) for pos, ratio, gc in zip(positions, cpg_ratios, gc_contents) if ratio > 0.6 and gc > 50]

# Print results
print(f"Promoter Region Length: {len(promoter_seq)} bases")
print(f"CpG Ratios (first 5): {cpg_ratios[:5]}")
print(f"GC Contents (first 5): {gc_contents[:5]}")
print(f"CpG Islands (position, ratio, GC%): {cpg_islands}")

# Plot CpG ratios and GC content
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax1.plot(positions, cpg_ratios, color="#1f77b4", label="CpG Ratio")
ax1.axhline(0.6, color="red", linestyle="--", label="CpG Island Threshold")
ax1.set_ylabel("Observed/Expected CpG Ratio")
ax1.set_title("CpG Island Prediction in Dp427m Promoter (NG_012232.1)")
ax1.legend()
ax2.plot(positions, gc_contents, color="#ff7f0e", label="GC Content (%)")
ax2.axhline(50, color="red", linestyle="--", label="GC Content Threshold")
ax2.set_xlabel("Position Relative to TSS (bp)")
ax2.set_ylabel("GC Content (%)")
ax2.set_title("GC Content in Dp427m Promoter")
ax2.legend()
plt.tight_layout()
plt.show()