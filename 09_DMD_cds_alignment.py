from Bio import SeqIO, pairwise2
import matplotlib.pyplot as plt

# Load GenBank files
dp427m_file = "DMD_sequence.gb"
dp427c_file = "DMD_Dp427c.gb"
try:
    dp427m_record = SeqIO.read(dp427m_file, "genbank")
    dp427c_record = SeqIO.read(dp427c_file, "genbank")
except FileNotFoundError as e:
    print(f"Error: {e}")
    exit()

# Extract CDS sequences
dp427m_cds = None
dp427c_cds = None
for feature in dp427m_record.features:
    if feature.type == "CDS" and "Dp427m" in feature.qualifiers.get("product", [""])[0]:
        dp427m_cds = feature.extract(dp427m_record.seq)
        break
for feature in dp427c_record.features:
    if feature.type == "CDS" and "Dp427c" in feature.qualifiers.get("product", [""])[0]:
        dp427c_cds = feature.extract(dp427c_record.seq)
        break

if dp427m_cds is None or dp427c_cds is None:
    print("Error: Could not find Dp427m or Dp427c CDS in GenBank files")
    exit()

# Perform pairwise alignment (global, match=2, mismatch=-1, gap open=-0.5, gap extend=-0.1)
alignments = pairwise2.align.globalms(dp427m_cds, dp427c_cds, 2, -1, -0.5, -0.1, score_only=True)

# Print results
print(f"Dp427m CDS Length: {len(dp427m_cds)} bases")
print(f"Dp427c CDS Length: {len(dp427c_cds)} bases")
print(f"Alignment Score: {alignments:.2f}")
print(f"Sequence Identity (%): {(alignments / (2 * min(len(dp427m_cds), len(dp427c_cds))) * 100):.2f}%")

# Plot alignment score
plt.figure(figsize=(8, 5))
plt.bar([1], [alignments], color="#1f77b4", label="Dp427m vs Dp427c")
plt.ylabel("Alignment Score")
plt.title("CDS Alignment Score (NG_012232.1 vs NM_000109.4)")
plt.xticks([1], ["Dp427m vs Dp427c"])
plt.legend()
plt.show()