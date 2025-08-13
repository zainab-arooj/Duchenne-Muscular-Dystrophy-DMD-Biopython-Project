from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt

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

# Apply ClinVar variant c.2971C>T (position 2971, 0-based)
variant_position = 2971
wild_type_cds = str(dp427m_cds)
if wild_type_cds[variant_position] != 'C':
    print(f"Warning: Position {variant_position} is {wild_type_cds[variant_position]}, not C. Variant may not be applicable.")
mutated_cds = wild_type_cds[:variant_position] + 'T' + wild_type_cds[variant_position + 1:]
mutated_cds_seq = Seq(mutated_cds)

# Translate wild-type and mutated CDS
wild_type_protein = dp427m_cds.transcribe().translate(to_stop=True)
mutated_protein = mutated_cds_seq.transcribe().translate(to_stop=True)

# Align proteins using PairwiseAligner
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.1
alignment_score = aligner.score(wild_type_protein, mutated_protein)
identity = (alignment_score / (2 * min(len(wild_type_protein), len(mutated_protein))) * 100)

# Print results
print(f"Dp427m CDS Length: {len(dp427m_cds)} bases")
print(f"Variant: c.2971C>T at position {variant_position} (codon {variant_position // 3 + 1})")
print(f"Wild-type Protein Length: {len(wild_type_protein)} amino acids")
print(f"Mutated Protein Length: {len(mutated_protein)} amino acids")
print(f"Alignment Score: {alignment_score:.2f}")
print(f"Sequence Identity (%): {identity:.2f}%")

# Plot alignment score
plt.figure(figsize=(8, 5))
plt.bar([1], [alignment_score], color="#1f77b4", label="Wild-type vs Mutated")
plt.ylabel("Alignment Score")
plt.title("Protein Alignment Score for c.2971C>T Variant (NG_012232.1)")
plt.xticks([1], ["Wild-type vs Mutated"])
plt.legend()
plt.show()