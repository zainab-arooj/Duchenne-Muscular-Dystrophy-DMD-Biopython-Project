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

# Define variants (position, ref_base, alt_base)
variants = [
    (2971, 'C', 'T'),  # c.2971C>T, p.Arg991Ter
    (3103, 'C', 'T')   # c.3103C>T, p.Arg1035Ter
]
wild_type_cds = str(dp427m_cds)
alignment_scores = []

# Apply variants and align
wild_type_protein = dp427m_cds.transcribe().translate(to_stop=True)
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.1

for pos, ref, alt in variants:
    if wild_type_cds[pos] != ref:
        print(f"Warning: Position {pos} is {wild_type_cds[pos]}, not {ref}. Variant may not be applicable.")
    mutated_cds = wild_type_cds[:pos] + alt + wild_type_cds[pos + 1:]
    mutated_cds_seq = Seq(mutated_cds)
    mutated_protein = mutated_cds_seq.transcribe().translate(to_stop=True)
    score = aligner.score(wild_type_protein, mutated_protein)
    alignment_scores.append(score)
    print(f"Variant at {pos} ({ref}>{alt}): Mutated Protein Length: {len(mutated_protein)} aa, Score: {score:.2f}")

# Calculate average score
avg_score = sum(alignment_scores) / len(alignment_scores) if alignment_scores else 0

# Print results
print(f"Dp427m CDS Length: {len(dp427m_cds)} bases")
print(f"Wild-type Protein Length: {len(wild_type_protein)} amino acids")
print(f"Average Alignment Score: {avg_score:.2f}")

# Plot alignment scores
plt.figure(figsize=(10, 5))
variant_labels = [f"Pos {pos} ({ref}>{alt})" for pos, ref, alt in variants]
plt.bar(variant_labels, alignment_scores, color="#1f77b4", label="Individual Scores")
plt.bar(len(variant_labels), avg_score, color="#ff7f0e", label="Average Score", width=0.5)
plt.ylabel("Alignment Score")
plt.title("Alignment Scores for Multiple Variants in Dp427m CDS (NG_012232.1)")
plt.legend()
plt.tight_layout()
plt.show()