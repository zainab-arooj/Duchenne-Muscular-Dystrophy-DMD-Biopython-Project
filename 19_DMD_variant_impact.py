from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import cyvcf2
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

# Load VCF file
vcf_file = "dmd_variants.vcf"
vcf = cyvcf2.VCF(vcf_file)
cds_start = 3114019  # Approximate CDS start (adjust based on your file)
cds_end = cds_start + 11058 - 1

# Apply variants and align
wild_type_cds = str(dp427m_cds)
wild_type_protein = dp427m_cds.transcribe().translate(to_stop=True)
alignment_scores = []
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.1

for variant in vcf:
    pos = variant.POS - cds_start  # Convert to 0-based CDS position
    if 0 <= pos < len(wild_type_cds) and variant.ALT and variant.FILTER[0] == "PASS":
        ref = variant.REF
        alt = variant.ALT[0]
        if wild_type_cds[pos] != ref:
            print(f"Warning: Position {pos} is {wild_type_cds[pos]}, not {ref}. Variant may not be applicable.")
        mutated_cds = wild_type_cds[:pos] + alt + wild_type_cds[pos + 1:]
        mutated_cds_seq = Seq(mutated_cds)
        mutated_protein = mutated_cds_seq.transcribe().translate(to_stop=True)
        score = aligner.score(wild_type_protein, mutated_protein)
        alignment_scores.append(score)
        print(f"Variant at {variant.POS} ({ref}>{alt}): Mutated Protein Length: {len(mutated_protein)} aa, Score: {score:.2f}")

# Calculate average score
avg_score = sum(alignment_scores) / len(alignment_scores) if alignment_scores else 0

# Print results
print(f"Dp427m CDS Length: {len(dp427m_cds)} bases")
print(f"Wild-type Protein Length: {len(wild_type_protein)} amino acids")
print(f"Average Alignment Score: {avg_score:.2f}")

# Plot alignment scores
plt.figure(figsize=(10, 5))
variant_positions = [str(v.POS) for v in vcf if v.POS - cds_start >= 0 and v.POS - cds_start < len(wild_type_cds) and v.FILTER[0] == "PASS"]
plt.bar(variant_positions, alignment_scores, color="#1f77b4", label="Individual Scores")
plt.bar(len(variant_positions), avg_score, color="#ff7f0e", label="Average Score", width=0.5)
plt.ylabel("Alignment Score")
plt.title("Alignment Scores for Variants in Dp427m CDS (NG_012232.1)")
plt.legend()
plt.tight_layout()
plt.show()