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

# Initialize metrics
wild_type_cds = str(dp427m_cds)
wild_type_protein = dp427m_cds.transcribe().translate(to_stop=True)
wild_type_score = len(wild_type_protein) * 2  # Max score for perfect alignment
variant_counts = {"stop-gained": 0, "frameshift": 0, "missense": 0}
alignment_scores = []
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.1

# Process variants
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
        info = dict(variant.INFO)
        clnvc = info.get("CLNVC", "other").lower()
        if "stop-gained" in clnvc:
            variant_counts["stop-gained"] += 1
        elif "frameshift" in clnvc:
            variant_counts["frameshift"] += 1
        elif "missense" in clnvc:
            variant_counts["missense"] += 1

# Calculate metrics
total_variants = sum(variant_counts.values())
frequencies = {k: v / total_variants * 100 if total_variants > 0 else 0 for k, v in variant_counts.items()}
score_reductions = [wild_type_score - s for s in alignment_scores]
avg_score_reduction = sum(score_reductions) / len(score_reductions) if score_reductions else 0
severity = {"high": variant_counts["stop-gained"] + variant_counts["frameshift"], "low": variant_counts["missense"]}

# Print results
print(f"Dp427m CDS Length: {len(dp427m_cds)} bases")
print(f"Wild-type Protein Length: {len(wild_type_protein)} amino acids")
print(f"Variant Counts: {variant_counts}")
print(f"Variant Frequencies (%): {frequencies}")
print(f"Average Alignment Score Reduction: {avg_score_reduction:.2f}")
print(f"Severity Distribution: {severity}")

# Plot metrics
metrics = ["Frequency (%)"] + list(frequencies.keys())
values = [frequencies[k] for k in frequencies] + [avg_score_reduction / 73.7]  # Normalize score reduction to % scale
plt.figure(figsize=(10, 5))
plt.bar(metrics, values, color="#1f77b4")
plt.ylabel("Value (%)")
plt.title("Statistical Analysis of Variants in Dp427m CDS (NG_012232.1)")
plt.tight_layout()
plt.show()