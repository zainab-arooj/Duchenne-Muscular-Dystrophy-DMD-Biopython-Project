import cyvcf2
import matplotlib.pyplot as plt

# Load GenBank file to get CDS coordinates
genbank_file = "DMD_sequence.gb"
with open(genbank_file, "r") as f:
    for line in f:
        if line.startswith("LOCUS"):
            total_length = int(line.split()[1])
            break
cds_start = 3114019  # Approximate CDS start on chrX for NG_012232.1 (adjust based on your file)
cds_end = cds_start + 11058 - 1  # CDS length 11,058 bp

# Load VCF file
vcf_file = "dmd_variants.vcf"
vcf = cyvcf2.VCF(vcf_file)

# Filter variants within CDS
variant_counts = {"stop-gained": 0, "frameshift": 0, "missense": 0, "other": 0}
for variant in vcf:
    pos = variant.POS
    if cds_start <= pos <= cds_end and variant.ALT and variant.FILTER[0] == "PASS":
        info = dict(variant.INFO)
        clnvc = info.get("CLNVC", "other").lower()
        if "stop-gained" in clnvc:
            variant_counts["stop-gained"] += 1
        elif "frameshift" in clnvc:
            variant_counts["frameshift"] += 1
        elif "missense" in clnvc:
            variant_counts["missense"] += 1
        else:
            variant_counts["other"] += 1

# Print results
print(f"Dp427m CDS Length: 11058 bases")
print(f"Variant Counts: {variant_counts}")

# Plot variant counts
plt.figure(figsize=(8, 5))
plt.bar(variant_counts.keys(), variant_counts.values(), color="#1f77b4")
plt.ylabel("Number of Variants")
plt.title("Variant Types in Dp427m CDS (NG_012232.1)")
plt.tight_layout()
plt.show()