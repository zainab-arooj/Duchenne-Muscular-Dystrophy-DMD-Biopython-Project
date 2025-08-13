import cyvcf2
import matplotlib.pyplot as plt

# Load annotated VCF file
vcf_file = "dmd_variants_annotated.vcf"
vcf = cyvcf2.VCF(vcf_file)

# Extract effect types
effect_counts = {"HIGH": 0, "MODERATE": 0, "LOW": 0, "MODIFIER": 0}
for variant in vcf:
    if variant.INFO.get("ANN"):
        ann = variant.INFO["ANN"][0].split("|")
        effect = ann[1]  # Effect impact (e.g., HIGH, MODERATE)
        if effect in effect_counts:
            effect_counts[effect] += 1

# Print results
print(f"Variant Effect Counts: {effect_counts}")

# Plot effect distribution
plt.figure(figsize=(8, 5))
plt.bar(effect_counts.keys(), effect_counts.values(), color="#1f77b4")
plt.ylabel("Number of Variants")
plt.title("Distribution of Variant Effects in Dp427m CDS (NG_012232.1)")
plt.tight_layout()
plt.show()