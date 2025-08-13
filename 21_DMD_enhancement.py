import matplotlib.pyplot as plt

# Define timeline data
days = list(range(1, 21))
tasks = [
    "Sequence Extraction", "Annotation Parsing", "Sequence Analysis", "Restriction Sites",
    "CDS Analysis", "Motif Analysis", "Promoter Analysis", "CDS Alignment",
    "Protein Analysis", "Secondary Structure", "Codon Usage", "CpG Islands",
    "6-mer Analysis", "Variant Simulation", "Real Variant", "Multi-Variant",
    "VCF Parsing", "Variant Impact", "Statistical Analysis"
]
milestones = [1, 5, 9, 12, 14, 15, 17, 18, 19, 20]  # Key days

# Plot timeline
plt.figure(figsize=(12, 6))
plt.plot(days, [0] * len(days), "o-", color="#1f77b4", label="Tasks")
for day, task in zip(days[:len(tasks)], tasks):
    plt.text(day, 0.1, task, rotation=45, ha="right", va="bottom", fontsize=8)
for milestone in milestones:
    plt.plot(milestone, 0, "ro", label="Milestone" if milestone == milestones[0] else "")
plt.xlabel("Day")
plt.ylabel("Progress")
plt.title("Timeline of DMD Biopython Project (NG_012232.1)")
plt.yticks([])  # Hide y-axis
plt.legend()
plt.tight_layout()
plt.show()

# Print confirmation
print("Timeline plot generated.")