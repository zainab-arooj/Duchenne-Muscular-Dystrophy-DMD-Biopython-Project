from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Restriction import EcoRI, BamHI, Analysis

# Load DMD sequence from FASTA
fasta_file = "DMD_sequence.fasta"
record = SeqIO.read(fasta_file, "fasta")
dmd_sequence = record.seq

# Calculate GC content
gc_content = gc_fraction(dmd_sequence) * 100
print(f"GC Content of DMD Gene: {gc_content:.2f}%")

# Restriction enzyme analysis (EcoRI: GAATTC, BamHI: GGATCC)
rb = Analysis([EcoRI, BamHI], dmd_sequence)
eco_sites = rb.with_sites()[EcoRI]
bam_sites = rb.with_sites()[BamHI]
print(f"EcoRI cut sites: {len(eco_sites)}")
print(f"BamHI cut sites: {len(bam_sites)}")

# Simple motif search (e.g., TATA box, a promoter motif: "TATAAA")
motif = "TATAAA"
motif_count = str(dmd_sequence).count(motif)
print(f"Occurrences of TATA box motif (TATAAA): {motif_count}")

# Optional: Print positions of first few TATA motifs
positions = [i for i in range(len(dmd_sequence)) if str(dmd_sequence[i:i+6]) == motif]
print(f"First 5 TATA box positions: {positions[:5]}")