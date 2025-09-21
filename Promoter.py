from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import os

# === File paths ===
genome_file = "genome.fasta"
gff_file = "annotation.gff"
output_file = "promoters_300bp_fullID.fasta"

# === Parameters ===
upstream_length = 300   # 300 bp upstream
min_length = 200        # minimum length of promoter to keep

# === Load genome into dictionary ===
print("Loading genome...")
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
print(f"✅ Loaded {len(genome_dict)} chromosomes/scaffolds.")

promoters = []
seen_genes = set()  # track duplicate genes

# === Counters for summary ===
total_genes = 0
skipped_short = 0

# === Parse GFF manually ===
print("Parsing GFF and extracting promoters...")
with open(gff_file) as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue

        chrom, source, feature_type, start, end, score, strand, phase, attributes = parts

        if feature_type.lower() != "gene":
            continue

        total_genes += 1
        start = int(start) - 1  # convert to 0-based
        end = int(end)

        # Extract gene ID from attributes
        gene_id_match = re.search(r"ID=([^;]+)", attributes)
        gene_id = gene_id_match.group(1) if gene_id_match else f"{chrom}_{start}_{end}"

        # Skip duplicate genes
        if gene_id in seen_genes:
            continue
        seen_genes.add(gene_id)

        # Extract promoter sequence
        if chrom not in genome_dict:
            print(f"⚠️ Warning: Chromosome {chrom} not found in genome. Skipping gene {gene_id}.")
            continue

        seq = genome_dict[chrom].seq
        if strand == "+":
            promoter_start = max(0, start - upstream_length)
            promoter_end = start
            promoter_seq = seq[promoter_start:promoter_end]
        else:
            promoter_start = end
            promoter_end = min(len(seq), end + upstream_length)
            promoter_seq = seq[promoter_start:promoter_end].reverse_complement()

        # Skip sequences shorter than min_length
        if len(promoter_seq) < min_length:
            skipped_short += 1
            continue

        # Add to list
        promoters.append(SeqRecord(
            promoter_seq,
            id=gene_id,
            description=f"{chrom}:{promoter_start}-{promoter_end} strand={strand}"
        ))

# === Write output FASTA ===
if promoters:
    SeqIO.write(promoters, output_file, "fasta")
    print(f"✅ Written {len(promoters)} promoters to {output_file}")
else:
    print("⚠️ No promoters to write. Check genome and GFF compatibility.")

# === Summary ===
print(f"Total genes scanned: {total_genes}")
print(f"Promoters extracted (≥ {min_length} bp): {len(promoters)}")
print(f"Sequences skipped for being too short: {skipped_short}")
