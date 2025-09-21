from Bio import SeqIO
import re
import pandas as pd

# === Files ===
input_fasta = "promoters_300bp_fullID.fasta"
output_excel = "promoters_with_YACCWACY.xlsx"

# === Motif (IUPAC) ===
iupac_dict = {
    "A":"A","C":"C","G":"G","T":"T",
    "R":"[AG]","Y":"[CT]","S":"[GC]","W":"[AT]",
    "K":"[GT]","M":"[AC]","B":"[CGT]","D":"[AGT]",
    "H":"[ACT]","V":"[ACG]","N":"[ACGT]"
}

motif = "YACCWACY"
motif_regex = "".join(iupac_dict[nuc] for nuc in motif)

# === Function to find motif matches ===
def find_motif(seq, pattern):
    return [m.start() for m in re.finditer(pattern, str(seq), flags=re.IGNORECASE)]

# === Scan sequences ===
data = []

for record in SeqIO.parse(input_fasta, "fasta"):
    seq = record.seq
    total_hits = 0

    # Sense strand
    total_hits += len(find_motif(seq, motif_regex))
    # Antisense strand
    total_hits += len(find_motif(seq.reverse_complement(), motif_regex))

    if total_hits > 0:
        data.append({
            "Sequence_ID": record.id,
            "Description": record.description,
            "Motif_Hits": total_hits,
            "Sequence": str(seq)
        })

# === Write to Excel ===
df = pd.DataFrame(data)
df.to_excel(output_excel, index=False)

print(f"✅ Total sequences scanned: {len(list(SeqIO.parse(input_fasta, 'fasta')))}")
print(f"✅ Sequences containing motif '{motif}': {len(df)}")
print(f"✅ Excel output saved to: {output_excel}")
