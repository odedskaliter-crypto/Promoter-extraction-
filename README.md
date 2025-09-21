# Petunia Promoter Extraction and Motif Analysis Scripts

This repository contains Python scripts used for analyzing **promoter sequences of genes** and scanning for a specific motif (`YACCWACY`) in petunia, as part of the manuscript by Oded Skaliter, Alexander Vainstein, and collaborators.

---

## **Scripts Overview**

### 1. `extract_promoters.py`
- **Purpose:** Extracts promoter sequences upstream of genes from a genome FASTA and a GFF annotation file.
- **Input:**
  - `genome.fasta` — genome sequences (downloaded from **Sol Genomics Network**, [https://solgenomics.net](https://solgenomics.net))
  - `annotation.gff` — gene annotation file (downloaded from **Sol Genomics Network**)
- **Output:**
  - `promoters_300bp_fullID.fasta` — promoter sequences of at least 200 bp
- **Parameters you can adjust:**
  - `upstream_length` — length upstream of gene start (default: 300 bp)
  - `min_length` — minimum promoter length to keep (default: 200 bp)

---

### 2. `scan_motif.py`
- **Purpose:** Scans promoter sequences for a specific motif (`YACCWACY`) on both sense and antisense strands.
- **Input:**
  - `promoters_300bp_fullID.fasta` — promoter sequences from previous script
- **Output:**
  - `promoters_with_YACCWACY.xlsx` — Excel file listing sequences with motif hits, number of hits, and sequence details

---

## **Requirements**

- Python 3.10 or higher  
- Python packages:
  - `biopython`
  - `pandas`
  - `openpyxl` (for writing Excel files)

Install dependencies with:

```bash
pip install biopython pandas openpyxl

Usage

Extract promoters:

python extract_promoters.py


Scan for motif:

python scan_motif.py

File Structure
EOBV_promoter_analysis/
│
├── extract_promoters.py           # Script to extract promoters from genome/GFF
├── scan_motif.py                 # Script to scan promoters for motif
├── genome.fasta                  # Genome file (input, from Sol Genomics)
├── annotation.gff                # GFF annotation file (input, from Sol Genomics)
├── promoters_300bp_fullID.fasta  # Output of promoter extraction
├── promoters_with_YACCWACY.xlsx  # Output of motif scanning
└── README.md                     # This file

Authors

Oded Skaliter

License

This code is provided for academic use. Please cite the corresponding manuscript if you use these scripts.
