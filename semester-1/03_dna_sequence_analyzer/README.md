# DNA Sequence Analyzer

## Problem Statement

DNA sequences consist of four nucleotides: Adenine (A), Thymine (T), Guanine (G), and Cytosine (C). Analyzing the composition of DNA sequences provides insights into gene structure, GC content (important for PCR primer design), and sequence characteristics.

**Objective:** Build a C program that reads DNA sequences from FASTA files and computes sequence statistics.

**Questions to Answer:**
1. What is the nucleotide composition (count of A, T, G, C)?
2. What is the GC content percentage?
3. What is the total sequence length?
4. What is the AT/GC ratio?

---

## Solution Approach

1. **File Reading:** Read FASTA format files
2. **Sequence Parsing:** Extract sequence (skip header lines starting with ">")
3. **Counting:** Count each nucleotide (A, T, G, C)
4. **Calculations:**
   - GC content = (G + C) / Total × 100
   - AT content = (A + T) / Total × 100
   - AT/GC ratio
5. **Output:** Display statistics summary

---

## Key Concepts

- **FASTA Format:** Standard format for DNA/protein sequences. Header starts with ">", followed by sequence.
- **GC Content:** Percentage of G and C bases. Important for melting temperature, genome stability.
- **Chargaff's Rules:** A pairs with T, G pairs with C (in double-stranded DNA)

---

## Files

```
03_dna_sequence_analyzer/
├── README.md
├── src/
│   └── dna_analyzer.c
├── data/
│   ├── sample1.fasta
│   └── sample2.fasta
└── output/
    └── (analysis results)
```

---

## How to Compile and Run

```bash
cd src
gcc -o dna_analyzer dna_analyzer.c
./dna_analyzer ../data/sample1.fasta
```

---

## Expected Output

```
=== DNA SEQUENCE ANALYSIS ===
Sequence: sample_gene
Length: 150 bp

Nucleotide Counts:
  A: 40 (26.67%)
  T: 38 (25.33%)
  G: 35 (23.33%)
  C: 37 (24.67%)

GC Content: 48.00%
AT Content: 52.00%
AT/GC Ratio: 1.08
```

---

## Skills Demonstrated

- C programming
- File handling
- String manipulation
- Cell and Molecular Biology concepts
