# Protein Amino Acid Composition Analyzer

## Problem Statement

Proteins are made up of 20 different amino acids. The amino acid composition of a protein influences its structure, function, and properties. Analyzing composition helps in understanding protein characteristics and classification.

**Objective:** Build a C program that analyzes protein sequences to determine amino acid frequency and composition.

**Questions to Answer:**
1. What is the frequency of each amino acid?
2. Which amino acids are most/least abundant?
3. What is the proportion of hydrophobic vs hydrophilic amino acids?
4. What is the estimated molecular weight?

---

## Solution Approach

1. **File Reading:** Read protein sequence from FASTA file
2. **Counting:** Count each of the 20 amino acids
3. **Classification:** Group amino acids by properties:
   - Hydrophobic: A, V, L, I, M, F, W, P
   - Hydrophilic: S, T, N, Q
   - Charged: D, E, K, R, H
   - Special: C, G, Y
4. **Calculations:**
   - Frequency percentage
   - Molecular weight estimation
5. **Output:** Detailed composition report

---

## Key Concepts

- **20 Standard Amino Acids:** Each with unique properties
- **Hydrophobicity:** Important for protein folding
- **Molecular Weight:** Sum of amino acid weights (average ~110 Da per residue)

---

## Files

```
04_protein_amino_acid_analyzer/
├── README.md
├── src/
│   └── protein_analyzer.c
├── data/
│   ├── insulin.fasta
│   └── hemoglobin.fasta
└── output/
    └── (analysis results)
```

---

## How to Compile and Run

```bash
cd src
gcc -o protein_analyzer protein_analyzer.c
./protein_analyzer ../data/insulin.fasta
```

---

## Expected Output

```
=== PROTEIN COMPOSITION ANALYSIS ===
Protein: Human Insulin
Length: 110 amino acids

Amino Acid Frequency:
  Ala (A):  8 ( 7.27%)
  Cys (C): 12 (10.91%)
  ...

Properties:
  Hydrophobic: 45.5%
  Hydrophilic: 30.2%
  Charged: 24.3%

Estimated Molecular Weight: 12100 Da
```

---

## Skills Demonstrated

- C programming
- File handling and parsing
- Arrays and structures
- Biochemistry concepts (protein chemistry)
