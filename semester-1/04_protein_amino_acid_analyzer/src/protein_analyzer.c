/*
 * Protein Amino Acid Composition Analyzer
 * Course: M.Sc. Medical Bioinformatics - Semester 1
 * Skills: C Programming, File Handling, Biochemistry
 *
 * Description: Analyzes protein sequences from FASTA files
 * Calculates amino acid composition and protein properties
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE 1000
#define MAX_SEQ 50000
#define MAX_HEADER 500
#define NUM_AMINO_ACIDS 20

/* Amino acid information structure */
typedef struct {
    char one_letter;
    char three_letter[4];
    char name[15];
    double mol_weight;  /* Molecular weight in Daltons */
    char property;      /* H=hydrophobic, P=polar, C=charged, S=special */
} AminoAcidInfo;

/* Initialize amino acid data */
AminoAcidInfo amino_acids[NUM_AMINO_ACIDS] = {
    {'A', "Ala", "Alanine",       89.1,  'H'},
    {'C', "Cys", "Cysteine",      121.2, 'S'},
    {'D', "Asp", "Aspartate",     133.1, 'C'},
    {'E', "Glu", "Glutamate",     147.1, 'C'},
    {'F', "Phe", "Phenylalanine", 165.2, 'H'},
    {'G', "Gly", "Glycine",       75.1,  'S'},
    {'H', "His", "Histidine",     155.2, 'C'},
    {'I', "Ile", "Isoleucine",    131.2, 'H'},
    {'K', "Lys", "Lysine",        146.2, 'C'},
    {'L', "Leu", "Leucine",       131.2, 'H'},
    {'M', "Met", "Methionine",    149.2, 'H'},
    {'N', "Asn", "Asparagine",    132.1, 'P'},
    {'P', "Pro", "Proline",       115.1, 'H'},
    {'Q', "Gln", "Glutamine",     146.2, 'P'},
    {'R', "Arg", "Arginine",      174.2, 'C'},
    {'S', "Ser", "Serine",        105.1, 'P'},
    {'T', "Thr", "Threonine",     119.1, 'P'},
    {'V', "Val", "Valine",        117.1, 'H'},
    {'W', "Trp", "Tryptophan",    204.2, 'H'},
    {'Y', "Tyr", "Tyrosine",      181.2, 'S'}
};

/* Protein sequence structure */
typedef struct {
    char header[MAX_HEADER];
    char sequence[MAX_SEQ];
    int length;
    int counts[NUM_AMINO_ACIDS];
    int unknown_count;
} ProteinSequence;

/* Function to initialize protein structure */
void init_protein(ProteinSequence *prot) {
    prot->header[0] = '\0';
    prot->sequence[0] = '\0';
    prot->length = 0;
    prot->unknown_count = 0;
    for (int i = 0; i < NUM_AMINO_ACIDS; i++) {
        prot->counts[i] = 0;
    }
}

/* Function to find amino acid index */
int find_amino_acid_index(char aa) {
    for (int i = 0; i < NUM_AMINO_ACIDS; i++) {
        if (amino_acids[i].one_letter == aa) {
            return i;
        }
    }
    return -1;  /* Not found */
}

/* Function to read FASTA file */
int read_fasta(const char *filename, ProteinSequence *prot) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return 0;
    }

    char line[MAX_LINE];
    int seq_index = 0;

    while (fgets(line, MAX_LINE, file) != NULL) {
        line[strcspn(line, "\n")] = '\0';

        if (line[0] == '>') {
            strcpy(prot->header, line + 1);
        } else {
            for (int i = 0; line[i] != '\0'; i++) {
                char aa = toupper(line[i]);
                if (aa >= 'A' && aa <= 'Z') {
                    prot->sequence[seq_index++] = aa;
                }
            }
        }
    }
    prot->sequence[seq_index] = '\0';
    prot->length = seq_index;

    fclose(file);
    return 1;
}

/* Function to count amino acids */
void count_amino_acids(ProteinSequence *prot) {
    for (int i = 0; i < prot->length; i++) {
        int idx = find_amino_acid_index(prot->sequence[i]);
        if (idx >= 0) {
            prot->counts[idx]++;
        } else {
            prot->unknown_count++;
        }
    }
}

/* Function to calculate molecular weight */
double calculate_molecular_weight(ProteinSequence *prot) {
    double weight = 0.0;
    for (int i = 0; i < NUM_AMINO_ACIDS; i++) {
        weight += prot->counts[i] * amino_acids[i].mol_weight;
    }
    /* Subtract water molecules lost in peptide bond formation */
    weight -= (prot->length - 1) * 18.0;
    return weight;
}

/* Function to display results */
void display_results(ProteinSequence *prot) {
    printf("\n");
    printf("============================================================\n");
    printf("         PROTEIN AMINO ACID COMPOSITION ANALYSIS           \n");
    printf("============================================================\n\n");

    printf("Protein: %s\n", prot->header);
    printf("Length: %d amino acids\n\n", prot->length);

    printf("--- Amino Acid Composition ---\n\n");
    printf("  AA   Name           Count    Percentage\n");
    printf("  ---  -------------  -----    ----------\n");

    double total = (double)prot->length;

    for (int i = 0; i < NUM_AMINO_ACIDS; i++) {
        if (prot->counts[i] > 0) {
            printf("  %c    %-13s  %5d    %6.2f%%\n",
                   amino_acids[i].one_letter,
                   amino_acids[i].name,
                   prot->counts[i],
                   (prot->counts[i] / total) * 100);
        }
    }

    if (prot->unknown_count > 0) {
        printf("  X    Unknown        %5d    %6.2f%%\n",
               prot->unknown_count,
               (prot->unknown_count / total) * 100);
    }

    /* Calculate property-based composition */
    printf("\n--- Property-Based Composition ---\n\n");

    int hydrophobic = 0, polar = 0, charged = 0, special = 0;

    for (int i = 0; i < NUM_AMINO_ACIDS; i++) {
        switch (amino_acids[i].property) {
            case 'H': hydrophobic += prot->counts[i]; break;
            case 'P': polar += prot->counts[i]; break;
            case 'C': charged += prot->counts[i]; break;
            case 'S': special += prot->counts[i]; break;
        }
    }

    printf("  Hydrophobic (A,F,I,L,M,P,V,W): %5d  (%5.2f%%)\n",
           hydrophobic, (hydrophobic / total) * 100);
    printf("  Polar       (N,Q,S,T):         %5d  (%5.2f%%)\n",
           polar, (polar / total) * 100);
    printf("  Charged     (D,E,H,K,R):       %5d  (%5.2f%%)\n",
           charged, (charged / total) * 100);
    printf("  Special     (C,G,Y):           %5d  (%5.2f%%)\n",
           special, (special / total) * 100);

    /* Most and least abundant */
    printf("\n--- Abundance Analysis ---\n\n");

    int max_idx = 0, min_idx = -1;
    int max_count = 0, min_count = prot->length + 1;

    for (int i = 0; i < NUM_AMINO_ACIDS; i++) {
        if (prot->counts[i] > max_count) {
            max_count = prot->counts[i];
            max_idx = i;
        }
        if (prot->counts[i] > 0 && prot->counts[i] < min_count) {
            min_count = prot->counts[i];
            min_idx = i;
        }
    }

    printf("  Most abundant:  %s (%c) - %d (%.2f%%)\n",
           amino_acids[max_idx].name,
           amino_acids[max_idx].one_letter,
           max_count,
           (max_count / total) * 100);

    if (min_idx >= 0) {
        printf("  Least abundant: %s (%c) - %d (%.2f%%)\n",
               amino_acids[min_idx].name,
               amino_acids[min_idx].one_letter,
               min_count,
               (min_count / total) * 100);
    }

    /* Molecular weight */
    printf("\n--- Physical Properties ---\n\n");
    double mol_weight = calculate_molecular_weight(prot);
    printf("  Estimated Molecular Weight: %.2f Da (%.2f kDa)\n",
           mol_weight, mol_weight / 1000);

    printf("\n============================================================\n");
}

/* Function to save results */
void save_results(ProteinSequence *prot, const char *output_file) {
    FILE *file = fopen(output_file, "w");
    if (file == NULL) {
        printf("Error: Cannot create output file\n");
        return;
    }

    fprintf(file, "Protein Composition Analysis Results\n");
    fprintf(file, "=====================================\n\n");
    fprintf(file, "Protein: %s\n", prot->header);
    fprintf(file, "Length: %d amino acids\n\n", prot->length);

    fprintf(file, "Amino Acid Counts:\n");
    for (int i = 0; i < NUM_AMINO_ACIDS; i++) {
        if (prot->counts[i] > 0) {
            fprintf(file, "%c (%s): %d\n",
                    amino_acids[i].one_letter,
                    amino_acids[i].three_letter,
                    prot->counts[i]);
        }
    }

    double mol_weight = calculate_molecular_weight(prot);
    fprintf(file, "\nMolecular Weight: %.2f Da\n", mol_weight);

    fclose(file);
    printf("\nResults saved to: %s\n", output_file);
}

/* Main function */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <fasta_file>\n", argv[0]);
        printf("Example: %s ../data/insulin.fasta\n", argv[0]);
        return 1;
    }

    ProteinSequence prot;
    init_protein(&prot);

    /* Read FASTA file */
    if (!read_fasta(argv[1], &prot)) {
        return 1;
    }

    /* Count amino acids */
    count_amino_acids(&prot);

    /* Display results */
    display_results(&prot);

    /* Save results */
    save_results(&prot, "../output/composition_results.txt");

    return 0;
}
