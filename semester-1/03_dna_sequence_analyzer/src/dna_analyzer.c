/*
 * DNA Sequence Analyzer
 * Course: M.Sc. Medical Bioinformatics - Semester 1
 * Skills: C Programming, File Handling, Cell and Molecular Biology
 *
 * Description: Analyzes DNA sequences from FASTA files
 * Calculates nucleotide composition, GC content, and sequence statistics
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE 1000
#define MAX_SEQ 100000
#define MAX_HEADER 500

/* Structure to hold sequence data */
typedef struct {
    char header[MAX_HEADER];
    char sequence[MAX_SEQ];
    int length;
    int count_A;
    int count_T;
    int count_G;
    int count_C;
    int count_N;  /* Unknown nucleotides */
} DNASequence;

/* Function to initialize sequence structure */
void init_sequence(DNASequence *seq) {
    seq->header[0] = '\0';
    seq->sequence[0] = '\0';
    seq->length = 0;
    seq->count_A = 0;
    seq->count_T = 0;
    seq->count_G = 0;
    seq->count_C = 0;
    seq->count_N = 0;
}

/* Function to read FASTA file */
int read_fasta(const char *filename, DNASequence *seq) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return 0;
    }

    char line[MAX_LINE];
    int seq_index = 0;

    while (fgets(line, MAX_LINE, file) != NULL) {
        /* Remove newline character */
        line[strcspn(line, "\n")] = '\0';

        if (line[0] == '>') {
            /* Header line */
            strcpy(seq->header, line + 1);  /* Skip '>' */
        } else {
            /* Sequence line */
            for (int i = 0; line[i] != '\0'; i++) {
                char nucleotide = toupper(line[i]);
                if (nucleotide == 'A' || nucleotide == 'T' ||
                    nucleotide == 'G' || nucleotide == 'C' ||
                    nucleotide == 'N') {
                    seq->sequence[seq_index++] = nucleotide;
                }
            }
        }
    }
    seq->sequence[seq_index] = '\0';
    seq->length = seq_index;

    fclose(file);
    return 1;
}

/* Function to count nucleotides */
void count_nucleotides(DNASequence *seq) {
    for (int i = 0; i < seq->length; i++) {
        switch (seq->sequence[i]) {
            case 'A': seq->count_A++; break;
            case 'T': seq->count_T++; break;
            case 'G': seq->count_G++; break;
            case 'C': seq->count_C++; break;
            default:  seq->count_N++; break;
        }
    }
}

/* Function to calculate and display statistics */
void display_statistics(DNASequence *seq) {
    printf("\n");
    printf("============================================\n");
    printf("       DNA SEQUENCE ANALYSIS REPORT        \n");
    printf("============================================\n\n");

    printf("Sequence Name: %s\n", seq->header);
    printf("Total Length: %d bp\n", seq->length);

    printf("\n--- Nucleotide Composition ---\n\n");

    double total = (double)seq->length;
    printf("  Adenine  (A): %5d  (%6.2f%%)\n",
           seq->count_A, (seq->count_A / total) * 100);
    printf("  Thymine  (T): %5d  (%6.2f%%)\n",
           seq->count_T, (seq->count_T / total) * 100);
    printf("  Guanine  (G): %5d  (%6.2f%%)\n",
           seq->count_G, (seq->count_G / total) * 100);
    printf("  Cytosine (C): %5d  (%6.2f%%)\n",
           seq->count_C, (seq->count_C / total) * 100);

    if (seq->count_N > 0) {
        printf("  Unknown  (N): %5d  (%6.2f%%)\n",
               seq->count_N, (seq->count_N / total) * 100);
    }

    printf("\n--- Derived Statistics ---\n\n");

    /* GC and AT content */
    int gc_count = seq->count_G + seq->count_C;
    int at_count = seq->count_A + seq->count_T;
    double gc_content = (gc_count / total) * 100;
    double at_content = (at_count / total) * 100;

    printf("  GC Content: %.2f%%\n", gc_content);
    printf("  AT Content: %.2f%%\n", at_content);

    /* AT/GC ratio */
    if (gc_count > 0) {
        double at_gc_ratio = (double)at_count / gc_count;
        printf("  AT/GC Ratio: %.2f\n", at_gc_ratio);
    }

    /* Purine/Pyrimidine ratio */
    int purines = seq->count_A + seq->count_G;    /* A and G */
    int pyrimidines = seq->count_T + seq->count_C; /* T and C */
    if (pyrimidines > 0) {
        printf("  Purine/Pyrimidine Ratio: %.2f\n",
               (double)purines / pyrimidines);
    }

    printf("\n--- Biological Interpretation ---\n\n");

    if (gc_content > 60) {
        printf("  High GC content - Thermophilic organism or stable region\n");
    } else if (gc_content < 40) {
        printf("  Low GC content - AT-rich region, possibly regulatory\n");
    } else {
        printf("  Moderate GC content - Typical for many organisms\n");
    }

    printf("\n============================================\n");
}

/* Function to save results to file */
void save_results(DNASequence *seq, const char *output_file) {
    FILE *file = fopen(output_file, "w");
    if (file == NULL) {
        printf("Error: Cannot create output file\n");
        return;
    }

    fprintf(file, "DNA Sequence Analysis Results\n");
    fprintf(file, "=============================\n\n");
    fprintf(file, "Sequence: %s\n", seq->header);
    fprintf(file, "Length: %d bp\n\n", seq->length);
    fprintf(file, "Nucleotide Counts:\n");
    fprintf(file, "A: %d\n", seq->count_A);
    fprintf(file, "T: %d\n", seq->count_T);
    fprintf(file, "G: %d\n", seq->count_G);
    fprintf(file, "C: %d\n", seq->count_C);
    fprintf(file, "\nGC Content: %.2f%%\n",
            ((double)(seq->count_G + seq->count_C) / seq->length) * 100);

    fclose(file);
    printf("\nResults saved to: %s\n", output_file);
}

/* Main function */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <fasta_file>\n", argv[0]);
        printf("Example: %s ../data/sample1.fasta\n", argv[0]);
        return 1;
    }

    DNASequence seq;
    init_sequence(&seq);

    /* Read FASTA file */
    if (!read_fasta(argv[1], &seq)) {
        return 1;
    }

    /* Count nucleotides */
    count_nucleotides(&seq);

    /* Display statistics */
    display_statistics(&seq);

    /* Save results */
    save_results(&seq, "../output/analysis_results.txt");

    return 0;
}
