# Gene Expression Statistical Analysis

## Problem Statement

Gene expression data helps researchers understand which genes are active under different conditions (e.g., healthy vs diseased tissue). However, raw expression data requires statistical analysis to identify significant differences.

**Objective:** Analyze gene expression data to identify differentially expressed genes between two conditions (e.g., normal vs tumor samples) using statistical methods.

**Questions to Answer:**
1. Which genes show significant difference in expression between conditions?
2. What is the distribution of gene expression values?
3. Are the observed differences statistically significant?

---

## Solution Approach

1. **Data Collection:** Use sample gene expression data (simulated or from NCBI GEO)
2. **Data Preprocessing:** Load and clean data in R
3. **Statistical Analysis:**
   - Calculate mean, median, standard deviation for each condition
   - Perform t-test to identify significant differences
   - Calculate fold change
4. **Visualization:**
   - Boxplots comparing conditions
   - Histogram of expression distribution
   - Volcano plot (fold change vs p-value)
5. **Interpretation:** Identify top differentially expressed genes

---

## Files

```
01_gene_expression_analysis/
├── README.md
├── src/
│   └── analysis.R
├── data/
│   └── gene_expression.csv
└── output/
    └── (generated plots and results)
```

---

## How to Run

```bash
cd src
Rscript analysis.R
```

Or open `analysis.R` in RStudio and run interactively.

---

## Expected Output

- Statistical summary table
- Boxplot comparing conditions
- List of significantly differentially expressed genes
- Volcano plot

---

## Skills Demonstrated

- R programming
- Statistical analysis (t-test, p-values)
- Data visualization (ggplot2)
- Cell and Molecular Biology concepts
