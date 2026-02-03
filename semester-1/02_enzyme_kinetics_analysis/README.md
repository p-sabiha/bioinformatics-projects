# Enzyme Kinetics Analysis

## Problem Statement

Enzymes are biological catalysts that speed up biochemical reactions. Understanding enzyme kinetics helps in drug development and understanding metabolic disorders. The Michaelis-Menten model describes how reaction velocity depends on substrate concentration.

**Objective:** Analyze enzyme kinetics data to calculate kinetic parameters (Km and Vmax) and compare different enzymes.

**Questions to Answer:**
1. What are the Km and Vmax values for each enzyme?
2. Which enzyme has higher substrate affinity (lower Km)?
3. Which enzyme has higher maximum velocity (Vmax)?

---

## Solution Approach

1. **Data Collection:** Use enzyme kinetics data (substrate concentration vs reaction velocity)
2. **Michaelis-Menten Model:** V = (Vmax × [S]) / (Km + [S])
3. **Analysis Methods:**
   - Plot substrate concentration vs velocity
   - Lineweaver-Burk plot (double reciprocal): 1/V vs 1/[S]
   - Calculate Km and Vmax from linear regression
4. **Comparison:** Compare kinetic parameters across enzymes
5. **Interpretation:** Relate parameters to enzyme efficiency

---

## Key Concepts

- **Km (Michaelis constant):** Substrate concentration at half-maximal velocity. Lower Km = higher affinity.
- **Vmax:** Maximum reaction velocity when enzyme is saturated with substrate.
- **Lineweaver-Burk Plot:** Linear transformation where slope = Km/Vmax, y-intercept = 1/Vmax

---

## Files

```
02_enzyme_kinetics_analysis/
├── README.md
├── src/
│   └── kinetics.R
├── data/
│   └── enzyme_data.csv
└── output/
    └── (generated plots and results)
```

---

## How to Run

```bash
cd src
Rscript kinetics.R
```

---

## Expected Output

- Michaelis-Menten curves
- Lineweaver-Burk plots
- Table of Km and Vmax values
- Enzyme comparison summary

---

## Skills Demonstrated

- R programming
- Curve fitting and regression
- Data visualization
- Biochemistry concepts (enzyme kinetics)
