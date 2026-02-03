# Enzyme Kinetics Analysis
# Course: M.Sc. Medical Bioinformatics - Semester 1
# Skills: R, Statistics, Biochemistry

# ============================================
# 1. SETUP AND DATA LOADING
# ============================================

# Load data
data <- read.csv("../data/enzyme_data.csv", header = TRUE)

cat("=== ENZYME KINETICS DATA ===\n")
print(data)

# Extract values
S <- data$substrate_conc  # Substrate concentration [S]
V_A <- data$enzyme_A_velocity
V_B <- data$enzyme_B_velocity
V_C <- data$enzyme_C_velocity

# ============================================
# 2. LINEWEAVER-BURK TRANSFORMATION
# ============================================

# Double reciprocal transformation: 1/V vs 1/[S]
inv_S <- 1 / S
inv_V_A <- 1 / V_A
inv_V_B <- 1 / V_B
inv_V_C <- 1 / V_C

# ============================================
# 3. LINEAR REGRESSION FOR KINETIC PARAMETERS
# ============================================

# Function to calculate Km and Vmax from Lineweaver-Burk
calculate_kinetics <- function(inv_S, inv_V, enzyme_name) {
  model <- lm(inv_V ~ inv_S)
  intercept <- coef(model)[1]  # 1/Vmax
  slope <- coef(model)[2]      # Km/Vmax

  Vmax <- 1 / intercept
  Km <- slope * Vmax

  cat(paste("\n=== ", enzyme_name, " ===\n", sep = ""))
  cat(paste("Vmax:", round(Vmax, 2), "units\n"))
  cat(paste("Km:", round(Km, 2), "mM\n"))
  cat(paste("R-squared:", round(summary(model)$r.squared, 4), "\n"))

  return(list(Vmax = Vmax, Km = Km, model = model))
}

# Calculate for each enzyme
kinetics_A <- calculate_kinetics(inv_S, inv_V_A, "Enzyme A")
kinetics_B <- calculate_kinetics(inv_S, inv_V_B, "Enzyme B")
kinetics_C <- calculate_kinetics(inv_S, inv_V_C, "Enzyme C")

# ============================================
# 4. SUMMARY TABLE
# ============================================

cat("\n=== KINETIC PARAMETERS SUMMARY ===\n")
summary_table <- data.frame(
  Enzyme = c("Enzyme A", "Enzyme B", "Enzyme C"),
  Vmax = c(kinetics_A$Vmax, kinetics_B$Vmax, kinetics_C$Vmax),
  Km = c(kinetics_A$Km, kinetics_B$Km, kinetics_C$Km)
)
summary_table$Vmax <- round(summary_table$Vmax, 2)
summary_table$Km <- round(summary_table$Km, 2)
summary_table$Efficiency <- round(summary_table$Vmax / summary_table$Km, 2)
print(summary_table)

# Interpretation
cat("\n=== INTERPRETATION ===\n")
best_affinity <- summary_table$Enzyme[which.min(summary_table$Km)]
best_velocity <- summary_table$Enzyme[which.max(summary_table$Vmax)]
most_efficient <- summary_table$Enzyme[which.max(summary_table$Efficiency)]

cat(paste("Highest substrate affinity (lowest Km):", best_affinity, "\n"))
cat(paste("Highest maximum velocity (Vmax):", best_velocity, "\n"))
cat(paste("Most efficient enzyme (Vmax/Km):", most_efficient, "\n"))

# ============================================
# 5. VISUALIZATION
# ============================================

pdf("../output/enzyme_kinetics_plots.pdf", width = 12, height = 10)

# Plot 1: Michaelis-Menten Curves
par(mfrow = c(2, 2))

plot(S, V_A, type = "b", pch = 19, col = "red",
     main = "Michaelis-Menten Curves",
     xlab = "Substrate Concentration [S] (mM)",
     ylab = "Reaction Velocity (V)",
     ylim = c(0, max(V_C) + 2))
points(S, V_B, type = "b", pch = 19, col = "blue")
points(S, V_C, type = "b", pch = 19, col = "green")
legend("bottomright",
       legend = c("Enzyme A", "Enzyme B", "Enzyme C"),
       col = c("red", "blue", "green"),
       pch = 19, lty = 1)

# Plot 2: Lineweaver-Burk Plot
plot(inv_S, inv_V_A, pch = 19, col = "red",
     main = "Lineweaver-Burk Plot (Double Reciprocal)",
     xlab = "1/[S] (1/mM)",
     ylab = "1/V",
     ylim = c(0, max(inv_V_B) + 0.1))
points(inv_S, inv_V_B, pch = 19, col = "blue")
points(inv_S, inv_V_C, pch = 19, col = "green")
abline(kinetics_A$model, col = "red", lty = 2)
abline(kinetics_B$model, col = "blue", lty = 2)
abline(kinetics_C$model, col = "green", lty = 2)
legend("topright",
       legend = c("Enzyme A", "Enzyme B", "Enzyme C"),
       col = c("red", "blue", "green"),
       pch = 19, lty = 2)

# Plot 3: Bar plot of Km values
barplot(summary_table$Km,
        names.arg = summary_table$Enzyme,
        main = "Km Comparison (Lower = Higher Affinity)",
        ylab = "Km (mM)",
        col = c("red", "blue", "green"))

# Plot 4: Bar plot of Vmax values
barplot(summary_table$Vmax,
        names.arg = summary_table$Enzyme,
        main = "Vmax Comparison",
        ylab = "Vmax (units)",
        col = c("red", "blue", "green"))

dev.off()

cat("\n=== PLOTS SAVED TO ../output/enzyme_kinetics_plots.pdf ===\n")

# ============================================
# 6. SAVE RESULTS
# ============================================

write.csv(summary_table, "../output/kinetic_parameters.csv", row.names = FALSE)
cat("=== RESULTS SAVED TO ../output/kinetic_parameters.csv ===\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
