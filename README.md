# Stabilizer Variable Test (SVT) - Simulation Code Repository

**Supplementary Materials for:**  
*Stabilizer Variables for Measurement Invariance–Induced Heterogeneity: Identification Theory and Testing in Multi-Group Models*

**Authors:** Salim Yilmaz¹, Erhan Cene²  
¹ Department of Health Management, Faculty of Health Sciences, Acıbadem Mehmet Ali Aydınlar University, İstanbul, Türkiye  
² Department of Statistics, Faculty of Science, Yıldız Technical University, İstanbul, Türkiye

**Journal:** *Mathematics* 
**Last Updated:** February 2026

---

## Overview

This repository contains complete R code for reproducing the Monte Carlo simulation study reported in Sections 4–5 of the main manuscript. The simulation evaluates the Stabilizer Variable Test (SVT) across **817,100 total replications** spanning four phases:

- **Phase 0:** Adaptive MI scoring validation (3,600 simulations)
- **Phase 1:** Core performance evaluation (800,000 simulations)
- **Phase 2:** Sensitivity analyses (9,900 simulations)
- **Phase 3:** Boundary condition testing (3,600 simulations)

All simulations validate the SVT's ability to detect and classify stabilizer variables under varying conditions of measurement invariance violations, group counts, sample sizes, and noise levels.

---

## Repository Structure

```
SVT_Simulation/
│
├── LICENSE                            # MIT License
├── README.md                          # This file
│
├── phase0_simulation.R                # Phase 0: MI scoring validation functions
├── phase1_simulation.R                # Phase 1: Core functions and DGPs
├── phase2_simulation.R                # Phase 2: Sensitivity analysis functions
├── phase3_simulation.R                # Phase 3: Boundary condition functions
│
├── run_phase0.R                       # Phase 0 execution script
├── run_phase1.R                       # Phase 1 execution script
├── run_phase2.R                       # Phase 2 execution script
├── run_phase3.R                       # Phase 3 execution script
│
└── SimulationRecords/                 # Output directory (created on first run)
    ├── phase0_results.rds
    ├── phase0_summary.rds
    ├── phase1_results_full.rds
    ├── phase1_results_full.csv
    ├── phase2A_bootstrap_convergence.rds
    ├── phase2B_noise_trajectory.rds
    ├── phase2C_mi_trajectory.rds
    └── phase3_boundary_conditions.rds
```

---

## System Requirements

### Software
- **R version:** 4.5.2 or higher (tested on R 4.5.2)
- **Operating System:** Windows, macOS, or Linux
- **Recommended RAM:** 16 GB minimum, 32 GB for optimal performance

### Required R Packages

```r
# Phase 0 (MI scoring validation)
install.packages("lavaan")

# Phases 1–3 (Monte Carlo simulations)
install.packages(c("parallel", "foreach", "doParallel", "boot", "dplyr"))
```

### Package Versions (as used in original simulations)

```r
R version 4.5.2 (2025-06-13)

attached packages:
- lavaan_0.6-19       # Phase 0 only
- parallel_4.5.2
- foreach_1.5.2
- doParallel_1.0.17
- boot_1.3-30
- dplyr_1.1.4
```

---

## Installation

1. **Clone or download this repository**
2. **Set working directory** to the repository folder
3. **Create output directory:**

```r
dir.create("SimulationRecords", showWarnings = FALSE)
```

4. **Verify package installation:**

```r
required_packages <- c("lavaan", "parallel", "foreach", "doParallel", "boot", "dplyr")
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}
```

---

## Running Simulations

### Quick Start (Small Test Run)

To verify setup without full simulation:

```r
# Test Phase 1 with reduced parameters
source("phase1_simulation.R")

test_result <- run_single_sim(
  scenario = "TypeAB",
  K = 5,
  n = 50,
  mi_severity = 0.45,
  seed = 9186
)

print(test_result)
```

Expected output: A list containing metrics like `delta_log`, `p_value`, `cv_0`, `cv_1`, etc.

---

### Full Simulation Runs

#### Phase 0: MI Scoring Validation (~4–8 hours with 14 cores)

```r
source("run_phase0.R")
```

**Parameters:**
- MI severity: 0.10, 0.20, 0.30, 0.45, 0.65, 0.90
- Sample sizes per group: 100, 200, 500
- Moderators per replication: 10 (5 non-invariant, 5 invariant)
- Replications: 200 per condition
- **Total: 3,600 simulations** (108,000 CFA model fits)

**Outputs:**
- `phase0_results.rds` — Full replication-level results
- `phase0_summary.rds` — Aggregated summary statistics

---

#### Phase 1: Core Performance (~18–22 hours with 14 cores)

```r
source("run_phase1.R")
```

**Parameters:**
- Scenarios: TypeA, TypeB, TypeAB, Null, Moderator
- K: 5, 6, 7, 8, 9, 10, 15, 20
- n: 50, 100, 200, 500, 1000
- MI: 0.20, 0.30, 0.45, 0.65
- Replications: 1,000 per condition
- Bootstrap iterations: B = 200
- **Total: 800,000 simulations**

**Outputs:**
- `phase1_results_full.rds` (compressed R data)
- `phase1_results_full.csv` (human-readable)

---

#### Phase 2: Sensitivity Analysis (~2–3 hours with 14 cores)

```r
source("run_phase2.R")
```

**Subphases:**
- **2A: Bootstrap convergence** — B ∈ {500, 1000, 2000}, TypeAB + Null scenarios, K=20, n=200, MI=0.45, 500 reps each (3,000 simulations)
- **2B: Noise robustness** — σ_ε from 0.20 to 0.70 in increments of 0.05 (11 levels), TypeAB, K=20, n=200, MI=0.45, B=1000, 300 reps each (3,300 simulations)
- **2C: MI trajectory** — MI from 0.15 to 0.70 in increments of 0.05 (12 levels), TypeAB, K=20, n=200, B=1000, 300 reps each (3,600 simulations)

**Total: 9,900 simulations**

**Outputs:**
- `phase2A_bootstrap_convergence.rds`
- `phase2B_noise_trajectory.rds`
- `phase2C_mi_trajectory.rds`

---

#### Phase 3: Boundary Conditions (~30–45 minutes with 14 cores)

```r
source("run_phase3.R")
```

**Configurations (18 total):**

| Category | K | n | MI | Description |
|----------|---|---|----|-------------|
| Extreme K (min) | 3 | 100, 200, 500 | 0.45 | Minimum group count |
| Extreme K (max) | 50 | 100, 200, 500 | 0.45 | Maximum group count |
| Extreme MI | 20 | 200 | 0.10, 0.90 | Weak and extreme violations |
| Extreme MI | 10, 5 | 200 | 0.90 | Severe violation, fewer groups |
| Minimal n | 10, 20 | 30 | 0.45, 0.65 | Small within-group samples |
| Worst case | 3 | 30 | 0.90 | Triple challenge |
| Edge cases | 50/3/50 | 30/500/500 | 0.45/0.90/0.10 | Mixed extremes |

- Replications: 200 per configuration
- Scenario: TypeAB only
- **Total: 3,600 simulations**

**Output:**
- `phase3_boundary_conditions.rds`

---

## Data Generating Processes (DGPs)

Five scenarios are implemented in `phase1_simulation.R`, corresponding to Equations (49)–(54) in the manuscript:

### Scenario 1 — Type A: Variance Purification

- Negative correlation between X and Z via opposite-sign latent confound U (ρ ≈ 0.15)
- Group-specific bias: a_k ~ N(0, MI · 0.6)
- **Mechanism:** Z absorbs heterogeneous group-specific confounds, compressing β̂_k distribution

### Scenario 2 — Type B: Directional Alignment

- Positive correlation between X and Z via same-sign U (ρ ≈ 0.35)
- Constant directional bias: c = 0.45 · MI + 0.15
- **Mechanism:** Z aligns slopes toward a common direction without compressing variance

### Scenario 3 — Type AB: Combined Mechanism

- Mixture weighting: λ_k ~ Beta(3, 7) between a_k and γ_k components
- Sign flips: 15% probability of negative γ_k
- **Mechanism:** Simultaneous variance compression and directional alignment

### Scenario 4 — Null: No Stabilizer

- Independent X and Z (no shared confound)
- Natural heterogeneity: β_k ~ N(−0.20, 0.10²)
- **Purpose:** Type I error control verification

### Scenario 5 — Moderator: Interaction Only

- X × Z interaction coefficient: δ_k = 0.25
- Violates structural independence (∂β/∂Z ≠ 0)
- **Purpose:** Distinguish stabilization from classical moderation

Full mathematical specifications are provided in Supplementary Material, Section S4.

---

## Inference Implementation

The Monte Carlo study implements the SVT dual-criterion framework (Section 4.1.4 of the manuscript):

**Criterion 1 — Permutation test for stabilization magnitude:**
For computational efficiency, the bootstrap criterion is implemented using a sign-flip permutation procedure. Centered group-level effects are randomly sign-flipped to construct a null distribution of Δℓ, which is asymptotically equivalent to the nonparametric bootstrap under the symmetric null hypothesis of no stabilization. The empirical application (not included in this repository) uses the group-level bootstrap variant via the `boot` package.

**Criterion 2 — Binomial test for directional consistency:**
Counts the number of groups where |β̂_k^(1) − β̄^(1)| < |β̂_k^(0) − β̄^(0)|, tested against Binomial(K, 0.5).

**Decision rule:** Both criteria must be satisfied at α = 0.05.

---

## Output Data Structure

Each simulation produces a row with the following key variables:

```r
# Scenario identifiers
scenario, K, n, mi_severity

# Baseline estimates
mean_0, sd_0, cv_0       # Cross-group statistics before stabilization

# Adjusted estimates
mean_1, sd_1, cv_1       # Cross-group statistics after stabilization

# Effect sizes
delta_log                 # Log-ratio of CVs: log(CV⁰) − log(CV¹) [primary metric]
pct_reduction             # Percentage CV reduction
var_reduction_pct         # Variance reduction percentage
orientation_share         # Type B contribution: 0 = pure variance, 1 = pure alignment

# Inference
p_value                   # Permutation test p-value
z_stat                    # Standardized effect size
ci_lower, ci_upper        # 95% percentile confidence interval

# Mechanism classification
bootstrap_pass            # Criterion 1 satisfied (p < 0.05 and delta_log > 0.05)
binom_variance_pass       # Binomial test for variance compression
binom_alignment_pass      # Binomial test for directional alignment
mechanism_detected        # TypeA, TypeB, TypeAB, Marginal, None
```

---

## Performance Metrics

### Expected Runtimes (14 cores, AMD Ryzen / Intel i7)

| Phase | Simulations | Estimated Time | Output Size (.rds) |
|-------|-------------|----------------|---------------------|
| Phase 0 | 3,600 | 4–8 hours | ~500 KB |
| Phase 1 | 800,000 | 18–22 hours | ~102 MB |
| Phase 2A | 3,000 | 15–20 min | ~387 KB |
| Phase 2B | 3,300 | 20–25 min | ~416 KB |
| Phase 2C | 3,600 | 25–30 min | ~454 KB |
| Phase 3 | 3,600 | 30–45 min | ~460 KB |

### Memory Usage

- **Peak RAM:** ~12–16 GB during parallel execution
- **Per-core:** ~1 GB average
- **Output storage:** ~105 MB total (compressed .rds)

---

## Adjusting Computational Resources

### Reduce Number of Cores

Edit the `n_cores` parameter in run scripts:

```r
# Default: 14 cores
n_cores <- 14

# Reduce to 4 cores (slower but less resource-intensive)
n_cores <- 4
```

### Reduce Number of Replications

For faster testing:

```r
# Full run: 1,000 reps per condition (Phase 1)
conditions <- expand.grid(..., rep = 1:1000, ...)

# Quick test: 100 reps per condition
conditions <- expand.grid(..., rep = 1:100, ...)
```

---

## Reproducibility

### Random Seed Management

All simulations use a **base seed of 9186** with deterministic offsets:

```r
# Phase 0: seed = 9186 + iteration_number
# Phase 1: seed = 9186 + iteration_number
# Phase 2A/2B: seed = 10000 + iteration_number
# Phase 2C: seed = 20000 + iteration_number
# Phase 3: seed = 30000 + iteration_number
```

This ensures:
- Full reproducibility with identical R version and packages
- Independence across phases
- Deterministic results for each condition

### Verification

To verify exact reproducibility:

```r
result1 <- run_single_sim("TypeAB", K=10, n=200, mi_severity=0.45, seed=9186)
result2 <- run_single_sim("TypeAB", K=10, n=200, mi_severity=0.45, seed=9186)
all.equal(result1, result2)  # Should return TRUE
```

---

## Data File Formats

### .rds Files (Recommended)
- Compressed R-native format; preserves exact data types
- Load with: `results <- readRDS("phase1_results_full.rds")`

### .csv Files (Optional)
- Human-readable; compatible with Python, Excel, etc.
- Load with: `results <- read.csv("phase1_results_full.csv")`

### Note on .RData Files
Do **not** use or share `.RData` or `.RDataTmp` files — they contain entire workspace environments, are not portable across systems, and are unnecessarily large.

---

## Troubleshooting

### Common Issues

**1. Parallel backend errors**

```r
library(parallel)
detectCores()        # Check available cores
cl <- makeCluster(4) # Reduce if needed
```

**2. Memory errors**

```r
# Reduce batch size or run sequentially
results <- lapply(1:nrow(conditions), function(i) {
  run_single_sim(...)
})
```

**3. Phase 0 — lavaan convergence warnings**

Phase 0 uses CFA model fitting which may produce convergence warnings for some conditions. Non-converging replications are automatically excluded (`tryCatch` with `NULL` return). This is expected behavior, particularly at extreme MI severity levels.

**4. Path issues**

```r
# Update working directory in run scripts before execution
setwd("/path/to/SVT_Simulation")
```

---

## Data Availability

- **Simulation Code:** This repository — [https://github.com/sy142/stabilizer-variable-simulations](https://github.com/sy142/stabilizer-variable-simulations)
- **Simulation Results & Datasets:** Figshare — [https://doi.org/10.6084/m9.figshare.30731633](https://doi.org/10.6084/m9.figshare.30731633)

---

## Citation

```bibtex
@article{yilmaz2026stabilizer,
  title   = {Stabilizer Variables for Measurement Invariance--Induced Heterogeneity:
             Identification Theory and Testing in Multi-Group Models},
  author  = {Yilmaz, Salim and Cene, Erhan},
  journal = {Mathematics},
  year    = {2026},
  note    = {Manuscript submitted for publication}
}
```

---

## License

This code is released under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## Contact

- **Salim Yilmaz** — salim.yilmaz@acibadem.edu.tr | ORCID: [0000-0003-2405-5084](https://orcid.org/0000-0003-2405-5084)
- **Erhan Cene** — ecene@yildiz.edu.tr | ORCID: [0000-0001-5336-6004](https://orcid.org/0000-0001-5336-6004)

---

## Acknowledgments

This manuscript is based on a master's thesis completed in the Department of Statistics at Yıldız Technical University, İstanbul, Türkiye. Article processing charges were funded by Acıbadem Mehmet Ali Aydınlar University.

---

**Repository Version:** 2.0.0  
**Last Updated:** February 2026  
**R Version Tested:** 4.5.2
