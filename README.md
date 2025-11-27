# Stabilizer Variable Test (SVT) - Simulation Code Repository

**Supplementary Materials for:**  
*Identification Theory and Testing of Stabilizer Variables in Multi-Group Structural Models*

**Authors:** This manuscript is currently under double-blind peer review. Information will be provided upon publication. 
 
**Journal:** Journal of the American Statistical Association (JASA)  
**Status:** Manuscript submitted for publication  
**Last Updated:** November 2025
**Expected Publication Year:** 2026 (pending acceptance)

---

## Overview

This repository contains complete R code for reproducing the Monte Carlo simulation study reported in Section 5 of the main manuscript. The simulation evaluates the Stabilizer Variable Test (SVT) across 813,500 datasets spanning three phases:

- **Phase 1:** Core performance evaluation (800,000 simulations)
- **Phase 2:** Sensitivity analyses (9,900 simulations)
- **Phase 3:** Boundary condition testing (3,600 simulations)

All simulations validate the SVT's ability to detect and classify stabilizer variables under varying conditions of measurement invariance violations, group counts, sample sizes, and noise levels.

---

## Repository Structure
```
SVT_Simulation/
|
+-- LICENSE                            # MIT License
+-- README.md                          # This file
|
+-- phase1_simulation.R                # Phase 1: Core functions and DGPs
+-- phase2_simulation.R                # Phase 2: Sensitivity analysis functions
+-- phase3_simulation.R                # Phase 3: Boundary condition functions
|
+-- run_phase1.R                       # Phase 1 execution script
+-- run_phase2.R                       # Phase 2 execution script
+-- run_phase3.R                       # Phase 3 execution script
|
+-- SimulationRecords/                 # Output directory (created on first run)
    +-- phase1_results_full.rds        # ~865 MB compressed
    +-- phase1_results_full.csv
    +-- phase2A_bootstrap_convergence.rds
    +-- phase2B_noise_trajectory.rds
    +-- phase2C_mi_trajectory.rds
    +-- phase3_boundary_conditions.rds
```

---

## System Requirements

### Software
- **R version:** 4.4.2 or higher (tested on R 4.4.2)
- **Platform:** x86_64-w64-mingw32/x64 (64-bit)
- **Operating System:** Windows, macOS, or Linux
- **Recommended RAM:** 16 GB minimum, 32 GB for optimal performance

### Required R Packages
```r
# Core packages
install.packages("parallel")
install.packages("foreach")
install.packages("doParallel")
install.packages("boot")
install.packages("dplyr")
```

### Package Versions (as used in original simulations)
```r
R version 4.4.2 (2024-10-31 ucrt) -- "Pile of Leaves"
Platform: x86_64-w64-mingw32/x64 (64-bit)

attached packages:
- parallel_4.4.2
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
required_packages <- c("parallel", "foreach", "doParallel", "boot", "dplyr")
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

#### Phase 1: Core Performance (~18-22 hours with 14 cores)
```r
source("run_phase1.R")
```

**Parameters:**
- Scenarios: TypeA, TypeB, TypeAB, Null, Moderator
- K: 5, 6, 7, 8, 9, 10, 15, 20
- n: 50, 100, 200, 500, 1000
- MI: 0.20, 0.30, 0.45, 0.65
- Replications: 1,000 per condition
- **Total: 800,000 simulations**

**Output:**
- `phase1_results_full.rds` (compressed R data, ~865 MB)
- `phase1_results_full.csv` (human-readable, larger file)

---

#### Phase 2: Sensitivity Analysis (~2-3 hours with 14 cores)
```r
source("run_phase2.R")
```

**Subphases:**
- **2A:** Bootstrap convergence (B = 500, 1000, 2000)
- **2B:** Noise robustness (sigma_epsilon = 0.20 to 0.70)
- **2C:** MI trajectory (MI = 0.15 to 0.70)

**Total: 9,900 simulations**

**Outputs:**
- `phase2A_bootstrap_convergence.rds`
- `phase2B_noise_trajectory.rds`
- `phase2C_mi_trajectory.rds`

---

#### Phase 3: Boundary Conditions (~30-45 minutes with 14 cores)
```r
source("run_phase3.R")
```

**Extreme configurations:**
- K in {3, 50}
- n in {30, 500}
- MI in {0.10, 0.90}
- 18 configurations x 200 reps = **3,600 simulations**

**Output:**
- `phase3_boundary_conditions.rds`

---

## Data File Formats

### .rds Files (Recommended)
- Compressed R-native format
- Preserves exact data types
- Smaller file size
- Load with: `results <- readRDS("phase1_results_full.rds")`

### .csv Files (Optional)
- Human-readable
- Compatible with Excel, Python, etc.
- Larger file size
- Load with: `results <- read.csv("phase1_results_full.csv")`

### WARNING: Avoid .RData Files
Do **NOT** use or share `.RData` or `.RDataTmp` files:
- They contain entire workspace environments
- Not reproducible across systems
- Unnecessarily large
- Delete them if created accidentally

---

## Adjusting Computational Resources

### Reduce Number of Cores

Edit the `n_cores` parameter in run scripts:
```r
# Default: 14 cores
results <- run_phase1_parallel(n_cores = 14, n_reps = 1000)

# Reduce to 4 cores (slower but less resource-intensive)
results <- run_phase1_parallel(n_cores = 4, n_reps = 1000)
```

### Reduce Number of Replications

For faster testing:
```r
# Full run: 1,000 reps per condition
results <- run_phase1_parallel(n_cores = 14, n_reps = 1000)

# Quick test: 100 reps per condition
results <- run_phase1_parallel(n_cores = 14, n_reps = 100)
```

---

## Reproducibility

### Random Seed Management

All simulations use a **base seed of 9186** with deterministic offsets:
```r
# Phase 1
seed = 9186 + iteration_number

# Phase 2
seed = 10000 + iteration_number  # 2A-2B
seed = 20000 + iteration_number  # 2C

# Phase 3
seed = 30000 + iteration_number

## Data Availability

- **Simulation Code:** This repository (https://github.com/sy142/stabilizer-variable-simulations)
- **Simulation Results:** Figshare (https://doi.org/10.6084/m9.figshare.30731633)

```

This ensures:
- Full reproducibility with identical R version and packages
- Independence across phases
- Deterministic results for each condition

### Verification

To verify exact reproducibility:
```r
# Run same condition twice
set.seed(9186)
result1 <- run_single_sim("TypeAB", K=10, n=200, mi_severity=0.45, seed=9186)

set.seed(9186)
result2 <- run_single_sim("TypeAB", K=10, n=200, mi_severity=0.45, seed=9186)

all.equal(result1, result2)  # Should return TRUE
```

---

## Data Generating Processes (DGPs)

### Type A: Variance Compression
- Negative correlation: X and Z via opposite-sign U
- Group-specific bias: b_k ~ N(0, (0.60 x MI)^2)
- **Mechanism:** Absorbs heterogeneous confounds

### Type B: Directional Alignment
- Positive correlation: X and Z via same-sign U
- Constant directional bias: gamma = 0.45 x MI + 0.15
- **Mechanism:** Aligns slopes toward common direction

### Type AB: Combined
- Mixture: Beta(3,7) weighting between b_k and gamma_k
- Sign flips: 15% probability of negative gamma_k
- **Mechanism:** Both variance compression and alignment

### Null: No Stabilizer
- Independent X and Z
- Natural heterogeneity: beta_k ~ N(-0.20, 0.10^2)
- **Purpose:** Type I error control

### Moderator: Interaction Only
- X x Z interaction: 0.25
- Violates structural independence
- **Purpose:** Distinguish from stabilization

Full mathematical specifications in **Supplementary Section S3**.

---

## Output Data Structure

Each simulation produces a row with the following key variables:
```r
# Scenario identifiers
scenario, K, n, mi_severity

# Baseline estimates
mean_0, sd_0, cv_0, beta_baseline (per group)

# Adjusted estimates
mean_1, sd_1, cv_1, beta_adjusted (per group)

# Effect sizes
delta_log          # Log-ratio of CVs (primary metric)
pct_reduction      # Percentage CV reduction
var_reduction_pct  # Variance reduction
orientation_share  # Alignment vs compression

# Inference
p_value            # Bootstrap test
z_stat             # Standardized effect
ci_lower, ci_upper # 95% CI

# Mechanism classification
bootstrap_pass
binom_variance_pass
binom_alignment_pass
mechanism_detected  # TypeA, TypeB, TypeAB, Marginal, None
```

---

## Performance Metrics

### Expected Runtimes (14 cores, AMD Ryzen/Intel i7)

| Phase    | Simulations | Time       | Output Size      |
|----------|-------------|------------|------------------|
| Phase 1  | 800,000     | 18-22 hrs  | ~865 MB (.rds)   |
| Phase 2A | 3,000       | 15-20 min  | ~30 MB           |
| Phase 2B | 3,300       | 20-25 min  | ~35 MB           |
| Phase 2C | 3,600       | 25-30 min  | ~40 MB           |
| Phase 3  | 3,600       | 30-45 min  | ~40 MB           |

### Memory Usage

- **Peak RAM:** ~12-16 GB during parallel execution
- **Per-core:** ~1 GB average
- **Output storage:** ~1 GB total (compressed .rds)

---

## Troubleshooting

### Common Issues

#### 1. Parallel backend errors
```r
# If parallel cluster fails to start
library(parallel)
detectCores()  # Check available cores

# Reduce cores if needed
cl <- makeCluster(4)  # Instead of 14
```

#### 2. Memory errors
```r
# Reduce batch size or run sequentially
results <- lapply(1:nrow(conditions), function(i) {
  run_single_sim(...)
})
```

#### 3. Missing packages
```r
# Reinstall all dependencies
install.packages(c("parallel", "foreach", "doParallel", "boot", "dplyr"))
```

#### 4. Path issues
```r
# Always set working directory explicitly
setwd("C:/path/to/SVT_Simulation")
getwd()  # Verify
```

#### 5. Large file handling
```r
# If .csv files are too large, use .rds only
saveRDS(results, "output.rds")  # Compressed
results <- readRDS("output.rds")  # Fast loading
```

---

## Citation

If you use this code, please cite:
```bibtex
@article{This manuscript is currently under double-blind peer review. Information will be provided upon publication. ,
  title={Identification Theory and Testing of Stabilizer Variables in Multi-Group Structural Models},
  author={This manuscript is currently under double-blind peer review. Information will be provided upon publication. },
  journal={Journal of the American Statistical Association},
  year={2026},
  note={Manuscript submitted for publication. Expected publication year: 2026 (pending acceptance). DOI to be assigned upon acceptance}
}
```

---

## License

This code is released under the **MIT License**. See [LICENSE](LICENSE) file for details.

---

## Contact

This manuscript is currently under double-blind peer review.

Contact information will be provided upon publication. For questions regarding the code or methodology, please open an issue in this repository.


---

## Acknowledgments

Institutional support for article processing charges has been secured.

---

**Repository Version:** 1.0.0  
**Last Updated:** November 2025  

**R Version Tested:** 4.4.2 (2024-10-31)



