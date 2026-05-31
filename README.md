# MatchedGFormula

Reference implementation and simulation code for the paper:

> **Scalable Counterfactual Risk Estimation for Rare Events in Longitudinal Data**
> Xiaohui Yin, Avijit Mitra, Ying Zhou, Kun Chen, Hong Yu.

The repository implements the Iterative Conditional Expectation (ICE) and
Non-Iterative Conditional Expectation (NICE) g-formula estimators of the
counterfactual risk under static treatment regimes, together with a
case–control subsampling and reweighting strategy that substantially reduces
the computational cost of bootstrap-based inference in rare-outcome
longitudinal data.

## Contents

```
.
├── 0_regularizedParameterSetup.R   # Simulation parameters + data-generating function
├── 1_regularizedTrueACE.R          # Monte-Carlo evaluation of true counterfactual risks
├── 2_regularizedDataGenerating.R   # Generate one simulated dataset (n = 30,000, K = 6)
├── 3_regularizedDataFitting.R      # Penalized (glmnet) g-formula run on one dataset
├── 4_dataFittingBootstrap.R        # Main bootstrap reproduction script (B = 100)
├── 5_dataFittingBootstrapSingleRun.R  # One-bootstrap-per-task variant for HPC arrays
├── 6_Analysis_Boots.R              # Aggregate bootstrap outputs into tables / plots
├── updatedMatchedGformula.R        # Core ICE / NICE implementation used in the paper
├── regularizedMatchedGformula.R    # Penalized variant (glmnet) used in script 3_
├── MatchedGFormula.Rproj           # RStudio project file
├── LICENSE                         # MIT
└── README.md
```

`updatedMatchedGformula.R` exposes four estimator entry points:

| Function                  | Form  | Subsampling      |
|---------------------------|-------|------------------|
| `gform_noniter_complete`  | NICE  | Complete data    |
| `gform_noniter_match`     | NICE  | Case–control (J) |
| `gform_iter_complete`     | ICE   | Complete data    |
| `gform_iter_match`        | ICE   | Case–control (J) |

`regularizedMatchedGformula.R` provides drop-in penalized counterparts that fit
each nuisance model with `glmnet`; it is used by `3_regularizedDataFitting.R`
to illustrate the discussion in the paper on penalized nuisance fits.

## Requirements

- R ≥ 4.4 (results in the paper were produced with R 4.4.1)
- CRAN packages: `data.table`, `dplyr`, `tidyr`, `ggplot2`, `fastDummies`,
  `glmnet`, `coda`

Install with:

```r
install.packages(c("data.table", "dplyr", "tidyr", "ggplot2",
                   "fastDummies", "glmnet", "coda"))
```

## Configuration

Every script reads its main configuration from environment variables so that
the same code can be run locally and on a SLURM array without edits:

| Variable               | Default | Meaning                                         |
|------------------------|--------:|-------------------------------------------------|
| `GFORM_N`              | `30000` | Number of individuals per simulated dataset     |
| `GFORM_K`              | `6`     | Number of follow-up time points                 |
| `GFORM_SETUP`          | `5`     | Simulation setup: 5 ≈ 1%, 6 ≈ 2-3%, 7 ≈ 5-8%   |
| `GFORM_J`              | `5`     | Case–control sampling ratio (controls per case) |
| `GFORM_NBOOT`          | `100`   | Bootstrap replicates per dataset                |
| `SLURM_ARRAY_TASK_ID`  | `1`     | Array index dispatched to                       |

The mapping between `GFORM_SETUP` and the paper's prevalence labels is fixed
by the intercept of the outcome model:

| `GFORM_SETUP` | Outcome intercept | Paper label       |
|---------------|-------------------|-------------------|
| `5`           | `-6`              | low prevalence (~1%)     |
| `6`           | `-5`              | medium prevalence (~2.5%) |
| `7`           | `-4`              | high prevalence (~6%)     |

## Reproducing the paper's simulation results

The simulation pipeline runs four phases in order. Working directory is the
repository root for all commands below.

### 1. Generate the 100 simulated datasets

```sh
# For each of setup 5, 6, 7 generate replicates 1..100
for setup in 5 6 7; do
  for i in $(seq 1 100); do
    GFORM_SETUP=$setup SLURM_ARRAY_TASK_ID=$i Rscript 2_regularizedDataGenerating.R
  done
done
```

This populates `data/Data_setup{5,6,7}_replicate{1..100}.rda`.

On a SLURM cluster the loops above can be replaced by an array job:

```sh
sbatch --array=1-100 --export=GFORM_SETUP=5 run_data_gen.sh
```

### 2. Compute the true counterfactual risks

`1_regularizedTrueACE.R` is sourced after `0_regularizedParameterSetup.R`; the
two scalars `true_risk0` and `true_risk1` (length-`K` vectors of cumulative
risks under always-treat and never-treat) are produced for the active setup.

To save them in the format consumed by `6_Analysis_Boots.R`:

```r
setup <- 5
source('0_regularizedParameterSetup.R')
source('1_regularizedTrueACE.R')
df_true <- rbind(data.frame(Time = 1:6, TrueRisk = true_risk0, treatment = 0),
                 data.frame(Time = 1:6, TrueRisk = true_risk1, treatment = 1))
dir.create('Boots_Results', showWarnings = FALSE)
write.table(df_true,
            file = sprintf('Boots_Results/n_30000_setup%d_true.txt', setup),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
```

### 3. Run the bootstrap g-formula estimators

Each SLURM array task fits a single (estimator, treatment) combination on a
single bootstrap of a single replicate. With 8 estimator × treatment
combinations and 100 replicates the array index runs `1..800` per setup:

```sh
# Reproduce paper Table 4 / Figure 3 (setup 5, J = 5)
for i in $(seq 1 800); do
  GFORM_SETUP=5 GFORM_J=5 SLURM_ARRAY_TASK_ID=$i Rscript 4_dataFittingBootstrap.R
done
```

The paper's Table 4 also reports `J = 10` and `J = 20`. **Those columns are
not in this repository.** To reproduce them, re-run the same loop with
`GFORM_J=10` and `GFORM_J=20`. Each run writes
`results_bootstrap/n_30000_setup<S>_J<J>_iter<i>_matched<m>_treat<t>_simu<r>.RDS`.

A finer-grained variant (`5_dataFittingBootstrapSingleRun.R`) emits one text
file per (bootstrap, dataset, estimator, treatment) tuple and is more suitable
for very wide HPC arrays (`SLURM_ARRAY_TASK_ID` ranges over `1..80000`).

### 4. Aggregate the bootstrap output

`6_Analysis_Boots.R` expects two tab-separated files per setup under
`Boots_Results/`:

- `n_<N>_setup<S>.txt` — long-format table with columns
  `Risk Time iterative matched treatment dataID b Computation_Data Computation_Fit`.
- `n_<N>_setup<S>_true.txt` — output of step 2.

The script consolidates the bootstrap means, standard deviations, and three
flavours of 90% confidence intervals (percentile, HPD, asymptotic) and renders
the boxplot figures included in the paper.

If you ran step 3, the `.RDS` files produced by `4_dataFittingBootstrap.R`
must first be combined into the expected tabular layout, e.g.:

```r
library(dplyr); library(purrr)
files <- list.files('results_bootstrap', full.names = TRUE)
rst <- map_dfr(files, function(f) {
  meta <- regmatches(f, regexec(
    'n_(\\d+)_setup(\\d+)_J(\\d+)_iter(\\d+)_matched(\\d+)_treat(\\d+)_simu(\\d+)', f))[[1]]
  obj <- readRDS(f)
  map_dfr(seq_along(obj$boot_results), function(b) {
    data.frame(Risk             = obj$boot_results[[b]]$risks,
               Time             = seq_along(obj$boot_results[[b]]$risks),
               iterative        = as.integer(meta[5]),
               matched          = as.integer(meta[6]),
               treatment        = as.integer(meta[7]),
               dataID           = as.integer(meta[8]),
               b                = b,
               Computation_Data = as.numeric(obj$boot_results[[b]]$time_data %||% NA),
               Computation_Fit  = as.numeric(obj$boot_results[[b]]$time_fit))
  })
})
dir.create('Boots_Results', showWarnings = FALSE)
write.table(rst, 'Boots_Results/n_30000_setup5.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')
```

Then:

```sh
GFORM_SETUP=5 Rscript 6_Analysis_Boots.R
```

## Verifying against the paper

With the aggregated `Boots_Results/n_30000_setup5.txt` in place, the
following snippet reproduces the J = 5 and Complete columns of Table 4 in the
paper to two decimal places:

```r
library(dplyr)
rst <- read.table('Boots_Results/n_30000_setup5.txt', header = TRUE)
rst %>%
  group_by(treatment, iterative, matched, Time, dataID) %>%
  summarise(Boots_mean = mean(Risk), Boots_sd = sd(Risk), .groups = 'drop') %>%
  group_by(treatment, iterative, matched, Time) %>%
  summarise(mean_mean = round(mean(Boots_mean) * 100, 2),
            mean_sd   = round(mean(Boots_sd)   * 100, 2),
            .groups   = 'drop')
```

The same `.txt` file also drives the timing tables (`Computation_Fit`,
`Computation_Data` columns).

## What is *not* in this repository

To keep the repository minimal, the following are intentionally omitted; see
the paper for details:

- **SuperLearner pipeline** used in Table 3 of the paper.
- **Real-world VHA analysis** (Section 4 of the paper). The veteran cohort is
  governed by a data-use agreement and cannot be redistributed.
- **Aggregated bootstrap output** (the `Boots_Results/` directory). Regenerate
  it by following the pipeline above.

## License

MIT. See `LICENSE`.
