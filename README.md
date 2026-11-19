# Trophic reorganization from predator stomach contents

This repository contains the R code and analytical outputs associated with the study *A framework for diagnosing trophic reorganization in marine systems using predator stomach contents*.

The study examines how predator diets change over time by combining three complementary measures: prey composition (Bray-Curtis dissimilarity and PERMANOVA), Shannon diversity (H′), and standardized Levins niche breadth (Bs). Together, these measures are used to distinguish four patterns of dietary change: **Stability, Substitution, Functional change, and Reorganization**.

The framework is applied to predator stomach-content data from the southern Gulf of St. Lawrence, with analyses conducted at three spatial scales: the sGulf as a whole, ecoregions, and survey strata.

## Analysis scripts

The analyses are organized into 12 R scripts, covering data preparation, statistical testing, robustness assessments, and the production of tables and figures.

| Script | Description |
|---|---|
| `01_Extract_and_Clean.R` | Extraction and cleaning of stomach-content records and prey taxonomy |
| `02_Complete_Preparation.R` | Preparation of prey groups, predator size classes, and analytical datasets |
| `03_Statistical_Pipeline.R` | Core statistical functions used throughout the analyses |
| `04_Run_Main_Analyses.R` | Main temporal comparisons across the three spatial scales |
| `05_Run_Grouping_Sensitivity.R` | Sensitivity analyses using alternative prey-grouping approaches |
| `06_Robustness_PostHoc.R` | Post-hoc robustness analyses (M0-M6) |
| `07_Robustness_Resampling.R` | Sampling-effort sensitivity, balanced resampling, null simulations, and calibration |
| `08_Intermediate_Outputs.R` | Preparation of intermediate summaries and figure data |
| `09_Tables.R` | Main and supplementary tables |
| `10_Figures.R` | Main manuscript figures |
| `11_Supplementary.R` | Supplementary robustness figures |
| `12_Other_Appendices.R` | Additional tables, figures, and exploratory analyses |

Additional functions shared across scripts are available in `R_helpers/`.

## Running the analyses

The project can be opened in RStudio using `Trophic_Reorganisation_StomachData.Rproj`.

Scripts are numbered in the order of the analytical workflow, from data preparation through the final tables and figures. `03_Statistical_Pipeline.R` contains functions called by subsequent scripts and does not need to be executed separately.

Some analyses, particularly resampling and permutation-based procedures, are computationally intensive. Saved intermediate results can be reused to avoid repeating these calculations, provided they correspond to the same data, code, and analytical settings.

## Data availability and reproducibility

The original stomach-content survey data were collected by Fisheries and Oceans Canada (DFO) and are not included in this repository because their distribution is subject to institutional restrictions. The initial data preparation steps also rely on DFO-specific resources, including spatial datasets and software tools.

The repository includes R scripts and analytical outputs organized by analysis type, including sensitivity assessments, robustness analyses, sampling-effort comparisons, calibration, and manuscript tables and figures.

These intermediate results support the reproduction of downstream analyses, tables, and figures without requiring access to the original stomach-content records. Reproducing the complete workflow from the initial data preparation stage, however, requires access to the underlying DFO datasets and associated resources.

Access to the original survey data may be requested through Jenni McDermid at **Jenni.McDermid@dfo-mpo.gc.ca**, subject to institutional approval.
