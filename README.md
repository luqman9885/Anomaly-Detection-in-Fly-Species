# Anomaly-Detection-in-Fly-Species
## 📁 Repository Structure & Files 
This repository is organised into two folders: 

- **Data**
  - The datasets consists of 14 fly species across 3 families. The variables are the elliptic Fourier Transform (EFT) coefficients extracted from two wing cell compartments - discal medial (dm) cell and second post-anterior radial (pa2r) cell. 
      - h10_dm_final_norm.csv : This file contains the normalised EFT coefficients (10 harmonics) extracted from the dm cell compartment. 
      - h10_pa2r_final_norm.csv : This file contains the normalised EFT coefficients (10 harmonics) extracted from the pa2r cell compartment.

- **Source Code** 
  - The anomaly detection was performed at two hierarchical levels - species level and family level.
      - Family-level.R : This R file contains the implemented codes for the family-level. 
      - Species-level.R : This R file contains the implemented codes for the species-level.

## 💻 Environment & Execution Guide

### R Environment
The analysis was conducted using **R software (version 4.4.0)**. The three main packages that were used:

- `tidyverse` (v2.0.0): Data wrangling.
- `MASS` (v7.3.60.2): Perform linear discriminant analysis (LDA).
- `mvtnorm` (v1.3.1): Fit multivariate normal and multivariate t-distributions.

To install the necessary packages, run the following command in your R console:
```r
install.packages(c("tidyverse", "MASS", "mvtnorm"))

