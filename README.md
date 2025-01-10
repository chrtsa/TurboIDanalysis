# TurboIDanalysis
A package that contains functions to create the ROC curves following the data analysis in Cho et al., 2020 (Nat Protoc). It calculates the TPR and FPR based on user-provided true- and false-positive lists, as well as the provided dataset with protein measurements. A previously calculated ratio of the protein measurements (experimental sample/omit-ligase control) is required.

## Features
- Function that flags each protein as either a true or false positive, based on the user-provided lists, and calculates the True Positive and False Positive Rates.
- Function that creates the ROC curve based on the TPR and FPR calculated by the previous function.
- Function that plots the TPR-FPR against the log2 fold change of the protein measurements, allowing for the calculation of the optimal cutoff value.

### Main functions
- flag() - Flags each protein as either a true or false positive, based on the user-provided lists, and calculates the True Positive and False Positive Rates.
- roc() - Creates the ROC curve based on the TPR and FPR calculated by the previous function.

## Installation

Install the package directly from GitHub using `devtools`:

```R
# Install devtools if not already installed
install.packages("devtools")

# Install package from GitHub
devtools::install_github("chrtsa/TurboIDanalysis")
```

## Usage
```R
library(TurboIDanalysis)
# Example: normalise the data and calculate the True Positive and False Positive Rates for all proteins detected in the experiment
# Load the necessary datasets
data <- "path/to/data.csv"
TP_NES <- "path/to/TP_NES.csv"
FP <- "path/to/FP.csv"
# Create new columns with the TPR and FPR of each protein, arranged by descending order according to their levels of enrichment in the experimenta sample compared to the omitted-ligase control
data <- flag(data = data, sample = "NES", TP = TP_NES, FP = FP, gene_name = "gene_name", ratio = "NES", TP_flag = "NES_TP")

# Create the ROC curve for each potential cutoff value
roc(data = data)

# Plot the TPR-FPR against the log2 fold change of the protein measurements to determine the optimal cutoff value
opt_cutoff(data = data)
```

### Reference
Cho, K.F., Branon, T.C., Udeshi, N.D. et al. Proximity labeling in mammalian cells with TurboID and split-TurboID. Nat Protoc 15, 3971–3999 (2020). https://doi.org/10.1038/s41596-020-0399-0