# iSensors <img src="https://img.shields.io/badge/R-package-blue" alt="R package" height="24">

**iSensors** is a package for analyzing signaling activity scores of gene panels based on single-cell RNA-seq data stored in Seurat objects. The package allows you to compute signaling activity scores for specified gene panels and store the results as new assays.

## Installation

You can install the package from GitHub using `devtools`:

```R
# Make sure you have devtools installed
install.packages("devtools")

# Install the iSensor package from main
devtools::install_github("MironovaLab/iSensors")
# or from other branch
devtools::install_github("MironovaLab/iSensors", ref = "devel")
```
In case dealing with a mistake like:
```
Failed to install 'unknown package' from GitHub
```
You need to go to GitHub Settings -> Developer Settings and generate your token.
Run
```
usethis::edit_r_environ()
```
and copy your token to the GITHUB_PAT= your token/

## Quick Start

Here's an example of how to load the package, check its version, explore documentation, and run a basic analysis:

```R
# Load the package
library(iSensors)

# Check installed version
packageVersion("iSensors")

# Access general help page
help("iSensors")

# Load test data (Seurat object)
testData <- readRDS("testData/testSeurData.rds")

# Load test panel with meta-panels
testPanel <- LoadSensors(setName = 'testPanelSet', species = 'AT', hormone = 'cyt', customPanels = TRUE,
                          randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), majortrend = TRUE),
                          metaPanels = list(
                            'meta1' = list('srcPanels' = c("AT_aux_cis_DR5_ARF1", "AT_aux_cistrans_DR5_ARF5_2_up"), rule = mean),
                            'meta2' = list('srcPanels' = c("AT_aux_cis_DR5_ARF1", "AT_aux_cistrans_DR5_ARF5_2_up"), rule = prod))
                          )

# Calculate signaling scores for selected panels
result <- CalcSensors(testData,
                      seurLayer = "data",
                      panelSet = testPanel,
                      signals = c("mean_normed", "median"))
```