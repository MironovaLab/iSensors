# iSensors

**iSensors** is a package for analyzing signaling activity scores of gene panels based on single-cell RNA-seq data stored in Seurat objects. The package allows you to compute signaling activity scores for specified gene panels and store the results as new assays.

## Installation

You can install the package from GitHub using `devtools`:

```R
# Make sure you have devtools installed
install.packages("devtools")

# Install the iSensor package
devtools::install_github("MironovaLab/iSensors")
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

## Running iSensors via default panels
```R
library(iSensors)
library(Seurat)
in_path <- "D://FILES/work/Sensor/ISensors_data/Auxin_transcriptoms/"
seurat_obj <- readRDS(paste0(in_path, "seurat_auxin_5bulkRNA-Seq.rds"))
seurat_obj <- iSensor_pipeline(seuratObject=seurat_obj, species = 'AT', hormone = 'aux', type = c('cis','trans'))
print(names(seurat_obj@assays))
seurat_obj@assays$iSensor_mean_normed@counts[1:10, 1:3]
```
