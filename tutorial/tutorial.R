devtools::load_all()
library(Seurat)

packageVersion("iSensors")

seurat_obj <- readRDS("tutorial/testData/testDataClean.rds")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n")

cat("\nAssays:\n")
print(Assays(seurat_obj))

cat("\nReductions:\n")
print(Reductions(seurat_obj))

cat("\nMetadata columns:\n")
print(colnames(seurat_obj@meta.data))

cat("\nVariable features (first 10):\n")
print(head(VariableFeatures(seurat_obj), 10))

cat("\nPreview of meta.data:\n")
print(head(seurat_obj@meta.data, 5))

default_assay <- DefaultAssay(seurat_obj)
cat("\nDefault assay:", default_assay, "\n")

cat("\nAvailable slots in default assay:\n")
print(slotNames(seurat_obj[[default_assay]]))

testPanel <- LoadSensors(setName = 'testPanelSet', species = 'AT', hormone = 'cyt', customPanels = FALSE,
                         randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), majortrend = TRUE),
                         metaPanels = list(
                           'meta1' = list('srcPanels' = c("AT_cyt_cis_ARR1_1", "AT_cyt_cistrans_ARR1_1down"), rule = mean),
                           'meta2' = list('srcPanels' = c("AT-cyt-cis-ARR1-1", "AT-cyt-cistrans-ARR1-1down"), rule = prod))
)

result <- CalcSensors(seurat_obj,
                      seurLayer = 'data',
                      panelSet = testPanel,
                      signals = c("mean_normed", "median"))
View(result)
rownames(result[['iSensors_median']])
