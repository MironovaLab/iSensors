devtools::load_all()
library(Seurat)
seurat_obj <- readRDS("D:/FILES/work/Sensor/ISensors_data/ra_integrated_Shahan_20250114.rds")
# 2. Размеры
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n")

# 3. Доступные слоты
cat("\nAssays:\n")
print(Assays(seurat_obj))

cat("\nReductions:\n")
print(Reductions(seurat_obj))

cat("\nMetadata columns:\n")
print(colnames(seurat_obj@meta.data))

cat("\nVariable features (first 10):\n")
print(head(VariableFeatures(seurat_obj), 10))

# 4. Посмотрим первые 5 строк метаданных
cat("\nPreview of meta.data:\n")
print(head(seurat_obj@meta.data, 5))

# 5. Посмотрим имена слоёв внутри assay "RNA" (или другого, если не RNA)
default_assay <- DefaultAssay(seurat_obj)
cat("\nDefault assay:", default_assay, "\n")

cat("\nAvailable slots in default assay:\n")
print(slotNames(seurat_obj[[default_assay]]))

# 2. Извлекаем counts-матрицу из RNA-ассая (Seurat v5)
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

# 3. Сэмплируем 1000 случайных клеток
set.seed(42)
cells_to_keep <- sample(colnames(counts_mat), 1000)

# 4. Обрезаем counts-матрицу (по клеткам и всем генам)
sub_counts <- counts_mat[, cells_to_keep]

# 5. Создаём новый Seurat-объект
sub_obj <- CreateSeuratObject(counts = sub_counts)

# 6. Препроцессинг: нормализация, HVG, скейлинг
sub_obj <- NormalizeData(sub_obj)
sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
sub_obj <- ScaleData(sub_obj)

# 7. PCA и UMAP
sub_obj <- RunPCA(sub_obj, features = VariableFeatures(sub_obj))
sub_obj <- RunUMAP(sub_obj, dims = 1:10)

saveRDS(object = sub_obj, file = paste0('data/', 'testSeurData.rds'))
