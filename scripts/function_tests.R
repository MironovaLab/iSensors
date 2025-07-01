devtools::load_all()
compute_signal <- function(exprData, detectors, method, normed, normBy) {
  if (normed) {
    if (method == "mean") {
      exprData <- if (normBy == "cols") exprData / colMeans(exprData) else exprData / rowMeans(exprData)
    } else if (method == "median") {
      median_func <- function(x) median(x[x != 0], na.rm = TRUE)
      exprData <- if (normBy == "cols") exprData / apply(exprData, 2, median_func) else exprData / apply(exprData, 1, median_func)
    }
  }
  
  matchedGenes <- intersect(rownames(exprData), detectors$Genes)
  filteredMatrix <- exprData[matchedGenes, , drop = FALSE]
  if (nrow(filteredMatrix) == 1) {
    filteredMatrix <- t(as.matrix(filteredMatrix))
  }
  
  if (method == "mean") {
    return(colMeans(filteredMatrix))
  } else {
    return(apply(filteredMatrix, 2, function(x) median(x[x != 0], na.rm = TRUE)))
  }
}
# Тестовая матрица экспрессии (5 генов, 3 образца)
exprData <- matrix(c(
  1, 2, 3,    # GeneA
  4, 5, 6,    # GeneB
  7, 8, 9,    # GeneC
  0, 0, 0,    # GeneD
  10, 11, 12  # GeneE
), nrow = 5, byrow = TRUE)

rownames(exprData) <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")
colnames(exprData) <- c("Sample1", "Sample2", "Sample3")

# Панель генов (только некоторые из них)
detectors <- data.frame(Genes = c("GeneA", "GeneC", "GeneE"))
# 1. Без нормировки, метод mean
compute_signal(exprData, detectors, method = "mean", normed = FALSE, normBy = "cols")

# 2. Без нормировки, метод median
compute_signal(exprData, detectors, method = "median", normed = FALSE, normBy = "cols")

# 3. С нормировкой по столбцам, метод mean
compute_signal(exprData, detectors, method = "mean", normed = TRUE, normBy = "cols")

# 4. С нормировкой по строкам, метод median
testRes <- compute_signal(exprData, detectors, method = "median", normed = TRUE, normBy = "rows")
str(testRes)
