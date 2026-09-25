#' Calculate Sensor Signals from Expression Data
#'
#' Computes sensor signals based on predefined gene panels from expression data.
#' Supports input as Seurat objects or expression matrices, and calculates the
#' mean or median expression of the panel genes in each cell.
#'
#' @param data A \code{Seurat} object or numeric expression matrix (genes x cells).
#' @param seurLayer Character. Layer of the default assay of a Seurat object to
#'   read expression from. Default \code{"data"} (log-normalised expression).
#' @param panelSet An \code{iSensorsPanelSet} object containing gene panels to calculate signals for,
#'   see \code{\link{LoadSensors}}.
#' @param signals Character vector. Types of signals to compute: \code{"mean"}
#'   (default) and/or \code{"median"}. The \code{"mean_normed"} and
#'   \code{"median_normed"} signals of earlier versions were removed in 1.3.0.
#'
#' @return If input is a Seurat object, returns the same Seurat object with added assays for each signal.
#' If input is a numeric matrix, returns an \code{iSensors} object with calculated signals.
#'
#' @details
#' For each gene panel and cell, the signal is the mean (or the median of the
#' non-zero values) of the expression of the panel genes. Genes with zero
#' variance across all cells are left out, and panels with fewer than 3 remaining
#' genes are skipped with a warning.
#'
#' Signals are stored in assays within Seurat objects named as \code{"iSensors_<signal>"},
#' e.g. \code{iSensors_mean}.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Using a Seurat object
#' seurat_obj <- CreateSeuratObject(counts = Read10X(data.dir = "path/to/data"))
#' seurat_obj <- NormalizeData(seurat_obj)
#' panelSet <- LoadSensors(setName = "ArabidopsisAuxin", species = "ATH", hormone = "aux")
#'
#' seurat_obj <- CalcSensors(seurat_obj, panelSet = panelSet, signals = c("mean", "median"))
#'
#' # Access sensor signal assay
#' head(GetAssayData(seurat_obj, assay = "iSensors_mean", layer = "counts"))
#'
#' # Using an expression matrix
#' expr_mat <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", layer = "data"))
#' iSensor_obj <- CalcSensors(expr_mat, panelSet = panelSet)
#' }
#'
#' @export
CalcSensors <- function(data, seurLayer = 'data', panelSet, signals = "mean") {
  
  # Allowed signals
  removed_signals <- c("mean_normed", "median_normed")
  if (any(signals %in% removed_signals)) {
    stop('The signals "mean_normed" and "median_normed" were removed in iSensors 1.3.0. ',
         'Use signals = "mean" or "median".')
  }
  allowed_signals <- c("mean", "median")
  if (!all(signals %in% allowed_signals)) {
    stop("Error: allowed signals are: ", paste(allowed_signals, collapse = ", "))
  }
  
  # Seurat object: signals are added as new assays
  if (inherits(data, "Seurat")) {
    exprData <- as.matrix(Seurat::GetAssayData(data, layer = seurLayer))
    exprData[is.nan(exprData)] <- 0
    iSensor_obj <- create_iSensors(data = exprData, panelSet = panelSet)
    
    for (signal in signals) {
      iSensor_obj <- iSensor_signal(iSensor_obj, transform = signal)
      assayName <- paste0("iSensors_", signal)
      newAssay <- Seurat::CreateAssayObject(counts = as(iSensor_obj$signals[[signal]], "dgCMatrix"))
      data[[assayName]] <- newAssay
    }
    
    return(data)
    
  } else if (is.matrix(data) || inherits(data, "dgCMatrix")) {
    iSensor_obj <- create_iSensors(data = data, panelSet = panelSet)
    for (signal in signals) {
      iSensor_obj <- iSensor_signal(iSensor_obj, transform = signal)
    }
    return(iSensor_obj)
    
  } else {
    stop("Error: `data` must be either a Seurat object or an expression matrix")
  }
}
