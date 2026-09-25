#' Calculate Sensor Signals from Expression Data
#'
#' Computes sensor signals based on predefined gene panels from expression data.
#' Supports input as Seurat objects or raw expression matrices, and calculates multiple
#' summary signals including mean, median, and their normalized versions.
#'
#' @param data A \code{Seurat} object or numeric expression matrix (genes x cells).
#' @param seurLayer Character. Layer of the default assay of a Seurat object to
#'   read expression from. Default \code{"data"} (log-normalised expression).
#' @param panelSet An \code{iSensorsPanelSet} object containing gene panels to calculate signals for,
#'   see \code{\link{LoadSensors}}.
#' @param signals Character vector. Types of signals to compute. Allowed values:
#'   \code{"mean"}, \code{"mean_normed"}, \code{"median"}, \code{"median_normed"}.
#'   Default is \code{"mean_normed"}.
#' @param normBy Character. For the \code{"_normed"} signals, divide expression by the
#'   mean or median of each cell (\code{"cols"}, default) or of each gene (\code{"rows"}).
#'
#' @return If input is a Seurat object, returns the same Seurat object with added assays for each signal.
#' If input is a numeric matrix, returns an \code{iSensors} object with calculated signals.
#'
#' @details
#' The function computes sensor signals for each gene panel by summarizing gene expression values
#' across panel genes in each cell. Genes with zero variance across all cells are left out.
#' Normalization is done by dividing gene expression by the mean or median expression per cell
#' or gene, depending on \code{normBy}.
#'
#' Signals are stored in assays within Seurat objects named as \code{"iSensors_<signal>"},
#' e.g. \code{iSensors_mean_normed}.
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
#' iSensor_obj <- CalcSensors(expr_mat, panelSet = panelSet, signals = "mean")
#' }
#'
#' @export
CalcSensors <- function(data, seurLayer = 'data', panelSet,
                        signals = c("mean_normed"), normBy = "cols") {
  
  # Allowed signals
  allowed_signals <- c("mean", "mean_normed", "median", "median_normed")
  if (!all(signals %in% allowed_signals)) {
    stop("Error: allowed signals are: ", paste(allowed_signals, collapse = ", "))
  }
  
  # Seurat object: signals are added as new assays
  if (inherits(data, "Seurat")) {
    exprData <- as.matrix(Seurat::GetAssayData(data, layer = seurLayer))
    exprData[is.nan(exprData)] <- 0
    iSensor_obj <- create_iSensors(data = exprData, panelSet = panelSet)
    
    for (signal in signals) {
      transform <- if (grepl("mean", signal)) "mean" else "median"
      normed <- grepl("normed", signal)
      iSensor_obj <- iSensor_signal(iSensor_obj, transform = transform, normed = normed, normBy = normBy)
      
      signalName <- signal
      assayName <- paste0("iSensors_", signalName)
      newAssay <- Seurat::CreateAssayObject(counts = as(iSensor_obj$signals[[signalName]], "dgCMatrix"))
      data[[assayName]] <- newAssay
    }
    
    return(data)
    
  } else if (is.matrix(data) || inherits(data, "dgCMatrix")) {
    iSensor_obj <- create_iSensors(data = data, panelSet = panelSet)
    for (signal in signals) {
      transform <- if (grepl("mean", signal)) "mean" else "median"
      normed <- grepl("normed", signal)
      iSensor_obj <- iSensor_signal(iSensor_obj, transform = transform, normed = normed, normBy = normBy)
    }
    return(iSensor_obj)
    
  } else {
    stop("Error: `data` must be either a Seurat object or an expression matrix")
  }
}
