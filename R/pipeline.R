#' Calculate Sensor Signals from Expression Data
#'
#' Computes sensor signals based on predefined gene panels from expression data.
#' Supports input as Seurat objects or raw expression matrices, and calculates multiple
#' summary signals including mean, median, and their normalized versions.
#'
#' @param data A \code{Seurat} object or numeric expression matrix (genes x samples).
#' @param seurLayer Character. Assay layer name to extract data from Seurat object (e.g., "RNA"). Default "RNA".
#' @param panelSet An \code{iSensorsPanelSet} object containing gene panels to calculate signals for.
#' @param signals Character vector. Types of signals to compute. Allowed values:
#'   \code{"mean"}, \code{"mean_normed"}, \code{"median"}, \code{"median_normed"}.
#'   Default is \code{"mean_normed"}.
#'
#' @return If input is a Seurat object, returns the same Seurat object with added assays for each signal.
#' If input is a numeric matrix, returns an \code{iSensors} object with calculated signals.
#'
#' @details
#' The function computes sensor signals for each gene panel by summarizing gene expression values
#' across panel genes in each sample or cell. Normalization is done by dividing gene expression by
#' the mean or median expression per sample or gene, depending on the signal type.
#'
#' Signals are stored in assays within Seurat objects named as \code{"iSensors_<signal>"},
#' e.g. \code{iSensors_mean_normed}.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' 
#' # Using Seurat object
#' seurat_obj <- Read10X(data.dir = "path/to/data")
#' seurat_obj <- CreateSeuratObject(counts = seurat_obj)
#' panelSet <- LoadSensors(setName = "ArabidopsisAuxin")
#' 
#' seurat_obj <- CalcSensors(seurat_obj, seurLayer = "RNA", panelSet = panelSet,
#'                           signals = c("mean_normed", "median"))
#' 
#' # Access sensor signal assay
#' head(seurat_obj@assays$iSensors_mean_normed@counts)
#' 
#' # Using raw matrix
#' expr_mat <- as.matrix(GetAssayData(seurat_obj, slot = "counts"))
#' iSensor_obj <- CalcSensors(expr_mat, panelSet = panelSet,
#'                            signals = c("mean", "median_normed"))
#' }
#'
#' @export
CalcSensors <- function(data, seurLayer = 'RNA', panelSet,
                        signals = c("mean_normed")) {
  
  # Проверка допустимых сигналов
  allowed_signals <- c("mean", "mean_normed", "median", "median_normed")
  if (!all(signals %in% allowed_signals)) {
    stop("Error: allowed signals are: ", paste(allowed_signals, collapse = ", "))
  }
  
  # Если объект Seurat
  if (inherits(data, "Seurat")) {
    exprData <- as.matrix(Seurat::GetAssayData(data, layer = seurLayer))
    exprData[is.nan(exprData)] <- 0
    iSensor_obj <- create_iSensors(data = exprData, panelSet = panelSet)
    
    for (signal in signals) {
      transform <- if (grepl("mean", signal)) "mean" else "median"
      normed <- grepl("normed", signal)
      iSensor_obj <- iSensor_signal(iSensor_obj, transform = transform, normed = normed)
      
      # signalName <- names(iSensor_obj$signals)[[length(iSensor_obj$signals)]]
      signalName <- signal
      assayName <- paste0("iSensors_", signalName)
      newAssay <- Seurat::CreateAssayObject(counts = as(iSensor_obj$signals[[signalName]], "dgCMatrix"))
      data[[assayName]] <- newAssay
    }
    
    # if (!is.null(defaultAssay)) {
    #   defaultAssayName <- paste0("iSensor_", defaultAssay)
    #   if (!defaultAssayName %in% names(data@assays)) {
    #     warning("Default assay name not found, skipping DefaultAssay assignment")
    #   } else {
    #     Seurat::DefaultAssay(data) <- defaultAssayName
    #     cat("Set DefaultAssay to", defaultAssayName, "\n")
    #   }
    # }
    
    return(data)
    
  } else if (is.matrix(data) || inherits(data, "dgCMatrix")) {
    iSensor_obj <- create_iSensors(data = data, panelSet = panelSet)
    for (signal in signals) {
      transform <- if (grepl("mean", signal)) "mean" else "median"
      normed <- grepl("normed", signal)
      iSensor_obj <- iSensor_signal(iSensor_obj, transform = transform, normed = normed)
    }
    return(iSensor_obj)
    
  } else {
    stop("Error: `data` must be either a Seurat object or an expression matrix")
  }
}