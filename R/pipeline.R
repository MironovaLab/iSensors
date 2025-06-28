#' @title iSensor Analysis Pipeline for Seurat Objects
#'
#' @description 
#' Main pipeline for analyzing single-cell RNA-seq data using iSensor framework.
#' Computes gene panel signatures and stores results as new assays in Seurat object.
#' Supports parallel processing, custom panels, and meta-panel creation.
#'
#' @param seuratObject A Seurat object (v3+) containing single-cell data
#' @param species Species for default panels ("AT", etc.)
#' @param hormone Hormone type for panel filtering (optional)
#' @param type Panel type for panel filtering (optional)
#' @param signals Vector of signals to compute: 'mean', 'mean_normed',
#'        'median', 'median_normed' (default: 'mean_normed')
#' @param seurLayer Seurat assay layer to use ("data", "counts", etc.)
#' @param usePanelPreset Custom panel set to add (overrides species/hormone/type)
#' @param presetPanelName Name for custom panel set (required if usePanelPreset specified)
#' @param randPanels Number of random panels to generate (default: 2)
#' @param randSize Size range for random panels (default: c(200, 500))
#' @param majortrend Whether to include all genes as "majortrend" panel (default: TRUE)
#' @param usePanels Specific panels to analyze (NULL uses all available)
#' @param defaultAssay Name of iSensor assay to set as default (optional)
#' @param useParallel Enable parallel processing (default: FALSE)
#' @param nCores Number of cores for parallel processing (NULL uses available-1)
#' @param metaPanels List of meta-panels to create (see Details)
#'
#' @details
#' Meta-panels should be specified as named lists with:
#' \itemize{
#'   \item \code{srcPanels}: Vector of source panel names
#'   \item \code{rule}: Combination function (e.g., \code{prod}, \code{sum})
#' }
#'
#' @return
#' Modified Seurat object with new assays:
#' \itemize{
#'   \item Each signal type becomes separate assay (names: 'iSensor_[signal]')
#'   \item Original data remains unchanged
#'   \item Meta-panels included if specified
#' }
#'
#' @examples
#' \dontrun{
#' # Basic analysis with default parameters
#' seurat_obj <- iSensor_pipeline(seurat_obj)
#'
#' # Advanced analysis with custom settings
#' seurat_obj <- iSensor_pipeline(
#'   seurat_obj,
#'   species = "human",
#'   signals = c("mean", "median_normed"),
#'   randPanels = 3,
#'   metaPanels = list(
#'     "MyMeta" = list(
#'       srcPanels = c("Panel1", "Panel2"),
#'       rule = function(x) sqrt(mean(x))
#'     )
#'   )
#' )
#' }
#'
#' @seealso
#' \code{\link{create_iSensor}}, \code{\link{iSensor_signal}}
#'
#' @importFrom Seurat GetAssayData CreateAssayObject DefaultAssay
#' @export
iSensor_pipeline <- function(seuratObject, species = NULL, hormone = NULL, type = NULL, signals = c('mean_normed'), seurLayer = "data",
                             usePanelPreset = NULL, presetPanelName = NULL,
                             randPanels = 2, randSize = c(200, 500), majortrend = TRUE,
                             usePanels = NULL, defaultAssay = NULL, useParallel = FALSE, nCores = NULL,
                             metaPanels = NULL) {
  if (class(seuratObject) != 'Seurat') {
    stop("Error: function must be applied to Seurat object")
  }
  
  expression_data <- as.matrix(GetAssayData(seuratObject, layer = seurLayer))
  expression_data[is.nan(expression_data)] <- 0
  iSensor_test <- create_iSensor(data = expression_data, species = species, hormone = hormone, type = type)
  if(!is.null(usePanelPreset)) {
    iSensor_test <- add_iSensor_panelSet(iSensor_test, panelsSet = usePanelPreset, setName = presetPanelName)
    iSensor_test <- set_iSensor_workPanel(iSensor_test, presetPanelName)
  }
  if(is.null(presetPanelName)) {
    presetPanelName <- names(iSensor_test$genePanelSets)[1]
  }
  cat(presetPanelName, ' panel preset will be used', '\n')
  iSensor_test <- set_iSensor_workPanel(iSensor_test, presetPanelName)
  iSensor_test <- add_iSensor_random_panel(iSensor_test, randNum=randPanels, randSize=randSize,
                                           majortrend = majortrend)
  
  for (signal in signals) {
    if (signal == 'mean') {
      iSensor_test <- iSensor_signal(iSensor_test, transform = 'mean', normed = FALSE,
                                     panels = usePanels, doParallel = useParallel, nCores = nCores,
                                     metaPanels = metaPanels)
    }
    if (signal == 'mean_normed') {
      iSensor_test <- iSensor_signal(iSensor_test, transform = 'mean', normed = TRUE,
                                     panels = usePanels, doParallel = useParallel, nCores = nCores,
                                     metaPanels = metaPanels)
    }
    if (signal == 'median') {
      iSensor_test <- iSensor_signal(iSensor_test, transform = 'median', normed = FALSE,
                                     panels = usePanels, doParallel = useParallel, nCores = nCores,
                                     metaPanels = metaPanels)
    }
    if (signal == 'median_normed') {
      iSensor_test <- iSensor_signal(iSensor_test, transform = 'median', normed = TRUE,
                                     panels = usePanels, doParallel = useParallel, nCores = nCores,
                                     metaPanels = metaPanels)
    }
  }
  
  for (assay in names(iSensor_test$Signals)) {
    seuratObject[[paste0('iSensor_', assay)]] <- CreateAssayObject(counts = as(t(iSensor_test$Signals[[assay]]), "dgCMatrix"))
  }
  
  if (!is.null(defaultAssay)) {
    DefaultAssay(seuratObject) <- paste0('iSensor_', defaultAssay)
    cat('Set DefaultAssay to', paste0('iSensor_', defaultAssay, '\n'))
  }
  
  return(seuratObject)
}