#' iSensor: Gene Panel Analysis for Single-Cell RNA-Seq Data
#'
#' The iSensor package provides tools for analyzing signaling pathway activity
#' in single-cell RNA-seq data using predefined gene panels. It calculates
#' pathway activity scores and integrates results with Seurat objects.
#'
#' @section Key Features:
#' \itemize{
#'   \item Integration with Seurat objects
#'   \item Multiple signal calculation methods (mean, median, normalized)
#'   \item Support for custom gene panels
#'   \item Background modeling with random panels
#'   \item Parallel processing support
#' }
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{iSensor_pipeline}}}{Main analysis pipeline}
#'   \item{\code{\link{create_iSensor}}}{Create iSensor analysis object}
#'   \item{\code{\link{set_iSensor_workPanel}}}{Set active gene panels}
#'   \item{\code{\link{add_defaultPanels}}}{Load default gene panels}
#'   \item{\code{\link{iSensor_signal}}}{Calculate signaling activity}
#' }
#'
#' @section Data Management:
#' \describe{
#'   \item{\code{\link{add_iSensor_data}}}{Add expression data}
#'   \item{\code{\link{add_iSensor_labels}}}{Add sample labels}
#'   \item{\code{\link{add_iSensor_panelSet}}}{Add panel sets}
#'   \item{\code{\link{read_genePanels}}}{Read panel files}
#' }
#'
#' @section Visualization:
#' \describe{
#'   \item{\code{\link{plot_iSensor_signals}}}{Plot pathway activity}
#' }
#'
#' @docType package
#' @name iSensor
NULL