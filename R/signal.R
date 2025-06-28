#' @title Calculate Panel Signals for iSensor Object
#'
#' @description 
#' Computes aggregated signals (mean/median) for each gene panel in an iSensor object.
#' Supports normalization options and parallel processing for large datasets.
#'
#' @param iSensor_obj An iSensor object created by `create_iSensor`
#' @param panels Character vector of panel names to analyze (NULL uses all panels)
#' @param transform Aggregation method: 'mean' or 'median'
#' @param normed Logical indicating whether to normalize signals (default: FALSE)
#' @param normBy Normalization direction: 'cols' (by samples) or 'rows' (by genes)
#' @param metaPanels List specifying meta-panels to create (see Details)
#' @param doParallel Logical indicating whether to use parallel processing (default: FALSE)
#' @param nCores Number of cores for parallel processing (NULL uses available-1)
#'
#' @details
#' The `metaPanels` parameter should be a named list where each element contains:
#' \itemize{
#'   \item `srcPanels`: Vector of source panel names
#'   \item `rule`: Function to combine signals (e.g., `prod`, `sum`)
#' }
#'
#' @return The modified iSensor object with computed signals stored in `Signals` slot
#'
#' @examples
#' \dontrun{
#' # Calculate mean signals
#' iSensor_obj <- iSensor_signal(iSensor_obj, transform = "mean")
#'
#' # Calculate normalized median signals with parallel processing
#' iSensor_obj <- iSensor_signal(iSensor_obj, 
#'                              transform = "median",
#'                              normed = TRUE,
#'                              doParallel = TRUE)
#'
#' # With custom meta-panel
#' iSensor_obj <- iSensor_signal(iSensor_obj,
#'                              transform = "mean",
#'                              metaPanels = list(
#'                                "MyMeta" = list(
#'                                  srcPanels = c("Panel1", "Panel2"),
#'                                  rule = prod
#'                                )
#'                              ))
#' }
#'
#' @export
iSensor_signal <- function(iSensor_obj, panels = NULL, transform, normed = FALSE, normBy = 'cols',
                           metaPanels = NULL, doParallel = FALSE, nCores = NULL) {
  if (is.null(iSensor_obj$exprData)) {
    stop("Error: iSensor expr data is empty")
  }
  
  if (is.null(panels)) {
    panelsNames <- names(iSensor_obj$genePanels)
  } else {
    panelsNames <- panels
  }
  
  # Mean or median, normed or not normed
  iSensor_mean_median <- function(expression_data, detectors, method, norm = normed, normOpt = normBy) {
    
    # отбираем строки в полной матрице с ненулевой дисперсией
    # non_zero_variance_rows <- apply(expression_data, 1, var) > 0
    # expression_data <- expression_data[non_zero_variance_rows, , drop = FALSE]
    
    if (normed) {
      if (method == "mean") {
        if (normOpt == 'cols') {
          expression_data <- expression_data / colMeans(expression_data)
        } else if (normOpt == 'rows') {
          expression_data <- expression_data / rowMeans(expression_data)
        } else {
          stop("Error: for iSensor signal norm 'normBy' can be 'cols' or 'rows'")
        }
      } else if (method == "median") {
        if (normOpt == 'cols') {
          # expression_data <- expression_data / apply(expression_data, 2, median)
          expression_data <- expression_data / apply(expression_data, 2, function(x) median(x[x != 0], na.rm = TRUE))
        } else if (normOpt == 'rows') {
          # expression_data <- expression_data / apply(expression_data, 1, median)
          expression_data <- expression_data / apply(expression_data, 1, function(x) median(x[x != 0], na.rm = TRUE))
        } else {
          stop("Error: for iSensor signal norm 'normBy' can be 'cols' or 'rows'")
        }
      } else {
        stop("Error: for iSensor signal 'transform' can be 'mean' or 'median'")
      }
    }
    
    # подсписок генов из списка
    matched_genes <- intersect(rownames(expression_data), detectors$Genes)
    geneLogicVec <- rownames(expression_data) %in% matched_genes
    filtered_detector_matrix <- expression_data[geneLogicVec, ]
    if(sum(geneLogicVec) == 1) {
      filtered_detector_matrix <- t(as.matrix(filtered_detector_matrix))
    }
    
    # Compute mean or median for each sample
    if (method == "mean") {
      aggregate_signal <- colMeans(filtered_detector_matrix)
    } else if (method == "median") {
      # aggregate_signal <- apply(filtered_detector_matrix, 2, median)
      aggregate_signal <- apply(filtered_detector_matrix, 2, function(x) median(x[x != 0], na.rm = TRUE))
    }
    
    return(aggregate_signal)
  }
  
  if (transform == 'mean' || transform == 'median') {
    if(is.null(iSensor_obj$metaData)) {
      non_zero_variance_rows <- apply(iSensor_obj$exprData, 1, var) > 0
      iSensor_obj$exprData <- iSensor_obj$exprData[non_zero_variance_rows, , drop = FALSE]
      iSensor_obj$metaData <- 'obtained'
    }
    # panelsNames <- names(iSensor_obj$genePanels) # delete ?
    if (doParallel) {
      library(future.apply)
      if (is.null(nCores)) {
        nCores <- availableCores() - 1
      } else {
        if (nCores >= availableCores()) {
          cat('nCores argument is larger than number of available Cores, set to `availableCores() - 1`', '\n')
          nCores <- availableCores() - 1
        }
      }
      
      plan(multisession, workers = nCores)  # Параллельный режим
      signalData <- suppressWarnings({future_lapply(panelsNames, function(x) {
        panelSignal <- iSensor_mean_median(
          expression_data = iSensor_obj$exprData,
          detectors = iSensor_obj$genePanels[[x]],
          method = transform, normOpt = normBy)
        return(panelSignal)
      })                          
      })
      # Если нужна матрица, а не список:
      signalData <- do.call(cbind, signalData)
    } else {
      pb <- txtProgressBar(min = 0, max = length(panelsNames), style = 3, width = 50)
      signalData <- lapply(panelsNames, function(x) {
        panelSignal <- iSensor_mean_median(expression_data = iSensor_obj$exprData,
                                           detectors = iSensor_obj$genePanels[[x]],
                                           method = transform, normOpt = normBy)
        setTxtProgressBar(pb, which(panelsNames == x))
        return(panelSignal)
      })
      # names(signalData) <- panelsNames
      # signalData <- as.data.frame(signalData)
      close(pb)
    }
    signalData <- as.data.frame(signalData)
    # signalData <- do.call(cbind, signalData)
    colnames(signalData) <- panelsNames
    
    if (!is.null(metaPanels)) {
      for (newPanel in names(metaPanels)) {
        metaPanel <- apply(signalData[, metaPanels[[newPanel]][['srcPanels']]], 1, metaPanels[[newPanel]][['rule']])
        signalData[[newPanel]] <- metaPanel
        signalData <- signalData[, c(newPanel, setdiff(names(signalData), newPanel))]
      }
    }
    
    if (normed) {
      transform_name <- paste0(transform, '_normed')
      if (normBy == 'rows') {
        transform_name <- paste0(transform, '_by_', normBy)
      }
    } else {
      transform_name <- transform
    }
    
    if (transform_name %in% names(iSensor_obj$Signals)) {
      iSensor_obj$Signals[[transform_name]] <- signalData
    } else {
      iSensor_obj$Signals <- c(iSensor_obj$Signals, setNames(list(signalData), transform_name))
    }
  }
  
  iSensor_obj$Signals <- iSensor_obj$Signals[order(names(iSensor_obj$Signals))]
  cat(transform_name, 'signal was calculated\n')
  
  return(iSensor_obj)
}
