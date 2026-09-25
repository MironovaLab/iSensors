#' @keywords internal
expand_panels_from_additional <- function(exprData, panelsList, additional) {
  if (is.null(additional)) return(panelsList)
  
  # --- RANDOM PANELS ---
  if (!is.null(additional$random)) {
    random <- additional$random
    
    if (!all(c("n", "sizes") %in% names(random))) {
      warning("random must contain 'n' and 'sizes'")
    } else {
      gene_pool <- rownames(exprData)
      n <- random$n
      sizes <- random$sizes
      if (length(sizes) != n) {
        stop("Length of 'sizes' must be equal to 'n'")
      }
      
      for (i in seq_len(n)) {
        size <- sizes[i]
        if (size > length(gene_pool)) {
          warning(paste0("Not enough genes to sample size ", size, "; skipping random", i))
          next
        }
        panelsList[[paste0("random", i)]] <- sample(gene_pool, size)
      }
      
      if (!is.null(random$majortrend) && isTRUE(random$majortrend)) {
        panelsList[["majortrend"]] <- gene_pool
      }
    }
  }
  
  return(panelsList)
}

#' @keywords internal
iSensor_signal <- function(iSensor_obj, transform = "mean") {
  # Checks
  if (is.null(iSensor_obj$exprData)) stop("Error: iSensor exprData is empty")
  if (is.null(iSensor_obj$panelSet)) stop("Error: iSensor panelSet is missing")

  if (!transform %in% c("mean", "median")) {
    stop("transform must be 'mean' or 'median'")
  }

  # Panels, including random and meta panels
  panelsList <- iSensor_obj$panelSet$panels
  additional <- iSensor_obj$panelSet$additional
  panelsList <- expand_panels_from_additional(iSensor_obj$exprData, panelsList, additional)
  
  panelNames <- names(panelsList)
  
  # Keep only genes with non-zero variance across cells (done once per object)
  if (is.null(iSensor_obj$metaData)) {
    nonZeroVarRows <- apply(iSensor_obj$exprData, 1, var) > 0
    iSensor_obj$exprData <- iSensor_obj$exprData[nonZeroVarRows, , drop = FALSE]
    iSensor_obj$metaData <- list(filtered = TRUE)
  }
  safe_colMeans <- function(mat) {
    if (inherits(mat, "dgCMatrix")) {
      library(Matrix)
      return(Matrix::colMeans(mat))
    } else {
      return(colMeans(mat))
    }
  }
  # Signal of one panel in every cell; NULL when too few panel genes are
  # detected (non-zero variance) to give a reliable score
  compute_signal <- function(exprData, detectors, method) {
    matchedGenes <- intersect(rownames(exprData), detectors)
    if (length(matchedGenes) < min_panel_genes) return(NULL)
    filteredMatrix <- exprData[matchedGenes, , drop = FALSE]

    if (method == "mean") {
      return(safe_colMeans(filteredMatrix))
    } else {
      return(apply(filteredMatrix, 2, function(x) median(x[x != 0], na.rm = TRUE)))
    }
  }
  
  # Signals of all panels
  n <- length(panelNames)
  pb <- txtProgressBar(min = 0, max = n, style = 3, width = 50)
  signalList <- lapply(panelNames, function(pName) {
    res <- compute_signal(iSensor_obj$exprData, panelsList[[pName]], transform)
    setTxtProgressBar(pb, which(panelNames == pName))
    return(res)
  })
  close(pb)
  
  names(signalList) <- panelNames
  skipped <- panelNames[vapply(signalList, is.null, logical(1))]
  if (length(skipped) > 0) {
    warning(sprintf("%d panel(s) skipped: fewer than %d of their genes are detected in the data: %s",
                    length(skipped), min_panel_genes, paste(skipped, collapse = ", ")),
            call. = FALSE)
  }
  signalList <- signalList[!vapply(signalList, is.null, logical(1))]
  if (length(signalList) == 0) {
    stop("No panel has at least ", min_panel_genes, " genes detected in the data.", call. = FALSE)
  }
  signalDF <- as.data.frame(do.call(rbind, signalList))
  signalDF <- as.matrix(signalDF)
  # Meta panels: combine source panel signals with their rule, added as new rows
  if (!is.null(additional$meta)) {
    for (metaName in names(additional$meta)) {
      metaInfo <- additional$meta[[metaName]]
      src <- metaInfo$srcPanels
      rule <- metaInfo$rule
      if (!all(src %in% rownames(signalDF))) {
        missing <- setdiff(src, rownames(signalDF))
        warning(paste("Meta panel", metaName, "skipped: missing source panels:", paste(missing, collapse = ", ")))
        next
      }
      
      new_row <- apply(signalDF[src, , drop = FALSE], 2, rule)
      signalDF <- rbind(signalDF, new_row)
      rownames(signalDF)[nrow(signalDF)] <- metaName
    }
  }
  
  # Store the signal in the object, named "mean" or "median"
  if (is.null(iSensor_obj$signals)) {
    iSensor_obj$signals <- list()
  }
  iSensor_obj$signals[[transform]] <- signalDF
  iSensor_obj$signals <- iSensor_obj$signals[order(names(iSensor_obj$signals))]

  message(transform, " signal was calculated")
  return(iSensor_obj)
}
