# Smallest number of genes a panel may have, when created and when scored
min_panel_genes <- 3L

#' @keywords internal
check_panel_size <- function(genes, panel_name) {
  n <- length(unique(genes))
  if (n < min_panel_genes) {
    stop(sprintf('Panel "%s" has %d gene(s); a panel needs at least %d genes.',
                 panel_name, n, min_panel_genes), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
check_iSensorsPanel_obj <- function(panel_obj) {
  
  if (!('genes' %in% names(panel_obj))) {
    stop('Panel object must have "genes" slot')
  }
  
  if (!('name' %in% names(panel_obj))) {
    stop('Panel object must have "name" slot')
  }
  
  check_iSensors_panel(panel_obj$genes)
  check_panel_size(panel_obj$genes, panel_obj$name)
  return(TRUE)
}

#' @keywords internal
check_iSensors_obj <- function(iSensor_obj, check_exprData = FALSE, check_panelSet = FALSE,
                               check_signals = FALSE, check_sampleLabels = FALSE) {
  # Class check
  if (!inherits(iSensor_obj, "iSensors")) {
    stop("Object is not of class 'iSensors'")
  }
  if (check_exprData) {
    if (is.null(iSensor_obj$exprData)) {
      stop("'exprData' field is NULL in iSensors object.")
    }
  }
  if (check_panelSet) {
    if (is.null(iSensor_obj$panelSet)) {
      stop("'panelSet' field is NULL in iSensors object.")
    }
  }
  
  return(TRUE)
}

#' @keywords internal
check_iSensors_panel <- function(panel, print_massage = FALSE) {
  if (!is.vector(panel)) {
    if (print_massage) {
      stop(print_massage)
    } else {
      stop('Panel slot "genes" must be a vector of gene names')
    }
  }
  return(TRUE)
}

#' @keywords internal
check_iSensorsPanelSet_obj <- function(iSensorsPanelSet_obj) {
  
  if (!is.character(iSensorsPanelSet_obj$panelSetName)) {
    stop('Name of a panel set must be character')
  }
  if (class(iSensorsPanelSet_obj) != 'iSensorsPanelSet') {
    stop('Object must be "iSensorsPanelSet" class')
  }
  if(!is.list(iSensorsPanelSet_obj$panels)) {
    stop('iSensorsPanelSet_obj$panels should be a list')
  }
  for (panel_ind in seq_along(iSensorsPanelSet_obj$panels)) {
    panel_name <- names(iSensorsPanelSet_obj$panels)[panel_ind]
    panel <- iSensorsPanelSet_obj$panels[[panel_ind]]
    str_print <- paste0('Panel', panel_name, 'must be a vector of gene names')
    check_iSensors_panel(panel, print_massage = str_print)
  }
  return(TRUE)
}

#' @keywords internal
check_random_info <- function(randomInfo) {
  if (!is.list(randomInfo)) {
    stop("`randomInfo` must be a list.")
  }
  
  required_fields <- c("n", "sizes", "majortrend")
  
  missing_fields <- setdiff(required_fields, names(randomInfo))
  if (length(missing_fields) > 0) {
    stop("`randomInfo` is missing required fields: ", paste(missing_fields, collapse = ", "))
  }
  
  if (!is.numeric(randomInfo$n) || length(randomInfo$n) != 1 || randomInfo$n <= 0) {
    stop("`randomInfo$n` must be a single positive number.")
  }
  
  if (!is.numeric(randomInfo$sizes) || length(randomInfo$sizes) != randomInfo$n) {
    stop("`randomInfo$sizes` must be a numeric vector of length equal to `randomInfo$n`.")
  }

  if (any(randomInfo$sizes < min_panel_genes)) {
    stop("`randomInfo$sizes` must be at least ", min_panel_genes, " genes each.")
  }
  
  if (!is.logical(randomInfo$majortrend) || length(randomInfo$majortrend) != 1) {
    stop("`randomInfo$majortrend` must be a single logical value (TRUE or FALSE).")
  }
  
  return(TRUE)
}

#' @keywords internal
check_metaPanels <- function(metaPanels) {
  if (!is.list(metaPanels)) {
    stop("metaPanels must be a list")
  }
  
  for (panel_name in names(metaPanels)) {
    panel <- metaPanels[[panel_name]]
    
    if (!is.list(panel)) {
      stop(paste0("Each meta panel must be a list. Problem with: ", panel_name))
    }
    
    # srcPanels and rule must be present
    if (!("srcPanels" %in% names(panel))) {
      stop(paste0("Missing 'srcPanels' in meta panel: ", panel_name))
    }
    if (!("rule" %in% names(panel))) {
      stop(paste0("Missing 'rule' in meta panel: ", panel_name))
    }
    
    # srcPanels must be a character vector
    if (!is.character(panel$srcPanels)) {
      stop(paste0("'srcPanels' must be a character vector in meta panel: ", panel_name))
    }

    # rule must be a function
    if (!is.function(panel$rule)) {
      stop(paste0("'rule' must be a function in meta panel: ", panel_name))
    }
  }
  
  return(TRUE)
}