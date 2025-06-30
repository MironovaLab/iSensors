check_iSensorsPanel_obj <- function(panel_obj) {
  
  if (!('genes' %in% names(panel_obj))) {
    stop('Panel object must have "genes" slot')
  }
  
  if (!('name' %in% names(panel_obj))) {
    stop('Panel object must have "name" slot')
  }
  
  check_iSensors_panel(panel_obj$genes)
  return(TRUE)
}

check_iSensors_obj <- function(iSensor_obj, check_exprData = TRUE, check_panelSet = TRUE,
                               check_signals = TRUE, check_sampleLabels = TRUE) {
  
  return(TRUE)
}

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
