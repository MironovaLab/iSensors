#' @keywords internal
create_iSensors <- function(data = NULL, panelSet = NULL, labels = NULL) {
  iSensor_obj <- list('exprData' = data,
                      'panelSet' = panelSet,
                      'signals' = NULL,
                      'sampleLabels' = labels,
                      'importance' = NULL,
                      'metaData' = NULL)
  class(iSensor_obj) <- 'iSensors'
  
  return(iSensor_obj)
}

#' @keywords internal
create_iSensorsPanelSet <- function(name, panelsList, additional = NULL) {
  iSensorsPanelSet_obj <- list('panelSetName' = name,
                               'panels' = panelsList,
                               'additional' = additional)
  class(iSensorsPanelSet_obj) <- 'iSensorsPanelSet'
  
  check_iSensorsPanelSet_obj(iSensorsPanelSet_obj)
  return(iSensorsPanelSet_obj)
}
