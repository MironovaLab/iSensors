#' @title Create iSensor Analysis Object
#' @seealso \code{\link{iSensor_pipeline}} for the complete analysis workflow
#' @export
create_iSensors <- function(data = NULL, panelSet = NULL, labels = NULL) {
  iSensor_obj <- list('exprData' = data,
                      'panelSet' = panelSet,
                      'signals' = NULL,
                      'sampleLabels' = labels,
                      'importance' = NULL,
                      'metaData' = NULL)
  class(iSensor_obj) <- 'iSensors'
  
  # iSensor_obj <- add_defaultPanels(iSensor_obj, species = species, hormone = hormone, type = type, format = 'rda')
  
  return(iSensor_obj)
}

create_iSensorsPanelSet <- function(name, panelsList) {
  iSensorsPanelSet_obj <- list('panelSetName' = name,
                               'panels' = panelsList)
  class(iSensorsPanelSet_obj) <- 'iSensorsPanelSet'
  
  check_iSensorsPanelSet_obj(iSensorsPanelSet_obj)
  return(iSensorsPanelSet_obj)
}

#' @title Add Default Gene Panels to iSensor Object
#' @seealso \code{\link{read_genePanels}} for panel loading implementation
#' @export
add_defaultPanels <- function(iSensor_obj, species = NULL, hormone = NULL, type = NULL, format = 'rda') {
  panelSet_name_gen <- function() {
    speciesName <- if(is.null(species)) 'allSpec' else species
    hormoneName <- if(is.null(hormone)) 'allHorm' else hormone
    typeName <- if(is.null(type)) 'allType' else type
    name <- paste0(c(speciesName, hormoneName, typeName), collapse='_')
    return(name)
  }
  
  panelsDir <- system.file("extdata/geneSensors", package = "iSensors")
  # panelsDir <- 'geneSensors/panels_rda/'
  
  panelsList <- read_genePanels(panelsDir = panelsDir,
                                species = species, hormone = hormone, type = type,
                                format = format)
  iSensor_obj <- add_iSensor_panelSet(iSensor_obj, panelsSet = panelsList, setName = panelSet_name_gen())
  
  return(iSensor_obj)
}