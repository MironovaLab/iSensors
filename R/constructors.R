#' @title Create iSensor Analysis Object
#' 
#' @description 
#' Initializes an iSensor analysis object that stores expression data and gene panels 
#' for subsequent processing in the iSensor pipeline.
#'
#' @param data Numeric matrix of expression data (genes in rows, samples in columns)
#' @param panels Optional list of custom gene panels (named character vectors)
#' @param labels Optional vector of sample labels or conditions
#' @param species Species specification for default panels ("AT")
#' @param hormone Optional hormone type specification for panel filtering
#' @param type Optional tissue/cell type specification for panel filtering 
#' @param format Format for default panels ("rda" or "list", default: "rda")
#'
#' @return An object of class "iSensor" containing:
#' \itemize{
#'   \item Input expression data
#'   \item Gene panels (custom and/or default)
#'   \item Placeholders for computed signals and importance scores
#' }
#'
#' @examples
#' # Basic usage with expression data only
#' expr_data <- matrix(rnorm(1000), nrow = 100)
#' sensor_obj <- create_iSensor(data = expr_data)
#'
#' @seealso \code{\link{iSensor_pipeline}} for the complete analysis workflow
#' @export
create_iSensor <- function(data = NULL, panels = NULL, labels = NULL,
                           species = NULL, hormone = NULL, type = NULL, format = 'rda') {
  iSensor_obj <- list('exprData' = data,
                      'genePanels' = panels,
                      'Signals' = NULL,
                      'sampleLabels' = labels,
                      'importance' = NULL,
                      'genePanelSets' = NULL,
                      'metaData' = NULL)
  class(iSensor_obj) <- "iSensor"
  
  iSensor_obj <- add_defaultPanels(iSensor_obj, species = species, hormone = hormone, type = type, format = 'rda')
  
  return(iSensor_obj)
}

#' @title Add Default Gene Panels to iSensor Object
#'
#' @description 
#' Loads and adds pre-defined gene panels to an iSensor object based on specified
#' biological filters (species, hormone, tissue type). Panels are stored in RDA format.
#'
#' @param iSensor_obj An iSensor object created by `create_iSensor`
#' @param species Optional species filter (e.g., "AT")
#' @param hormone Optional hormone type filter (e.g., "aux", "cyt")
#' @param type Optional type filter (e.g., "cis", "trans")
#' @param format File format for panels ("rda" or "txt", default: "rda")
#'
#' @return The modified iSensor object with default panels added to `genePanelSets`
#'
#' @examples
#' \dontrun{
#' # Add all default panels
#' iSensor_obj <- add_defaultPanels(iSensor_obj)
#'
#' # Add species-specific panels
#' iSensor_obj <- add_defaultPanels(iSensor_obj, species = "AT")
#'
#' # Add filtered panels
#' iSensor_obj <- add_defaultPanels(iSensor_obj, 
#'                                species = "AT",
#'                                hormone = "aux",
#'                                type = "trans")
#' }
#'
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
  
  panelsDir <- system.file("extdata/geneSensors", package = "iSensor")
  # panelsDir <- 'geneSensors/panels_rda/'
  
  panelsList <- read_genePanels(panelsDir = panelsDir,
                                species = species, hormone = hormone, type = type,
                                format = format)
  iSensor_obj <- add_iSensor_panelSet(iSensor_obj, panelsSet = panelsList, setName = panelSet_name_gen())
  
  return(iSensor_obj)
}