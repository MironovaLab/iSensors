#' @title Load Gene Panel Files
#' 
#' @description 
#' Reads gene panel files from a directory and filters them based on biological criteria
#' (species, hormone, type). Supports both RDA and TXT file formats.
#'
#' @param panelsDir Directory containing panel files
#' @param species Optional species filter (e.g., "AT")
#' @param hormone Optional hormone type filter (e.g., "aux", "cyt")
#' @param type Optional type filter (e.g., "cis", "trans")
#' @param panelsFiles Optional vector of specific files to read (NULL reads all matching files)
#' @param format File format ("rda" or "txt", default: "rda")
#'
#' @return A named list of gene panels, where each element is a data frame with a "Genes" column
#'
#' @examples
#' \dontrun{
#' # Read all RDA panels
#' panels <- read_genePanels("path/to/panels")
#'
#' # Read filtered panels
#' panels <- read_genePanels("path/to/panels",
#'                          species = "AT",
#'                          hormone = "aux")
#' }
#'
#' @export
read_genePanels <- function(panelsDir, species = NULL, hormone = NULL, type = NULL, panelsFiles = NULL, format = 'rda') {
  
  if (is.null(panelsFiles)) {
    panelsFiles <- list.files(path = panelsDir, pattern = paste0("*", format), full.names = FALSE)
  }
  panelFilesToGet <- sapply(strsplit(panelsFiles, '_'), function(x) {
    sceckSpecies <- if(is.null(species)) TRUE else any(species %in% x)
    checkHormone <- if(is.null(hormone)) TRUE else any(hormone %in% x)
    checkType <- if(is.null(type)) TRUE else any(type %in% x)
    
    return (all(sceckSpecies, checkHormone, checkType))
  })
  panelsFiles <- panelsFiles[panelFilesToGet]
  
  if (format == 'rda') {
    # cat('Format: ', format, '\n')
    panelEnv <- new.env()
    for (panelFile in panelsFiles) {
      tmpDir <- file.path(panelsDir, panelFile)
      load(tmpDir, envir = panelEnv)
    }
    panelsNames <- ls(envir = panelEnv)
    # cat('Panels: ', panelsNames, '\n')
    panelsList <- lapply(panelsNames, function(x) {
      tmpPanel <- get(x, envir = panelEnv)[['genes']]
      # cat('Panel current:', x, '\n')
      # cat('Panel current length: ', length(tmpPanel),'\n')
      tmpPanel <- data.frame(tmpPanel)
      colnames(tmpPanel) <- c('Genes')
      # cat('Panel info: ', dim(tmpPanel), '\n')
      return (tmpPanel)
    })
    names(panelsList) <- panelsNames
    # cat('Panels list info:', class(panelsList), length(panelsList), '\n')
    rm(panelEnv)
    return(panelsList)
  }
  
  if (format == 'txt') {
    panelsList <- lapply(panelsFiles, function(x) {
      readDir <- file.path(panelsDir, x)
      readData <- data.frame(read.delim(readDir, header = FALSE)[,1])
      colnames(readData) <- c('Genes')
      # cat(x, 'was read', '\n')
      return(readData)
    })
    names(panelsList) <- sub("\\.txt$", "", panelsFiles)
    
    return(panelsList)
  }
}