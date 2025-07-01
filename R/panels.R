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

load_panel_from_rda <- function(file_path) {
  panelEnv <- new.env()
  load(file_path, envir = panelEnv)
  panelName <- ls(envir = panelEnv)
  if (length(panelName) != 1) {
    warning("Expected one object per panel file, got: ", paste(panelName, collapse = ", "))
  }
  panel <- get(panelName[[1]], envir = panelEnv)
  panel$name <- panelName
  check_iSensorsPanel_obj(panel)
  return(panel)
}

load_panel_files <- function(panelsDir, files) {
  full_paths <- file.path(panelsDir, files)
  panel_objs <- lapply(full_paths, load_panel_from_rda)
  return(panel_objs)
}

filter_and_load_panels <- function(panelsDir, species = NULL, hormone = NULL, type = NULL, filter = TRUE) {
  
  panelsFiles <- list.files(path = panelsDir, pattern = '\\.rda$', full.names = FALSE)
  
  if (filter) {
    name_delim <- '_'
    firstThreeParts <- lapply(strsplit(panelsFiles, name_delim), function(x) head(x, 3))
    all_species  <- unique(sapply(firstThreeParts, function(x) x[1]))
    all_hormones <- unique(sapply(firstThreeParts, function(x) x[2]))
    all_types    <- unique(sapply(firstThreeParts, function(x) x[3]))
    
    if (!is.null(species) && !any(species %in% all_species)) {
      stop('No such "species" found in files')
    }
    if (!is.null(hormone) && !any(hormone %in% all_hormones)) {
      stop('No such "hormone" found in files')
    }
    if (!is.null(type) && !any(type %in% all_types)) {
      stop('No such "type" found in files')
    }
    
    panelFilesToGet <- sapply(strsplit(panelsFiles, name_delim), function(x) {
      checkSpecies <- is.null(species) || any(species %in% x)
      checkHormone <- is.null(hormone) || any(hormone %in% x)
      checkType <- is.null(type) || any(type %in% x)
      all(c(checkSpecies, checkHormone, checkType))
    })
    panelsFiles <- panelsFiles[panelFilesToGet]
  }
  
  panel_objs <- load_panel_files(panelsDir = panelsDir, files = panelsFiles)
  panel_names <- sapply(panel_objs, function(p) p$name)
  panelsList <- setNames(lapply(panel_objs, function(p) p$genes), panel_names)
  return(panelsList)
}

load_sensors <- function(setName, species = NULL, hormone = NULL, type = NULL,
                         defaultPanels = TRUE, customPanels = FALSE,
                         random = TRUE, randomInfo = NULL, metaPanels = NULL) {
  resultList <- list()
  if (defaultPanels) {
    defaultDir <- system.file("extdata/geneSensors", package = "iSensors")
    defaultPanelsList <- filter_and_load_panels(defaultDir, species, hormone, type)
    resultList <- c(resultList, defaultPanelsList)
  }
  if (customPanels) {
    customDir <- 'iSensors/'
    customDir_path <- file.path(getwd(), customDir)
    if (dir.exists(customDir_path)) {
      customPanelsList <- load_panels_from_dir(customDir_path, filter = FALSE)
      resultList <- c(resultList, customPanelsList)
    } else {
      message('No "iSensors" folder is found in the working directory, no custom panels will be used')
    }
  }
  
  # Создание additional
  additional <- list()
  
  if (isTRUE(random)) {
    # Пользователь запросил рандомы по умолчанию
    additional$random <- list(n = 2, sizes = c(200, 500), majortrend = TRUE)
  }
  
  if (!is.null(randomInfo)) {
    # Пользователь указал собственные настройки
    check_random_info(randomInfo)
    additional$random <- randomInfo
  }
  
  if (!is.null(metaPanels)) {
    # check_metaPanels(metaPanels)
    additional$meta <- metaPanels
  }
  
  panelSet <- create_iSensorsPanelSet(name = setName, panelsList = resultList, additional = additional)
  return(panelSet)
}
