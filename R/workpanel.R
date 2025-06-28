#' @title Add Random Gene Panels to iSensor Object
#'
#' @description 
#' Generates and adds random gene panels to an iSensor object. Can optionally include
#' a "majortrend" panel containing all genes. Each random panel samples genes without replacement.
#'
#' @param iSensor_obj An iSensor object created by `create_iSensor`
#' @param randNum Number of random panels to generate (must match length of `randSize`)
#' @param randSize Vector specifying sizes for each random panel
#' @param majortrend Logical indicating whether to include a panel with all genes (default: TRUE)
#'
#' @return The modified iSensor object with added random panels
#'
#' @examples
#' # Create iSensor object with 1000 genes
#' iSensor_obj <- create_iSensor(data = matrix(rnorm(10000), nrow = 1000))
#' 
#' # Add 3 random panels of different sizes
#' iSensor_obj <- add_iSensor_random_panel(iSensor_obj, 
#'                                       randNum = 3,
#'                                       randSize = c(50, 100, 200))
#'
#' @export
add_iSensor_random_panel <- function(iSensor_obj, randNum, randSize, majortrend = TRUE) { 
  if (is.null(iSensor_obj$exprData)) {
    stop("Error: iSensor expr data is empty")
  }
  
  if (length(randSize) != randNum) {
    stop("Error: number of random panels 'randNum' does not match length of 'randSize'")
  }
  
  featuresAll <- rownames(iSensor_obj$exprData)
  majorPanel <- as.data.frame(featuresAll)
  colnames(majorPanel) <- c('Genes')
  randomPanels <- list('majortrend' = majorPanel)
  randomLst <- lapply(randSize, function(x) {
    randomPanel <- sample(featuresAll, size = x, replace = FALSE)
    randomPanel <- as.data.frame(randomPanel)
    colnames(randomPanel) <- c('Genes')
    return(randomPanel)
  })
  names(randomLst) <- paste0("random", 1:randNum, sep = "")
  randomPanels <- c(randomPanels, randomLst)
  iSensor_obj <- add_iSensor_workPanel(iSensor_obj, randomPanels)
  
  return (iSensor_obj)
}