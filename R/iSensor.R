#' Execute iSensor Analysis Pipeline on Seurat Object
#' 
#' Performs comprehensive gene panel analysis of single-cell RNA-seq data stored in a Seurat object,
#' calculating signaling activity scores for specified gene panels and storing results as new assays.
#'
#' @import Seurat
#'
#' @param seuratObject A Seurat object containing single-cell expression data
#' @param presetPanels Character vector specifying which predefined gene panels to use.
#'        Default is 'Auxin'. Multiple panels can be specified (e.g., c('Auxin', 'Cytokinin'))
#' @param signals Character vector of signal calculation methods to apply. Options include:
#'        \itemize{
#'          \item 'mean' - Simple mean expression
#'          \item 'mean_normed' - Mean expression with normalization (default)
#'          \item 'median' - Median expression
#'          \item 'median_normed' - Median expression with normalization
#'        }
#' @param seurLayer Character specifying which Seurat assay layer to use ('data', 'counts', or 'scale.data').
#'        Default is "data"
#' @param randPanels Integer specifying number of random gene panels to generate for background comparison.
#'        Default is 2
#' @param randSize Integer vector of length 2 specifying size range (min,max) for random panels.
#'        Default is c(200, 500)
#' @param majortrend Logical indicating whether to maintain major expression trends in random panels.
#'        Default is TRUE
#' @param usePanels Optional character vector of specific panel names to analyze. If NULL, uses all available panels.
#'        Default is NULL
#' @param defaultAssay Character specifying which iSensor assay to set as default in output Seurat object.
#'        If NULL, keeps current default assay. Default is NULL
#' @param useParallel Logical indicating whether to use parallel processing. Default is FALSE
#' @param nCores Integer specifying number of cores for parallel processing if useParallel=TRUE.
#'        If NULL, uses available cores - 1. Default is NULL
#' @param metaPanels List defining composite meta-panels and their calculation rules. Each element should be:
#'        \itemize{
#'          \item A list with 'srcPanels' (character vector of source panels)
#'          \item 'rule' (function to combine panels, e.g., prod for product)
#'        }
#'        Default example combines two auxin panels using product
#'
#' @return Modified Seurat object with new assays containing iSensor signaling activity scores.
#'         Each assay is named as 'iSensor_[signal_type]' (e.g., 'iSensor_mean_normed')
#'
#' @examples
#' \dontrun{
#' # Basic analysis with defaults
#' seurat_obj <- iSensor_pipeline(seurat_obj)
#'
#' # Advanced analysis with multiple panels and signals
#' seurat_obj <- iSensor_pipeline(
#'     seurat_obj,
#'     presetPanels = c('Auxin', 'Cytokinin'),
#'     signals = c('mean_normed', 'median')
#'     metaPanels = list(
#'         'MyMeta' = list(
#'             srcPanels = c('Panel1', 'Panel2'),
#'             rule = mean
#'         )
#'     )
#' )
#' }
#'
#' @section Analysis Workflow:
#' The pipeline performs these key steps:
#' \enumerate{
#'   \item Extracts expression matrix from Seurat object
#'   \item Initializes iSensor analysis object
#'   \item Loads specified gene panels
#'   \item Generates random background panels (if enabled)
#'   \item Calculates all requested signal types
#'   \item Stores results as new assays in Seurat object
#' }
#'
#' @section Technical Details:
#' Signal calculation methods:
#' \describe{
#'   \item{mean/median}{Simple average or median expression of panel genes}
#'   \item{*_normed}{Normalized by background (random panels or specified method)}
#' }
#'
#' @seealso
#' \code{\link{create_iSensor}} for object initialization,
#' \code{\link{iSensor_signal}} for signal calculation details,
#' \code{\link{set_iSensor_workPanel}} for panel management
#'
#' @export
iSensor_pipeline <- function(seuratObject, presetPanels = 'Auxin', signals = c('mean_normed'), seurLayer = "data",
                             randPanels = 2, randSize = c(200, 500), majortrend = TRUE,
                             usePanels = NULL, defaultAssay = NULL, useParallel = FALSE, nCores = NULL,
                             metaPanels = list('R2D2' = list('srcPanels' = c('AT_aux_trans_A-ARF', 'AT_aux_trans_IAA'),
                                                             'rule' = prod))) {
    if (class(seuratObject) != 'Seurat') {
        stop("Error: function must be applied to Seurat object")
    }

    expression_data <- as.matrix(GetAssayData(seuratObject, layer = seurLayer))
    expression_data[is.nan(expression_data)] <- 0
    iSensor_test <- create_iSensor(data = expression_data)
    iSensor_test <- set_iSensor_workPanel(iSensor_test, presetPanels)
    iSensor_test <- add_iSensor_random_panel(iSensor_test, randNum=randPanels, randSize=randSize,
                                         majortrend = majortrend)

    for (signal in signals) {
        if (signal == 'mean') {
            iSensor_test <- iSensor_signal(iSensor_test, transform = 'mean', normed = FALSE,
                                           panels = usePanels, doParallel = useParallel, nCores = nCores,
                                           metaPanels = metaPanels)
        }
        if (signal == 'mean_normed') {
            iSensor_test <- iSensor_signal(iSensor_test, transform = 'mean', normed = TRUE,
                                           panels = usePanels, doParallel = useParallel, nCores = nCores,
                                           metaPanels = metaPanels)
        }
        if (signal == 'median') {
            iSensor_test <- iSensor_signal(iSensor_test, transform = 'median', normed = FALSE,
                                           panels = usePanels, doParallel = useParallel, nCores = nCores,
                                           metaPanels = metaPanels)
        }
        if (signal == 'median_normed') {
            iSensor_test <- iSensor_signal(iSensor_test, transform = 'median', normed = TRUE,
                                           panels = usePanels, doParallel = useParallel, nCores = nCores,
                                           metaPanels = metaPanels)
        }
    }

    for (assay in names(iSensor_test$Signals)) {
        seuratObject[[paste0('iSensor_', assay)]] <- CreateAssayObject(counts = as(t(iSensor_test$Signals[[assay]]), "dgCMatrix"))
    }

    if (!is.null(defaultAssay)) {
        DefaultAssay(seuratObject) <- paste0('iSensor_', defaultAssay)
        cat('Set DefaultAssay to', paste0('iSensor_', defaultAssay, '\n'))
    }

    return(seuratObject)
}

#' Initialize an iSensor Analysis Object
#'
#' Creates a new iSensor object for gene panel analysis, serving as the foundation
#' for all subsequent iSensor pipeline operations. The object stores expression data,
#' gene panels, and analysis results in a structured format.
#'
#' @param data A numeric matrix of gene expression data where:
#'        \itemize{
#'          \item Rows represent genes (features)
#'          \item Columns represent samples/cells
#'          \item Row names should be gene identifiers
#'          \item Column names should be sample/cell identifiers
#'        }
#'        If NULL, creates an empty object (default: NULL)
#' @param panels Optional named list of gene panels where each element is:
#'        \itemize{
#'          \item A character vector of gene names
#'          \item Named with the panel identifier
#'        }
#'        (default: NULL)
#' @param labels Optional vector of sample/cell labels corresponding to columns
#'        in the expression matrix (default: NULL)
#'
#' @return An object of class `iSensor` containing:
#' \describe{
#'   \item{exprData}{Input expression matrix}
#'   \item{genePanels}{List of gene panels}
#'   \item{Signals}{Placeholder for calculated signaling activity (initialized NULL)}
#'   \item{sampleLabels}{Sample/cell labels}
#'   \item{importance}{Placeholder for feature importance scores (initialized NULL)}
#'   \item{genePanelSets}{Placeholder for panel sets (initialized NULL)}
#' }
#'
#' @examples
#' # Create empty iSensor object
#' iSensor_obj <- create_iSensor()
#'
#' # Initialize with expression data
#' expr_data <- matrix(rnorm(1000), nrow=100, 
#'                     dimnames=list(paste0("gene", 1:100), paste0("cell", 1:10)))
#' iSensor_obj <- create_iSensor(data = expr_data)
#'
#' # With custom panels
#' custom_panels <- list(
#'   "Pathway_A" = c("gene1", "gene5", "gene20"),
#'   "Pathway_B" = c("gene3", "gene8", "gene15")
#' )
#' iSensor_obj <- create_iSensor(data = expr_data, panels = custom_panels)
#'
#' @section Object Structure:
#' The iSensor object uses S3 class system and automatically loads default panels
#' via \code{\link{add_defaultPanels}}. Key features:
#' \itemize{
#'   \item Maintains compatibility with Seurat objects via \code{\link{iSensor_pipeline}}
#'   \item Designed for both single-cell and bulk RNA-seq analysis
#'   \item Extensible structure for additional analysis components
#' }
#'
#' @seealso
#' \code{\link{add_defaultPanels}} for loading default gene panels,
#' \code{\link{iSensor_pipeline}} for Seurat integration,
#' \code{\link{set_iSensor_workPanel}} for panel management
#'
#' @export
create_iSensor <- function(data = NULL, panels = NULL, labels = NULL) {
    iSensor_obj <- list('exprData' = data,
                        'genePanels' = panels,
                        'Signals' = NULL,
                        'sampleLabels' = labels,
                        'importance' = NULL,
                        'genePanelSets' = NULL)
    class(iSensor_obj) <- "iSensor"
    
    iSensor_obj <- add_defaultPanels(iSensor_obj)
    
    return(iSensor_obj)
}

#' Set Working iSensor Gene Panel
#'
#' Configures the iSensor object to use a specified gene panel set for analysis.
#'
#' @param iSensor_obj An iSensor object containing gene panel sets.
#' @param panelSet Character string specifying the name of the gene panel set to use. 
#'        Default is 'Auxin'.
#'
#' @return Updated iSensor object with the selected gene panel set assigned.
#'
#' @examples
#' \dontrun{
#' # Example of setting a working panel
#' iSensor_obj <- set_iSensor_workPanel(iSensor_obj, panelSet = 'Cytokinin')
#' }
#' 
#' @section Details:
#' This function checks if the specified panel set exists in the iSensor object. 
#' If it does, the function assigns the corresponding gene panels to the iSensor object 
#' for subsequent analysis.
#'
#' @seealso
#' \code{\link{iSensor_pipeline}}, \code{\link{create_iSensor}}, 
#' \code{\link{add_iSensor_random_panel}}, \code{\link{iSensor_signal}}
#'
#' @export
set_iSensor_workPanel <- function(iSensor_obj, panelSet = 'Auxin') {
    if (!(panelSet %in% names(iSensor_obj$genePanelSets))) {
        stop("Error: no such panel in iSensor_obj$genePanelSets")
    }

    getPanel <- iSensor_obj$genePanelSets[[panelSet]]
    iSensor_obj$genePanels <- getPanel

    return(iSensor_obj)
}

#' Add Default Gene Panels to iSensor Object
#'
#' Loads default gene panels from the package's extdata directory and adds them 
#' to the specified iSensor object.
#'
#' @param iSensor_obj An iSensor object to which the default gene panels will be added.
#'
#' @return Updated iSensor object with default gene panels included.
#'
#' @examples
#' \dontrun{
#' # Example of adding default panels to an iSensor object
#' iSensor_obj <- add_defaultPanels(iSensor_obj)
#' }
#' 
#' @section Details:
#' This function retrieves gene panel directories from the package's extdata folder, 
#' reads the panel files, and adds them to the iSensor object as new panel sets. 
#' Each folder in the directory corresponds to a different hormone panel.
#'
#' @seealso
#' \code{\link{iSensor_pipeline}}, \code{\link{set_iSensor_workPanel}}, 
#' \code{\link{add_iSensor_panelSet}}, \code{\link{read_genePanels}}
#'
#' @export
add_defaultPanels <- function(iSensor_obj) {
    panelsDir <- system.file("extdata/geneSensors", package = "iSensor")
    # panelsDir <- 'geneSensors/'
    folders <- list.dirs(path = panelsDir, full.names = FALSE, recursive = FALSE)
    
    for (horm in folders) {
        tmpDir <- file.path(panelsDir, horm)
        panel <- read_genePanels(panelsDir = tmpDir,
                                 panelsFiles = list.files(path = tmpDir, full.names = FALSE),
                                )
        iSensor_obj <- add_iSensor_panelSet(iSensor_obj, panelsSet=panel, setName = horm)
    }

    return(iSensor_obj)
}

#' Add Expression Data to iSensor Object
#'
#' Assigns expression data to the specified iSensor object for further analysis.
#'
#' @param iSensor_obj An iSensor object to which the expression data will be added.
#' @param data A matrix or data frame containing the expression data to be assigned.
#'
#' @return Updated iSensor object with the new expression data included.
#'
#' @examples
#' \dontrun{
#' # Example of adding expression data to an iSensor object
#' iSensor_obj <- add_iSensor_data(iSensor_obj, expression_data_matrix)
#' }
#' 
#' @section Details:
#' This function updates the expression data slot of the iSensor object with the 
#' provided data, allowing for subsequent analysis using the iSensor framework.
#'
#' @seealso
#' \code{\link{iSensor_pipeline}}, \code{\link{set_iSensor_workPanel}}, 
#' \code{\link{add_defaultPanels}}
#'
#' @export
add_iSensor_data <- function(iSensor_obj, data) {
    iSensor_obj$exprData <- data
    return(iSensor_obj)
}

#' Add Sample Labels to iSensor Object
#'
#' Assigns sample labels to the specified iSensor object, linking them to the expression data.
#'
#' @param iSensor_obj An iSensor object to which the sample labels will be added.
#' @param labels A character vector containing labels for each sample in the expression data.
#'
#' @return Updated iSensor object with the new sample labels included.
#'
#' @examples
#' \dontrun{
#' # Example of adding sample labels to an iSensor object
#' iSensor_obj <- add_iSensor_labels(iSensor_obj, c("Sample1", "Sample2", "Sample3"))
#' }
#' 
#' @section Details:
#' This function checks if the length of the provided labels matches the number of samples 
#' in the expression data. If the lengths match, the labels are assigned to the iSensor object, 
#' allowing for easier identification of samples during analysis.
#'
#' @seealso
#' \code{\link{iSensor_pipeline}}, \code{\link{add_iSensor_data}}, 
#' \code{\link{set_iSensor_workPanel}}
#'
#' @export
add_iSensor_labels <- function(iSensor_obj, labels) {
    if (length(colnames(iSensor_obj$exprData)) != length(labels)) {
        stop("Error: length of 'labels' does not match number of samples in iSensor expr data")
    }
    iSensor_obj$sampleLabels <- labels

    return(iSensor_obj)
}

#' Add Gene Panels to iSensor Object
#'
#' Incorporates additional gene panels into the specified iSensor object, updating existing panels 
#' and adding new ones as necessary.
#'
#' @param iSensor_obj An iSensor object to which the gene panels will be added.
#' @param panels A list of gene panels to be added or updated in the iSensor object.
#' @param drop Logical indicating whether to drop existing panels that are not included in the new panels. Default is FALSE.
#'
#' @return Updated iSensor object with the new and updated gene panels included.
#'
#' @examples
#' \dontrun{
#' # Example of adding gene panels to an iSensor object
#' new_panels <- list(panel1 = gene_data1, panel2 = gene_data2)
#' iSensor_obj <- add_iSensor_workPanel(iSensor_obj, new_panels)
#' }
#' 
#' @section Details:
#' This function checks for existing gene panels in the iSensor object and updates them 
#' with the new panels provided. If a panel already exists, it will be replaced; otherwise, 
#' it will be added to the iSensor object. The `drop` parameter controls whether panels 
#' not included in the new list should be removed.
#'
#' @seealso
#' \code{\link{iSensor_pipeline}}, \code{\link{set_iSensor_workPanel}}, 
#' \code{\link{add_defaultPanels}}
#'
#' @export
add_iSensor_workPanel <- function(iSensor_obj, panels, drop = FALSE) {
    havePanel <- names(panels) %in% names(iSensor_obj$genePanels)
    oldPanel <- panels[havePanel]
    newPanel <- panels[!havePanel]
    iSensor_obj$genePanels[names(oldPanel)] <- oldPanel
    iSensor_obj$genePanels <- c(iSensor_obj$genePanels, newPanel)
    return(iSensor_obj)
}

#' Add Gene Panel Set to iSensor Object
#'
#' Adds a new set of gene panels to the specified iSensor object, allowing for organized 
#' management of multiple panel sets.
#'
#' @param iSensor_obj An iSensor object to which the gene panel set will be added.
#' @param panelsSet A list of gene panels to be included in the new panel set.
#' @param setName Optional character string specifying the name of the new panel set. 
#'        If NULL, a default name will be generated. Default is NULL.
#'
#' @return Updated iSensor object with the new gene panel set included.
#'
#' @examples
#' \dontrun{
#' # Example of adding a gene panel set to an iSensor object
#' new_panel_set <- list(panel1 = gene_data1, panel2 = gene_data2)
#' iSensor_obj <- add_iSensor_panelSet(iSensor_obj, new_panel_set, setName = "CustomPanelSet")
#' }
#' 
#' @section Details:
#' This function allows users to add a new set of gene panels to the iSensor object. 
#' If a set name is not provided, a default name is generated based on the current number 
#' of panel sets. The function ensures that the new panel set is properly integrated into 
#' the existing structure of the iSensor object.
#'
#' @seealso
#' \code{\link{iSensor_pipeline}}, \code{\link{add_iSensor_workPanel}}, 
#' \code{\link{add_defaultPanels}}
#'
#' @export
add_iSensor_panelSet <- function(iSensor_obj, panelsSet, setName = NULL) {
    if (is.null(setName)) {
        setName <- paste0('panelSet_', length(iSensor_obj$genePanelSets) + 1)
    }
    if(length(panelsSet) > 1) {
        # panelsSet <- setNames(panelsSet, setName)
        panelsSet_tmp <- panelsSet
        panelsSet <- list()
        panelsSet[[setName]] <- panelsSet_tmp
    }
    iSensor_obj$genePanelSets <- c(iSensor_obj$genePanelSets, panelsSet)

    return(iSensor_obj)
}

#' Read Gene Panels from Files
#'
#' Reads gene panel data from specified files in a directory and returns a list of data frames 
#' containing gene information.
#'
#' @param panelsDir Character string specifying the directory containing the gene panel files.
#' @param panelsFiles Character vector of file names to be read from the specified directory.
#' @param panelsNames Optional character vector of names to assign to the gene panels. 
#'        If NULL, names will be derived from the file names. Default is NULL.
#'
#' @return A list of data frames, each containing a column of gene names corresponding to 
#'         the specified panel files.
#'
#' @examples
#' \dontrun{
#' # Example of reading gene panels from a directory
#' gene_panels <- read_genePanels("path/to/panels", c("panel1.txt", "panel2.txt"))
#' }
#' 
#' @section Details:
#' This function reads gene panel files from the specified directory, creating a data frame 
#' for each file that contains the gene names. If panel names are provided, they will be used 
#' as the names of the list elements; otherwise, names will be generated from the file names 
#' by removing the ".txt" extension.
#'
#' @seealso
#' \code{\link{add_defaultPanels}}, \code{\link{add_iSensor_panelSet}}, 
#' \code{\link{set_iSensor_workPanel}}
#'
#' @export
read_genePanels <- function(panelsDir, panelsFiles, panelsNames = NULL) {
    iSensor_panel <- lapply(panelsFiles, function(x) {
        readDir <- file.path(panelsDir, x)
        readData <- data.frame(read.delim(readDir, header = FALSE)[,1])
        colnames(readData) <- c('Genes')
        # cat(x, 'was read', '\n')
        return(readData)
    })
    if (!is.null(panelsNames)) {
        names(iSensor_panel) <- panelsNames
    } else {
        names(iSensor_panel) <- sub("\\.txt$", "", panelsFiles)
    }
    return(iSensor_panel)
}

#' Add Random Gene Panels to iSensor Object
#'
#' Generates random gene panels based on the expression data in the iSensor object and 
#' adds them to the specified iSensor object.
#'
#' @param iSensor_obj An iSensor object to which the random gene panels will be added.
#' @param randNum Integer specifying the number of random panels to generate.
#' @param randSize Numeric vector specifying the sizes of the random panels. The length 
#'        of this vector must match `randNum`.
#' @param majortrend Logical indicating whether to include a major trend panel based on 
#'        all features. Default is TRUE.
#'
#' @return Updated iSensor object with new random gene panels included.
#'
#' @examples
#' \dontrun{
#' # Example of adding random gene panels to an iSensor object
#' iSensor_obj <- add_iSensor_random_panel(iSensor_obj, randNum = 3, randSize = c(100, 200, 300))
#' }
#' 
#' @section Details:
#' This function checks if the expression data in the iSensor object is available. It then 
#' generates random gene panels of specified sizes and adds them to the iSensor object. 
#' A major trend panel, containing all features, is also included by default.
#'
#' @seealso
#' \code{\link{iSensor_pipeline}}, \code{\link{add_iSensor_workPanel}}, 
#' \code{\link{add_defaultPanels}}
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

#' Calculate iSensor Signals for Gene Panels
#'
#' Computes signaling activity scores for specified gene panels in the iSensor object, 
#' using either mean or median transformations, with optional normalization.
#'
#' @param iSensor_obj An iSensor object containing expression data and gene panels.
#' @param panels Optional character vector of specific gene panel names to use. 
#'        If NULL, all available panels will be used.
#' @param transform Character string specifying the transformation method to use. 
#'        Options are 'mean' or 'median'.
#' @param normed Logical indicating whether to normalize the signals. Default is FALSE.
#' @param normBy Character string specifying the normalization method. Options are 'cols' or 'rows'. 
#'        Default is 'cols'.
#' @param metaPanels Optional list defining meta-panels and their calculation rules.
#' @param doParallel Logical indicating whether to use parallel processing. Default is FALSE.
#' @param nCores Number of cores for parallel processing if doParallel=TRUE. Default is NULL 
#'        (autodetect).
#'
#' @return Updated iSensor object with calculated signals added to the Signals slot.
#'
#' @examples
#' \dontrun{
#' # Example of calculating signals for an iSensor object
#' iSensor_obj <- iSensor_signal(iSensor_obj, panels = c("panel1", "panel2"), 
#'                                transform = "mean", normed = TRUE)
#' }
#' 
#' @section Details:
#' This function checks if the expression data in the iSensor object is available. It then 
#' computes the specified signal (mean or median) for each gene panel, with optional normalization 
#' based on the specified method. The results are stored in the Signals slot of the iSensor object.
#' If meta-panels are provided, they will be calculated based on the specified rules.
#'
#' @seealso
#' \code{\link{iSensor_pipeline}}, \code{\link{add_iSensor_workPanel}}, 
#' \code{\link{add_defaultPanels}}, \code{\link{add_iSensor_random_panel}}
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
        non_zero_variance_rows <- apply(expression_data, 1, var) > 0
        expression_data <- expression_data[non_zero_variance_rows, , drop = FALSE] 
        
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