#' @title iSensor Analysis Pipeline for Seurat Objects
#'
#' @description 
#' Main pipeline for analyzing single-cell RNA-seq data using iSensor framework.
#' Computes gene panel signatures and stores results as new assays in Seurat object.
#' Supports parallel processing, custom panels, and meta-panel creation.
#'
#' @param seuratObject A Seurat object (v3+) containing single-cell data
#' @param species Species for default panels ("AT", etc.)
#' @param hormone Hormone type for panel filtering (optional)
#' @param type Panel type for panel filtering (optional)
#' @param signals Vector of signals to compute: 'mean', 'mean_normed',
#'        'median', 'median_normed' (default: 'mean_normed')
#' @param seurLayer Seurat assay layer to use ("data", "counts", etc.)
#' @param usePanelPreset Custom panel set to add (overrides species/hormone/type)
#' @param presetPanelName Name for custom panel set (required if usePanelPreset specified)
#' @param randPanels Number of random panels to generate (default: 2)
#' @param randSize Size range for random panels (default: c(200, 500))
#' @param majortrend Whether to include all genes as "majortrend" panel (default: TRUE)
#' @param usePanels Specific panels to analyze (NULL uses all available)
#' @param defaultAssay Name of iSensor assay to set as default (optional)
#' @param useParallel Enable parallel processing (default: FALSE)
#' @param nCores Number of cores for parallel processing (NULL uses available-1)
#' @param metaPanels List of meta-panels to create (see Details)
#'
#' @details
#' Meta-panels should be specified as named lists with:
#' \itemize{
#'   \item \code{srcPanels}: Vector of source panel names
#'   \item \code{rule}: Combination function (e.g., \code{prod}, \code{sum})
#' }
#'
#' @return
#' Modified Seurat object with new assays:
#' \itemize{
#'   \item Each signal type becomes separate assay (names: 'iSensor_[signal]')
#'   \item Original data remains unchanged
#'   \item Meta-panels included if specified
#' }
#'
#' @examples
#' \dontrun{
#' # Basic analysis with default parameters
#' seurat_obj <- iSensor_pipeline(seurat_obj)
#'
#' # Advanced analysis with custom settings
#' seurat_obj <- iSensor_pipeline(
#'   seurat_obj,
#'   species = "human",
#'   signals = c("mean", "median_normed"),
#'   randPanels = 3,
#'   metaPanels = list(
#'     "MyMeta" = list(
#'       srcPanels = c("Panel1", "Panel2"),
#'       rule = function(x) sqrt(mean(x))
#'     )
#'   )
#' )
#' }
#'
#' @seealso
#' \code{\link{create_iSensor}}, \code{\link{iSensor_signal}}
#'
#' @importFrom Seurat GetAssayData CreateAssayObject DefaultAssay
#' @export
iSensor_pipeline <- function(seuratObject, species = NULL, hormone = NULL, type = NULL, signals = c('mean_normed'), seurLayer = "data",
                             usePanelPreset = NULL, presetPanelName = NULL,
                             randPanels = 2, randSize = c(200, 500), majortrend = TRUE,
                             usePanels = NULL, defaultAssay = NULL, useParallel = FALSE, nCores = NULL,
                             metaPanels = NULL) {
    if (class(seuratObject) != 'Seurat') {
        stop("Error: function must be applied to Seurat object")
    }

    expression_data <- as.matrix(GetAssayData(seuratObject, layer = seurLayer))
    expression_data[is.nan(expression_data)] <- 0
    iSensor_test <- create_iSensor(data = expression_data, species = species, hormone = hormone, type = type)
    if(!is.null(usePanelPreset)) {
        iSensor_test <- add_iSensor_panelSet(iSensor_test, panelsSet = usePanelPreset, setName = presetPanelName)
        iSensor_test <- set_iSensor_workPanel(iSensor_test, presetPanelName)
    }
    if(is.null(presetPanelName)) {
        presetPanelName <- names(iSensor_test$genePanelSets)[1]
    }
    cat(presetPanelName, ' panel preset will be used', '\n')
    iSensor_test <- set_iSensor_workPanel(iSensor_test, presetPanelName)
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

#' @title Set Active Gene Panels for iSensor Analysis
#' 
#' @description 
#' Selects a predefined set of gene panels to use for subsequent analysis steps
#' in the iSensor pipeline. The specified panel set must exist in the iSensor object.
#'
#' @param iSensor_obj An iSensor analysis object created by \code{create_iSensor}
#' @param panelSet Name of the panel set to activate (default: 'allSpec_allHorm_allType')
#'
#' @return The modified iSensor object with updated active gene panels
#'
#' @examples
#' \dontrun{
#' # After creating an iSensor object:
#' iSensor_obj <- set_iSensor_workPanel(iSensor_obj, panelSet = "my_panels")
#' }
#'
#' @export
set_iSensor_workPanel <- function(iSensor_obj, panelSet = 'allSpec_allHorm_allType') {
    if (!(panelSet %in% names(iSensor_obj$genePanelSets))) {
        stop("Error: no such panel in iSensor_obj$genePanelSets")
    }

    getPanel <- iSensor_obj$genePanelSets[[panelSet]]
    iSensor_obj$genePanels <- getPanel

    return(iSensor_obj)
}

#' @title Add Expression Data to iSensor Object
#'
#' @description 
#' Updates or replaces the expression data matrix in an existing iSensor analysis object.
#' The input matrix should have genes as rows and samples as columns.
#'
#' @param iSensor_obj An iSensor object created by \code{create_iSensor}
#' @param data Numeric matrix of expression values (genes x samples)
#'
#' @return The modified iSensor object with updated expression data
#'
#' @examples
#' # Create initial object
#' iSensor_obj <- create_iSensor(data = matrix(rnorm(100), nrow = 10))
#' 
#' # Add new data
#' new_data <- matrix(rnorm(200), nrow = 10)
#' iSensor_obj <- add_iSensor_data(iSensor_obj, new_data)
#'
#' @export
add_iSensor_data <- function(iSensor_obj, data) {
    iSensor_obj$exprData <- data
    return(iSensor_obj)
}

#' @title Add Sample Labels to iSensor Object
#'
#' @description 
#' Adds or updates sample labels in an iSensor object. The labels must match
#' the number of samples (columns) in the expression data matrix.
#'
#' @param iSensor_obj An iSensor object created by \code{create_iSensor}
#' @param labels Character vector of sample labels (length must match number of samples)
#'
#' @return The modified iSensor object with updated sample labels
#'
#' @examples
#' # Create iSensor object
#' iSensor_obj <- create_iSensor(data = matrix(rnorm(100), nrow = 10))
#' 
#' # Add sample labels
#' sample_labels <- paste0("Sample", 1:10)
#' iSensor_obj <- add_iSensor_labels(iSensor_obj, sample_labels)
#'
#' @export
add_iSensor_labels <- function(iSensor_obj, labels) {
    if (length(colnames(iSensor_obj$exprData)) != length(labels)) {
        stop("Error: length of 'labels' does not match number of samples in iSensor expr data")
    }
    iSensor_obj$sampleLabels <- labels

    return(iSensor_obj)
}

#' @title Add Gene Panels to iSensor Object
#'
#' @description 
#' Adds new gene panels or updates existing panels in an iSensor object. 
#' Existing panels with matching names will be updated unless `drop = TRUE`.
#'
#' @param iSensor_obj An iSensor object created by `create_iSensor`
#' @param panels Named list of gene panels to add/update (each element should be a character vector of genes)
#' @param drop Logical indicating whether to replace existing panels (FALSE) or keep both (TRUE) (default: FALSE)
#'
#' @return The modified iSensor object with updated gene panels
#'
#' @examples
#' # Create iSensor object
#' iSensor_obj <- create_iSensor(data = matrix(rnorm(100), nrow = 10))
#' 
#' # Add new gene panel
#' new_panel <- list("MyPanel" = c("Gene1", "Gene2", "Gene3"))
#' iSensor_obj <- add_iSensor_workPanel(iSensor_obj, new_panel)
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

#' @title Add Panel Set to iSensor Object
#' 
#' @description 
#' Adds a new set of gene panels to an iSensor object. The panel set can be either
#' a single panel or multiple panels combined into a named set.
#'
#' @param iSensor_obj An iSensor object (created by `create_iSensor`)
#' @param panelsSet List of gene panels to add (either a single panel or multiple panels)
#' @param setName Optional name for the panel set (autogenerated if NULL)
#'
#' @return The modified iSensor object with the new panel set added to `genePanelSets`
#'
#' @examples
#' # Create iSensor object
#' iSensor_obj <- create_iSensor(data = matrix(rnorm(1000), nrow = 100))
#' 
#' # Add a new panel set
#' new_panels <- list(
#'   "CellCycle" = c("CDK1", "CDK2", "CCNB1"),
#'   "Apoptosis" = c("BAX", "BCL2", "CASP3")
#' )
#' iSensor_obj <- add_iSensor_panelSet(iSensor_obj, new_panels, setName = "MyPanels")
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

#' @title Calculate Panel Signals for iSensor Object
#'
#' @description 
#' Computes aggregated signals (mean/median) for each gene panel in an iSensor object.
#' Supports normalization options and parallel processing for large datasets.
#'
#' @param iSensor_obj An iSensor object created by `create_iSensor`
#' @param panels Character vector of panel names to analyze (NULL uses all panels)
#' @param transform Aggregation method: 'mean' or 'median'
#' @param normed Logical indicating whether to normalize signals (default: FALSE)
#' @param normBy Normalization direction: 'cols' (by samples) or 'rows' (by genes)
#' @param metaPanels List specifying meta-panels to create (see Details)
#' @param doParallel Logical indicating whether to use parallel processing (default: FALSE)
#' @param nCores Number of cores for parallel processing (NULL uses available-1)
#'
#' @details
#' The `metaPanels` parameter should be a named list where each element contains:
#' \itemize{
#'   \item `srcPanels`: Vector of source panel names
#'   \item `rule`: Function to combine signals (e.g., `prod`, `sum`)
#' }
#'
#' @return The modified iSensor object with computed signals stored in `Signals` slot
#'
#' @examples
#' \dontrun{
#' # Calculate mean signals
#' iSensor_obj <- iSensor_signal(iSensor_obj, transform = "mean")
#'
#' # Calculate normalized median signals with parallel processing
#' iSensor_obj <- iSensor_signal(iSensor_obj, 
#'                              transform = "median",
#'                              normed = TRUE,
#'                              doParallel = TRUE)
#'
#' # With custom meta-panel
#' iSensor_obj <- iSensor_signal(iSensor_obj,
#'                              transform = "mean",
#'                              metaPanels = list(
#'                                "MyMeta" = list(
#'                                  srcPanels = c("Panel1", "Panel2"),
#'                                  rule = prod
#'                                )
#'                              ))
#' }
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
        # non_zero_variance_rows <- apply(expression_data, 1, var) > 0
        # expression_data <- expression_data[non_zero_variance_rows, , drop = FALSE]
        
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
        if(is.null(iSensor_obj$metaData)) {
            non_zero_variance_rows <- apply(iSensor_obj$exprData, 1, var) > 0
            iSensor_obj$exprData <- iSensor_obj$exprData[non_zero_variance_rows, , drop = FALSE]
            iSensor_obj$metaData <- 'obtained'
        }
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
