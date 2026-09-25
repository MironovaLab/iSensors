#' @title The iSensorsCisTransPanelCreate function for creating CisTrans panels based on motif recognition in promoters.
#'
#' @description 
#' The iSensorsCisTransPanelCreate function performs the recognition of binding sites in promoters using positional weight matrices and generates a GenePanel object for a CisTrans type panel. 
#' Format for CisTrans panels is compatible with iSensors R package. 
#' The function receives (1) a set of promoters in FASTA format, (2) a positional probability matrix in Homer format, (3) a list of RNA-seq experiments with logFC and adjusted p-values for each gene (optional) and (4) the the corresponding trivial gene names in txt format (optional).
#' The function generates a list of three items. 
#' The first is a list genes containing gene IDs. 
#' The second is a data frame gene_metadata containing gene IDs, recognized sites, coordinates of each site in the genome, location (forward or reverse strand) and distance of each site relative to the TSS, and the corresponding short and full gene names. 
#' The third is a data frame panel_metadata containing the information about the species, gene panel type, gene panel description, creation date. 
#' The function writes the panel as an rda file in the user_panels directory (and creates it if this directory doesn't exist).
#' 
#' @param panel_name The name of the panel to create. The function will create a GenePanel object with this name and save it as rda file in the user_panels directory.
#' @param species Name of species for which the trans panel is created.
#' @param promoters_set Path to a FASTA file containing promoter sequences. 
#' @param ppm Path to a positional probability matrix file in Homer format. 
#' @param deg_list Path to a file containing logFC and adjusted p-value information obtained from a set of transcriptomic experiments that determine the differential expression status of a gene.
#' @param min_dataset_number The minimum number of transcriptomics experiments in which a gene is differentially expressed to be considered differentially expressed in the cis-trans panel. For example, if min_dataset_number = 1, then a gene must be differentially expressed in at least one experiment to be considered differentially expressed in the CisTrans panel.
#' @param trivial_names_file List of trivial gene names as a .txt file. (optional)
#' @param panel_type Type of panel. Could be cis, UP and DOWN.
#' @param transcriptomes_info 	Information about the used transcriptomes, written by the user.
#' 
#' @examples
#' \dontrun{
#' # Arabidopsis thaliana, cis panel.
#' 
#' iSensorsCisTransPanelCreate(panel_name = 'cis_panel',  
#' species = 'Arabidopsis thaliana',  
#' promoters_set = 'At_TAIR10_promoters.fas', 
#' ppm = 'ARF1.txt', panel_type = 'cis',  
#' trivial_names_file = 'at_trivial.txt')
#' 
#' # Arabidopsis thaliana, UP cis-trans panel.
#' 
#' iSensorsCisTransPanelCreate(panel_name = 'cis-trans_panel',  
#' species = 'Arabidopsis thaliana', 
#' promoters_set = 'At_TAIR10_promoters.fas', 
#' ppm = 'ARF1.txt', 
#' deg_list = 'auxin_degs.txt', 
#' min_dataset_number = 1, 
#' panel_type = 'UP', 
#' transcriptomes_info = 'Auxin 1h, auxin 4h')
#' }
#'
#' @import Biostrings
#' @import universalmotif
#' @import dplyr
#' @import magrittr
#' @import stringr
#' @import purrr
#' @export

iSensorsCisTransPanelCreate <- function(panel_name, species, promoters_set, ppm, deg_list, min_dataset_number, trivial_names_file, panel_type, transcriptomes_info)
{
  
  if(missing(promoters_set))
  {
    stop("Error: The promoters set FASTA file must be loaded")
  }
  
  if(missing(panel_name))
  {
    stop("Error: panel name must be specified")
  }
  
  if(missing(panel_type))
  {
    stop("Error: panel type must be specified")
  }
  
  promoters <- readDNAStringSet(promoters_set)
  
  header_check <- names(promoters)[1] %>% str_split(., '-') %>% unlist() %>% length()
  
  if (header_check != 5) {
    stop("Error: invalid header format in promoters set FASTA file")
  }
  
  rm(header_check)
  
  if(missing(ppm))
  {
    stop("Error: The positional probability matrix must be loaded")
  }
  
  ppm_name <- ppm
  ppm <- ppm
  ppm <- read.table(ppm, skip = 1, col.names = c('A', 'C', 'G', 'T'))
  
  if(ncol(ppm) != 4)
  {
    stop("Error: Invalid positional probability matrix format")
  }
  
  motif_length <- nrow(ppm)
  ppm <- log2((ppm + 0.00000001) / 0.25)
  ppm <- as.matrix(ppm) %>% t()
  ppm <- create_motif(ppm, type = 'PWM')
  
  if(!missing(deg_list)) #new
  {
    deg_list <- deg_list
    
    if(missing(min_dataset_number))
    {
      stop('Error: enter the minimum number of datasets in which the gene should be differentially expressed')
    }
    
    
    min_dataset_number <- min_dataset_number
    
    if(class(deg_list) == 'data.frame')
    {
      transcriptomes <- deg_list
    }
    
    else if(class(deg_list) == 'character')
    {
      transcriptomes <- read.table(deg_list, skip = 1)
    }
    
    number_of_datasets <- (ncol(transcriptomes) - 1)/2
    
    columns <- c()
    columns <- append(columns, 'gene')
    logFC <- c()
    padj <- c()
    
    for(i in 1:number_of_datasets)
    {
      logFC_single <- paste('logFC',i, sep = '_')
      padj_single <- paste('padj',i, sep = '_')
      columns <- append(columns, logFC_single)
      columns <- append(columns, padj_single)
      logFC <- append(logFC, logFC_single)
      padj <- append(padj, padj_single)
    }
    
    colnames(transcriptomes) <- columns
    
    padjes <- transcriptomes %>% select(gene, starts_with("padj"))
    FCs <- transcriptomes %>% select(gene, starts_with("logFC"))
    
    padjes[is.na(padjes)] <- 1
    FCs[is.na(FCs)] <- 0
    
    # Функция для определения статуса гена
    get_deg_status <- function(row)
    {
      count_de <- sum(row < 0.05)
      
      if (count_de >= min_dataset_number) 
      {
        return("DEG")
      } 
      else 
      {
        return("NON")
      }
    }
    
    padjes$status <- apply(padjes[, -1], 1, get_deg_status)
    FCs$status <- padjes$status
    
    
    get_expression_status <- function(row) {
      status <- row["status"]
      
      if (status == "NON") 
      {
        return("NON")
      } 
      
      else if (status == "DEG") 
      {
        expression_values <- as.numeric(row, na.rm=TRUE) %>% .[!is.na(.)]
        
        if (any(is.na(expression_values))) 
        {
          print(paste("NA found in row:", row["gene"]))
          print(expression_values)
          return("NON")
        }
        
        if (all(expression_values > 0)) 
        {
          return("UP")
        } 
        else if (all(expression_values < 0)) 
        {
          return("DOWN")
        } 
        else 
        {
          return("NON")
        }
      }
    }
    
    FCs$fc <- suppressWarnings(apply(FCs, 1, get_expression_status))
    
    transcriptomes$deg_status <- FCs$fc
  }
  
  score_unweighted <- motif_pvalue(ppm, pvalue = 0.0001, k = motif_length)
  scan <- scan_sequences(ppm, promoters, threshold.type = 'logodds.abs', threshold = score_unweighted, RC = TRUE,
                         no.overlaps = TRUE)
  scan <- as.data.frame(scan)
  motif_scanning_results <- data.frame(row.names = c(1:nrow(scan)))
  motif_scanning_results <- cbind(motif_scanning_results, 
                                  str_split_fixed(scan$sequence,"-", 5),
                                  scan$strand,
                                  scan$match,
                                  scan$start,
                                  scan$stop)
  
  colnames(motif_scanning_results) <- c("gene", "orientation", "chromosome", "promoter_sequence_start", "promoter_sequence_end", "strand", "match",
                                        "motif_start", "motif_end")
  
  motif_scanning_results$strand <- chartr("+", "1", motif_scanning_results$strand)
  motif_scanning_results$strand <- chartr("-", "0", motif_scanning_results$strand)
  
  to_tss <- c()
  absolute_motif_start <- c()
  absolute_motif_end <- c()
  to_tss_motif_start <- c()
  to_tss_motif_end <- c()
  
  i <- 1
  while(i <= nrow(scan))
  {
    if(as.numeric(motif_scanning_results$orientation[i]) == 1 & as.numeric(motif_scanning_results$strand[i] == 1))
    {single_absolute_motif_start <- as.numeric(motif_scanning_results$promoter_sequence_start[i]) + as.numeric(motif_scanning_results$motif_start[i]) - 1
    single_absolute_motif_end <- as.numeric(motif_scanning_results$promoter_sequence_start[i]) + as.numeric(motif_scanning_results$motif_end[i])
    single_tss_motif_start <- as.numeric(motif_scanning_results$promoter_sequence_end[i]) - as.numeric(single_absolute_motif_start)
    single_tss_motif_end <- as.numeric(motif_scanning_results$promoter_sequence_end[i]) - as.numeric(single_absolute_motif_end) + 1
    absolute_motif_start <- append(absolute_motif_start, single_absolute_motif_start)
    absolute_motif_end <- append(absolute_motif_end, single_absolute_motif_end)
    to_tss_motif_start <- append(to_tss_motif_start, single_tss_motif_start)
    to_tss_motif_end <- append(to_tss_motif_end, single_tss_motif_end)} # done
    
    else if(as.numeric(motif_scanning_results$orientation[i]) == 1 & as.numeric(motif_scanning_results$strand[i] == 0))
    {single_absolute_motif_start <- as.numeric(motif_scanning_results$promoter_sequence_start[i]) + as.numeric(motif_scanning_results$motif_end[i]) - 1
    single_absolute_motif_end <- as.numeric(motif_scanning_results$promoter_sequence_start[i]) + as.numeric(motif_scanning_results$motif_start[i])
    single_tss_motif_start <- as.numeric(motif_scanning_results$promoter_sequence_end[i]) - as.numeric(single_absolute_motif_end) 
    single_tss_motif_end <- as.numeric(motif_scanning_results$promoter_sequence_end[i]) - as.numeric(single_absolute_motif_start) + 1
    absolute_motif_start <- append(absolute_motif_start, single_absolute_motif_start)
    absolute_motif_end <- append(absolute_motif_end, single_absolute_motif_end)
    to_tss_motif_start <- append(to_tss_motif_start, single_tss_motif_start)
    to_tss_motif_end <- append(to_tss_motif_end, single_tss_motif_end)} # done
    
    else if(as.numeric(motif_scanning_results$orientation[i]) == 0 & as.numeric(motif_scanning_results$strand[i] == 1))
    {single_absolute_motif_start <- as.numeric(motif_scanning_results$promoter_sequence_start[i]) + as.numeric(motif_scanning_results$motif_start[i]) - 1
    single_absolute_motif_end <- as.numeric(motif_scanning_results$promoter_sequence_start[i]) + as.numeric(motif_scanning_results$motif_end[i])
    single_tss_motif_start <- as.numeric(single_absolute_motif_start) - as.numeric(motif_scanning_results$promoter_sequence_start[i]) + 1
    single_tss_motif_end <- as.numeric(single_absolute_motif_end) - as.numeric(motif_scanning_results$promoter_sequence_start[i]) 
    absolute_motif_start <- append(absolute_motif_start, single_absolute_motif_start)
    absolute_motif_end <- append(absolute_motif_end, single_absolute_motif_end)
    to_tss_motif_start <- append(to_tss_motif_start, single_tss_motif_start)
    to_tss_motif_end <- append(to_tss_motif_end, single_tss_motif_end)} # done
    
    else if(as.numeric(motif_scanning_results$orientation[i]) == 0 & as.numeric(motif_scanning_results$strand[i] == 0))
    {single_absolute_motif_start <- as.numeric(motif_scanning_results$promoter_sequence_start[i]) + as.numeric(motif_scanning_results$motif_end[i]) - 1
    single_absolute_motif_end <- as.numeric(motif_scanning_results$promoter_sequence_start[i]) + as.numeric(motif_scanning_results$motif_start[i])
    single_tss_motif_start <- as.numeric(single_absolute_motif_start) - as.numeric(motif_scanning_results$promoter_sequence_start[i])
    single_tss_motif_end <- as.numeric(single_absolute_motif_end) - as.numeric(motif_scanning_results$promoter_sequence_start[i]) + 1
    absolute_motif_start <- append(absolute_motif_start, single_absolute_motif_start)
    absolute_motif_end <- append(absolute_motif_end, single_absolute_motif_end)
    to_tss_motif_start <- append(to_tss_motif_start, single_tss_motif_start)
    to_tss_motif_end <- append(to_tss_motif_end, single_tss_motif_end)} # done
    
    i <- i + 1
  }
  
  motif_scanning_results <- cbind(motif_scanning_results, 
                                  absolute_motif_start,
                                  absolute_motif_end,
                                  to_tss_motif_start,
                                  to_tss_motif_end)
  
  
  
  motif_scanning_results <- motif_scanning_results[
    order(motif_scanning_results$chromosome, motif_scanning_results$absolute_motif_start, -motif_scanning_results$to_tss_motif_start),
  ]
  
  i <- 1
  repeater <- c()
  while(i <= nrow(motif_scanning_results))
  {value <- as.numeric(motif_scanning_results$absolute_motif_start[i + 1]) - as.numeric(motif_scanning_results$absolute_motif_start[i])
  repeater <- append(repeater, value)
  i <- i + 1}
  
  motif_scanning_results <- cbind(motif_scanning_results, repeater)
  motif_scanning_results <- subset(motif_scanning_results, motif_scanning_results$repeater > 0)
  motif_scanning_results$repeater <- NULL
  
  i <- 1
  to_tss <- c()
  while (i <= nrow(motif_scanning_results))
  {
    to_tss_start <- motif_scanning_results$to_tss_motif_start[i]
    to_tss_end <- motif_scanning_results$to_tss_motif_end[i]
    single_to_tss <- as.numeric(max(to_tss_start,to_tss_end)) - 
      ((as.numeric(max(to_tss_start,to_tss_end)) - as.numeric(min(to_tss_start,to_tss_end))) / 2)
    single_to_tss <- as.numeric(as.integer(single_to_tss))
    to_tss <- append(to_tss,single_to_tss)
    i <- i + 1
  }
  motif_scanning_results <- cbind(motif_scanning_results, to_tss)
  
  output <- data.frame(motif_scanning_results$gene, motif_scanning_results$chromosome, 
                       motif_scanning_results$absolute_motif_start,
                       motif_scanning_results$absolute_motif_end,
                       motif_scanning_results$strand,
                       motif_scanning_results$match,
                       motif_scanning_results$to_tss)
  
  colnames(output) <- c('gene', 'chromosome', 'cis-element_start', 'cis-element_end', 'strand', 'site', 'to_tss')
  
  if(!missing(trivial_names_file))
  {
    if(class(trivial_names_file) == 'data.frame')
    {
      trivial_names <- trivial_names_file
    }
    
    else if(class(trivial_names_file) == 'character')
    {
      tryCatch(
        {
          trivial_names <- suppressWarnings(read.table(trivial_names_file, header = T, sep = '\t', fill = T))
        },
        error = function(e) {
          message("Error while uploading trivial names list")
          stop()
        }
      )
    }
    
  }
  
  if(missing(trivial_names_file))
  {
    
    gene_name <- rep('none', nrow(output))
    gene_trivial_name <- rep('none', nrow(output))
    trivial_names <- data.frame(output$gene, gene_name, gene_trivial_name)
  }
  
  colnames(trivial_names) <- c('gene', 'gene_name', 'gene_full_name')
  
  output <- suppressWarnings(left_join(output, trivial_names, by = 'gene'))
  
  output[is.na(output)] <- 'none'
  
  
  if(!missing(deg_list))
  {
    out <- inner_join(output, transcriptomes, by = 'gene')
    out <- subset(out, out$deg_status != 'NON')
    
    if(panel_type == 'cis')
    {
      stop("Error: you can't create a cis panel with loaded lists of differentially expressed genes")
    }
    
    
    else if(panel_type == 'UP')
    {
      out <- subset(out, out$deg_status == 'UP')
    }
    
    else if(panel_type == 'DOWN')
    {
      out <- subset(out, out$deg_status == 'DOWN')
    }
    
    output <- data.frame(out$gene,
                         out$chromosome, 
                         out$'cis-element_start',
                         out$'cis-element_end',
                         out$strand,
                         out$site,
                         out$to_tss,
                         out$gene_name,
                         out$gene_full_name)
    
    
    colnames(output) <- c('GeneID', 'Chromosome', 'SiteStart', 'SiteEnd', 
                          'Strand', 'Site', 'ToTSS','GeneName', 'GeneFullName')
    
  }
  
  if(missing(deg_list))
  {
    colnames(output) <- c('GeneID', 'Chromosome', 'SiteStart', 'SiteEnd', 
                          'Strand', 'Site', 'ToTSS','GeneName', 'GeneFullName')
  }
  
  
  
  # creating of panel metadata
  
  date_of_analysis <- Sys.Date()
  
  promoter_length <- promoters[1] %>% nchar()
  
  species_name <- species
  
  if(!missing(deg_list))
  {
    transcriptomes_info <- transcriptomes_info
    transcriptomes_info <- paste(transcriptomes_info, collapse = ",")
    panel_type <- 'CisTrans'
  }
  
  if(missing(deg_list))
  {
    transcriptomes_info <- 'None'
    panel_type <- 'CisTrans'
  }
  
  panel_metadata <- data.frame(species_name, promoter_length, ppm_name, panel_type, transcriptomes_info, date_of_analysis)
  colnames(panel_metadata) <- c('Species', 'PromoterLength', 'MotifModelName', 
                                'PanelType', 'TranscriptomesExperimentInfo', 'DateCreated')
  
  panels <- output$'GeneID' %>% unique()
  
  output <- unique(output)
  
  panels_object <- list(
    genes = panels,
    genes_metadata = output,
    panel_metadata = panel_metadata
  )
  
  class(panels_object) <- 'GenePanel'
  
  if (!dir.exists('iSensors')) {
    dir.create('iSensors') }
  
  assign(panel_name, panels_object, envir = parent.frame())
  save(list = panel_name, file = paste0("iSensors/", panel_name, ".rda"))
  rm(list = panel_name, envir = parent.frame())
}


#' @title The iSensorsTransPanelCreate function generates a GenePanel object for a Trans type gene panel compatible with iSensors R package.
#'
#' @description 
#' The function receives a list of gene IDs as input (either as a vector or as a txt file), and the the corresponding trivial gene names in txt format (optional).
#' The function generates a list of three items. The first is a list genes containing gene IDs. 
#' The second is a data frame gene_metadata containing gene IDs and the corresponding short and full gene names. 
#' The third is a data frame panel_metadata containing the information about the species, gene panel type, gene panel description, creation date. 
#' The function writes the panel to an object in the current environment and saves it as an rda file in the current working directory.
#' 
#' @param panel_name The name of the panel to create. The function will create a GenePanel object with this name and save it as rda file in the user_panels directory.
#' @param species Name of species for which the trans panel is created.
#' @param gene_list Vector of genes. Can be supplied as a vector c('AT1G01010', 'AT1G01030', 'AT1G01040'), or as a .txt file.
#' @param trivial_names_file List of trivial gene names as a .txt file. (optional)
#' @param panel_description Description of the panel
#' 
#' @examples
#' \dontrun{
#' # Arabidopsis thaliana, input as vector.
#' iSensorsTransPanelCreate(panel_name = 'trans_panel', 
#' species = 'Arabidopsis thaliana', 
#' gene_list = c('AT1G01010', 'AT1G01030', 'AT1G01040'), 
#' panel_description = 'This is an example of trans panel')
#' 
#' # Arabidopsis thaliana, input as txt-file.
#' 
#' iSensorsTransPanelCreate(panel_name = 'trans_panel', 
#' species = 'Arabidopsis thaliana', 
#' gene_list = 'gene_list.txt', 
#' panel_description = 'This is an example of trans panel')
#' 
#' # Arabidopsis thaliana, input as txt-file, with trivial names list.
#' 
#' iSensorsTransPanelCreate(panel_name = 'trans_panel', 
#' species = 'Arabidopsis thaliana', 
#' gene_list = 'gene_list.txt', 
#' trivial_names_file = 'Arabidopsis_trivial_names_example.txt', 
#' panel_description = 'This is an example of trans panel')
#' }
#'
#' @import Biostrings
#' @import dplyr
#' @import magrittr
#' @import stringr
#' @export

iSensorsTransPanelCreate <- function(panel_name, species, gene_list, trivial_names_file, panel_description)
{
  
  if(missing(panel_name))
  {
    stop("Error: panel name must be specified")
  }
  
  if(class(gene_list) != 'character')
  {
    stop("Invalid input file format")
  }
  
  
  if(length(gene_list) == 1)
  {
    tryCatch(
      {
        suppressWarnings(genes_metadata <- read.table(gene_list) %>% as.data.frame())
      },
      error = function(e) {
        message("Invalid input file format")
        stop()
      }
    )
    
    if(ncol(genes_metadata) > 1)
    {
      stop("Invalid input file format")
    }
    
    colnames(genes_metadata) <- 'gene'
  }
  
  if(length(gene_list) > 1)
  {
    genes_metadata <- as.data.frame(gene_list)
    colnames(genes_metadata) <- 'gene'
  }
  
  if(!missing(trivial_names_file))
  {
    if(class(trivial_names_file) == 'data.frame')
    {
      trivial_names <- trivial_names_file
    }
    
    else if(class(trivial_names_file) == 'character')
    {
      tryCatch(
        {
          
          trivial_names <- suppressWarnings(read.table(trivial_names_file, header = T, sep = '\t', fill = T))
        },
        error = function(e) {
          message("Error while uploading trivial names list")
          stop()
        }
      )
    }
    
  }
  
  if(missing(trivial_names_file))
  {

    gene_name <- rep('none', nrow(genes_metadata))
    gene_trivial_name <- rep('none', nrow(genes_metadata))
    trivial_names <- data.frame(genes_metadata, gene_name, gene_trivial_name)
  }

  colnames(trivial_names) <- c('gene', 'gene_name', 'gene_trivial_name')
  
  genes_metadata <- left_join(genes_metadata, trivial_names, by = 'gene')
  genes_metadata[is.na(genes_metadata)] <- 'none'
  genes_metadata <- unique(genes_metadata)
  colnames(genes_metadata) <- c('GeneID', 'GeneName', 'GeneFullName')
  genes <- genes_metadata$GeneID %>% unique()
  date_of_analysis <- Sys.Date()
  panel_type <- 'Trans'
  panel_description <- panel_description
  
  if(missing(species))
  {
    species_name <- 'none' 
  }
  
  if(!missing(species))
  {
    species_name <- species 
  }
  
  panel_metadata <- data.frame(species_name, panel_type, panel_description, date_of_analysis)
  colnames(panel_metadata) <- c('Species', 'PanelType', 'PanelDescription', 'DateCreated')
  
  panels_object <- list(
    genes = genes,
    genes_metadata = genes_metadata,
    panel_metadata = panel_metadata
  )
  
  class(panels_object) <- 'GenePanel'
  
  
  if (!dir.exists('iSensors')) {
    dir.create('iSensors') }
  
  assign(panel_name, panels_object, envir = parent.frame())
  save(list = panel_name, file = paste0("iSensors/", panel_name, ".rda"))
  rm(list = panel_name, envir = parent.frame())

}
