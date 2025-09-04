testList <- list(1, 2, 3, 4)
names(testList) <- c('a', 'b', 'c', 'd')
for (i in seq_along(testList)) {
  name <- names(testList)[i]
  value <- testList[[i]]
  cat("Name:", name, "Value:", value, "\n")
}

testList <- c(FALSE, FALSE, FALSE)
all(!testList)

a <- list(species = "human")
b <- list(hormone = "estrogen")

combined <- c(a, b)
combined

my_function <- function(method = c("mean", "median")) {
  method <- match.arg(method)
  print(paste("Selected method:", method))
}
my_function(c('mean', 'median'))

packageVersion('iSensors')

load_panel_from_rda_all <- function(file_path) {
  panelEnv <- new.env()
  load(file_path, envir = panelEnv)
  panelName <- ls(envir = panelEnv)
  if (length(panelName) != 1) {
    warning("Expected one object per panel file, got: ", paste(panelName, collapse = ", "))
  }
  panel <- get(panelName, envir = panelEnv)
  panel$name <- panelName
  check_iSensorsPanel_obj(panel)
  return(panel)
}

pamel_PAT <- load_panel_from_rda_all(paste0(panelsDir_tmp, '/AT_aux_trans_EffluxInflux.rda'))
length(pamel_PAT)
View(pamel_PAT[[3]])
write.table(
  pamel_PAT[[1]],
  file = "C:/Users/user/Desktop/ИЦиГ/Sensor proj/PAT_1_slot.txt",
  sep = "\t",          # Разделитель — табуляция
  row.names = FALSE,   # Без номеров строк
  quote = FALSE        # Без кавычек вокруг строк
)

panelsDir_tmp <- system.file("extdata/geneSensors", package = "iSensors")
panelsFiles_tmp <- list.files(path = panelsDir_tmp, pattern = '\\.rda$', full.names = FALSE)
panels_tmp <- load_panel_files(panelsDir = panelsDir_tmp, files = panelsFiles_tmp)

str(panels_tmp)

testPanel <- LoadSensors(setName = 'testPanelSet', species = 'AT',
                         random = FALSE)

names(testPanel$panels)
testPanel$panels[['AT_aux_trans_EffluxInflux']]
length(testPanel$panels[['AT_aux_trans_EffluxInflux']])

library(readxl)
file_path <- "C:/Users/user/Desktop/ИЦиГ/Sensor proj/2025-08-21_Supplementary_Tables.xlsx"
df <- read_excel(path = file_path, sheet = "S3", skip = 1)
View(df)
colnames(df)

panel_type <- tolower(df[["Panel type"]])
panel_type[panel_type == "reg"] <- "cistrans"
panel_name <- gsub("-", "_", df[["Panel name"]])
panel_name[panel_name == "PAT"] <- "EffluxInflux"
names_list <- paste(df[["Species"]], "aux", panel_type, panel_name, sep = "_")
names_list
df$package_name <- names_list
df
df[["Number of genes"]] <- sapply(df[["package_name"]], function(pkg) {
  # Проверяем, что панель существует в testPanel$panels
  if (!is.null(testPanel$panels[[pkg]])) {
    length(testPanel$panels[[pkg]])
  } else {
    NA  # Если панели нет, ставим NA
  }
})
df
write.table(
  df,
  file = "C:/Users/user/Desktop/ИЦиГ/Sensor proj/2025-08-21_Supplementary_Tables_S3.txt",
  sep = "\t",          # Разделитель — табуляция
  row.names = FALSE,   # Без номеров строк
  quote = FALSE        # Без кавычек вокруг строк
)
