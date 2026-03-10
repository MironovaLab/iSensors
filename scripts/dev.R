devtools::load_all()
library(iSensors)
packageVersion('iSensors')
help("iSensors")
help('LoadSensors')
help('CalcSensors')
library(Seurat)
exists('check_iSensorsPanelSet_obj')
help(LoadSensors)

testData <- readRDS('testData/testSeurData.rds')
testMatr <- GetAssayData(object = testData, assay = "RNA", layer = "counts")
testMatr

help(create_iSensor)
testSens <- create_iSensor()
class(testSens)

available_panels <- c("random12", "random2", "majortrend1", "panelX")
metaPanelsTest <- list(
  meta1 = list(srcPanels = c("random1", "random2"), rule = mean),
  meta2 = list(srcPanels = c("majortrend", "panelX"), rule = median)
)
check_metaPanels(metaPanelsTest, available_panels)


testPanel <- load_panel_from_rda('inst/extdata/geneSensors/AT_aux_cis_DR5_ARF1.rda')
print(testPanel$genes)
View(testPanel)
names(testPanel)
is.vector(testPanel$genes)
is.character('asff')
str(testPanel)

testPanel <- load_sensors(setName = 'testPanelSet', customPanels = TRUE,
                          randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), majortrend = FALSE),
                          metaPanels = list(
                            'meta1' = list('srcPanels' = c("random1", "random2"), rule = mean),
                            'meta2' = list('srcPanels' = c("majortrend", "panelX"), rule = median)))
is.vector(testPanel[['panels']][['AT_aux_cistrans_IR8_ARF5syn_down']])
str(testPanel)

testPanel <- LoadSensors(setName = 'testPanelSet', species = 'AT', hormone = 'cyt', customPanels = TRUE,
                          randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), majortrend = TRUE),
                          metaPanels = list(
                            'meta1' = list('srcPanels' = c("AT_aux_cis_DR5_ARF1", "AT_aux_cistrans_DR5_ARF5_2_up"), rule = mean),
                            'meta2' = list('srcPanels' = c("AT_aux_cis_DR5_ARF1", "AT_aux_cistrans_DR5_ARF5_2_up"), rule = prod))
                          )
View(testPanel)
testSens <- create_iSensors(data = testMatr, panelSet = testPanel)
str(testSens)
testSens <- iSensor_signal(iSensor_obj = testSens, transform = 'mean', normed = TRUE)
testSens <- iSensor_signal(iSensor_obj = testSens, transform = 'median', normed = FALSE)
View(testSens$signals$mean_normed)

result <- CalcSensors(testData,
                      seurLayer = "data",
                      panelSet = testPanel,
                      signals = c("mean_normed", "median"))
result <- CalcSensors(testMatr,
                      panelSet = testPanel,
                      signals = c("mean_normed", "median"))
View(result)
slotNames(result)
View(result@assays$iSensors_median$counts)

save_panel_to_rda <- function(panel, new_name, folder_path) {
  # Проверка, что panel - валидный объект панели
  check_iSensorsPanel_obj(panel)
  
  # Создаем временное окружение
  temp_env <- new.env()
  
  # Используем panel$name как имя переменной в .rda файле
  panel$name <- new_name
  variable_name <- panel$name
  
  # Сохраняем панель в окружение под нужным именем
  assign(new_name, panel, envir = temp_env)
  
  # Сохраняем в файл
  save(list = new_name, 
       file = paste0(folder_path, new_name, '.rda'),
       envir = temp_env)
  
  invisible(TRUE)
}

create_panel_folder <- function(folder_path) {
  # Проверяем, существует ли уже папка
  if (!dir.exists(folder_path)) {
    # Создаем папку (recursive = TRUE создает все промежуточные папки)
    dir.create(folder_path, recursive = TRUE)
    message("Папка создана: ", folder_path)
    return(TRUE)
  } else {
    message("Папка уже существует: ", folder_path)
    return(FALSE)
  }
}

change_panel_names <- function() {
  create_panel_folder('geneSensors_new/')
  panels_info <- ListSensorPanels()
  
  for (panel_file in panels_info$default) {
    panel_file_no_extn <- gsub("\\.rda$", "", panel_file)
    panel <- InspectSensorPanel(panel_file, verbose = FALSE)
    
    if (panel_file_no_extn != panel$name) {
      message('Renaming ', panel$name, ' with ', panel_file_no_extn)
      panel$name <- panel_file_no_extn
    }
    if (grepl("cistrans", panel$name)) {
      panel$name <- gsub("cistrans", "reg", panel$name)
    }
    save_panel_to_rda(panel, panel$name, 'geneSensors_new/')
    
  }
}
change_panel_names()

system.file("extdata/geneSensors", package = "iSensors")