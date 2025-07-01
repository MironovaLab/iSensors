devtools::load_all()
library(Seurat)
exists('check_iSensorsPanelSet_obj')

testData <- readRDS('data/testSeurData.rds')
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

testPanel <- load_sensors(setName = 'testPanelSet', species = 'AT', hormone = 'cyt', customPanels = TRUE,
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
View(testSens$signals)
