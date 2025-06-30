devtools::load_all()
library(Seurat)
exists('check_iSensorsPanelSet_obj')

help(create_iSensor)
testSens <- create_iSensor()
class(testSens)


testPanel <- load_panel_from_rda('inst/extdata/geneSensors/AT_aux_cis_DR5_ARF1.rda')
print(testPanel$genes)
View(testPanel)
names(testPanel)
is.vector(testPanel$genes)
is.character('asff')
str(testPanel)

load_sensors()
