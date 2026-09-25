# Empty temporary directory to work in (the panel functions write iSensors/)
make_dir <- function() {
  d <- tempfile("isensors-"); dir.create(d); d
}

test_that("all shipped panels have at least 3 genes", {
  for (f in ListSensorPanels()$default) {
    panel <- InspectSensorPanel(f, verbose = FALSE)
    expect_gte(length(unique(panel$genes)), 3)
  }
})

test_that("trans panels with fewer than 3 genes cannot be created", {
  d <- make_dir(); old <- setwd(d); on.exit(setwd(old))
  expect_error(
    iSensorsTransPanelCreate(panel_name = "tooSmall", species = "Arabidopsis thaliana",
                             gene_list = c("AT1G01010", "AT1G01030"),
                             panel_description = "two genes"),
    "has 2 gene\\(s\\); a panel needs at least 3"
  )
  expect_false(file.exists(file.path("iSensors", "tooSmall.rda")))

  genes <- c("AT1G01010", "AT1G01030", "AT1G01040")
  threeGenes <- "my own object"
  panel <- iSensorsTransPanelCreate(panel_name = "threeGenes", species = "Arabidopsis thaliana",
                                    gene_list = genes, panel_description = "three genes")
  expect_true(file.exists(file.path("iSensors", "threeGenes.rda")))
  expect_setequal(panel$genes, genes)               # returned invisibly
  expect_identical(threeGenes, "my own object")     # caller's object untouched
  expect_setequal(InspectSensorPanel("threeGenes.rda", verbose = FALSE)$genes, genes)
})

test_that("custom panels with fewer than 3 genes are refused when loading", {
  d <- make_dir(); old <- setwd(d); on.exit(setwd(old))
  dir.create("iSensors")
  small <- list(genes = c("AT1G01010", "AT1G01030"))
  save(small, file = file.path("iSensors", "small.rda"))
  expect_error(LoadSensors(setName = "t", defaultPanels = FALSE, customPanels = TRUE,
                           random = FALSE),
               "small.*at least 3")
})

test_that("random panels must have at least 3 genes", {
  expect_error(LoadSensors(setName = "t", species = "ATH", hormone = "aux",
                           randomInfo = list(n = 2, sizes = c(2, 100), majortrend = TRUE)),
               "at least 3")
})
