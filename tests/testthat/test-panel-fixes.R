test_that("Arabidopsis trans ARF panel contains all 23 ARFs, including ARF1", {
  panel <- InspectSensorPanel("ATH-aux-trans-ARF.rda", verbose = FALSE)
  expect_length(panel$genes, 23)
  expect_true("AT1G59750" %in% panel$genes)  # ARF1
  expect_false(anyDuplicated(panel$genes) > 0)
})

test_that("activator ARF panel stays the five class A ARFs", {
  panel <- InspectSensorPanel("ATH-aux-trans-A-ARF.rda", verbose = FALSE)
  expect_setequal(panel$genes, c("AT1G19850", "AT1G30330", "AT5G20730",
                                 "AT5G37020", "AT1G19220"))
})

test_that("documented filter codes select panels", {
  panelSet <- LoadSensors(setName = "test", species = "ATH", hormone = "aux",
                          random = FALSE)
  expect_true("ATH-aux-trans-ARF" %in% names(panelSet$panels))
  expect_true(all(startsWith(names(panelSet$panels), "ATH-aux-")))
  expect_error(LoadSensors(setName = "test", species = "AT"), "species")
})
