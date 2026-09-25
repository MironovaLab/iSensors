# Expression matrix: 8 genes x 4 cells; G8 has zero variance
expr <- matrix(c(0, 1, 2, 3,
                 1, 1, 0, 2,
                 2, 0, 0, 4,
                 0, 3, 1, 1,
                 1, 2, 3, 0,
                 4, 0, 2, 2,
                 0, 0, 1, 5,
                 5, 5, 5, 5),
               nrow = 8, byrow = TRUE,
               dimnames = list(paste0("G", 1:8), paste0("cell", 1:4)))

panel_set <- function(panels) {
  iSensors:::create_iSensorsPanelSet(name = "test", panelsList = panels)
}
five <- paste0("G", 1:5)

test_that("mean signal is the mean of the panel genes per cell", {
  res <- CalcSensors(expr, panelSet = panel_set(list(p = five)))
  expect_equal(unname(res$signals$mean["p", ]), unname(colMeans(expr[five, ])))
})

test_that("median signal ignores zeros", {
  res <- CalcSensors(expr, panelSet = panel_set(list(p = five)), signals = "median")
  expected <- apply(expr[five, ], 2, function(x) median(x[x != 0]))
  expect_equal(unname(res$signals$median["p", ]), unname(expected))
})

test_that("genes with zero variance do not count towards the panel", {
  # six genes listed, but G8 is constant: five genes are scored
  res <- CalcSensors(expr, panelSet = panel_set(list(p = c(five, "G8"))))
  expect_equal(unname(res$signals$mean["p", ]), unname(colMeans(expr[five, ])))
})

test_that("panels with fewer than 3 detected genes are skipped with a warning", {
  panels <- list(ok = five,
                 three = c("G1", "G2", "G3"),
                 small = c("G1", "G2"),
                 constant = c("G1", "G2", "G8"),                 # 2 after G8 is left out
                 absent = c("G1", "X1", "X2", "X3", "X4"))
  expect_warning(res <- CalcSensors(expr, panelSet = panel_set(panels)),
                 "3 panel\\(s\\) skipped.*small, constant, absent")
  expect_equal(rownames(res$signals$mean), c("ok", "three"))
})

test_that("an error is raised when no panel can be scored", {
  expect_error(suppressWarnings(CalcSensors(expr, panelSet = panel_set(list(s = c("G1", "G2"))))),
               "No panel has at least 3 genes")
})

test_that("mean is the default signal", {
  res <- CalcSensors(expr, panelSet = panel_set(list(p = five)))
  expect_named(res$signals, "mean")
})

test_that("removed normed signals give an informative error", {
  ps <- panel_set(list(p = five))
  expect_error(CalcSensors(expr, panelSet = ps, signals = "mean_normed"), "removed in iSensors 1.3.0")
  expect_error(CalcSensors(expr, panelSet = ps, signals = "median_normed"), "removed in iSensors 1.3.0")
  expect_error(CalcSensors(expr, panelSet = ps, signals = "sum"), "allowed signals")
})

test_that("Seurat objects get one assay per signal", {
  skip_if_not_installed("Seurat")
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = as(expr * 10, "dgCMatrix")))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- CalcSensors(obj, panelSet = panel_set(list(p = five)), signals = c("mean", "median"))
  expect_true(all(c("iSensors_mean", "iSensors_median") %in% names(obj@assays)))
  expect_false("iSensors_mean_normed" %in% names(obj@assays))
})
