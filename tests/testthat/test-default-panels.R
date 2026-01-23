test_that("Default gene sensor panels are valid", {
  
  panels_info <- ListSensorPanels()
  
  # 1. Default panels exist
  expect_true(length(panels_info$default) > 0)
  
  for (panel_file in panels_info$default) {
    
    panel <- InspectSensorPanel(panel_file, verbose = FALSE)
    
    # 2. Panel is a list
    expect_true(is.list(panel))
    
    # 3. Mandatory fields exist
    expect_true("genes" %in% names(panel))
    expect_true("name" %in% names(panel))
    
    # 4. Genes field
    expect_true(is.character(panel$genes))
    expect_true(length(panel$genes) > 0)
    
    # 5. Name consistency
    expect_true(is.character(panel$name))
    expect_equal(length(panel$name), 1)
  }
})
