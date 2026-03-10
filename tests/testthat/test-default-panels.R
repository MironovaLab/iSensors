test_that("Default gene sensor panels are valid", {
  
  panels_info <- ListSensorPanels()
  
  # 1. Default panels exist
  expect_true(length(panels_info$default) > 0)
  
  mismatched_panels <- c()
  
  for (panel_file in panels_info$default) {
    
    panel_file_no_extn <- gsub("\\.rda$", "", panel_file)
    panel <- InspectSensorPanel(panel_file, verbose = FALSE)
    
    if (panel_file_no_extn != panel$name) {
      mismatched_panels <- c(mismatched_panels, panel_file)
    }
    
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
  
  expect_equal(length(mismatched_panels), 0, 
               info = paste("Diffrent panel filename and panel$name in panels:", 
                            paste(mismatched_panels, collapse = ", ")))
})
