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
