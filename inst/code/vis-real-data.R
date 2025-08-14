### Visualize results from real data benchmarks
### LK 2025

library("tidyverse")
library("scales")
save <- TRUE

### List files
fin <- "./inst/results/datasets"
files <- list.files(fin, pattern = "*-sumtab.rds", full.names = TRUE, recursive = TRUE)

### Read files
res <- lapply(files, \(x) {
  sumtab <- readRDS(x)
  knitr::kable(sumtab, format = "latex", booktabs = TRUE, digits = 2)
})
names(res) <- files

res
