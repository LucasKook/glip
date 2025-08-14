### Visualize results from real data benchmarks
### LK 2025

library("tidyverse")
library("scales")
save <- TRUE

### List files
fin <- "./inst/results/datasets/admg-dmax-6"
files <- list.files(fin, pattern = "all-tab.rds", full.names = TRUE, recursive = TRUE)

res <- lapply(files, readRDS)
names(res) <- files
res

### ### Read files
### res <- lapply(files, \(x) {
###   sumtab <- readRDS(x)
###   knitr::kable(sumtab, format = "latex", booktabs = TRUE, digits = 2)
### })
### names(res) <- files
###
### res
