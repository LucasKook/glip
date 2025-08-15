### Visualize results from real data benchmarks
### LK 2025

library("tidyverse")
library("scales")
library("knitr")
save <- TRUE

### List files
fin <- "./inst/results/datasets/2025-08-15"
files <- list.files(fin, pattern = "all-tab.rds", full.names = TRUE, recursive = TRUE)

to_table <- function(x, y, digits = 2) {
  del <- paste0("%.", digits, "f")
  paste0(sprintf(del, x), " (", sprintf(del, y), ")")
}

res <- lapply(files, \(x) {
  sumtab <- readRDS(x)
  texout <- sumtab |>
    mutate(
      SHD = to_table(SHD, sdSHD),
      SEP = to_table(SEP, sdSEP),
      PREC = to_table(PREC, sdPREC),
      REC = to_table(REC, sdREC)
    ) |>
    mutate(Method = method, Dataset = toupper(str_extract(x, "alarm|asia|child|sachs|hepar2"))) |>
    select(Dataset, Method, SHD, SEP, PREC, REC)
}) |> bind_rows()
res

knitr::kable(res, format = "latex", booktabs = TRUE, digits = 2, align = "lrrrrr") |>
  kableExtra::collapse_rows(columns = 1, latex_hline = "major")

### ### Read files
### res <- lapply(files, \(x) {
###   sumtab <- readRDS(x)
###   knitr::kable(sumtab, format = "latex", booktabs = TRUE, digits = 2)
### })
### names(res) <- files
###
### res
