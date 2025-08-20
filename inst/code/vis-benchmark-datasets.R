### Visualize results from real data benchmarks
### LK 2025

library("tidyverse")
library("scales")
library("knitr")
save <- TRUE

### List files
fin <- paste0("./inst/results/datasets/2025-08-18/d8ms", 1:2)

to_table <- function(x, y, digits = 2) {
  del <- paste0("%.", digits, "f")
  paste0(sprintf(del, x), " (", sprintf(del, y), ")")
}

read_results <- function(dir) {
  files <- list.files(dir, pattern = "all-tab.rds", full.names = TRUE, recursive = TRUE)
  lapply(files, \(x) {
    sumtab <- readRDS(x)
    texout <- sumtab |>
      mutate(
        SHD = to_table(SHD, sdSHD),
        SEP = to_table(SEP, sdSEP),
        FSEP = to_table(FSEP, sdFSEP),
        PREC = to_table(PREC, sdPREC),
        REC = to_table(REC, sdREC)
      ) |>
      mutate(Method = method, Dataset = toupper(str_extract(x, "alarm|asia|child|sachs|hepar2"))) |>
      select(Dataset, Method, SHD, SEP, FSEP, PREC, REC)
  }) |> bind_rows()
}

res1 <- read_results(fin[1]) |>
  rename(`1-SEP` = SEP, SEP = FSEP, FDR = PREC, FNR = REC)
res2 <- read_results(fin[2]) |>
  rename(`2-SEP` = SEP, SEP = FSEP, FDR = PREC, FNR = REC)
res <- full_join(res1, res2) |>
  arrange(Dataset, Method) |>
  select(Dataset, Method, SHD, `1-SEP`, `2-SEP`, SEP, FDR, FNR)

out <- knitr::kable(res, format = "latex", booktabs = TRUE, digits = 2, align = "lrrrrrrr") |>
  kableExtra::collapse_rows(columns = 1, latex_hline = "major")
write_lines(out, file.path(fin[1], "table.tex"))
