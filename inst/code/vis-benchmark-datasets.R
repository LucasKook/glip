### Visualize results from real data benchmarks
### LK 2025

library("tidyverse")
library("scales")
library("knitr")
save <- TRUE

### List files
which <- "d8"
fin <- paste0("./inst/results/datasets/weak/", which, "ms", c(1, 2, 4))

to_table <- function(x, y, digits = 2) {
  del <- paste0("%.", digits, "f")
  paste0(sprintf(del, x), " (", sprintf(del, y), ")")
}

read_results <- function(dir) {
  files <- list.files(dir, pattern = "all.rds", full.names = TRUE, recursive = TRUE)
  res <- data.frame(file = files) |>
    mutate(
      res = map(files, readRDS),
      dataset = c(str_match(file, "alarm|asia|child|hepar2|sachs"))
    ) |>
    unnest(res)
  sumtab <- res |>
    group_by(dataset, method) |>
    summarize(
      SHD = mean(shd * d^2, na.rm = TRUE),
      SEP = mean(sep, na.rm = TRUE),
      FSEP = mean(full_sep, na.rm = TRUE),
      HF1 = mean(head_f1, na.rm = TRUE),
      TF1 = mean(tail_f1, na.rm = TRUE),
      sdSHD = sd(shd * d^2), sdSEP = sd(sep), sdFSEP = sd(full_sep),
      sdHF1 = sd(head_f1, na.rm = TRUE),
      sdTF1 = sd(tail_f1, na.rm = TRUE)
    ) |>
    dplyr::mutate_all(~ dplyr::case_when(is.nan(.x) ~ NA, !is.nan(.x) ~ .x)) |>
    ungroup()
  sumtab |>
    mutate(
      SHD = to_table(SHD, sdSHD),
      SEP = to_table(SEP, sdSEP),
      FSEP = to_table(FSEP, sdFSEP),
      HF1 = to_table(HF1, sdHF1),
      TF1 = to_table(TF1, sdTF1)
    ) |>
    mutate(Method = method, Dataset = toupper(dataset)) |>
    select(Dataset, Method, SHD, SEP, FSEP, HF1, TF1)
}

res1 <- read_results(fin[1]) |>
  rename(`1-SEP` = SEP, SEP = FSEP, `F1 (HEAD)` = HF1, `F1 (TAIL)` = TF1) |>
  mutate(Method = ifelse(Method == "GLIP", "GLIP (k=1)", Method))
res2 <- read_results(fin[2]) |>
  rename(`2-SEP` = SEP, SEP = FSEP, `F1 (HEAD)` = HF1, `F1 (TAIL)` = TF1) |>
  mutate(Method = ifelse(Method == "GLIP", "GLIP (k=2)", Method))
res3 <- read_results(fin[2]) |>
  rename(`d-2-SEP` = SEP, SEP = FSEP, `F1 (HEAD)` = HF1, `F1 (TAIL)` = TF1) |>
  mutate(Method = ifelse(Method == "GLIP", "GLIP (k=d-2)", Method))
res <- full_join(res1, res2) |>
  full_join(res3) |>
  arrange(Dataset, Method) |>
  select(Dataset, Method, SHD, `1-SEP`, `2-SEP`, SEP, `F1 (HEAD)`, `F1 (TAIL)`)

out <- knitr::kable(res, format = "latex", booktabs = TRUE, digits = 2, align = "lrrrrrrr") |>
  kableExtra::collapse_rows(columns = 1, latex_hline = "major")
if (save) {
  odir <- "./inst/tables/"
  if (!dir.exists(odir)) {
    dir.create(odir, recursive = TRUE)
  }
  write_lines(out, file.path(odir, paste0(which, ".tex")))
}
