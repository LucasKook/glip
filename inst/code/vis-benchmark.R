### Visualize results from benchmark
### LK 2025

library("tidyverse")
save <- TRUE

### List files
fin <- "./inst/results/benchmark/2025-07-29/oracle"
fout <- str_replace(fin, "results", "figures")
if (!dir.exists(fout)) {
  dir.create(fout, recursive = TRUE)
}
files <- list.files(fin, pattern = "*.rds", full.names = TRUE)

### Read files
res <- tibble(file = files) |>
  mutate(data = map(file, ~ readRDS(.x))) |>
  unnest(data)

### Timings
timings <- res |>
  group_by(method, n, d, ms, mode) |>
  mutate(time = as.numeric(time)) |>
  summarize(q50 = median(time), q25 = quantile(time, 0.25), q75 = quantile(time, 0.75))

p1 <- ggplot(timings, aes(x = ordered(d), y = q50, color = method, ymin = q25, ymax = q75)) +
  geom_pointrange(position = position_dodge(0.8)) +
  facet_grid(~mode, labeller = as_labeller(c("dag" = "DAG", "admg" = "ADMG"))) +
  theme_bw() +
  labs(y = "median runtime in seconds", x = "number of nodes") +
  scale_y_log10() +
  theme(text = element_text(size = 13.5), legend.position = "top")

### Performance
pdat <- res |>
  pivot_longer(
    cols = c("shd", "sep", "tail_prec", "tail_rec", "tail_f1", "head_prec", "head_rec", "head_f1"),
    names_to = "metric",
    values_to = "value"
  )

lbs <- c(
  "dag" = "DAG", "admg" = "ADMG",
  "shd" = "SHD", "sep" = "SEP",
  "tail_prec" = "1 - Precision (tail)",
  "tail_rec" = "1 - Recall (tail)",
  "tail_f1" = "1 - F1 (tail)",
  "head_prec" = "1 - Precision (head)",
  "head_rec" = "1 - Recall (head)",
  "head_f1" = "1 - F1 (head)"
)

p2 <- ggplot(
  pdat |> mutate(metric = factor(metric, levels = names(lbs))),
  aes(x = ordered(d), y = value, color = method)
) +
  stat_summary(position = position_dodge(0.8), fun.data = "mean_se") +
  facet_grid(mode ~ metric, labeller = as_labeller(lbs)) +
  theme_bw() +
  labs(y = "score", x = "number of nodes") +
  theme(text = element_text(size = 13.5), legend.position = "top")

print(p1)
print(p2)

if (save) {
  ggsave(file.path(fout, "timings.pdf"), p1, height = 4, width = 6)
  ggsave(file.path(fout, "performance.pdf"), p2, height = 6, width = 18)
}
