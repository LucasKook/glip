### Visualize results from benchmark
### LK 2025

library("tidyverse")
library("scales")
save <- TRUE
max_time <- 1800

### List files
fin <- "./inst/results/benchmark/2025-08-12/fixed"
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
  mutate(time = as.numeric(time))

p1 <- ggplot(timings, aes(x = time, color = method)) +
  stat_ecdf(pad = FALSE) +
  theme_bw() +
  facet_wrap(mode ~ d, labeller = label_both) +
  labs(x = "runtime in seconds", y = "relative rank") +
  scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  theme(text = element_text(size = 13.5), legend.position = "top") +
  scale_color_brewer(palette = "Dark2", labels = c("asp" = "ASP", "glip" = "GLIP")) +
  coord_flip(xlim = c(min(timings$time) * 0.99, min(max_time, max(timings$time))))

### Performance
pdat <- res |>
  pivot_longer(
    cols = c(
      "shd", "sep", "input_sep",
      "tail_prec", "tail_rec", "tail_f1", "head_prec", "head_rec", "head_f1"
    ),
    names_to = "metric",
    values_to = "value"
  )

lbs <- c(
  "dag" = "DAG", "admg" = "ADMG",
  "shd" = "SHD", "sep" = "SEP",
  # "input_sep" = "iSEP",
  "tail_prec" = "1 - Precision (tail)",
  "tail_rec" = "1 - Recall (tail)",
  "tail_f1" = "1 - F1 (tail)",
  "head_prec" = "1 - Precision (head)",
  "head_rec" = "1 - Recall (head)",
  "head_f1" = "1 - F1 (head)"
)

p2 <- ggplot(
  pdat |>
    filter(metric != "input_sep") |>
    mutate(metric = factor(metric, levels = names(lbs))),
  aes(x = ordered(d), y = value, color = method)
) +
  stat_summary(position = position_dodge(0.8), fun.data = "mean_se") +
  facet_grid(mode ~ metric, labeller = as_labeller(lbs)) +
  theme_bw() +
  labs(y = "score", x = "number of nodes") +
  theme(text = element_text(size = 13.5), legend.position = "top")

p3 <- res |>
  mutate(mode = toupper(mode)) |>
  pivot_longer(c("sep", "input_sep")) |>
  filter(method == "GLIP") |>
  ggplot(aes(x = d, y = value, color = name)) +
  stat_summary(geom = "line") +
  stat_summary() +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(mode ~ method + n) +
  theme_bw() +
  labs(x = "number of nodes", y = "SEP", color = element_blank()) +
  scale_color_brewer(palette = "Dark2", labels = c("input_sep" = "input", "sep" = "learned")) +
  theme(text = element_text(size = 13.5), legend.position = "top")
p3

if (save) {
  ggsave(file.path(fout, "timings.pdf"), p1, height = 6.5, width = 8)
  ggsave(file.path(fout, "performance.pdf"), p2, height = 6, width = 18)
  ggsave(file.path(fout, "separation.pdf"), p3, height = 5.5, width = 8)
}
