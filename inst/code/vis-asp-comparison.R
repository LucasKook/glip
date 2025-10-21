### Timing comparison with Hyttinen et al 2014 ASP
### LK 2025

library("tidyverse")
library("scales")
save <- TRUE
max_time <- 600
walltime <- 600

### List files
fin <- "./inst/results/asp-comparison/weak-large-n-large-fx-sparser"
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
  pivot_longer(c("glip", "asp"), names_to = "method", values_to = "time") |>
  mutate(time = as.numeric(time), time = ifelse(
    is.infinite(time), walltime, time
  ))

p1 <- ggplot(timings, aes(x = time, color = method)) +
  stat_ecdf(pad = FALSE) +
  theme_bw() +
  facet_wrap(mode ~ d, labeller = label_both) +
  labs(x = "runtime in seconds", y = "relative rank") +
  scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  theme(text = element_text(size = 13.5), legend.position = "top") +
  scale_color_brewer(palette = "Dark2", labels = c("asp" = "ASP", "glip" = "GLIP")) +
  coord_flip(xlim = c(min(timings$time) * 0.99, max_time)) +
  geom_vline(xintercept = walltime, linetype = 3, color = "darkred")
p1

if (save) {
  ggsave(file.path(fout, "timings.pdf"), p1, height = 6.5, width = 9.5)
}
