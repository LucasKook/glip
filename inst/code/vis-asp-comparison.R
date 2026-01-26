### Timing comparison with Hyttinen et al 2014 ASP
### LK 2025

library("tidyverse")
library("scales")
save <- TRUE
max_time <- 600
walltime <- 600

### List files
fin <- "./inst/results/asp-comparison/full"
nrow <- c(1, 2)[1 + str_detect(fin, "weak")]
nm <- c("full", "weak")[1 + str_detect(fin, "weak")]
fout <- str_replace(fin, "results", "figures")
if (!dir.exists(fout)) {
  dir.create(fout, recursive = TRUE)
}
files <- list.files(fin, pattern = "*.rds", full.names = TRUE)

### Read files
res <- tibble(file = files) |>
  mutate(data = map(file, ~ readRDS(.x))) |>
  unnest(data) |>
  mutate(mode = toupper(mode))

### Timings
timings <- res |>
  pivot_longer(c("glip", "asp"), names_to = "method", values_to = "time") |>
  mutate(time = as.numeric(time), time = ifelse(
    is.infinite(time), walltime, time
  ), method = toupper(method))

p1 <- ggplot(timings, aes(x = time, linetype = method, color = ordered(d))) +
  stat_ecdf(pad = FALSE) +
  theme_bw() +
  facet_wrap(~mode, labeller = label_both) +
  labs(x = "runtime in seconds", y = "relative rank", color = "d") +
  scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  theme(text = element_text(size = 13.5), legend.position = "top") +
  guides(color = guide_legend(nrow = nrow)) +
  coord_flip(xlim = c(min(timings$time) * 0.99, max_time)) +
  geom_vline(xintercept = walltime, linetype = 3, color = "darkred", size = 1) +
  scale_color_viridis_d()
p1

p3 <- ggplot(timings |> filter(d >= 6), aes(x = time, linetype = method, color = ordered(d))) +
  stat_ecdf(pad = FALSE) +
  theme_bw() +
  facet_wrap(~mode, labeller = label_both) +
  labs(x = "runtime in seconds", y = "relative rank", color = "d") +
  scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  theme(text = element_text(size = 13.5), legend.position = "top") +
  guides(color = guide_legend(nrow = nrow)) +
  coord_flip(xlim = c(min(timings$time) * 0.99, max_time)) +
  geom_vline(xintercept = walltime, linetype = 3, color = "darkred", size = 1) +
  scale_color_viridis_d()
p3

rel_time <- res |> mutate(rel = as.numeric(glip / ifelse(is.infinite(asp), walltime, asp)))
p2 <- ggplot(rel_time, aes(x = rel, color = ordered(d))) +
  geom_vline(xintercept = 1, linetype = 3, color = "darkred", size = 1) +
  stat_ecdf(pad = FALSE) +
  theme_bw() +
  facet_wrap(~mode, labeller = label_both) +
  labs(x = "relative runtime GLIP/ASP", y = "relative rank", color = "d") +
  scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  theme(text = element_text(size = 13.5), legend.position = "top") +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_brewer(palette = "Dark2", labels = c("asp" = "ASP", "glip" = "GLIP")) +
  coord_flip() +
  scale_color_viridis_d() +
  annotate(x = 100, y = 0.5, label = "ASP faster", size = 4, color = "gray", geom = "text") +
  annotate(x = 1 / 100, y = 0.5, label = "GLIP faster", size = 4, color = "gray", geom = "text") +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
p2

p4 <- ggplot(rel_time |> filter(d >= 6), aes(x = rel, color = ordered(d))) +
  geom_vline(xintercept = 1, linetype = 3, color = "darkred", size = 1) +
  stat_ecdf(pad = FALSE) +
  theme_bw() +
  facet_wrap(~mode, labeller = label_both) +
  labs(x = "relative runtime GLIP/ASP", y = "relative rank", color = "d") +
  scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  theme(text = element_text(size = 13.5), legend.position = "top") +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_brewer(palette = "Dark2", labels = c("asp" = "ASP", "glip" = "GLIP")) +
  coord_flip() +
  scale_color_viridis_d() +
  annotate(x = 100, y = 0.5, label = "ASP faster", size = 4, color = "gray", geom = "text") +
  annotate(x = 1 / 100, y = 0.5, label = "GLIP faster", size = 4, color = "gray", geom = "text") +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
p4

if (save) {
  ggsave(file.path(fout, paste0("timings-", nm, ".pdf")), p1, height = 4, width = 7.5, scale = 0.93)
  ggsave(file.path(fout, paste0("rel-timings-", nm, ".pdf")), p2, height = 4, width = 7.5, scale = 0.93)
  ggsave(file.path(fout, paste0("timings-", nm, "-large-d.pdf")), p4, height = 4, width = 7.5, scale = 0.93)
  ggsave(file.path(fout, paste0("rel-timings-", nm, "-large-d.pdf")), p4, height = 4, width = 7.5, scale = 0.93)
}
