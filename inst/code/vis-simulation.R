### Visualize results from simulation
### LK 2025

library("tidyverse")
library("scales")
library("knitr")
library("kableExtra")
save <- TRUE
max_time <- 600
walltime <- 300

folders <- c("full", "weak", "k2", "full-hpc")[4]

lapply(folders, \(which) {
  ### List files
  fin <- paste0("./inst/results/benchmark/", which)
  fout <- str_replace(fin, "results", "figures")
  if (!dir.exists(fout)) {
    dir.create(fout, recursive = TRUE)
  }
  files <- list.files(fin, pattern = "*.rds", full.names = TRUE)

  ### Read files
  res <<- tibble(file = files) |>
    mutate(data = map(file, ~ readRDS(.x))) |>
    unnest(data)

  ### Timings
  timings <- res |>
    group_by(method, n, d, ms, mode) |>
    mutate(time = as.numeric(time))

  tt <- timings |>
    filter(method == "GLIP") |>
    mutate(mode = toupper(mode)) |>
    summarize(frac_optimal = sprintf("%.3f", mean(time < walltime))) |>
    ungroup() |>
    select(n, mode, d, frac_optimal) |>
    pivot_wider(names_from = "d", values_from = "frac_optimal")
  print(tt)
  tex <<- tt |>
    kable(
      format = "latex", booktabs = TRUE,
      align = paste0("l", paste0(rep("r", ncol(tt) - 1), collapse = ""))
    ) |>
    add_header_above(c(" " = 2, "d" = ncol(tt) - 2))

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
    "tail_f1" = "1 - F1 (tail)",
    "head_f1" = "1 - F1 (head)"
  )

  lapply(unique(pdat$n), \(tn) {
    p2 <<- ggplot(
      pdat |>
        filter(metric %in% names(lbs), n == tn) |>
        mutate(metric = factor(metric, levels = names(lbs))),
      aes(x = ordered(d), y = value, color = method, shape = method)
    ) +
      stat_summary(position = position_dodge(0.8), fun.data = "mean_se", size = rel(0.3)) +
      facet_grid(mode ~ metric, labeller = as_labeller(lbs)) +
      theme_bw() +
      labs(y = "score", x = "number of nodes") +
      theme(
        text = element_text(size = 13.5), legend.position = "top",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
      )

    p3 <<- res |>
      mutate(mode = toupper(mode)) |>
      pivot_longer(c("sep", "input_sep")) |>
      filter(method %in% c("GLIP", "PC", "FCI")[1]) |>
      ggplot(aes(x = d, y = value, color = name)) +
      stat_summary(geom = "line") +
      stat_summary() +
      geom_abline(slope = 1, intercept = 0) +
      facet_grid(mode ~ method + n) +
      theme_bw() +
      labs(x = "number of nodes", y = "SEP", color = element_blank()) +
      scale_color_brewer(palette = "Dark2", labels = c("input_sep" = "input", "sep" = "learned")) +
      theme(text = element_text(size = 13.5), legend.position = "top")

    if (save) {
      ggsave(file.path(fout, paste0("n-", tn, "_timings-", which, ".pdf")), p1, height = 6.5, width = 8, bg = "transparent")
      ggsave(file.path(fout, paste0("n-", tn, "_performance-", which, ".pdf")), p2, height = 5.5, width = 9, bg = "transparent")
      ggsave(file.path(fout, paste0("n-", tn, "_separation-", which, ".pdf")), p3, height = 5.5, width = 8, bg = "transparent")
      save_kable(tex, file.path(fout, "tab-completion.tex"))
    }
  })
})
