### Figures and tables of the rerun (inst/code/rerun) from the collected CSVs,
### adapted from inst/code/vis-*.R: both GLIP encodings are shown as methods
### (GLIP-minlength, GLIP-cuts); certified = solver status OPTIMAL (the paper
### used runtime < time limit); plus a comparison of the two encodings.
### Usage (repository root):
###   Rscript inst/code/rerun/vis.R <out>/collected <figures_dir> <tables_dir> [encodings]
### encodings: "minlength,cuts" (default, both as methods) or one encoding,
### e.g. "minlength" to reproduce the figures and tables of the paper (that
### encoding labelled "GLIP", the other one left out, no encoding comparison)
suppressPackageStartupMessages({
  library("tidyverse")
  library("scales")
  library("knitr")
  library("kableExtra")
})
args <- commandArgs(trailingOnly = TRUE)
inp <- if (length(args) > 0) args[1] else "inst/results/rerun-2026-11/raw/rerun-2026-11/collected"
fig <- if (length(args) > 1) args[2] else "inst/results/rerun-2026-11/figures"
tabd <- if (length(args) > 2) args[3] else "inst/results/rerun-2026-11/tables"
encs <- if (length(args) > 3) strsplit(args[4], ",")[[1]] else c("minlength", "cuts")
stopifnot(all(encs %in% c("minlength", "cuts")))
single <- length(encs) == 1
glip_lab <- function(m) ifelse(single & m == paste0("GLIP-", encs[1]), "GLIP", m)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
dir.create(tabd, recursive = TRUE, showWarnings = FALSE)
rd <- function(e) {
  f <- file.path(inp, paste0(e, ".csv"))
  if (!file.exists(f)) {
    return(NULL)
  }
  x <- read_csv(f, show_col_types = FALSE, guess_max = 1e5)
  if ("encoding" %in% names(x)) x <- filter(x, encoding %in% encs)
  if ("method" %in% names(x)) {
    x <- filter(x, !grepl("^GLIP-", method) | method %in% paste0("GLIP-", encs))
  }
  x
}
theme_glip <- theme_bw() + theme(text = element_text(size = 13.5), legend.position = "top")
log_x <- scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))
cols <- c("GLIP" = "#1b9e77", "GLIP-minlength" = "#1b9e77", "GLIP-cuts" = "#d95f02", "ASP" = "#7570b3",
  "PC" = "#e7298a", "HPC" = "#66a61e", "FCI" = "#e6ab02", "NOTEARS" = "#a6761d",
  "R2SORT" = "#666666")
save_tex <- function(x, file) writeLines(as.character(x), file.path(tabd, file))
certified <- function(status) status %in% c("OPTIMAL", "USER_OBJ_LIMIT")

### Encoding comparison (min-length vs cuts) for any experiment
encoding_table <- function(res, by, file) {
  if (single) {
    return(invisible(NULL))
  }
  g <- res |>
    filter(method %in% c("GLIP-minlength", "GLIP-cuts")) |>
    group_by(across(all_of(by)), method) |>
    summarize(
      n = n(), certified = mean(certified(status)),
      median_runtime = median(runtime, na.rm = TRUE),
      mean_runtime = mean(runtime, na.rm = TRUE),
      mean_gap = mean(ifelse(certified(status), NA, mipgap), na.rm = TRUE),
      errors = sum(!is.na(error)), .groups = "drop"
    )
  write_csv(g, file.path(tabd, paste0(file, ".csv")))
  ### Paired: equal certified objectives, objective of time-limited fits
  w <- res |>
    filter(method %in% c("GLIP-minlength", "GLIP-cuts")) |>
    select(task, all_of(by), method, status, objval, runtime) |>
    pivot_wider(names_from = method, values_from = c(status, objval, runtime))
  p <- w |>
    group_by(across(all_of(by))) |>
    summarize(
      pairs = n(),
      both_certified = sum(certified(`status_GLIP-minlength`) & certified(`status_GLIP-cuts`)),
      equal_certified_obj = sum(certified(`status_GLIP-minlength`) & certified(`status_GLIP-cuts`) &
        abs(`objval_GLIP-minlength` - `objval_GLIP-cuts`) < 1e-6),
      cuts_faster = sum(`runtime_GLIP-cuts` < `runtime_GLIP-minlength`, na.rm = TRUE),
      cuts_lower_obj = sum(`objval_GLIP-cuts` < `objval_GLIP-minlength` - 1e-6, na.rm = TRUE),
      minlength_lower_obj = sum(`objval_GLIP-minlength` < `objval_GLIP-cuts` - 1e-6, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(p, file.path(tabd, paste0(file, "-paired.csv")))
  if (any(p$both_certified != p$equal_certified_obj)) {
    warning("certified objectives differ between the encodings in ", file)
  }
  invisible(list(summary = g, paired = p))
}

completion_table <- function(res, rows, file) {
  tt <- res |>
    filter(method %in% c("GLIP-minlength", "GLIP-cuts")) |>
    mutate(mode = toupper(mode)) |>
    group_by(across(all_of(rows)), method, d) |>
    summarize(frac = sprintf("%.3f", mean(certified(status))), .groups = "drop") |>
    mutate(method = glip_lab(method)) |>
    pivot_wider(names_from = "d", values_from = "frac")
  tex <- kable(tt, format = "latex", booktabs = TRUE,
    align = paste0(strrep("l", length(rows) + 1), strrep("r", ncol(tt) - length(rows) - 1))
  ) |> add_header_above(c(" " = length(rows) + 1, "d" = ncol(tt) - length(rows) - 1))
  save_tex(tex, file)
  tt
}

lbs <- c(
  "dag" = "DAG", "admg" = "ADMG", "chain" = "CG", "shd" = "SHD", "sep" = "SEP",
  "tail_f1" = "1 - F1 (tail)", "head_f1" = "1 - F1 (head)"
)

### A: ASP timing as in the paper: GLIP time = Gurobi runtime vs clingo's time
### (wall time of the clingo call). The wall time of the whole GLIP fit (model
### build incl. reading the cached matrix, complete warm start, solve; for
### cuts also the Python process) is reported alongside in asp-ratio.csv
asp <- rd("asp")
if (!is.null(asp)) {
  wt <- max(asp$walltime)
  asp <- asp |> mutate(setting = str_extract(task, "full|weak"), mode = toupper(mode),
    asp = ifelse(is.infinite(asp) | is.na(asp), wt, asp))
  tim <- bind_rows(
    asp |> transmute(setting, mode, d, iter, method = glip_lab(paste0("GLIP-", encoding)), time = runtime),
    asp |> distinct(setting, mode, d, iter, asp) |> transmute(setting, mode, d, iter, method = "ASP", time = asp)
  )
  for (s in unique(tim$setting)) {
    td <- filter(tim, setting == s)
    p1 <- ggplot(td, aes(x = time, linetype = method, color = ordered(d))) +
      stat_ecdf(pad = FALSE) + facet_wrap(~mode) + theme_glip + log_x +
      labs(x = "runtime in seconds", y = "relative rank", color = "d") +
      coord_flip() + geom_vline(xintercept = wt, linetype = 3, color = "darkred") +
      scale_color_viridis_d()
    rel <- asp |> filter(setting == s) |> mutate(rel = runtime / asp)
    p2 <- ggplot(rel, aes(x = rel, color = ordered(d))) +
      geom_vline(xintercept = 1, linetype = 3, color = "darkred") +
      stat_ecdf(pad = FALSE) + (if (single) facet_wrap(~mode) else facet_grid(encoding ~ mode)) +
      theme_glip + log_x +
      labs(x = "relative runtime GLIP/ASP", y = "relative rank", color = "d") +
      coord_flip() + scale_color_viridis_d()
    ggsave(file.path(fig, paste0("asp-timings-", s, ".pdf")), p1, height = 4.5, width = 8)
    ggsave(file.path(fig, paste0("asp-rel-timings-", s, ".pdf")), p2, height = 6, width = 8)
  }
  asp_res <- asp |> transmute(task, setting, mode, d, method = paste0("GLIP-", encoding),
    status, objval, runtime, mipgap, error = NA)
  encoding_table(asp_res, c("setting", "mode", "d"), "asp-encodings")
  write_csv(asp |> group_by(setting, mode, d, encoding) |>
    summarize(median_ratio = median(runtime / asp), glip_faster = mean(runtime < asp),
      median_ratio_wall = median(wall / asp), glip_faster_wall = mean(wall < asp),
      median_runtime = median(runtime), median_wall = median(wall), median_asp = median(asp),
      asp_at_limit = mean(asp >= wt), .groups = "drop"),
  file.path(tabd, "asp-ratio.csv"))
}

### B: simulations and D: chain graphs
for (e in c("simulation", "chain")) {
  res <- rd(e)
  if (is.null(res)) next
  res <- res |> mutate(setting = str_extract(task, "full|weak"))
  for (s in unique(res$setting)) {
    rs <- filter(res, setting == s)
    completion_table(rs, c("n", "mode", "ms")[c(TRUE, e == "simulation", s == "weak")],
      paste0(e, "-", s, "-completion.tex"))
    rs <- mutate(rs, method = glip_lab(method))
    p1 <- ggplot(rs, aes(x = time, color = method)) + stat_ecdf(pad = FALSE) +
      facet_wrap(mode ~ d, labeller = label_both) + theme_glip + log_x +
      labs(x = "runtime in seconds", y = "relative rank") + coord_flip() +
      scale_color_manual(values = cols)
    ggsave(file.path(fig, paste0(e, "-", s, "-timings.pdf")), p1, height = 6.5, width = 8)
    pdat <- rs |> pivot_longer(c("shd", "sep", "tail_f1", "head_f1"), names_to = "metric")
    tms_vals <- if (s == "weak") unique(rs$ms) else NA
    for (tn in unique(rs$n)) for (tms in tms_vals) {
      dd <- filter(pdat, n == tn, is.na(tms) | ms == tms)
      p2 <- ggplot(dd |> mutate(metric = factor(metric, levels = names(lbs))),
        aes(x = ordered(d), y = value, color = method, shape = method)) +
        stat_summary(position = position_dodge(0.8), fun.data = "mean_se", size = rel(0.3)) +
        facet_grid(mode ~ metric, labeller = as_labeller(lbs)) + theme_glip +
        labs(y = "score", x = "number of nodes") + scale_color_manual(values = cols) +
        guides(color = guide_legend(nrow = 1))
      sfx <- paste0("n-", tn, if (!is.na(tms)) paste0("_ms-", tms))
      ggsave(file.path(fig, paste0(e, "-", s, "-", sfx, "-performance.pdf")), p2, height = 4.5, width = 10)
      p3 <- rs |> filter(n == tn, is.na(tms) | ms == tms,
        method %in% c("GLIP", "GLIP-minlength", "GLIP-cuts", "PC", "FCI")) |>
        pivot_longer(c("sep", "input_sep")) |>
        ggplot(aes(x = d, y = value, color = name, linetype = method)) +
        stat_summary(geom = "line", fun = mean) + facet_wrap(~mode) + theme_glip +
        labs(x = "number of nodes", y = "SEP", color = NULL) +
        scale_color_brewer(palette = "Dark2", labels = c("input_sep" = "input", "sep" = "learned"))
      ggsave(file.path(fig, paste0(e, "-", s, "-", sfx, "-separation.pdf")), p3, height = 4, width = 7)
    }
  }
  encoding_table(res, c("setting", "mode", "n", "ms", "d"), paste0(e, "-encodings"))
}

### C: benchmarks (Tables 4 and 6)
bench <- rd("benchmark")
if (!is.null(bench)) {
  pm <- function(x, y) sprintf("%.2f (%.2f)", x, y)
  for (dm in unique(bench$d_max)) {
    b <- filter(bench, d_max == dm) |> mutate(k = ifelse(ms == d - 2, "d-2", as.character(ms)))
    glip <- b |> filter(grepl("^GLIP", method)) |> mutate(method = paste0(glip_lab(method), " (k=", k, ")"))
    comp <- b |> filter(!grepl("^GLIP", method), ms == d - 2)
    tab <- bind_rows(glip, comp) |>
      group_by(dataset, method) |>
      summarize(SHD = pm(mean(shd * d^2, na.rm = TRUE), sd(shd * d^2, na.rm = TRUE)),
        `1-SEP` = pm(mean(sep_k1, na.rm = TRUE), sd(sep_k1, na.rm = TRUE)),
        `2-SEP` = pm(mean(sep_k2, na.rm = TRUE), sd(sep_k2, na.rm = TRUE)),
        SEP = pm(mean(full_sep, na.rm = TRUE), sd(full_sep, na.rm = TRUE)),
        `F1 (HEAD)` = pm(mean(head_f1, na.rm = TRUE), sd(head_f1, na.rm = TRUE)),
        `F1 (TAIL)` = pm(mean(tail_f1, na.rm = TRUE), sd(tail_f1, na.rm = TRUE)),
        certified = mean(certified(status)), n = n(), .groups = "drop")
    write_csv(tab, file.path(tabd, paste0("benchmark-d", dm, ".csv")))
    tex <- kable(tab |> mutate(dataset = toupper(dataset)) |> select(-certified, -n),
      format = "latex", booktabs = TRUE) |>
      collapse_rows(columns = 1, latex_hline = "major")
    save_tex(tex, paste0("benchmark-d", dm, ".tex"))
  }
  encoding_table(bench, c("d_max", "dataset", "ms"), "benchmark-encodings")
}
cat("figures in", fig, "\ntables in", tabd, "\n")
