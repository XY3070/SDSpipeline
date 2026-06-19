#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(args) {
  opts <- list(
    title = "chr1 SDS compare",
    width = 15,
    height = 5.2,
    dpi = 400,
    point_cex = 0.28
  )

  opts[["panel"]] <- character()
  opts[["window"]] <- character()

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected argument: %s", arg), call. = FALSE)
    }
    key <- substring(arg, 3)
    if (key %in% c("panel", "window")) {
      if (i == length(args)) stop(sprintf("Missing value for argument: %s", arg), call. = FALSE)
      opts[[key]] <- c(opts[[key]], args[[i + 1]])
      i <- i + 2
      next
    }
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", arg), call. = FALSE)
    opts[[key]] <- args[[i + 1]]
    i <- i + 2
  }

  required <- c("output-prefix")
  missing <- required[!vapply(required, function(x) !is.null(opts[[x]]), logical(1))]
  if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(sprintf("--%s", missing), collapse = ", ")), call. = FALSE)
  }
  if (length(opts[["panel"]]) != 3) {
    stop("Exactly three --panel LABEL=PATH arguments are required.", call. = FALSE)
  }
  if (length(opts[["window"]]) == 0) {
    stop("At least one --window LABEL=START-END argument is required.", call. = FALSE)
  }

  opts
}

parse_panel <- function(spec) {
  parts <- strsplit(spec, "=", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop(sprintf("Invalid --panel spec: %s", spec), call. = FALSE)
  list(label = parts[[1]], path = parts[[2]])
}

parse_window <- function(spec) {
  parts <- strsplit(spec, "=", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop(sprintf("Invalid --window spec: %s", spec), call. = FALSE)
  span <- strsplit(parts[[2]], "-", fixed = TRUE)[[1]]
  if (length(span) != 2) stop(sprintf("Invalid window span: %s", spec), call. = FALSE)
  list(label = parts[[1]], start = as.integer(span[[1]]), end = as.integer(span[[2]]))
}

load_summary_value <- function(path, target_key) {
  if (!file.exists(path)) return(NA_real_)
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  hit <- dt[["value"]][dt[["key"]] == target_key]
  if (length(hit) == 0) return(NA_real_)
  suppressWarnings(as.numeric(hit[[1]]))
}

load_panel_data <- function(spec) {
  panel <- parse_panel(spec)
  dt <- fread(panel$path, sep = "\t", header = TRUE, showProgress = TRUE)
  dt <- dt[is_common_variant == 1 & passes_plot_filter == 1]
  dt[, POS := suppressWarnings(as.integer(POS))]
  dt[, neg_log10_p := suppressWarnings(as.numeric(neg_log10_p))]
  dt <- dt[!is.na(POS) & !is.na(neg_log10_p) & is.finite(neg_log10_p)]
  dt[, panel := panel$label]
  summary_path <- sub("\\.normalized\\.tsv$", ".summary.tsv", panel$path)
  bonf <- load_summary_value(summary_path, "bonferroni_threshold")
  list(label = panel$label, data = dt, bonf = bonf)
}

draw_compare <- function(panels, windows, output_prefix, title, width, height, dpi, point_cex) {
  all_y <- unlist(lapply(panels, function(x) x$data$neg_log10_p))
  ymax <- max(8, ceiling(max(all_y, na.rm = TRUE) * 1.05))

  xmins <- vapply(panels, function(x) min(x$data$POS, na.rm = TRUE), numeric(1))
  xmaxs <- vapply(panels, function(x) max(x$data$POS, na.rm = TRUE), numeric(1))
  xmin <- min(xmins)
  xmax <- max(xmaxs)

  draw_once <- function(device_open) {
    device_open()
    on.exit(dev.off(), add = TRUE)
    par(mfrow = c(1, 3), mar = c(4.5, 4.8, 3.2, 1.2), oma = c(0, 0, 3, 0))

    for (panel in panels) {
      dt <- panel$data
      plot(
        dt$POS,
        dt$neg_log10_p,
        pch = 16,
        cex = point_cex,
        col = grDevices::rgb(78 / 255, 121 / 255, 167 / 255, alpha = 0.65),
        xlab = "chr1 position (bp)",
        ylab = expression(-log[10](italic(p))),
        main = panel$label,
        xlim = c(xmin, xmax),
        ylim = c(0, ymax),
        xaxt = "n"
      )
      axis(1, at = pretty(c(xmin, xmax), n = 5), labels = pretty(c(xmin, xmax), n = 5), cex.axis = 0.8)

      for (window in windows) {
        rect(
          xleft = window$start,
          ybottom = 0,
          xright = window$end,
          ytop = ymax,
          col = grDevices::rgb(214 / 255, 39 / 255, 40 / 255, alpha = 0.07),
          border = NA
        )
        text(
          x = (window$start + window$end) / 2,
          y = ymax * 0.97,
          labels = window$label,
          cex = 0.65,
          col = "#7f1d1d"
        )
      }

      points(
        dt$POS,
        dt$neg_log10_p,
        pch = 16,
        cex = point_cex,
        col = grDevices::rgb(78 / 255, 121 / 255, 167 / 255, alpha = 0.65)
      )

      if (is.finite(panel$bonf) && panel$bonf > 0) {
        abline(h = -log10(panel$bonf), col = "#D62728", lty = 2, lwd = 1.1)
      }
    }

    mtext(title, outer = TRUE, cex = 1.25, font = 2)
  }

  png_path <- sprintf("%s.png", output_prefix)
  pdf_path <- sprintf("%s.pdf", output_prefix)

  draw_once(function() {
    png(png_path, width = width, height = height, units = "in", res = dpi, bg = "white")
  })
  draw_once(function() {
    cairo_pdf(pdf_path, width = width, height = height, bg = "white", fallback_resolution = dpi)
  })
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  panels <- lapply(opts[["panel"]], load_panel_data)
  windows <- lapply(opts[["window"]], parse_window)
  output_prefix <- opts[["output-prefix"]]
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  draw_compare(
    panels = panels,
    windows = windows,
    output_prefix = output_prefix,
    title = opts[["title"]],
    width = as.numeric(opts[["width"]]),
    height = as.numeric(opts[["height"]]),
    dpi = as.integer(opts[["dpi"]]),
    point_cex = as.numeric(opts[["point_cex"]])
  )
}

main()
