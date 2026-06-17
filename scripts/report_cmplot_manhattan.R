#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CMplot)
  library(data.table)
})

parse_args <- function(args) {
  opts <- list(
    title = "NCN common-variant standardized SDS Manhattan plot",
    width = 16,
    height = 7,
    dpi = 400,
    "point-cex" = 0.28,
    "maf-threshold" = 0.01,
    "plot-p-threshold" = 1
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected argument: %s", arg), call. = FALSE)
    }
    key <- substring(arg, 3)
    if (i == length(args)) {
      stop(sprintf("Missing value for argument: %s", arg), call. = FALSE)
    }
    opts[[key]] <- args[[i + 1]]
    i <- i + 2
  }

  required <- c("input-normalized-tsv", "output-prefix")
  missing <- required[!vapply(required, function(x) !is.null(opts[[x]]), logical(1))]
  if (length(missing) > 0) {
    stop(
      sprintf(
        "Missing required arguments: %s",
        paste(sprintf("--%s", missing), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  opts
}

load_bonferroni_threshold <- function(summary_path, fallback_count) {
  if (!is.null(summary_path) && file.exists(summary_path)) {
    summary_dt <- fread(summary_path, sep = "\t", header = TRUE, showProgress = FALSE)
    hit <- summary_dt[key == "bonferroni_threshold", value]
    if (length(hit) > 0) {
      threshold <- suppressWarnings(as.numeric(hit[[1]]))
      if (!is.na(threshold) && is.finite(threshold) && threshold > 0) {
        return(threshold)
      }
    }
  }

  if (fallback_count <= 0) {
    stop("Unable to determine Bonferroni threshold because there are no plottable variants.", call. = FALSE)
  }
  0.05 / fallback_count
}

load_plot_data <- function(normalized_path, maf_threshold, plot_p_threshold) {
  dt <- fread(
    normalized_path,
    sep = "\t",
    header = TRUE,
    select = c("ID", "chr", "pos", "DAF", "MAF", "p_bothside", "is_common_variant"),
    showProgress = TRUE
  )

  # 统计过滤前的位点数
  total_variants_before_filter <- nrow(dt)

  dt <- dt[is_common_variant == 1]
  dt[, chr := suppressWarnings(as.integer(chr))]
  dt[, pos := suppressWarnings(as.integer(pos))]
  dt[, DAF := suppressWarnings(as.numeric(DAF))]
  dt[, MAF := suppressWarnings(as.numeric(MAF))]
  dt[, p_bothside := suppressWarnings(as.numeric(p_bothside))]
  dt <- dt[
    !is.na(chr) &
      chr >= 1L &
      chr <= 22L &
      !is.na(pos) &
      pos > 0L &
      !is.na(DAF) &
      is.finite(DAF) &
      !is.na(p_bothside) &
      is.finite(p_bothside) &
      p_bothside > 0 &
      p_bothside <= 1 &
      !is.na(MAF) &
      is.finite(MAF)
  ]

  # 统计基础过滤后的位点数
  total_variants_after_basic_qc <- nrow(dt)

  # 联合频率过滤: common-variant 掩码 + DAF 在 (maf_threshold, 1-maf_threshold) + MAF > maf_threshold
  dt <- dt[DAF > maf_threshold & DAF < (1 - maf_threshold) & MAF > maf_threshold]
  total_variants_after_frequency_filter <- nrow(dt)

  # p值过滤
  dt <- dt[p_bothside < plot_p_threshold]
  total_variants_after_p_filter <- nrow(dt)

  if (nrow(dt) == 0) {
    stop("No variants passed DAF/MAF and p-value filters.", call. = FALSE)
  }

  dt[p_bothside < 1e-300, p_bothside := 1e-300]
  setorder(dt, chr, pos)

  list(
    data = data.frame(
      SNP = dt$ID,
      CHR = dt$chr,
      BP = dt$pos,
      P = dt$p_bothside,
      stringsAsFactors = FALSE
    ),
    stats = list(
      before_filter = total_variants_before_filter,
      after_basic_qc = total_variants_after_basic_qc,
      after_frequency_filter = total_variants_after_frequency_filter,
      after_p_filter = total_variants_after_p_filter,
      final_plotted = nrow(dt)
    )
  )
}

plot_cmplot_manhattan <- function(plot_df, threshold, output_prefix, title, width, height, dpi, point_cex) {
  ymax <- max(
    ceiling(max(-log10(plot_df$P), na.rm = TRUE) * 1.05),
    ceiling(-log10(threshold) * 1.05),
    8
  )

  colors <- c("#4E79A7", "#9C755F")
  chr_labels <- c("chr1", as.character(2:22))

  render_once <- function(device_open) {
    device_open()
    on.exit(dev.off(), add = TRUE)
    CMplot(
      plot_df,
      plot.type = "m",
      LOG10 = TRUE,
      file.output = FALSE,
      col = colors,
      cex = point_cex,
      pch = 19,
      points.alpha = 80L,
      threshold = threshold,
      threshold.col = "#D62728",
      threshold.lwd = 1.1,
      threshold.lty = 2,
      amplify = FALSE,
      chr.labels = chr_labels,
      chr.den.col = NULL,
      axis.cex = 0.8,
      axis.lwd = 1.1,
      lab.cex = 1.0,
      lab.font = 2,
      main = title,
      main.cex = 1.2,
      main.font = 2,
      ylab = expression(-log[10](italic(p))),
      band = 0.6,
      ylim = c(0, ymax),
      mar = c(4.5, 5.5, 3.5, 1.5),
      verbose = TRUE
    )
  }

  png_path <- sprintf("%s.png", output_prefix)

  render_once(function() {
    png(png_path, width = width, height = height, units = "in", res = dpi, bg = "white")
  })

  invisible(list(png = png_path))
}

main <- function() {
  Sys.setenv(XDG_CACHE_HOME = "/tmp", HOME = Sys.getenv("HOME", unset = "/tmp"))

  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  normalized_path <- normalizePath(opts[["input-normalized-tsv"]], mustWork = TRUE)
  output_prefix <- opts[["output-prefix"]]
  summary_path <- opts[["input-summary-tsv"]]
  bonferroni_arg <- opts[["bonferroni-threshold"]]
  title <- opts[["title"]]
  width <- as.numeric(opts[["width"]])
  height <- as.numeric(opts[["height"]])
  dpi <- as.integer(opts[["dpi"]])
  point_cex <- as.numeric(opts[["point-cex"]])
  maf_threshold <- as.numeric(opts[["maf-threshold"]])
  plot_p_threshold <- as.numeric(opts[["plot-p-threshold"]])

  if (is.na(width) || width <= 0 || is.na(height) || height <= 0) {
    stop("Plot width and height must be positive numbers.", call. = FALSE)
  }
  if (is.na(dpi) || dpi <= 0) {
    stop("DPI must be a positive integer.", call. = FALSE)
  }
  if (is.na(point_cex) || point_cex <= 0) {
    stop("point-cex must be a positive number.", call. = FALSE)
  }
  if (is.na(maf_threshold) || maf_threshold <= 0 || maf_threshold >= 0.5) {
    stop("maf-threshold must be between 0 and 0.5.", call. = FALSE)
  }
  if (is.na(plot_p_threshold) || plot_p_threshold <= 0 || plot_p_threshold > 1) {
    stop("plot-p-threshold must be between 0 and 1.", call. = FALSE)
  }

  out_dir <- dirname(output_prefix)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # 加载数据并获取统计信息
  result <- load_plot_data(normalized_path, maf_threshold, plot_p_threshold)
  plot_df <- result$data
  stats <- result$stats

  # 输出统计信息
  cat("=== Filtering Statistics ===\n")
  cat(sprintf("Total variants in input file:\t%d\n", stats$before_filter))
  cat(sprintf("After basic QC (chr 1-22, valid p-value):\t%d\n", stats$after_basic_qc))
  cat(sprintf("After DAF/MAF filter (DAF > %.4g, DAF < %.4g, MAF > %.4g):\t%d\n", maf_threshold, 1 - maf_threshold, maf_threshold, stats$after_frequency_filter))
  cat(sprintf("After p-value filter (<%.4g):\t%d\n", plot_p_threshold, stats$after_p_filter))
  cat(sprintf("Final variants to plot:\t%d\n", stats$final_plotted))
  cat(sprintf("Frequency filtering removed:\t%d variants\n", stats$after_basic_qc - stats$after_frequency_filter))
  cat(sprintf("P-value filtering removed:\t%d variants\n", stats$after_frequency_filter - stats$after_p_filter))
  cat("============================\n")

  bonferroni_threshold <- NA_real_
  if (!is.null(bonferroni_arg)) {
    parsed <- suppressWarnings(as.numeric(bonferroni_arg))
    if (length(parsed) > 0) {
      bonferroni_threshold <- parsed[[1]]
    }
  }
  if (!is.finite(bonferroni_threshold) || bonferroni_threshold <= 0) {
    bonferroni_threshold <- load_bonferroni_threshold(summary_path, nrow(plot_df))
  }

  outputs <- plot_cmplot_manhattan(
    plot_df = plot_df,
    threshold = bonferroni_threshold,
    output_prefix = output_prefix,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    point_cex = point_cex
  )

  cat(sprintf("input\t%s\n", normalized_path))
  cat(sprintf("variants_before_filter\t%d\n", stats$before_filter))
  cat(sprintf("variants_after_frequency_filter\t%d\n", stats$after_frequency_filter))
  cat(sprintf("variants_after_p_filter\t%d\n", stats$after_p_filter))
  cat(sprintf("variants_plotted\t%d\n", stats$final_plotted))
  cat(sprintf("maf_threshold\t%.10g\n", maf_threshold))
  cat(sprintf("plot_p_threshold\t%.10g\n", plot_p_threshold))
  cat(sprintf("bonferroni_threshold\t%.10g\n", bonferroni_threshold))
  cat(sprintf("png\t%s\n", outputs$png))
}

main()
