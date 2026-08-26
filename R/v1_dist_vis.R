#' @keywords internal
vis_dist_impl <- function(
  .data,
  mode = c("nearest", "all"),
  xval = c("norm_dist", "dist", "sim"),
  facet = "<auto>",
  binwidth = NULL,
  title = NULL,
  dir = c("h", "v"),
  ...
) {
  checkmate::assert_data_frame(.data)
  mode <- match.arg(mode)
  xval <- match.arg(xval)
  dir <- match.arg(dir)
  checkmate::assert_number(
    binwidth,
    lower = 0,
    finite = TRUE,
    null.ok = TRUE
  )
  if (!is.null(binwidth) && binwidth == 0) {
    cli::cli_abort("{.arg binwidth} must be greater than 0.")
  }
  checkmate::assert_string(title, null.ok = TRUE)
  checkmate::assert_subset(
    c("imd_receptor_id_1", "imd_receptor_id_2", "dist", xval),
    names(.data)
  )

  if (!is.null(facet)) {
    checkmate::assert_string(facet)
  }
  if (identical(facet, "<auto>")) {
    facet <- infer_dist_group_col(names(.data))
  } else if (!is.null(facet)) {
    checkmate::assert_choice(facet, names(.data))
  }

  if (mode == "nearest") {
    edge_group_cols <- setdiff(
      names(.data),
      c(
        "imd_receptor_id_1",
        "imd_receptor_id_2",
        "seq_len",
        "dist",
        "norm_dist",
        "sim"
      )
    )
    group_cols <- c(edge_group_cols, ".idist_receptor")
    endpoint_cols <- c(
      edge_group_cols,
      .idist_receptor = "imd_receptor_id_1",
      .idist_value = xval,
      .idist_raw_dist = "dist"
    )
    endpoint_1 <- .data |>
      dplyr::select(dplyr::all_of(endpoint_cols))
    endpoint_cols[".idist_receptor"] <- "imd_receptor_id_2"
    endpoint_2 <- .data |>
      dplyr::select(dplyr::all_of(endpoint_cols))

    nearest_data <- dplyr::union_all(endpoint_1, endpoint_2) |>
      dplyr::filter(.data$.idist_raw_dist > 0)
    plot_data <- if (xval == "sim") {
      nearest_data |>
        dplyr::summarise(
          .idist_value = max(.data$.idist_value, na.rm = TRUE),
          .by = dplyr::all_of(group_cols)
        )
    } else {
      nearest_data |>
        dplyr::summarise(
          .idist_value = min(.data$.idist_value, na.rm = TRUE),
          .by = dplyr::all_of(group_cols)
        )
    }
    plot_data <- dplyr::collect(plot_data)
  } else {
    plot_cols <- c(facet, .idist_value = xval)
    plot_data <- .data |>
      dplyr::select(dplyr::all_of(plot_cols)) |>
      dplyr::collect()
  }

  if (is.null(binwidth)) {
    binwidth <- if (xval == "dist") 1 else 0.02
  }
  boundary <- if (xval == "dist") -0.5 else 0
  if (is.null(title)) {
    title <- if (mode == "nearest") {
      "Nearest-neighbor distance distribution"
    } else {
      "Pairwise distance distribution"
    }
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$.idist_value)
  ) +
    ggplot2::geom_histogram(
      binwidth = binwidth,
      boundary = boundary,
      colour = "white",
      na.rm = TRUE
    ) +
    ggplot2::labs(
      title = title,
      x = xval,
      y = "Count"
    )

  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste0("~", facet)),
      dir = dir,
      scales = "free_y"
    )
  }

  if (requireNamespace("ggthemes", quietly = TRUE)) {
    p <- p + ggthemes::theme_few()
  }

  p
}


register_immunarch_visualisation(
  vis_dist_impl,
  family = "dist",
  name = NULL
)
