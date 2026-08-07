#' @keywords internal
add_clonality_facets <- function(p, facet, dir) {
  if (is.null(facet)) {
    return(p)
  }

  if (length(facet) == 1) {
    p + ggplot2::facet_wrap(
      stats::as.formula(paste0("~", facet)),
      dir = dir,
      scales = "free_y"
    )
  } else {
    p + ggplot2::facet_grid(
      stats::as.formula(paste(facet[1], "~", facet[2])),
      scales = "free_y"
    )
  }
}


#' @keywords internal
add_clonality_colour_scale <- function(p, values, aesthetic, name) {
  k <- length(unique(values))

  if (requireNamespace("ggsci", quietly = TRUE) && k <= 11) {
    if (aesthetic == "colour") {
      p + ggsci::scale_color_locuszoom(name = name)
    } else {
      p + ggsci::scale_fill_locuszoom(name = name)
    }
  } else if (aesthetic == "colour") {
    p + ggplot2::scale_colour_viridis_d(option = "H", name = name)
  } else {
    p + ggplot2::scale_fill_viridis_d(option = "H", name = name)
  }
}


#' @keywords internal
vis_airr_clonality_line_impl <- function(
    .data,
    yval = immundata::imd_schema("proportion"),
    color = immundata::imd_schema("repertoire"),
    facet = NULL,
    title = "Clonal rank-abundance",
    log = TRUE,
    dir = c("h", "v"),
    ...) {
  checkmate::assert_data_frame(.data)
  checkmate::assert_string(yval)
  checkmate::assert_choice(
    yval,
    c(
      immundata::imd_schema("count"),
      immundata::imd_schema("proportion")
    )
  )
  checkmate::assert_string(color)
  checkmate::assert_string(title)
  checkmate::assert_logical(log, len = 1, any.missing = FALSE)
  checkmate::assert_subset(
    c("index", yval, color, immundata::imd_schema("repertoire")),
    names(.data)
  )

  dir <- match.arg(dir)

  if (!is.null(facet)) {
    checkmate::assert_character(
      facet,
      any.missing = FALSE,
      min.len = 1,
      max.len = 2
    )
    checkmate::assert_subset(facet, names(.data))
  }

  plot_data <- .data
  plot_data[[color]] <- as.factor(plot_data[[color]])
  repertoire_col <- immundata::imd_schema("repertoire")

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$index,
      y = .data[[yval]],
      colour = .data[[color]],
      group = .data[[repertoire_col]]
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::labs(
      title = title,
      x = "Receptor rank",
      y = if (yval == immundata::imd_schema("proportion")) {
        "Receptor proportion"
      } else {
        "Receptor count"
      }
    )

  if (isTRUE(log)) {
    p <- p + ggplot2::scale_y_log10(
      labels = if (yval == immundata::imd_schema("proportion")) {
        scales::label_percent()
      } else {
        scales::label_number()
      }
    )
  } else if (yval == immundata::imd_schema("proportion")) {
    p <- p + ggplot2::scale_y_continuous(labels = scales::label_percent())
  }

  p <- add_clonality_facets(p, facet, dir)
  p <- add_clonality_colour_scale(p, plot_data[[color]], "colour", color)

  if (requireNamespace("ggthemes", quietly = TRUE)) {
    p <- p + ggthemes::theme_few()
  }

  p
}

register_immunarch_visualisation(
  vis_airr_clonality_line_impl,
  "airr_clonality",
  "line"
)


#' @keywords internal
make_clonality_space_plot <- function(bin_col, title_default, fill_label) {
  checkmate::assert_string(bin_col)
  checkmate::assert_string(title_default)
  checkmate::assert_string(fill_label)

  space_plot <- function(
      .data,
      xval = immundata::imd_schema("repertoire"),
      yval = "occupied_prop",
      fill = bin_col,
      facet = NULL,
      title = title_default,
      dir = c("h", "v"),
      ...) {
    checkmate::assert_data_frame(.data)
    checkmate::assert_string(xval)
    checkmate::assert_string(yval)
    checkmate::assert_string(fill)
    checkmate::assert_string(title)
    checkmate::assert_subset(c(xval, yval, fill), names(.data))

    dir <- match.arg(dir)

    if (!is.null(facet)) {
      checkmate::assert_character(
        facet,
        any.missing = FALSE,
        min.len = 1,
        max.len = 2
      )
      checkmate::assert_subset(facet, names(.data))
    }

    plot_data <- .data
    plot_data <- plot_data[!is.na(plot_data[[fill]]), , drop = FALSE]
    plot_data[[xval]] <- as.factor(plot_data[[xval]])
    if (!is.factor(plot_data[[fill]])) {
      fill_values <- unique(plot_data[[fill]])
      if (is.numeric(fill_values)) {
        fill_values <- sort(fill_values)
      }
      plot_data[[fill]] <- factor(plot_data[[fill]], levels = fill_values)
    }

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data[[xval]],
        y = .data[[yval]],
        fill = .data[[fill]]
      )
    ) +
      ggplot2::geom_col(
        position = "stack",
        colour = "grey30",
        na.rm = TRUE
      ) +
      ggplot2::labs(
        title = title,
        x = "Repertoire",
        y = "Occupied repertoire space",
        fill = if (fill == bin_col) fill_label else fill
      ) +
      ggplot2::scale_y_continuous(labels = scales::label_percent())

    p <- add_clonality_facets(p, facet, dir)
    p <- add_clonality_colour_scale(
      p,
      plot_data[[fill]],
      "fill",
      if (fill == bin_col) fill_label else fill
    )

    if (requireNamespace("ggthemes", quietly = TRUE)) {
      p <- p + ggthemes::theme_few()
    }

    p
  }

  space_plot
}


#' @keywords internal
vis_airr_clonality_rank_impl <- make_clonality_space_plot(
  bin_col = "clonal_rank_bin",
  title_default = "Clonal space by receptor rank",
  fill_label = "Rank bin"
)

register_immunarch_visualisation(
  vis_airr_clonality_rank_impl,
  "airr_clonality",
  "rank"
)


#' @keywords internal
vis_airr_clonality_prop_impl <- make_clonality_space_plot(
  bin_col = "clonal_prop_bin",
  title_default = "Clonal space by receptor proportion",
  fill_label = "Proportion bin"
)

register_immunarch_visualisation(
  vis_airr_clonality_prop_impl,
  "airr_clonality",
  "prop"
)
