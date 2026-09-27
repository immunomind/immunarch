#' @title
#' Visualise public receptors between two repertoires
#'
#' @description
#' Visualise receptor frequencies shared between two repertoires.
#'
#' Receptor proportions are shown on log-scaled axes. Point fill shows the
#' log proportion ratio between the x and y repertoires.
#'
#' @param .data An [immundata::ImmunData] object with aggregated repertoires.
#' @param x Repertoire for the x-axis. Either a numeric repertoire index,
#'   repertoire ID, or a repertoire schema value.
#' @param y Repertoire for the y-axis. Either a numeric repertoire index,
#'   repertoire ID, or a repertoire schema value.
#' @param density Logical. If `TRUE`, show marginal distributions of unique
#'   and public receptor proportions.
#' @param lm Logical. If `TRUE`, add a linear model fit and adjusted R-squared.
#' @param ... Must be empty.
#'
#' @return A ggplot2 or patchwork object.
#'
#' @examples
#' \dontrun{
#' idata |> vis()
#' idata |> vis(x = 1, y = 2)
#' idata |> vis(x = "Blood", y = "Tumor")
#' idata |> vis(x = 1, y = 2, density = FALSE)
#' }
#'
#' @method vis ImmunData
#' @export
vis.ImmunData <- function(
  .data,
  x = 1,
  y = 2,
  density = TRUE,
  lm = TRUE,
  ...
) {
  idata <- .data
  checkmate::assert_r6(idata, "ImmunData")
  checkmate::assert_flag(density)
  checkmate::assert_flag(lm)
  rlang::check_dots_empty()

  repertoires <- idata$repertoires
  if (is.null(repertoires) || is.null(idata$schema_repertoire)) {
    cli::cli_abort(
      "Repertoire aggregation is required. Run {.code agg_repertoires()} first."
    )
  }

  repertoire_col <- immundata::imd_schema("repertoire")
  receptor_col <- immundata::imd_schema("receptor")
  prop_col <- immundata::imd_schema("proportion")

  required_cols <- c(
    repertoire_col,
    receptor_col,
    prop_col
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(idata$annotations)
  )

  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      "Required annotation column(s) {missing_cols} are missing from {.field idata$annotations}."
    )
  }

  # -------------------------------------------------------------------------
  # Repertoires
  # -------------------------------------------------------------------------

  if (nrow(repertoires) < 2L) {
    cli::cli_abort(
      "At least two repertoires are required for this visualisation."
    )
  }

  missing_schema <- setdiff(
    idata$schema_repertoire,
    names(repertoires)
  )

  if (length(missing_schema) > 0L) {
    cli::cli_abort(
      "Repertoire schema column(s) {missing_schema} are missing from {.field idata$repertoires}."
    )
  }

  resolve_repertoire <- function(value, arg) {
    if (is.numeric(value)) {
      checkmate::assert_int(
        value,
        lower = 1L,
        upper = nrow(repertoires)
      )

      hits <- as.integer(value)
    } else {
      checkmate::assert_string(value)

      hits <- which(
        as.character(repertoires[[repertoire_col]]) == value
      )

      for (schema_col in idata$schema_repertoire) {
        hits <- union(
          hits,
          which(as.character(repertoires[[schema_col]]) == value)
        )
      }

      if (length(hits) == 0L) {
        cli::cli_abort(
          "Could not find repertoire {.val {value}} supplied to {.arg {arg}}."
        )
      }

      if (length(hits) > 1L) {
        cli::cli_abort(
          "Repertoire {.val {value}} supplied to {.arg {arg}} is ambiguous."
        )
      }
    }

    schema_values <- vapply(
      idata$schema_repertoire,
      function(column) {
        value <- as.character(repertoires[[column]][hits])
        if (column == immundata::imd_schema("manifest_filename")) {
          value <- basename(value)
        }
        value
      },
      character(1)
    )

    label <- if (length(schema_values) == 1L) {
      unname(schema_values)
    } else {
      paste(names(schema_values), schema_values, sep = "=", collapse = ", ")
    }

    list(
      id = repertoires[[repertoire_col]][hits],
      label = label
    )
  }

  x_rep <- resolve_repertoire(x, "x")
  y_rep <- resolve_repertoire(y, "y")

  if (identical(x_rep$id, y_rep$id)) {
    cli::cli_abort(
      "{.arg x} and {.arg y} must refer to different repertoires."
    )
  }

  # -------------------------------------------------------------------------
  # Aggregate chain rows in duckplyr before collecting receptor proportions.
  # -------------------------------------------------------------------------

  plot_data <- idata$annotations |>
    dplyr::filter(
      !!rlang::sym(repertoire_col) %in% c(x_rep$id, y_rep$id),
      !!rlang::sym(prop_col) > 0
    ) |>
    dplyr::summarise(
      proportion = max(!!rlang::sym(prop_col), na.rm = TRUE),
      .by = dplyr::all_of(c(receptor_col, repertoire_col))
    ) |>
    dplyr::collect() |>
    dplyr::rename(
      repertoire = dplyr::all_of(repertoire_col),
      receptor = dplyr::all_of(receptor_col)
    )

  # imd_proportion is already calculated per repertoire.
  plot_data <- plot_data |>
    dplyr::filter(
      is.finite(proportion)
    ) |>
    dplyr::mutate(
      axis = ifelse(
        repertoire == x_rep$id,
        "x",
        "y"
      )
    ) |>
    dplyr::select(
      receptor,
      axis,
      proportion
    ) |>
    tidyr::pivot_wider(
      names_from = axis,
      values_from = proportion
    )

  if (!"x" %in% names(plot_data)) {
    cli::cli_abort(
      "Repertoire {.val {x_rep$label}} contains no receptor proportions."
    )
  }

  if (!"y" %in% names(plot_data)) {
    cli::cli_abort(
      "Repertoire {.val {y_rep$label}} contains no receptor proportions."
    )
  }

  # -------------------------------------------------------------------------
  # Public receptors
  # -------------------------------------------------------------------------

  plot_public <- plot_data |>
    dplyr::filter(
      !is.na(x),
      !is.na(y)
    )

  if (nrow(plot_public) == 0L) {
    cli::cli_abort(
      "Repertoires {.val {x_rep$label}} and {.val {y_rep$label}} have no public receptors."
    )
  }

  plot_public <- plot_public |>
    dplyr::mutate(
      enrichment = log(x / y)
    )

  # Use identical limits for the central and marginal plots.
  limits <- range(
    c(plot_data$x, plot_data$y),
    finite = TRUE
  )

  if (limits[1] == limits[2]) {
    limits <- limits * c(0.9, 1.1)
  }

  # -------------------------------------------------------------------------
  # Main scatter
  # -------------------------------------------------------------------------

  points <- ggplot2::ggplot(
    plot_public,
    ggplot2::aes(
      x = x,
      y = y
    )
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = enrichment
      ),
      shape = 21,
      size = 3,
      alpha = 0.8
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      name = "Log ratio (x / y)"
    ) +
    ggplot2::scale_x_log10(
      limits = limits
    ) +
    ggplot2::scale_y_log10(
      limits = limits
    ) +
    ggplot2::labs(
      x = x_rep$label,
      y = y_rep$label
    ) +
    ggthemes::theme_few() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )

  # No coord_equal() here:
  # a fixed aspect ratio prevents the central panel from filling its
  # patchwork cell and causes the marginal density plots to misalign.

  # -------------------------------------------------------------------------
  # Linear model
  # -------------------------------------------------------------------------

  if (
    lm &&
      nrow(plot_public) >= 3L &&
      length(unique(plot_public$x)) >= 2L &&
      length(unique(plot_public$y)) >= 2L
  ) {
    model <- stats::lm(
      log10(y) ~ log10(x),
      data = plot_public
    )

    adj_r2 <- summary(model)$adj.r.squared

    points <- points +
      ggplot2::geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        fullrange = TRUE,
        linewidth = 0.5
      )

    if (is.finite(adj_r2)) {
      points <- points +
        ggplot2::annotate(
          "text",
          x = limits[1],
          y = limits[2],
          label = sprintf(
            "R\u00b2(adj.) = %.2f",
            adj_r2
          ),
          hjust = 0,
          vjust = 1
        )
    }
  }

  if (!density) {
    return(points)
  }

  # -------------------------------------------------------------------------
  # Marginal densities
  # -------------------------------------------------------------------------

  x_density <- plot_data |>
    dplyr::filter(
      !is.na(x),
      x > 0
    ) |>
    dplyr::mutate(
      type = ifelse(
        is.na(y),
        "unique",
        "public"
      )
    )

  y_density <- plot_data |>
    dplyr::filter(
      !is.na(y),
      y > 0
    ) |>
    dplyr::mutate(
      type = ifelse(
        is.na(x),
        "unique",
        "public"
      )
    )

  top_plot <- ggplot2::ggplot(
    x_density,
    ggplot2::aes(
      x = x,
      fill = type
    )
  ) +
    ggplot2::geom_density(
      alpha = 0.3,
      colour = "grey25",
      na.rm = TRUE
    ) +
    ggplot2::scale_x_log10(
      limits = limits
    ) +
    ggplot2::scale_fill_brewer(
      palette = "Set2",
      name = "Receptors"
    ) +
    ggplot2::theme_void()

  right_plot <- ggplot2::ggplot(
    y_density,
    ggplot2::aes(
      x = y,
      fill = type
    )
  ) +
    ggplot2::geom_density(
      alpha = 0.3,
      colour = "grey25",
      na.rm = TRUE
    ) +
    ggplot2::scale_x_log10(
      limits = limits
    ) +
    ggplot2::scale_fill_brewer(
      palette = "Set2",
      guide = "none"
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_void()

  patchwork::wrap_plots(
    top_plot,
    patchwork::plot_spacer(),
    points,
    right_plot,
    ncol = 2,
    widths = c(4, 1),
    heights = c(1, 4),
    guides = "collect"
  ) &
    ggplot2::theme(
      legend.position = "right",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.just = "left",
      legend.title = ggplot2::element_text(hjust = 0)
    )
}
