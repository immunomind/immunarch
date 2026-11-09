#' @keywords internal
public_fill_scale <- function(n, name) {
  labels <- function(x) stringr::str_wrap(x, width = 30)
  if (n <= 7L) {
    ggsci::scale_fill_locuszoom(name = name, labels = labels)
  } else {
    ggplot2::scale_fill_viridis_d(option = "H", name = name, labels = labels)
  }
}


#' @keywords internal
vis_annotate_public_impl <- function(
    .data,
    title = NULL) {
  checkmate::assert_r6(.data, "ImmunData")
  checkmate::assert_string(title, null.ok = TRUE)
  idata <- .data
  annotations <- idata$annotations
  if (!"imd_public_incidence" %in% names(annotations)) {
    cli::cli_abort("Run {.fn annotate_public} before plotting receptor incidence.")
  }
  stratified <- !is.null(idata$strata)
  if (stratified) {
    strata <- dplyr::collect(idata$strata)
    columns <- paste0("imd_public_incidence_strata_", strata$imd_strata_id)
    labels <- strata$strata_name
  } else {
    columns <- "imd_public_incidence"
    labels <- "All repertoires"
  }
  checkmate::assert_subset(c("imd_receptor_id", columns), names(annotations))
  # Receptor-level statistics repeat across cell and chain rows. Count each
  # identity once per stratum, collecting only the small incidence histograms.
  histograms <- lapply(seq_along(columns), function(i) {
    annotations |>
      dplyr::select(imd_receptor_id, incidence = dplyr::all_of(columns[i])) |>
      dplyr::distinct() |>
      dplyr::filter(!is.na(incidence), incidence > 0) |>
      dplyr::count(incidence, name = "receptors") |>
      dplyr::collect() |>
      tibble::as_tibble() |>
      dplyr::mutate(stratum = labels[i])
  })
  plot_data <- dplyr::bind_rows(histograms)
  if (!nrow(plot_data)) cli::cli_abort("No receptor incidence observations to plot.")
  plot_data$stratum <- factor(plot_data$stratum, levels = sort(unique(labels)))
  totals <- plot_data |>
    dplyr::summarise(receptors = sum(receptors), .by = incidence)
  maximum <- max(plot_data$incidence)
  p <- ggplot2::ggplot(plot_data,
    ggplot2::aes(x = incidence, y = receptors, fill = stratum)) +
    ggplot2::geom_col(width = .7, colour = "grey30", linewidth = .3, alpha = .6) +
    ggplot2::geom_text(data = totals, inherit.aes = FALSE,
      ggplot2::aes(x = incidence, y = receptors, label = receptors),
      vjust = -.6, size = 3.2, colour = "grey30") +
    ggplot2::scale_x_continuous(
      breaks = if (maximum <= 20L) seq_len(maximum) else scales::breaks_pretty(),
      expand = ggplot2::expansion(add = .6)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .1))) +
    public_fill_scale(nlevels(plot_data$stratum), "Stratum") +
    ggplot2::guides(fill = if (stratified) ggplot2::guide_legend(ncol = 1) else "none") +
    ggplot2::labs(
      title = title,
      x = if (stratified) "Number of repertoires within stratum" else "Number of repertoires",
      y = "Receptor counts") +
    ggthemes::theme_few() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold"))
  p
}


#' @keywords internal
summarise_public_selection <- function(observations, reference) {
  checkmate::assert_r6(reference, "ImmunData")
  if (is.null(reference$schema_repertoire)) {
    cli::cli_abort("The reference must have defined repertoires.")
  }
  schema <- reference$schema_repertoire
  columns <- unique(c("imd_receptor_id", schema, "imd_orig_proportion"))
  checkmate::assert_subset(columns, names(observations))
  repertoires <- dplyr::collect(reference$repertoires)
  # Match repertoire metadata, not IDs rebuilt by ImmunData filtering.
  # Receptor IDs are used only within the selected result, never across objects.
  occurrences <- duckplyr::as_duckdb_tibble(observations) |>
    dplyr::select(dplyr::all_of(columns)) |>
    dplyr::distinct()
  n_selected <- occurrences |>
    dplyr::distinct(imd_receptor_id) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::collect() |>
    dplyr::pull(n)
  if (n_selected == 0L) cli::cli_abort("No selected receptors to plot.")
  population <- duckplyr::as_duckdb_tibble(repertoires)
  abundances <- occurrences |>
    dplyr::inner_join(
      population |> dplyr::select(dplyr::all_of(unique(c(schema, "imd_repertoire_id")))),
      by = schema
    ) |>
    dplyr::summarise(
      n_detected = dplyr::n(), abundance = sum(imd_orig_proportion),
      .by = imd_repertoire_id
    )
  population |>
    dplyr::left_join(abundances, by = "imd_repertoire_id") |>
    dplyr::mutate(
      n_detected = dplyr::coalesce(n_detected, 0L),
      abundance = dplyr::coalesce(abundance, 0),
      n_selected = !!n_selected,
      coverage = n_detected / !!n_selected
    ) |>
    dplyr::collect() |>
    tibble::as_tibble()
}


#' @keywords internal
vis_vdj_public_impl <- function(
    .data,
    reference = NULL,
    xval = "imd_repertoire_id",
    fill = NULL,
    title = "Selected immune receptors across repertoires") {
  checkmate::assert_data_frame(.data)
  if (is.null(reference)) {
    cli::cli_abort("Supply {.arg reference}: the original ImmunData passed to {.fn vdj_public}.")
  }
  plot_data <- summarise_public_selection(.data, reference)
  checkmate::assert_string(xval)
  checkmate::assert_string(fill, null.ok = TRUE)
  checkmate::assert_string(title, null.ok = TRUE)
  if (is.null(fill) && "strata_name" %in% names(plot_data)) fill <- "strata_name"
  checkmate::assert_subset(c(xval, fill), names(plot_data))
  if (anyNA(plot_data[c(xval, fill)])) cli::cli_abort("Plotting metadata must not contain missing values.")
  if (anyDuplicated(plot_data[[xval]])) cli::cli_abort("Each x value must identify one repertoire.")
  order_columns <- c(fill, xval)
  plot_data <- plot_data[do.call(order, plot_data[order_columns]), ]
  plot_data[[xval]] <- factor(plot_data[[xval]], levels = plot_data[[xval]])
  plot_data$public_group <- factor(if (is.null(fill)) "Selected receptors" else plot_data[[fill]])
  plot_data$coverage_label <- paste0(plot_data$n_detected, "/", plot_data$n_selected)
  plot_data$abundance_label <- scales::percent(plot_data$abundance, accuracy = .01)
  close_marks <- abs(plot_data$coverage - plot_data$abundance) < .06
  plot_data$abundance_vjust <- ifelse(close_marks & plot_data$abundance > .06, 1.5, -.6)
  plot_data$coverage_vjust <- ifelse(close_marks & plot_data$abundance <= .06, -2.4, -1)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = !!rlang::sym(xval), fill = public_group)) +
    ggplot2::geom_segment(ggplot2::aes(xend = !!rlang::sym(xval), y = 0, yend = coverage),
      colour = "grey80", linewidth = .5, show.legend = FALSE) +
    ggplot2::geom_col(ggplot2::aes(y = abundance), width = .65,
      colour = "grey30", linewidth = .3, alpha = .6) +
    ggplot2::geom_text(ggplot2::aes(y = abundance, label = abundance_label, vjust = abundance_vjust),
      size = 3.2, colour = "grey30", show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(y = coverage, shape = "Coverage"),
      colour = "grey30", size = 2.8, stroke = .4, show.legend = TRUE) +
    ggplot2::geom_text(ggplot2::aes(y = coverage, label = coverage_label, vjust = coverage_vjust),
      size = 3.2, colour = "grey30", show.legend = FALSE) +
    ggplot2::scale_shape_manual(name = "Measure",
      values = c("Abundance" = 22, "Coverage" = 21),
      limits = c("Abundance", "Coverage"), drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2),
      labels = scales::label_percent(), expand = ggplot2::expansion(mult = c(0, .08))) +
    ggplot2::labs(title = title, x = xval, y = "Percentage") +
    ggplot2::guides(
      fill = if (is.null(fill)) "none" else ggplot2::guide_legend(order = 1, ncol = 1,
        override.aes = list(shape = NA, alpha = 1)),
      shape = ggplot2::guide_legend(order = 2, ncol = 1,
        override.aes = list(fill = "grey70", colour = "grey30", size = 2.8))) +
    ggthemes::theme_few() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold"))
  fill_name <- if (identical(fill, "strata_name")) "Stratum" else fill
  p + public_fill_scale(nlevels(plot_data$public_group), fill_name)
}

register_immunarch_visualisation(vis_vdj_public_impl, "vdj", "public")


#' @title
#' Visualise public receptor selection, incidence, or pairwise frequencies
#'
#' @description
#' Visualise selected receptor coverage and abundance with `plot = "public"`,
#' or receptor frequencies shared between two repertoires with the default
#' `plot = "pairwise"`. For [annotate_public()] results, `plot = "incidence"`
#' shows receptor counts by the number of repertoires in which they occur.
#'
#' In pairwise plots, receptor proportions are shown on log-scaled axes.
#' Point fill shows the log proportion ratio between the x and y repertoires.
#'
#' @param .data An [immundata::ImmunData] object with aggregated repertoires.
#' @param x Repertoire for the x-axis. Either a numeric repertoire index,
#'   repertoire ID, or a repertoire schema value.
#' @param y Repertoire for the y-axis. Either a numeric repertoire index,
#'   repertoire ID, or a repertoire schema value.
#' @param density Logical. If `TRUE`, show marginal distributions of unique
#'   and public receptor proportions.
#' @param lm Logical. If `TRUE`, add a linear model fit and adjusted R-squared.
#' @param plot One of `"pairwise"` (the default), `"public"`, or `"incidence"`.
#' @param reference For `plot = "public"`, the original `ImmunData` supplied
#'   to [vdj_public()]. It defines the complete repertoire axis, including
#'   repertoires without selected receptors. Required for this plot.
#' @param xval For `plot = "public"`, the repertoire ID or metadata column
#'   used for x-axis labels. Each value must identify one repertoire.
#' @param fill For `plot = "public"`, a repertoire metadata column used for
#'   colour. `NULL` uses `strata_name` when strata are defined; otherwise the
#'   plot uses a single colour.
#' @param title For `plot = "public"` or `"incidence"`, the plot title, or
#'   `NULL` to omit it. Each view has its own default title.
#' @param ... Must be empty.
#' @details
#' The incidence view stacks within-stratum receptor counts when strata are
#' defined. An identity observed in several strata contributes once to each
#' stratum's distribution; stack totals are not globally distinct receptor
#' counts. Without strata, one unstacked global distribution is shown and
#' every identity contributes once. This view uses statistics from the last
#' [annotate_public()] call; rerun annotation after changing the population.
#' Both publicness views use vertical legends on the right and wrap long
#' colour labels for display without changing the underlying names.
#'
#' @return A ggplot2 or patchwork object.
#'
#' @examples
#' \dontrun{
#' idata |> vis()
#' idata |> vis(x = 1, y = 2)
#' idata |> vis(x = "Blood", y = "Tumor")
#' idata |> vis(x = 1, y = 2, density = FALSE)
#' selected <- vdj_public(idata)
#' vis(selected, plot = "public", reference = idata)
#' annotate_public(idata) |> vis(plot = "incidence")
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
  plot = c("pairwise", "public", "incidence"),
  reference = NULL,
  xval = "imd_repertoire_id",
  fill = NULL,
  title = "Selected immune receptors across repertoires",
  ...
) {
  idata <- .data
  checkmate::assert_r6(idata, "ImmunData")
  checkmate::assert_flag(density)
  checkmate::assert_flag(lm)
  rlang::check_dots_empty()
  plot <- match.arg(plot)
  if (plot == "incidence") {
    if (missing(title)) {
      title <- if (is.null(idata$strata)) "Immune receptor sharing across repertoires" else
        "Immune receptor sharing within strata"
    }
    return(vis_annotate_public_impl(idata, title = title))
  }
  if (plot == "public") {
    return(vis_vdj_public_impl(idata$annotations, reference = reference,
                              xval = xval, fill = fill, title = title))
  }

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
