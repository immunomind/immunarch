#
# 1) Look at this for examples: https://chatgpt.com/c/68eacc4b-785c-832a-88f5-b94086ad7f26
#
# 2) Register in the registry. Registry should check if all visualisations are assigned to existing (!) classes.
# vis_airr_diversity_chao1 <- ...
# It is probably impossible to do...
# Should use https://stat.ethz.ch/R-manual/R-devel/RHOME/library/base/html/S3method.html
#
# 3) I need registry to control the class system in one place. One for assigning classes following a specific schema, the other one is for parsing classes.
# I probably need to write a short manual on the logic of this. And on the structure of the code. And on the phiolosophy / architectural design decisions. Same for immundata.
#
# 4) Key question: can I use different names? https://stackoverflow.com/questions/61482561/whats-the-preferred-means-for-defining-an-s3-method-in-an-r-package-without-int/61483612#61483612
# https://vctrs.r-lib.org/reference/s3_register.html
#
# 5) [!!!] I still need some name for Roxygen to @inheritParams. There is no way around that.
# opt.1 - dynamically register S3 stuff, one big vis() or some weird names for docs
# opt.2 - pre-defined classes, vis() is constructed via @inheritParams
# Or...
# - https://github.com/r-lib/roxygen2/issues/1159
# - https://github.com/rstudio/renv/blob/4a8bcb4605f085fbea5f29a76ad9a291ac2bd363/R/roxygen.R#L2-L26

#' @keywords internal
make_fixed_col_plot <- function(x_col, y_col, title, xlab, ylab) {
  checkmate::assert_string(x_col)
  checkmate::assert_string(y_col)
  checkmate::assert_string(title)
  checkmate::assert_string(xlab)
  checkmate::assert_string(ylab)

  col_plot <- function(.data, fill = immundata::imd_schema("repertoire"), facet = NULL, dir = c("h", "v"), ...) {
    checkmate::assert_data_frame(.data)
    checkmate::assert_subset(c(x_col, y_col), choices = names(.data))

    dir <- match.arg(dir)

    # factors in case the fill is continuous
    if (!is.null(fill)) {
      checkmate::assert_string(fill)
      checkmate::assert_choice(fill, names(.data))
    } else {
      fill <- immundata::imd_schema("repertoire")
    }
    .data[[fill]] <- as.factor(.data[[fill]])

    if (!is.null(facet)) {
      checkmate::assert_character(facet, any.missing = FALSE, min.len = 1, max.len = 2)
      checkmate::assert_subset(facet, names(.data))
    }

    # core plot
    p <- ggplot2::ggplot(.data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill]]))

    if (fill == immundata::imd_schema("repertoire")) {
      p <- p +
        ggplot2::geom_col(
          na.rm = TRUE,
          position = ggplot2::position_dodge(),
          colour = "grey30"
        )
    } else {
      len_vals <- sort(unique(.data[[x_col]]))
      p <- p +
        ggplot2::geom_vline(xintercept = len_vals, colour = "grey90", linewidth = 0.25) +
        ggplot2::geom_boxplot(
          mapping = ggplot2::aes(group = interaction(.data[[x_col]], .data[[fill]], drop = TRUE)),
          na.rm = TRUE,
          colour = "grey30"
        ) +
        ggplot2::geom_violin(
          mapping = ggplot2::aes(group = interaction(.data[[x_col]], .data[[fill]], drop = TRUE)),
          na.rm = TRUE,
          colour = "grey30",
          alpha = 0.3
        )
    }

    # one string -> one facet
    # vector of two strings -> square facet
    if (!is.null(facet)) {
      if (length(facet) == 1) {
        p <- p + ggplot2::facet_wrap(dir = dir, stats::as.formula(paste0("~", facet)), scales = "free_y")
      } else {
        p <- p + ggplot2::facet_grid(stats::as.formula(paste(facet[1], "~", facet[2])),
          scales = "free_y"
        )
      }
    }

    # titles
    p <- p +
      ggplot2::ggtitle(title) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)

    # auto select the palette
    if (!is.null(fill)) {
      k <- length(unique(.data[[fill]]))
      if (requireNamespace("ggsci", quietly = TRUE) && k <= 11) {
        p <- p + ggsci::scale_fill_locuszoom(name = fill)
      } else {
        p <- p + ggplot2::scale_fill_viridis_d(option = "H")
      }
    }

    # light theme
    if (requireNamespace("ggthemes", quietly = TRUE)) {
      p <- p + ggthemes::theme_few()
    }

    p
  }

  col_plot
}


#' @keywords internal
make_dynam_col_plot <- function(y_default, title_default, position = c("dodge", "stack")) {
  checkmate::assert_string(y_default)
  checkmate::assert_string(title_default)
  position <- match.arg(position)

  col_plot <- function(.data, xval = immundata::imd_schema("repertoire"), yval = y_default, fill = NULL, facet = NULL, title = title_default, dir = c("h", "v"), ...) {
    checkmate::assert_data_frame(.data)
    checkmate::assert_subset(c(xval, yval), choices = names(.data))

    dir <- match.arg(dir)

    .data[[xval]] <- as.factor(.data[[xval]])

    # factors in case the fill is continuous
    if (!is.null(fill)) {
      checkmate::assert_string(fill)
      checkmate::assert_choice(fill, names(.data))
      .data[[fill]] <- as.factor(.data[[fill]])
    }

    if (!is.null(facet)) {
      checkmate::assert_character(facet, any.missing = FALSE, min.len = 1, max.len = 2)
      checkmate::assert_subset(facet, names(.data))
    }

    # check how many observation per group so we could select either bar or box plot
    grp_vars <- c(xval, if (!is.null(fill)) fill else NULL)
    n_by_grp <- .data |>
      dplyr::group_by(dplyr::across(all_of(grp_vars))) |>
      dplyr::summarise(n = sum(!is.na(.data[[yval]])), .groups = "drop")
    all_single <- all(n_by_grp$n == 1)

    # core plot
    p <- ggplot2::ggplot(.data, ggplot2::aes(x = .data[[xval]], y = .data[[yval]]))

    if (all_single) {
      # if (xval == immundata::imd_schema("repertoire") || (!is.null(fill) && fill == immundata::imd_schema("repertoire"))) {
      p <- p +
        ggplot2::geom_col(
          mapping = if (!is.null(fill)) ggplot2::aes(fill = .data[[fill]]) else NULL,
          na.rm = TRUE,
          position = position,
          colour = "grey30"
        )
    } else {
      p <- p +
        ggplot2::geom_boxplot(
          mapping = if (!is.null(fill)) ggplot2::aes(group = interaction(.data[[xval]], .data[[fill]], drop = TRUE), fill = .data[[fill]]) else ggplot2::aes(fill = .data[[xval]]),
          na.rm = TRUE,
          colour = "grey30"
        ) +
        ggplot2::geom_violin(
          mapping = if (!is.null(fill)) ggplot2::aes(group = interaction(.data[[xval]], .data[[fill]], drop = TRUE), fill = .data[[fill]]) else ggplot2::aes(fill = .data[[xval]]),
          colour = "grey30",
          alpha = 0.3,
          position = if (!is.null(fill)) ggplot2::position_dodge(width = 0.75) else "identity"
        )
    }


    # one string -> one facet
    # vector of two strings -> square facet
    if (!is.null(facet)) {
      if (length(facet) == 1) {
        p <- p + ggplot2::facet_wrap(dir = dir, stats::as.formula(paste0("~", facet)), scales = "free_y")
      } else {
        p <- p + ggplot2::facet_grid(stats::as.formula(paste(facet[1], "~", facet[2])),
          scales = "free_y"
        )
      }
    }

    # titles
    p <- p +
      ggplot2::ggtitle(title)

    # auto select the palette
    if (!is.null(fill)) {
      k <- length(unique(.data[[fill]]))
    } else {
      k <- length(unique(.data[[xval]]))
    }
    if (requireNamespace("ggsci", quietly = TRUE) && k <= 11) {
      p <- p + ggsci::scale_fill_locuszoom(name = fill)
    } else {
      p <- p + ggplot2::scale_fill_viridis_d(option = "H")
    }

    # light theme
    if (requireNamespace("ggthemes", quietly = TRUE)) {
      p <- p + ggthemes::theme_few() +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank()
        ) +
        ggplot2::scale_x_discrete(guide = guide_axis(angle = 45))
    }

    p
  }

  col_plot
}


#' @keywords internal
make_dotplot <- function(title_default, fill_default, size_default) {
  checkmate::assert_string(title_default)
  checkmate::assert_string(fill_default)
  checkmate::assert_string(size_default)

  dotplot_matrix <- function(.data,
                             size_max_mm = 7,
                             row = NULL, # string: row id (e.g., "v_call")
                             col = NULL, # string or character vector: column id(s) (e.g., "Cluster" or c("Cluster","Tissue"))
                             value = NULL, # string: numeric column for values; if NULL and "n" exists, uses "n"; else counts
                             col_sep = " | ", # used when col has length > 1
                             row_order = NULL,
                             col_order = NULL) {
    # ---------- input checks ----------
    if (is.matrix(.data)) {
      checkmate::assert_matrix(.data, mode = "numeric", any.missing = TRUE, min.rows = 1, min.cols = 1)
      mat <- .data
    } else {
      checkmate::assert_data_frame(.data, min.rows = 1, min.cols = 3)

      # Infer sensible defaults if not provided
      rep_col <- immundata::imd_schema("repertoire")

      if (is.null(row)) {
        # default row key = first non-repertoire column
        row <- setdiff(names(.data), rep_col)[1]
      }
      checkmate::assert_string(row)
      checkmate::assert_choice(row, names(.data))

      if (is.null(col)) {
        # default columns = repertoire
        col <- rep_col
      }
      checkmate::assert_character(col, any.missing = FALSE, min.len = 1)
      checkmate::assert_subset(col, names(.data))

      # value: prefer provided; else "n"; else count per cell
      if (!is.null(value)) {
        checkmate::assert_string(value)
        checkmate::assert_choice(value, names(.data))
        checkmate::assert_numeric(.data[[value]], any.missing = TRUE)
      } else if ("n" %in% names(.data)) {
        value <- "n"
        checkmate::assert_numeric(.data[[value]], any.missing = TRUE)
      } else {
        value <- NULL # we'll compute counts
      }

      df <- .data
      df[[row]] <- as.character(df[[row]])
      if (length(col) > 1L) {
        df <- tidyr::unite(df, "__col__", tidyselect::all_of(col), sep = col_sep, remove = FALSE)
        col_key <- "__col__"
      } else {
        col_key <- col
      }
      df[[col_key]] <- as.character(df[[col_key]])

      # Aggregate to one value per (row, col)
      if (is.null(value)) {
        df_sum <- df |>
          dplyr::group_by(.data[[row]], .data[[col_key]]) |>
          dplyr::summarise(.val = dplyr::n(), .groups = "drop")
      } else {
        df_sum <- df |>
          dplyr::group_by(.data[[row]], .data[[col_key]]) |>
          dplyr::summarise(.val = sum(.data[[value]], na.rm = TRUE), .groups = "drop")
      }

      # Pivot wider → numeric matrix
      wide <- tidyr::pivot_wider(
        df_sum,
        names_from  = tidyselect::all_of(col_key),
        values_from = .val,
        values_fill = 0
      )
      rn <- wide[[row]]
      checkmate::assert_atomic_vector(rn, all.missing = FALSE, min.len = 1)
      mat_df <- as.data.frame(wide, check.names = FALSE)
      mat_df[[row]] <- NULL
      mat <- as.matrix(mat_df)
      rownames(mat) <- rn
    }

    checkmate::assert_number(size_max_mm, lower = 0, finite = TRUE)

    #
    # Plotting
    #
    rn <- rownames(mat)
    if (is.null(rn)) rn <- as.character(seq_len(nrow(mat)))
    cn <- colnames(mat)
    if (is.null(cn)) cn <- as.character(seq_len(ncol(mat)))

    df <- as.data.frame(mat, check.names = FALSE)
    df$row <- rn
    df <- tidyr::pivot_longer(df, -row, names_to = "col", values_to = "value")

    # TODO: maybe make size relative to the source repertoire?
    rng <- range(df$value, na.rm = TRUE)
    df$size <- if (is.finite(diff(rng)) && diff(rng) > 0) (df$value - rng[1]) / diff(rng) else 0

    if (is.null(row_order)) row_order <- rn
    if (is.null(col_order)) col_order <- cn
    df$row <- factor(df$row, levels = rev(row_order))
    df$col <- factor(df$col, levels = col_order)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row)) +
      ggplot2::geom_point(
        ggplot2::aes(size = size, fill = value),
        colour = "grey30",
        stroke = 0.3,
        shape = 21, na.rm = TRUE
      ) +
      ggplot2::scale_size_area(
        max_size = size_max_mm,
        guide = ggplot2::guide_legend(title = size_default)
      ) +
      ggplot2::scale_fill_gradientn(
        colors = c("#0F172A", "#243B55", "#FFE3A3", "#FFC107"),
        na.value = "grey90",
        guide = ggplot2::guide_colorbar(title = fill_default)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = title_default) +
      ggplot2::coord_cartesian(clip = "off") +
      ggthemes::theme_few(base_size = 12) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )

    p
  }

  dotplot_matrix
}
