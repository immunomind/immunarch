if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "value", "variable", "Clonotype", "Count", "Sample", "Clone",
    "Value", "Value.mean", "Value.min", "Value.max", "quantile",
    "Kmer", "Group", "Size", "Val", "Grouping.var", "Xrep", "Yrep",
    "Clone.num", "Q", "Clonotypes", "Sample_subj", "Seq.count",
    "Overlap", "head", "Mean", "MeanVal", "MinVal", "MaxVal",
    "Q1", "Q2", "Type", "Length", "Gene", "Freq", "Sequence",
    "AA", "Clones", "Source.gr", "Target.gr", "Samples", "Samples.y",
    "CDR3.aa", "p.adj", "group1", "group2", "y.coord", "..p.adj..", ".SD",
    "name", "label", "."
  ))
}


##### Utility functions #####


.rem_legend <- function(.p) {
  .p + theme(legend.position = "none")
}


.rem_xlab <- function(.p) {
  .p + theme(axis.title.x = element_blank())
}


.rem_ylab <- function(.p) {
  .p + theme(axis.title.y = element_blank())
}


.colourblind_vector <- function() {
  c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6")
}


.colourblind_discrete <- function(.n, .colour = FALSE) {
  cs <- .colourblind_vector()
  if (.colour) {
    scale_colour_manual(values = colorRampPalette(cs)(.n))
  } else {
    scale_fill_manual(values = colorRampPalette(cs)(.n))
  }
}


.tweak_fill <- function(.n) {
  palette_name <- ""
  if (.n == 1) {
    palette_name <- "Set2"
  } else if (.n == 2) {
    palette_name <- "Set1"
  }
  # else if (.n < 4) { palette_name = "YlGnBu" }
  # else if (.n < 6) {  palette_name = "RdBu" }
  else if (.n < 12) {
    palette_name <- "Spectral"
  } else {
    return(scale_fill_hue())
  }

  scale_fill_brewer(palette = palette_name)
}

.tweak_col <- function(.n) {
  palette_name <- ""
  if (.n == 1) {
    palette_name <- "Set2"
  } else if (.n == 2) {
    palette_name <- "Set1"
  }
  # else if (.n < 4) { palette_name = "YlGnBu" }
  # else if (.n < 6) {  palette_name = "RdBu" }
  else if (.n < 12) {
    palette_name <- "Spectral"
  } else {
    return(scale_colour_hue())
  }

  scale_color_brewer(palette = palette_name)
}

theme_cleveland2 <- function(rotate = TRUE) {
  if (rotate) {
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(
        colour = "grey70",
        linetype = "dashed"
      )
    )
  } else {
    theme(
      panel.grid.major.x = element_line(
        colour = "grey70",
        linetype = "dashed"
      ), panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  }
}


##### The one and only - the ultimate vis() function #####


#' One function to visualise them all
#'
#' @concept vis
#'
#' @name vis
#'
#' @import ggplot2
#' @importFrom factoextra fviz_cluster fviz_dend fviz_pca_ind
#' @importFrom grDevices colorRampPalette
#'
#' @description Output from every function in immunarch can be visualised with a
#' single function - \code{vis}. The \code{vis} automatically detects
#' the type of the data and draws a proper visualisation. For example, output
#' from the \code{repOverlap} function will be identified as repertoire overlap values
#' and respective visualisation will be chosen without any additional arguments.
#' See "Details" for the list of available visualisations.
#'
#' @param .data Pass the output from any immunarch analysis tool to \code{vis()}.
#' @param ... Any other arguments, see the "Details" section for specific visualisation functions.
#'
#' @details
#' List of available visualisations for different kinds of data.
#'
#' Basic analysis:
#'
#' - Exploratory analysis results (from \link{repExplore}) - see \link{vis.immunr_exp_vol};
#'
#' - Clonality statistics (from \link{repClonality}) - see \link{vis.immunr_homeo}.
#'
#' Overlaps and public clonotypes:
#'
#' - Overlaps (from \link{repOverlap}) using heatmaps, circos plots, polar area plots - see \link{vis.immunr_ov_matrix};
#'
#' -  Overlap clustering (from \link{repOverlapAnalysis}) - see \link{vis.immunr_hclust};
#'
#' - Repertoire incremental overlaps (from \link{repOverlap}) - see \link{vis.immunr_inc_overlap};
#'
#' - Public repertoire abundance (from \link{pubRep}) - vis \link{vis.immunr_public_repertoire}.
#'
#' Gene usage:
#'
#' - Gene usage statistics (from \link{geneUsage}) using bar plots, box plots - see \link{vis.immunr_gene_usage};
#'
#' - Gene usage distances (from \link{geneUsageAnalysis}) using heatmaps, circos plots, polar area plots - see \link{vis.immunr_ov_matrix};
#'
#' - Gene usage clustering (from \link{geneUsageAnalysis}) - see \link{vis.immunr_hclust}.
#'
#' Diversity estimation:
#'
#' - Diversity estimations (from \link{repDiversity}) - see \link{vis.immunr_chao1}.
#'
#' Advanced analysis:
#'
#' - Repertoire dynamics (from \link{trackClonotypes}) - see \link{vis.immunr_dynamics};
#'
#' - Sequence logo plots of amino acid distributions (from \link{kmer_profile}) - see \link{vis_seqlogo};
#'
#' - Kmers distributions (from \link{getKmers}) - see \link{vis.immunr_kmer_table};
#'
#' - Mutation networks (from mutationNetwork) - Work In Progress on vis.immunr_mutation_network;
#'
#' - CDR3 amino acid properties, e.g., biophysical (from cdrProp) - Work In Progress on vis.immunr_cdr_prop.
#'
#' Additionaly, we provide a wrapper functions for visualisations of common data types:
#'
#' - Any data frames or matrices using heatmaps - see \link{vis_heatmap} and \link{vis_heatmap2};
#'
#' - Any data frames or matrices using circos plots - see \link{vis_circos}.
#'
#' @return
#' A ggplot2, pheatmap or circlize object.
#'
#' @seealso \link{fixVis} for precise manipulation of plots.
#'
#' @examples
#' # Load the test data
#' data(immdata)
#'
#' # Compute and visualise:
#' ov <- repOverlap(immdata$data)
#' vis(ov)
#'
#' gu <- geneUsage(immdata$data)
#' vis(gu)
#'
#' dv <- repDiversity(immdata$data)
#' vis(dv)
#' @export vis
vis <- function(.data, ...) {
  UseMethod("vis")
}


##### Overlap & heatmaps, circos plots and polar area plots #####


#' Repertoire overlap and gene usage visualisations
#'
#' @concept overlap
#'
#' @name vis.immunr_ov_matrix
#'
#' @aliases vis.immunr_ov_matrix vis.immunr_gu_matrix
#'
#' @description Visualise matrices with overlap values or gene usage distances among samples.
#' For details see links below.
#'
#' @param .data Output from \link{repOverlap} or \link{geneUsageAnalysis}.
#'
#' @param .plot A string specifying the plot type:
#'
#' - "heatmap" for heatmaps using \link{vis_heatmap};
#'
#' - "heatmap2" for heatmaps using \link{vis_heatmap2};
#'
#' - "circos" for circos plots using \link{vis_circos};
#'
#' @param ... Other arguments are passed through to the underlying plotting function:
#'
#' - "heatmap" - passes arguments to \link{vis_heatmap};
#'
#' - "heatmap2" - passes arguments to \link{vis_heatmap2} and \link{heatmap} from the "pheatmap" package;
#'
#' - "circos" - passes arguments to \link{vis_circos} and \link{chordDiagram} from the "circlize" package;
#'
#' @return
#' A ggplot2, pheatmap or circlize object.
#'
#' @examples
#' data(immdata)
#' ov <- repOverlap(immdata$data)
#' vis(ov)
#' vis(ov, "heatmap")
#' vis(ov, "heatmap2")
#' vis(ov, "circos")
#' @export
vis.immunr_ov_matrix <- function(.data, .plot = c("heatmap", "heatmap2", "circos"), ...) {
  args <- list(...)
  if (!(".title" %in% names(args))) {
    args$.title <- "Repertoire overlap"
  }

  .plot <- .plot[1]
  heatfun <- vis_heatmap
  if (.plot == "heatmap2") {
    heatfun <- vis_heatmap2
  } else if (.plot == "circos") {
    heatfun <- vis_circos
  }
  do.call(heatfun, c(list(.data = .data), args))
}

#' @export
vis.immunr_gu_matrix <- function(.data, .plot = c("heatmap", "heatmap2", "circos"), ...) {
  args <- list(...)
  if (!(".title" %in% names(args))) {
    args$.title <- "Gene usage"
  }

  .plot <- .plot[1]
  heatfun <- vis_heatmap
  if (.plot == "heatmap2") {
    heatfun <- vis_heatmap2
  } else if (.plot == "circos") {
    heatfun <- vis_circos
  }
  do.call(heatfun, c(list(.data = .data), args))
}


#' Visualisation of matrices and data frames using ggplo2-based heatmaps
#'
#' @concept vis
#'
#' @name vis_heatmap
#'
#' @aliases vis_heatmap
#'
#' @description Fast and easy visualisations of matrices or data frames
#' with functions based on the ggplot2 package.
#'
#' @param .data Input object: a matrix or a data frame.
#'
#' If matrix: column names and row names (if presented) will be used as names for labs.
#'
#' If data frame: the first column will be used for row names and removed from the data.
#' Other columns will be used for values in the heatmap.
#'
#' @param .text If TRUE then plot values in the heatmap cells. If FALSE do not plot values,
#' just plot coloured cells instead.
#'
#' @param .scientific If TRUE then use the scientific notation for numbers (e.g., "2.0e+2").
#'
#' @param .signif.digits Number of significant digits to display on plot.
#'
#' @param .text.size Size of text in the cells of heatmap.
#'
#' @param .axis.text.size Size of text on the axis labels.
#'
#' @param .labs A character vector of length two with names for x-axis and y-axis, respectively.
#'
#' @param .title The The text for the plot's title.
#'
#' @param .leg.title The The text for the plots's legend. Provide NULL to remove the legend's title completely.
#'
#' @param .legend If TRUE then displays a legend, otherwise removes legend from the plot.
#'
#' @param .na.value Replace NA values with this value. By default they remain NA.
#'
#' @param .transpose Logical. If TRUE then switch rows and columns.
#'
#' @param ... Other passed arguments.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \link{vis}, \link{repOverlap}.
#'
#' @examples
#' data(immdata)
#' ov <- repOverlap(immdata$data)
#' vis_heatmap(ov)
#' gu <- geneUsage(immdata$data, "hs.trbj")
#' vis_heatmap(gu)
#' @export
vis_heatmap <- function(.data, .text = TRUE, .scientific = FALSE, .signif.digits = 2,
                        .text.size = 4, .axis.text.size = NULL,
                        .labs = c("Sample", "Sample"), .title = "Overlap",
                        .leg.title = "Overlap values", .legend = TRUE,
                        .na.value = NA, .transpose = FALSE, ...) {
  if (has_class(.data, "data.frame")) {
    names <- .data[[1]]
    .data <- as.matrix(.data[, -1])
    row.names(.data) <- names

    if (.transpose) {
      .data <- t(.data)
    }
  } else if (is.null(dim(.data))) {
    .data <- as.matrix(.data)

    if (.transpose) {
      .data <- t(.data)
    }
  }

  if (is.null(colnames(.data))) {
    colnames(.data) <- paste0("C", seq_len(ncol(.data)))
  }

  if (is.null(row.names(.data))) {
    row.names(.data) <- paste0("C", seq_len(nrow(.data)))
  }

  .data[is.na(.data)] <- .na.value

  tmp <- as.data.frame(.data)
  tmp$name <- row.names(.data)

  m <- reshape2::melt(tmp, id.var = c("name"))
  m[, 1] <- factor(m[, 1], levels = rev(rownames(.data)))
  m[, 2] <- factor(m[, 2], levels = colnames(.data))

  format <- if (.scientific) "g" else "fg"
  m$label <- formatC(m$value, format = format, digits = .signif.digits)

  p <- ggplot(m, aes(x = variable, y = name, fill = value))
  p <- p + geom_tile(aes(fill = value), colour = "white")

  if (.text) {
    p <- p + geom_text(aes(label = label), size = .text.size)
  }

  p <- p + scale_fill_distiller(palette = "RdBu")

  p <- p + ggtitle(.title) +
    guides(fill = guide_colourbar(title = .leg.title)) +
    xlab(.labs[1]) + ylab(.labs[2]) + coord_fixed() +
    theme_linedraw() + theme(
      axis.text.x =
        element_text(angle = 90, vjust = .5, size = .axis.text.size)
    ) + theme(
      axis.text.y =
        element_text(size = .axis.text.size)
    ) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))

  if (!.legend) {
    p <- .rem_legend(p)
  }

  if (is.null(.labs)) {
    p <- .rem_xlab(.rem_ylab(p))
  }

  p
}


#' Visualisation of matrices using pheatmap-based heatmaps
#'
#' @concept vis
#'
#' @importFrom pheatmap pheatmap
#'
#' @name vis_heatmap2
#'
#' @description Visualise matrices with the functions based on the \link[pheatmap]{pheatmap}
#' package with minimum amount of arguments.
#'
#' @param .data Input matrix. Column names and row names (if presented) will be used as names for labs.
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as age, serostatus or hla.
#'
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' @param .title The text for the plot's title (same as the "main" argument in \link[pheatmap]{pheatmap}).
#'
#' @param .color A vector specifying the colors (same as the "color" argument in \link[pheatmap]{pheatmap}).
#' Pass NA to use the default pheatmap colors.
#'
#' @param ... Other arguments for the \link[pheatmap]{pheatmap} function.
#'
#' @return
#' A pheatmap object.
#'
#' @seealso \link{vis}, \link{repOverlap}
#'
#' @examples
#' data(immdata)
#' ov <- repOverlap(immdata$data)
#' vis_heatmap2(ov)
#' @export
vis_heatmap2 <- function(.data, .meta = NA, .by = NA, .title = NA, .color = colorRampPalette(c("#67001f", "#d6604d", "#f7f7f7", "#4393c3", "#053061"))(1024), ...) {
  args <- list(...)
  if (!is.na(.by)[1]) {
    if (!is.na(.meta)[1]) {
      args[["annotation_col"]] <- .meta %>%
        tibble::column_to_rownames(var = "Sample") %>%
        dplyr::select(tidyselect::any_of(.by))
    }
  }

  args[["mat"]] <- .data
  args[["main"]] <- .title
  if (!is.na(.color)[1]) {
    args[["color"]] <- .color
  }
  do.call(pheatmap, args)
}


#' Visualisation of matrices using circos plots
#'
#' @concept vis
#'
#' @importFrom circlize chordDiagram
#'
#' @name vis_circos
#'
#' @description Visualise matrices with the \link{chordDiagram} function
#' from the circlize package.
#'
#' @param .data Input matrix.
#'
#' @param .title The The text for the title of the plot.
#'
#' @param ... Other arguments passed to \link{chordDiagram} from the 'circlize' package.
#'
#' @return
#' A circlize object.
#'
#' @seealso \link{vis}, \link{repOverlap}.
#'
#' @examples
#' data(immdata)
#' ov <- repOverlap(immdata$data)
#' vis(ov, .plot = "circos")
#' @export
vis_circos <- function(.data, .title = NULL, ...) {
  if (has_class(.data, "tibble")) {
    .data <- as.data.frame(.data)
  }
  chordDiagram(.data, ...)
}


# vis_radar <- function(.data, .by = NA, .meta = NA,
#                       .ncol = NA, .which = NA,
#                       .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .expand = c(.25, 0),
#                       .title = NA, .subtitle = NULL,
#                       .legend = FALSE, .leg.title = NULL, ...) {
#   stop("vis_radar() is under complety re-working. Please use other functions such as vis_heatmap.
# Sorry for the inconvenience! Please contact us if you need vis_radar() ASAP.")
#
#   .data <- melt(.data)
#   colnames(.data) <- c("Source", "Target", "Value")
#   .data$Value[is.na(.data$Value)] <- 0
#
#   group_res <- process_metadata_arguments(.data, .by, .meta, .data.sample.col = "Source")
#   group_res2 <- process_metadata_arguments(.data, .by, .meta, .data.sample.col = "Target")
#   .data$Source.gr <- group_res$group_column
#   .data$Target.gr <- group_res2$group_column
#
#   .data_proc <- .data %>%
#     dplyr::group_by(Source.gr, Target.gr) %>%
#     dplyr::summarise(
#       Value = mean(Value, na.rm = TRUE),
#       Value.min = quantile(Value, .errorbars[1], na.rm = TRUE),
#       Value.max = quantile(Value, .errorbars[2], na.rm = TRUE)
#     )
#   .data_proc <- .data_proc[!is.na(.data_proc$Source.gr), ]
#
#   step <- length(unique(.data_proc$Source.gr))
#
#   if (!is.na(.which[1])) {
#     .data_proc <- .data_proc[.data_proc$Source.gr %in% .which, ]
#   }
#
#   if (is.na(.leg.title)) {
#     if (group_res$is_grouped) {
#       .leg.title <- "Group"
#     } else {
#       .leg.title <- "Sample"
#     }
#   }
#
#   if (!group_res$is_grouped) {
#     .errorbars.off <- TRUE
#   }
#
#   if (is.na(.ncol)) {
#     .ncol <- round(length(unique(.data_proc$Source.gr))**.5)
#   }
#
#   .data_proc$Value[.data_proc$Source.gr == .data_proc$Target.gr] <- 0
#   .data_proc$Value.min[.data_proc$Source.gr == .data_proc$Target.gr] <- NA
#   .data_proc$Value.max[.data_proc$Source.gr == .data_proc$Target.gr] <- NA
#
#   ps <- lapply(seq(1, nrow(.data_proc), step), function(l) {
#     p <- ggplot(.data_proc[l:(l + step - 1), ], aes(x = Target.gr, y = Value, fill = Target.gr)) +
#       geom_bar(colour = "black", stat = "identity")
#
#     if (!.errorbars.off) {
#       p <- p +
#         geom_errorbar(aes(ymin = Value.min, ymax = Value.max), width = 0.2, position = position_dodge(.9))
#     }
#
#     p <- p +
#       coord_polar() +
#       scale_y_continuous(expand = .expand) +
#       theme_linedraw() +
#       .colourblind_discrete(step)
#
#     if (!.legend) {
#       p <- p + guides(fill = "none")
#     } else {
#       p <- p + guides(fill = guide_legend(title = .leg.title))
#     }
#
#     p <- p +
#       labs(title = .data_proc$Source.gr[l], subtitle = .subtitle)
#
#     p
#   })
#   do.call(grid.arrange, c(ps, ncol = .ncol, top = .title))
# }


##### Complex overlaps - top overlap & public repertoires #####


#' Visualise incremental overlaps
#'
#' @concept overlap
#'
#' @importFrom dplyr rename
#'
#' @name vis.immunr_inc_overlap
#'
#' @param .data Output from the \link{repOverlap} function that uses "top" methods.
#'
#' @param .target Index of a repertoire to plot. Omitted if .grid is TRUE.
#'
#' @param .grid Logical. If TRUE then plot all similarities in a grid.
#'
#' @param .ncol Numeric. Number of columns in the resulting grid.
#'
#' @param ... Not used here.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \link{repOverlap}
#'
#' @examples
#' data(immdata)
#' tmp <- repOverlap(immdata$data[1:4], "inc+overlap", .verbose.inc = FALSE, .verbose = FALSE)
#' vis(tmp, .target = 1)
#' vis(tmp, .grid = TRUE)
#' @export
vis.immunr_inc_overlap <- function(.data, .target = 1, .grid = FALSE, .ncol = 2, ...) {
  data_is_bootstrapped <- !is.null(attr(.data, "bootstrap"))

  if (.grid) {
    if (data_is_bootstrapped) {
      sample_names <- colnames(.data[[1]][[1]])
    } else {
      sample_names <- colnames(.data[[1]])
    }

    p_list <- lapply(seq_along(sample_names), function(i_name) {
      p <- vis(.data, .target = i_name) + ggtitle(sample_names[i_name]) + theme_pubr(legend = "right") + theme_cleveland2()
      p
    })

    return(do.call(wrap_plots, c(p_list, list(ncol = .ncol))))
  } else {
    if (data_is_bootstrapped) {
      sample_names <- colnames(.data[[1]][[1]])

      .data <- lapply(.data, function(mat_list) {
        lapply(mat_list, function(mat) {
          replace(mat, lower.tri(mat, TRUE), NA)
        })
      })
      .data <- reshape2::melt(.data, na.rm = TRUE)
      names(.data) <- c("Sample_subj", "Sample", "Overlap", "N", "Seq.count")

      filtered_data <- .data[.data$Sample_subj == sample_names[.target], ]
      filtered_data$Seq.count <- as.numeric(filtered_data$Seq.count)

      filtered_data <- filtered_data %>%
        dplyr::group_by(Sample_subj, Sample, Seq.count, Overlap) %>%
        dplyr::summarise(
          Value = mean(Overlap, na.rm = TRUE),
          Value.min = quantile(Overlap, 0.025, na.rm = TRUE),
          Value.max = quantile(Overlap, 0.975, na.rm = TRUE)
        )
    } else {
      sample_names <- colnames(.data[[1]])

      .data <- lapply(.data, function(mat) replace(mat, lower.tri(mat, TRUE), NA))
      .data <- reshape2::melt(.data, na.rm = TRUE)
      names(.data) <- c("Sample_subj", "Sample", "Overlap", "Seq.count")

      filtered_data <- .data[.data$Sample_subj == sample_names[.target], ]
      filtered_data$Seq.count <- as.numeric(filtered_data$Seq.count)
      filtered_data <- filtered_data %>%
        dplyr::rename(Value = Overlap)
    }

    p <- ggplot(aes(x = Seq.count, y = Value, col = Sample, group = Sample), data = filtered_data) +
      geom_line()

    if (data_is_bootstrapped) {
      p <- p +
        geom_ribbon(aes(group = Sample, ymin = Value.min, ymax = Value.max, fill = Sample), alpha = 0.15)
    }

    p <- p +
      theme_pubr(legend = "right") + theme_cleveland2() +
      ggtitle("Similarity of repertoires") +
      xlab("Sequence count") + ylab("Similarity")

    return(p)
  }
}


#' Public repertoire visualisation
#'
#' @concept pubrep
#'
#' @name vis.immunr_public_repertoire
#'
#' @param .data Public repertoire, an output from \link{pubRep}.
#' @param .plot A string specifying the plot type:
#'
#' - "freq" for visualisation of the distribution of occurrences of clonotypes
#' and their frequencies using \link{vis_public_frequencies}.
#'
#' - "clonotypes" for visualisation of public clonotype frequenciy correlations between pairs of
#' samples using \link{vis_public_clonotypes}
#' @param ... Further arguments passed \link{vis_public_frequencies} or \link{vis_public_clonotypes},
#' depending on the ".plot" argument.
#'
#' @return
#' A ggplot2 object.
#'
#' @examples
#' data(immdata)
#' immdata$data <- lapply(immdata$data, head, 300)
#' pr <- pubRep(immdata$data, .verbose = FALSE)
#' vis(pr, "freq")
#' vis(pr, "freq", .type = "none")
#'
#' vis(pr, "clonotypes", 1, 2)
#' @export
vis.immunr_public_repertoire <- function(.data, .plot = c("freq", "clonotypes"), ...) {
  .plot <- .plot[1]
  if (.plot == "freq") {
    vis_public_frequencies(.data, ...)
  } else if (.plot == "clonotypes") {
    vis_public_clonotypes(.data, ...)
  } else {
    stop('Error: unknown argument for .plot, please provide one of the following: "freq", "clonotypes" or "???"')
  }
}


#' Visualise sharing of clonotypes among samples
#'
#' @concept pubrep
#'
#' @importFrom UpSetR upset fromExpression
#'
#' @name vis.immunr_public_statistics
#'
#' @description Visualise public clonotype frequencies.
#'
#' @param .data Public repertoire - an output from the \link{pubRep} function.
#'
#' @param ... Other arguments passsed directly to \link{upset}.
#'
#' @return
#' A ggplot2 object.
#'
#' @examples
#' data(immdata)
#' immdata$data <- lapply(immdata$data, head, 2000)
#' pr <- pubRep(immdata$data, .verbose = FALSE)
#' pubRepStatistics(pr) %>% vis()
#' @export
vis.immunr_public_statistics <- function(.data, ...) {
  upsetr_data <- as.list(.data$Count)
  names(upsetr_data) <- .data$Group
  upset(fromExpression(upsetr_data), ...)
}


#' Public repertoire visualisation
#'
#' @concept pubrep
#'
#' @name vis_public_frequencies
#'
#' @description Visualise public clonotype frequencies.
#'
#' @param .data Public repertoire - an output from the \link{pubRep} function.
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as age, serostatus or hla.
#' @param .type Character. Either "boxplot" for plotting distributions of frequencies,
#' "none" for plotting everything, or "mean" for plotting average values only.
#'
#' @return
#' A ggplot2 object.
#'
#' @examples
#' data(immdata)
#' immdata$data <- lapply(immdata$data, head, 500)
#' pr <- pubRep(immdata$data, .verbose = FALSE)
#' vis(pr, "freq", .type = "boxplot")
#' vis(pr, "freq", .type = "none")
#' vis(pr, "freq", .type = "mean")
#' vis(pr, "freq", .by = "Status", .meta = immdata$meta)
vis_public_frequencies <- function(.data, .by = NA, .meta = NA,
                                   .type = c("boxplot", "none", "mean")) {
  .type <- .type[1]

  # ToDo: make it work with V genes
  melted_pr <- reshape2::melt(.data, id.vars = colnames(.data)[1:(match("Samples", colnames(.data)))])

  colnames(melted_pr)[1] <- "Sequence"
  colnames(melted_pr)[ncol(melted_pr) - 1] <- "Sample"
  group_res <- process_metadata_arguments(melted_pr, .by, .meta)
  group_column <- group_res$name
  melted_pr$Group <- group_res$group_column

  if (.type == "boxplot") {
    p <- ggplot() +
      geom_boxplot(aes(y = value, x = as.factor(Samples), col = Group), data = melted_pr) +
      scale_y_log10() +
      ggtitle("Distribution of public clonotype frequencies") +
      xlab("Number of shared samples") +
      ylab("Public clonotype frequencies")
  } else if (.type == "mean") {
    melted_pr <- melted_pr %>%
      group_by(Sequence, Samples) %>%
      summarise(Value = mean(value, na.rm = TRUE)) %>%
      as.data.table()

    p <- ggplot() +
      geom_jitter(aes(x = Value, y = Samples), data = melted_pr) +
      scale_x_log10() +
      ggtitle("Mean frequencies of public clonotypes") +
      ylab("Number of shared samples") +
      xlab("Public clonotype frequencies")
  } else {
    # ToDo: add the .verbose argument for this
    message("Warning! Visualising ", sum(!is.na(melted_pr$value)), ' points. Too many points may take a while to visualise depending on your hardware. We highly recommend you to use other types: "mean" or "boxplot"')
    p <- ggplot() +
      geom_jitter(aes(x = value, y = Samples, col = Group), data = melted_pr) +
      scale_x_log10() +
      ggtitle("Distribution of public clonotype frequencies") +
      ylab("Number of shared samples") +
      xlab("Public clonotype frequencies")
  }

  # tmp = as_tibble(data.frame(N_samples = .data$Samples,
  #                            Value = rowMeans(public_matrix(.data), na.rm = TRUE)))

  # ggplot(aes(x=Value, y=N_samples, color=N_samples), data=tmp) + geom_jitter() +
  #   scale_x_log10() +
  #   xlab("Public clonotype frequency") + ylab("Number of samples") +
  #   ggtitle("Public clonotype frequencies") +
  #   theme_bw()

  p + theme_pubr(legend = "right") + theme_cleveland2()
}



#' Visualisation of public clonotypes
#'
#' @concept pubrep
#'
#' @name vis_public_clonotypes
#'
#' @importFrom grid rectGrob gpar
#' @importFrom stats lm
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @description Visualise correlation of public clonotype frequencies in pairs of repertoires.
#'
#' @param .data Public repertoire data - an output from the \link{pubRep} function.
#'
#' @param .x.rep Either indices of samples or character vector of sample names
#' for the x-axis. Must be of the same length as ".y.rep".
#'
#' @param .y.rep Either indices of samples or character vector of sample names
#' for the y-axis. Must be of the same length as ".x.rep".
#'
#' @param .title The text for the title of the plot.
#'
#' @param .ncol An integer number of columns to print in the grid of pairs of repertoires.
#'
#' @param .point.size.modif An integer value that is a modifier of the point size.
#' The larger the number, the larger the points.
#'
#' @param .cut.axes If TRUE then axes limits become shorter.
#'
#' @param .density If TRUE then displays density plot for distributions of clonotypes
#' for each sample. If FALSE then removes density plot from the visualisation.
#'
#' @param .lm If TRUE then fit a linear model and displays an R adjusted coefficient
#' that shows how similar samples are in terms of shared clonotypes.
#'
#' @param .radj.size An integer value, that defines the size of the The text
#' for the R adjusted coefficient.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \link{pubRep}, \link{vis.immunr_public_repertoire}
#'
#' @examples
#' data(immdata)
#' pr <- pubRep(immdata$data, .verbose = FALSE)
#' vis(pr, "clonotypes", 1, 2)
vis_public_clonotypes <- function(.data, .x.rep = NA, .y.rep = NA,
                                  .title = NA, .ncol = 3,
                                  .point.size.modif = 1, .cut.axes = TRUE,
                                  .density = TRUE, .lm = TRUE, .radj.size = 3.5) {
  .shared.rep <- .data

  mat <- public_matrix(.shared.rep)

  if (is.na(.x.rep) && is.na(.y.rep)) {
    ps <- list()
    for (i in 1:ncol(mat)) {
      for (j in 1:ncol(mat)) {
        ps <- c(ps, list(vis_public_clonotypes(.shared.rep, i, j, "",
          .point.size.modif = .point.size.modif,
          .cut.axes = .cut.axes, .density = .density, .lm = .lm, .radj.size = .radj.size
        )))
      }
    }
    do.call(wrap_plots, c(ps, ncol = .ncol)) + plot_annotation(title = .title)
  } else if (is.na(.x.rep)) {
    ps <- lapply(1:ncol(mat), function(i) {
      vis_public_clonotypes(.shared.rep, i, .y.rep, "",
        .point.size.modif = .point.size.modif,
        .cut.axes = .cut.axes, .density = .density, .lm = .lm
      )
    })
    do.call(wrap_plots, c(ps, ncol = .ncol)) + plot_annotation(title = .title)
  } else if (is.na(.y.rep)) {
    ps <- lapply(1:ncol(mat), function(j) {
      vis_public_clonotypes(.shared.rep, .x.rep, j, "",
        .point.size.modif = .point.size.modif,
        .cut.axes = .cut.axes, .density = .density, .lm = .lm
      )
    })
    do.call(wrap_plots, c(ps, ncol = .ncol)) + plot_annotation(title = .title)
  } else {
    if (!is.character(.x.rep)) {
      .x.rep <- colnames(mat)[.x.rep]
    }
    if (!is.character(.y.rep)) {
      .y.rep <- colnames(mat)[.y.rep]
    }

    if (.x.rep == .y.rep) {
      return(rectGrob(gp = gpar(col = "white")))
    }

    df <- data.frame(cbind(mat[, .x.rep], mat[, .y.rep]))
    df[, 1] <- df[, 1] / sum(df[, 1], na.rm = TRUE)
    df[, 2] <- df[, 2] / sum(df[, 2], na.rm = TRUE)
    df_full <- df
    df <- df[!is.na(df[, 1]) & !is.na(df[, 2]), ]
    freq <- log10(sqrt(as.numeric(df[, 1]) * df[, 2])) / 2
    names(df) <- c("Xrep", "Yrep")
    names(df_full) <- c("Xrep", "Yrep")

    pnt.cols <- log(df[, 1] / df[, 2])
    suppressWarnings(pnt.cols[pnt.cols > 0] <- pnt.cols[pnt.cols > 0] / max(pnt.cols[pnt.cols > 0]))
    suppressWarnings(pnt.cols[pnt.cols < 0] <- -pnt.cols[pnt.cols < 0] / min(pnt.cols[pnt.cols < 0]))

    if (.cut.axes) {
      mat.lims <- c(min(as.matrix(df_full), na.rm = TRUE), max(as.matrix(df_full), na.rm = TRUE))
    } else {
      mat.lims <- c(min(as.matrix(df_full), na.rm = TRUE), 1)
    }

    empty <- ggplot() +
      geom_point(aes(1, 1), colour = "white") +
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      )

    min_df <- min(floor(log10(min(df_full[, 1], na.rm = TRUE))), floor(log10(min(df_full[, 2], na.rm = TRUE))))
    max_df <- max(trunc(log10(max(df_full[, 1], na.rm = TRUE))), trunc(log10(max(df_full[, 2], na.rm = TRUE))))
    breaks_values <- 10^seq(min_df, 1)
    breaks_labels <- format(log10(breaks_values), scientific = FALSE)

    grey_col <- "#CCCCCC"

    points <- ggplot() +
      geom_point(aes(x = Xrep, y = Yrep, size = freq, fill = pnt.cols), data = df, shape = 21) +
      scale_radius(range = c(.point.size.modif, .point.size.modif * 6)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      theme_linedraw() +
      scale_fill_distiller(palette = "RdBu") +
      scale_x_log10(breaks = breaks_values, labels = breaks_labels, lim = mat.lims, expand = c(.015, .015)) +
      scale_y_log10(breaks = breaks_values, labels = breaks_labels, lim = mat.lims, expand = c(.015, .015)) +
      theme(legend.position = "none") +
      xlab(.x.rep) +
      ylab(.y.rep)
    if (!is.na(.title)) {
      points <- points + ggtitle(.title)
    }
    if (.lm) {
      adj.R.sq <- summary(lm(Yrep ~ Xrep, df))$adj.

      points <- points +
        geom_smooth(aes(x = Xrep, y = Yrep), method = "lm", data = df, fullrange = TRUE, colour = "grey20", size = .5) +
        geom_text(aes(
          x = max(df_full, na.rm = TRUE) / 4,
          y = min(df_full, na.rm = TRUE),
          label = paste0("R^2(adj.) = ", as.character(round(adj.R.sq, 2)))
        ), size = .radj.size)
      # ggtitle(paste0("R^2(adj.) = ", as.character(round(adj.R.sq, 2))))
    }

    if (.density) {
      df2 <- data.frame(Clonotype = df_full[!is.na(df_full[, 1]) & is.na(df_full[, 2]), 1], Type = "unique", stringsAsFactors = FALSE)
      df2 <- rbind(df2, data.frame(Clonotype = df_full[!is.na(df_full[, 1]) & !is.na(df_full[, 2]), 1], Type = "public", stringsAsFactors = FALSE))
      top_plot <- ggplot() +
        geom_density(aes(x = Clonotype, fill = Type), colour = "grey25", data = df2, alpha = .3) +
        scale_x_log10(breaks = 10^(seq(min_df, 0)), lim = mat.lims, expand = c(.12, .015)) +
        theme_bw() +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"
        ) +
        scale_fill_manual(values = colorRampPalette(c(.colourblind_vector()[5], grey_col))(2))

      df2 <- data.frame(Clonotype = df_full[!is.na(df_full[, 2]) & is.na(df_full[, 1]), 2], Type = "unique", stringsAsFactors = FALSE)
      df2 <- rbind(df2, data.frame(Clonotype = df_full[!is.na(df_full[, 2]) & !is.na(df_full[, 1]), 2], Type = "public", stringsAsFactors = FALSE))
      right_plot <- ggplot() +
        geom_density(aes(x = Clonotype, fill = Type), colour = "grey25", data = df2, alpha = .3) +
        scale_x_log10(breaks = 10^(seq(min_df, 0)), lim = mat.lims, expand = c(.12, .015)) +
        coord_flip() +
        theme_bw() +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"
        ) +
        scale_fill_manual(values = colorRampPalette(c(.colourblind_vector()[1], grey_col))(2))

      wrap_plots(top_plot, empty, points, right_plot, ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
    } else {
      points
    }
  }
}


##### Gene usage & histogram plot, boxplot #####


#' Histograms and boxplots (general case / gene usage)
#'
#' @concept gene_usage
#'
#' @name vis.immunr_gene_usage
#'
#' @description Visualise distributions of genes using heatmaps or other plots.
#'
#' @param .data Output from the \link{geneUsage} function.
#'
#' @param .plot String specifying the plot type:
#'
#' - "hist" for histograms using \link{vis_hist};
#'
#' - "heatmap" for heatmaps using \link{vis_heatmap};
#'
#' - "heatmap2" for heatmaps using \link{vis_heatmap2};
#'
#' - "circos" for circos plots using \link{vis_circos}.
#'
#' @param ... Other arguments passed to corresponding functions depending on the plot type:
#'
#' - "hist" - passes arguments to \link{vis_hist};
#'
#' - "box" - passes arguments to \link{vis_box};
#'
#' - "heatmap" - passes arguments to \link{vis_heatmap};
#'
#' - "heatmap2" - passes arguments to \link{vis_heatmap2} and \link{heatmap} from the "pheatmap" package;
#'
#' - "circos" - passes arguments to \link{vis_circos} and \link{chordDiagram} from the "circlize" package.
#'
#' @return
#' A ggplot2 object, pheatmap or circlize object.
#'
#' @examples
#' data(immdata)
#'
#' gu <- geneUsage(immdata$data[[1]])
#' vis(gu)
#'
#' gu <- geneUsage(immdata$data)
#' vis(gu, .by = "Status", .meta = immdata$meta)
#' vis(gu, "box", .by = "Status", .meta = immdata$meta)
#' @seealso \link{geneUsage}
#'
#' @export
vis.immunr_gene_usage <- function(.data, .plot = c("hist", "box", "heatmap", "heatmap2", "circos"), ...) {
  .plot <- .plot[1]
  if (.plot == "hist") {
    vis_hist(.data, ...)
  } else if (.plot == "box") {
    vis_box(.data,
      .melt = TRUE, .grouping.var = "Gene",
      .labs = c("Gene", "Usage"), .title = "Gene usage",
      .subtitle = NULL, ...
    )
  } else if (.plot == "heatmap") {
    row.names(.data) <- .data[[1]]
    .data <- t(as.matrix(.data[2:ncol(.data)]))
    vis_heatmap(.data, ...)
  } else if (.plot == "heatmap2") {
    row.names(.data) <- .data[[1]]
    .data <- t(as.matrix(.data[2:ncol(.data)]))
    vis_heatmap2(.data, ...)
  } else if (.plot == "circos") {
    row.names(.data) <- .data[[1]]
    .data <- t(as.matrix(.data[2:ncol(.data)]))
    vis_circos(.data, ...)
  } else {
    stop("Error: Unknown value of the .plot parameter. Please provide one of the following: 'hist', 'box', 'heatmap', 'heatmap2', 'circos'.")
  }
}


#' Visualisation of distributions using histograms
#'
#' @concept vis
#'
#' @name vis_hist
#'
#' @description Visualisation of distributions using ggplot2-based histograms.
#'
#' @param .data Input matrix or data frame.
#'
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as age, serostatus or hla.
#'
#' @param .title The text for the title of the plot.
#'
#' @param .ncol A number of columns to display. Provide NA (by default) if you want the function
#' to automatically detect the optimal number of columns.
#'
#' @param .points A logical value defining whether points will be visualised or not.
#'
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#'
#' @param .coord.flip If TRUE then swap x- and y-axes.
#'
#' @param .grid If TRUE then plot separate visualisations for each sample.
#'
#' @param .labs A character vector of length two with names for x-axis and y-axis, respectively.
#'
#' @param .melt If TRUE then apply \link{melt} to the ".data" before plotting.
#' In this case ".data" is supposed to be a data frame with the first character column reserved
#' for names of genes and other numeric columns reserved to counts or frequencies of genes.
#' Each numeric column should be associated with a specific repertoire sample.
#'
#' @param .legend If TRUE then plots the legend. If FALSE removes the legend from the plot.
#' If NA automatically detects the best way to display legend.
#'
#' @param .add.layer Addditional ggplot2 layers, that added to each plot in the output plot or grid of plots.
#'
#' @param ... Is not used here.
#'
#' @details
#' If data is grouped, then statistical tests for comparing means of groups will be performed, unless \code{.test = FALSE} is supplied.
#' In case there are only two groups, the Wilcoxon rank sum test (https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) is performed
#' (R function \code{\link{wilcox.test}} with an argument \code{exact = FALSE}) for testing if there is a difference in mean rank values between two groups.
#' In case there more than two groups, the Kruskal-Wallis test (https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) is performed (R function \code{\link{kruskal.test}}), that is equivalent to ANOVA for ranks and it tests whether samples from different groups originated from the same distribution.
#' A significant Kruskal-Wallis test indicates that at least one sample stochastically dominates one other sample.
#' Adjusted for multiple comparisons P-values are plotted on the top of groups.
#' P-value adjusting is done using the Holm method (https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method) (also known as Holm-Bonferroni correction).
#' You can execute the command \code{?p.adjust} in the R console to see more.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \link{vis.immunr_gene_usage}, \link{geneUsage}
#'
#' @examples
#' data(immdata)
#' imm_gu <- geneUsage(immdata$data[[1]])
#' vis(imm_gu,
#'   .plot = "hist", .add.layer =
#'     theme(axis.text.x = element_text(angle = 75, vjust = 1))
#' )
#' imm_gu <- geneUsage(immdata$data[1:4])
#' vis(imm_gu,
#'   .plot = "hist", .grid = TRUE, .add.layer =
#'     theme(axis.text.x = element_text(angle = 75, vjust = 1))
#' )
#' @export
vis_hist <- function(.data, .by = NA, .meta = NA, .title = "Gene usage", .ncol = NA,
                     .points = TRUE, .test = TRUE, .coord.flip = FALSE,
                     .grid = FALSE, .labs = c("Gene", NA),
                     .melt = TRUE, .legend = NA, .add.layer = NULL, ...) {
  res <- .data

  if (.melt) {
    res <- reshape2::melt(res)
    res <- res[1:nrow(res), ]
    if (ncol(.data) == 2) {
      res[[2]] <- "Data"
    }
    colnames(res) <- c("Gene", "Sample", "Freq")
  }

  if (is.na(.labs[2])) {
    .labs[2] <- "Frequency"
    if (sum(res$Freq > 1, na.rm = TRUE) > 0) {
      .labs[2] <- "Count"
    }
  }

  if (length(unique(res$Sample)) == 1) {
    p <- ggplot() +
      geom_bar(aes(x = Gene, y = Freq, fill = Freq), data = res, stat = "identity", colour = "black")

    expand_vec <- c(.02, 0)

    p <- p + theme_pubr() + rotate_x_text() + theme_cleveland2() +
      labs(x = .labs[1], y = .labs[2], title = .title, fill = .labs[2]) +
      scale_fill_distiller(palette = "RdBu")

    if (.coord.flip) {
      p <- p + coord_flip() + scale_y_continuous(expand = c(.005, 0))
    } else {
      p <- p + scale_y_continuous(expand = c(.02, 0))
    }

    if (!is.na(.legend)) {
      if (!.legend) {
        p <- p + guides(fill = "none")
      }
    } else {
      p <- p + guides(fill = "none")
    }

    return(p + .add.layer)
  } else {
    if (.grid) {
      if (is.na(.legend)) {
        .legend <- FALSE
      }

      if (is.na(.ncol)) {
        .ncol <- round(length(unique(res$Sample))^.5)
      }

      if (!is.na(.by[1])) {
        group_res <- process_metadata_arguments(res, .by, .meta, "Sample")
        group_column <- group_res$name
        res$Group <- group_res$group_column

        res <- split(res, res$Group)
        ps <- list()

        for (i in seq_along(res)) {
          res[[i]]$Value <- res[[i]]$Freq

          ps[[i]] <- vis_bar(
            .data = res[[i]], .by = .by, .meta = .meta,
            .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .stack = FALSE,
            .points = .points, .test = .test, .signif.label.size = 3.5, .errorbar.width = 0.45,
            .defgroupby = "Sample", .grouping.var = "Gene",
            .labs = .labs, .title = names(res)[i], .subtitle = NULL,
            .legend = FALSE, .leg.title = NA
          ) +
            .add.layer
        }

        p <- do.call(wrap_plots, c(ps, ncol = .ncol))
        p <- p + plot_annotation(title = .title)
        return(p)
      } else {
        res <- split(res, res$Sample)
        ps <- list()
        for (i in seq_along(res)) {
          ps[[i]] <- vis_hist(res[[i]],
            .title = names(res)[i], .ncol = NA, .coord.flip = .coord.flip,
            .grid = FALSE, .labs = c("Gene", NA),
            .melt = FALSE, .legend = .legend, .add.layer = .add.layer, ...
          )
        }

        p <- do.call(wrap_plots, c(ps, ncol = .ncol))
        p <- p + plot_annotation(title = .title)
        return(p)
      }
    } else {
      res$Value <- res$Freq

      p <- vis_bar(
        .data = res, .by = .by, .meta = .meta,
        .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .stack = FALSE,
        .points = .points, .test = .test, .signif.label.size = 3.5, .errorbar.width = 0.45,
        .defgroupby = "Sample", .grouping.var = "Gene",
        .labs = .labs,
        .title = .title, .subtitle = NULL,
        .legend = TRUE, .leg.title = NA
      )

      return(p + .add.layer)
    }
  }
}


#' Flexible box-plots for visualisation of distributions
#'
#' @concept vis
#'
#' @name vis_box
#'
#' @description Visualisation of distributions using ggplot2-based boxplots.
#'
#' @param .data Input matrix or data frame.
#'
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as age, serostatus or hla.
#' @param .title The text for the title of the plot.
#' @param .labs Character vector of length two with names for x-axis and y-axis, respectively.
#' @param .melt If TRUE then apply \link{melt} to the ".data" before plotting.
#' In this case ".data" is supposed to be a data frame with the first character column reserved
#' for names of genes and other numeric columns reserved to counts or frequencies of genes.
#' Each numeric column should be associated with a specific repertoire sample.
#'
#' @param .points A logical value defining whether points will be visualised or not.
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#' @param .signif.label.size An integer value defining the size of text for p-value.
#' @param .defgroupby A name for the column with sample names.
#' @param .grouping.var A name for the column to group by.
#' @param .subtitle The The text for the plot's subtitle.
#' @param .leg.title The The text for the plots's legend. Provide NULL to remove the legend's title completely.
#' @param .legend If TRUE then displays a legend, otherwise removes legend from the plot.
#' @param .legend.pos Positions of the legend: either "top", "bottom", "left" or "right".
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \link{vis.immunr_gene_usage}, \link{geneUsage}
#'
#' @examples
#' vis_box(data.frame(Sample = sample(c("A", "B", "C"), 100, TRUE), Value = rnorm(100)), .melt = FALSE)
#' @export
vis_box <- function(.data, .by = NA, .meta = NA, .melt = TRUE,
                    .points = TRUE, .test = TRUE, .signif.label.size = 3.5, .defgroupby = "Sample", .grouping.var = "Group",
                    .labs = c("X", "Y"), .title = "Boxplot (.title argument)",
                    .subtitle = "Subtitle (.subtitle argument)",
                    .legend = NA, .leg.title = "Legend (.leg.title argument)", .legend.pos = "right") {
  if (.melt) {
    res <- reshape2::melt(.data)
    res <- res[1:nrow(res), ]
    if (ncol(.data) == 2) {
      res[[2]] <- "Data"
    }
    colnames(res) <- c(.grouping.var, "Sample", "Value")
    .data <- res
  }

  group_res <- process_metadata_arguments(.data, .by, .meta, .defgroupby)
  group_column <- group_res$name
  .data$Group <- group_res$group_column

  if (!(.grouping.var %in% names(.data))) {
    stop("Error: no column '", .grouping.var, "' in the input data frame.")
  }

  .data$Grouping.var <- .data[[.grouping.var]]

  .data_proc <- .data %>%
    dplyr::group_by(Grouping.var, Group)

  if (is.na(.labs[1])) {
    if (group_res$is_grouped) {
      .labs[1] <- "Group"
    } else {
      .labs[1] <- "Sample"
    }
  }

  if (is.na(.leg.title)) {
    if (group_res$is_grouped) {
      .leg.title <- "Group"
    } else {
      .leg.title <- "Sample"
    }
  }

  if (group_res$is_grouped) {
    # p = ggplot(aes(x = Grouping.var, y = Value, fill = Group, color=Group, group = Group), data = .data) +
    #     geom_boxplot(position = "dodge", col = "black")
    p <- ggplot(aes(x = Grouping.var, y = Value, fill = Group), data = .data) +
      geom_boxplot(position = "dodge", col = "black")

    if (is.na(.legend)) {
      if (.grouping.var == "Group") {
        # p = p + guides(fill=FALSE)
        .legend.pos <- "none"
      }
    } else if (.legend == FALSE) {
      # p = p + guides(fill=FALSE)
      .legend.pos <- "none"
    }

    if (.points) {
      p <- p +
        geom_point(color = "black", position = position_jitterdodge(0.05), size = 1)
    }

    if (.test) {
      if (.grouping.var == "Group") {
        comparisons <- list()
        for (i in 1:(nrow(.data_proc) - 1)) {
          for (j in (i + 1):nrow(.data_proc)) {
            # if ((.data$Grouping.var[i] == .data$Grouping.var[i]) && (.data$Group[i] != .data$Group[j])) {
            comparisons <- c(comparisons, list(c(i, j)))
            # }
          }
        }

        p_df <- compare_means(Value ~ Group, .data, comparisons = comparisons, p.adjust.method = "holm")

        y_max <- max(.data$Value)
        p.value.y.coord <- rep(y_max, nrow(p_df))
        step.increase <- (1:nrow(p_df)) * (y_max / 10)
        p.value.y.coord <- p.value.y.coord + step.increase

        p_df <- p_df %>%
          mutate(
            y.coord = p.value.y.coord,
            p.adj = format.pval(p.adj, digits = 1)
          )

        p <- p + geom_signif(
          data = p_df,
          aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord),
          manual = TRUE, tip_length = 0.03, size = .5, inherit.aes = FALSE
        )
      } else {
        # Seems fine...
        # p_df = compare_means(Value ~ Group, group.by = "Grouping.var", method = "kruskal.test", .data, p.adjust.method = "holm")
        # print(p_df)

        p <- p +
          stat_compare_means(aes(label = ..p.adj..),
            bracket.size = .5, size = .signif.label.size,
            label.y = max(.data$Value, na.rm = TRUE) * 1.07
          )
      }
    }

    p <- p + .tweak_fill(length(unique(.data$Group)))
  } else {
    p <- ggplot() +
      geom_boxplot(aes(x = Grouping.var, y = Value, fill = Sample), position = "dodge", data = .data, col = "black")

    if (is.na(.legend)) {
      # p = p + guides(fill=FALSE)
      .legend.pos <- "none"
    }

    p <- p + .colourblind_discrete(length(unique(.data$Sample)))
  }

  p + theme_pubr(legend = .legend.pos) +
    labs(x = .labs[1], y = .labs[2], title = .title, subtitle = .subtitle, fill = group_column) +
    rotate_x_text() + theme_cleveland2()
}


##### Clustering #####


#' Visualisation of hierarchical clustering
#'
#' @concept post_analysis
#'
#' @description
#' Visualisation of the results of hierarchical clustering.
#' For other clustering visualisations see \link{vis.immunr_kmeans}.
#'
#' @aliases vis.immunr_hclust
#'
#' @param .data Clustering results from \link{repOverlapAnalysis} or \link{geneUsageAnalysis}.
#' @param .rect Passed to \link{fviz_dend} - whether to add a rectangle around groups.
#' @param .plot A character vector of length one or two specifying which plots to visualise.
#' If "clust" then plot only the clustering. If "best" then plot the number of optimal clusters.
#' If both then plot both.
#' @param ... Not used here.
#'
#' @return
#' Ggplot2 objects inside the patchwork container.
#'
#' @seealso \link{vis}, \link{repOverlapAnalysis}, \link{geneUsageAnalysis}
#'
#' @examples
#' data(immdata)
#' ov <- repOverlap(immdata$data)
#' repOverlapAnalysis(ov, "mds+hclust") %>% vis()
#' @export
vis.immunr_hclust <- function(.data, .rect = FALSE, .plot = c("clust", "best"), ...) {
  p1 <- NULL
  if ("clust" %in% .plot) {
    p1 <- fviz_dend(.data[[1]], main = "Hierarchical clustering", rect = .rect)
  }

  p2 <- NULL
  if ("best" %in% .plot) {
    p2 <- .data[[2]]
  }

  if (is.null(p1) | is.null(p2)) {
    if (is.null(p1)) {
      p2
    } else {
      p1
    }
  } else {
    p <- p1 + p2
    p
  }
}


#' Visualisation of K-means and DBSCAN clustering
#'
#' @concept post_analysis
#'
#' @description
#' Visualisation of the results of K-means and DBSCAN clustering.
#' For hierarhical clustering visualisations see \link{vis.immunr_hclust}.
#'
#' @aliases vis.immunr_kmeans vis.immunr_dbscan
#'
#' @param .data Clustering results from \link{repOverlapAnalysis} or \link{geneUsageAnalysis}.
#' @param .point If TRUE then plot sample points. Passed to \link{fviz_cluster}.
#' @param .text If TRUE then plot text labels. Passed to \link{fviz_cluster}.
#' @param .ellipse If TRUE then plot ellipses around all samples. Passed to "ellipse" from \link{fviz_cluster}.
#' @param .point.size Size of points, passed to "pointsize" from \link{fviz_cluster}.
#' @param .text.size Size of text labels, passed to labelsize from \link{fviz_cluster}.
#' @param .plot A character vector of length one or two specifying which plots to visualise.
#' If "clust" then plot only the clustering. If "best" then plot the number of optimal clusters.
#' If both then plot both.
#' @param ... Not used here.
#'
#' @return
#' Ggplot2 objects inside the pathwork container.
#'
#' @seealso \link{vis}, \link{repOverlapAnalysis}, \link{geneUsageAnalysis}
#'
#' @examples
#' data(immdata)
#' ov <- repOverlap(immdata$data)
#' repOverlapAnalysis(ov, "mds+kmeans") %>% vis()
#' @export
vis.immunr_kmeans <- function(.data, .point = TRUE, .text = TRUE, .ellipse = TRUE,
                              .point.size = 2, .text.size = 10, .plot = c("clust", "best"),
                              ...) {
  p1 <- NULL
  if ("clust" %in% .plot) {
    p1 <- fviz_cluster(.data[[1]],
      data = .data[[3]], main = "K-means clustering", geom = c("point", "text")[c(.point, .text)],
      show.legend.text = FALSE, show.clust.cent = FALSE, repel = TRUE, ellipse = .ellipse, shape = 16,
      pointsize = .point.size, labelsize = .text.size, label.rectangle = TRUE
    ) +
      guides(fill = "none", shape = "none") +
      theme_linedraw()
  }

  p2 <- NULL
  if ("best" %in% .plot) {
    p2 <- .data[[2]]
  }

  if (is.null(p1) | is.null(p2)) {
    if (is.null(p1)) {
      p2
    } else {
      p1
    }
    # p = ifelse(is.null(p1), p2, p1) # doesnt work for some reason
  } else {
    p1 + p2
  }
}

#' @export
vis.immunr_dbscan <- function(.data, .point = TRUE, .text = TRUE, .ellipse = TRUE,
                              .point.size = 2, .text.size = 10, .plot = c("clust", "best"), ...) {
  fviz_cluster(.data[[1]],
    data = .data[[2]], main = "DBSCAN clustering", geom = c("point", "text")[c(.point, .text)],
    show.legend.text = FALSE, show.clust.cent = FALSE, repel = TRUE, ellipse = .ellipse, shape = 16,
    pointsize = .point.size, labelsize = .text.size, label.rectangle = TRUE
  ) +
    guides(fill = "none", shape = "none") +
    theme_linedraw()
}


##### Dimension reduction #####


#' PCA / MDS / tSNE visualisation (mainly overlap / gene usage)
#'
#' @concept post_analysis
#'
#' @importFrom ggpubr ggscatter
#'
#' @aliases vis.immunr_mds vis.immunr_pca vis.immunr_tsne
#'
#' @param .data Output from analysis functions such as \link{geneUsageAnalysis} or
#' \link{immunr_pca}, \link{immunr_mds} or \link{immunr_tsne}.
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as age, serostatus or hla.
#' @param .point Logical. If TRUE then plot points corresponding to objects.
#' @param .text Logical. If TRUE then plot sample names.
#' @param .ellipse Logical. If TRUE then plot ellipses around clusters of grouped samples.
#' @param .point.size Numeric. A size of points to plot.
#' @param .text.size Numeric. A size of sample names' labels.
#' @param ... Not used here.
#'
#' @return
#' A ggplot2 object.
#'
#' @details
#' Other visualisation methods:
#'
#' - PCA - \link{vis.immunr_pca}
#'
#' - MDS - \link{vis.immunr_mds}
#'
#' - tSNE - \link{vis.immunr_tsne}
#'
#' @examples
#' data(immdata)
#' ov <- repOverlap(immdata$data)
#' repOverlapAnalysis(ov, "mds") %>% vis()
#' @export
vis.immunr_mds <- function(.data, .by = NA, .meta = NA,
                           .point = TRUE, .text = TRUE, .ellipse = TRUE,
                           .point.size = 2, .text.size = 4, ...) {
  if (!.point & !.text) {
    stop("Error: Please provide at least one of the arguments: .point and .text")
  }
  group_res <- process_metadata_arguments(data.frame(Sample = row.names(.data$x), stringsAsFactors = FALSE), .by, .meta)
  group_column <- group_res$name
  if (!group_res$is_grouped) {
    .ellipse <- FALSE
  }

  fviz_pca_ind(.data,
    habillage = group_res$group_column, geom = c("point", "text")[c(.point, .text)],
    repel = TRUE, addEllipses = .ellipse, mean.point = FALSE, pointshape = 16,
    pointsize = .point.size, labelsize = .text.size, label.rectangle = TRUE, show.legend.text = FALSE
  ) +
    theme_linedraw() +
    guides(fill = "none", shape = "none") +
    labs(x = "Dim1", y = "Dim2", title = "Multidimensional scaling", col = group_column)
}

#' @export
vis.immunr_pca <- function(.data, .by = NA, .meta = NA,
                           .point = TRUE, .text = TRUE, .ellipse = TRUE,
                           .point.size = 2, .text.size = 4, ...) {
  if (!.point & !.text) {
    stop("Error: Please provide at least one of the arguments: .point and .text")
  }
  group_res <- process_metadata_arguments(data.frame(Sample = row.names(.data$x), stringsAsFactors = FALSE), .by, .meta)
  group_column <- group_res$name
  if (!group_res$is_grouped) {
    .ellipse <- FALSE
  }

  fviz_pca_ind(.data,
    habillage = group_res$group_column, geom = c("point", "text")[c(.point, .text)],
    repel = TRUE, addEllipses = .ellipse, mean.point = FALSE, pointshape = 16,
    pointsize = .point.size, labelsize = .text.size, label.rectangle = TRUE, show.legend.text = FALSE
  ) +
    theme_linedraw() +
    guides(fill = "none", shape = "none") +
    labs(title = "Principal component analysis", col = group_column)
}

#' @export
vis.immunr_tsne <- function(.data, .by = NA, .meta = NA,
                            .point = TRUE, .text = TRUE, .ellipse = TRUE,
                            .point.size = 2, .text.size = 4, ...) {
  .data <- data.frame(.data)
  colnames(.data) <- c("Dim1", "Dim2")
  .data$Sample <- row.names(.data)
  group_res <- process_metadata_arguments(.data, .by, .meta)
  group_column <- group_res$name
  .data$Group <- group_res$group_column

  if (!group_res$is_grouped) {
    .ellipse <- FALSE
  }

  ggscatter(
    data = .data, x = "Dim1", y = "Dim2", color = "Group", ellipse = .ellipse, size = .point.size,
    point = .point, label = ifelse(.text, "Sample", NULL), repel = TRUE, label.rectangle = TRUE, show.legend.text = FALSE
  ) +
    theme_linedraw() +
    guides(fill = "none", shape = "none") +
    labs(title = "t-Distributed Stochastic Neighbor Embedding", col = group_column)
}

#
# ToDo: ideally, geneUsageAnalysis returns raw objects, and vis.immunr_mds transforms it into the fviz-compatible objects
#


##### Clonality analysis #####


vis_bar_stacked <- function(.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .stack = NA,
                            .grouping.var = NA,
                            .labs = c(NA, "Y"),
                            .title = "Barplot (.title argument)", .subtitle = "Subtitle (.subtitle argument)",
                            .legend = NA, .leg.title = NA) {
  group_res <- process_metadata_arguments(.data, .by, .meta)
  .data$Group <- group_res$group_column

  if (!(.grouping.var %in% names(.data))) {
    stop("Error: no column '", .grouping.var, "' in the input data frame.")
  }
  .data$Grouping.var <- .data[[.grouping.var]]

  position_bar <- "dodge"
  if (is.na(.stack)) {
    if (group_res$is_grouped) {
      if (.errorbars.off) {
        position_bar <- "stack"
      } else {
        position_bar <- "dodge"
      }
    } else {
      position_bar <- "stack"
    }
  } else if (.stack) {
    if (group_res$is_grouped) {
      if (.errorbars.off) {
        position_bar <- "stack"
      } else {
        warning("Warning: Provide .errorbars.off = TRUE in order to display stacked barplots.")
        position_bar <- "dodge"
      }
    } else {
      position_bar <- "stack"
    }
  }

  if (is.na(.labs[1])) {
    if (group_res$is_grouped) {
      .labs[1] <- "Group"
    } else {
      .labs[1] <- "Sample"
    }
  }

  if (group_res$is_grouped) {
    perc <- .data %>%
      dplyr::group_by(Grouping.var, Group) %>%
      dplyr::summarise(
        Value.mean = mean(Value, na.rm = TRUE),
        Value.min = quantile(Value, .errorbars[1], na.rm = TRUE),
        Value.max = quantile(Value, .errorbars[2], na.rm = TRUE)
      )

    p <- ggplot() +
      geom_bar(aes(x = Group, y = Value.mean, fill = Grouping.var), data = perc, colour = "black", stat = "identity", position = position_bar)

    if (!.errorbars.off) {
      p <- p + geom_errorbar(aes(x = Group, fill = Grouping.var, ymin = Value.min, ymax = Value.max),
        data = perc, colour = "black", width = .2, position = position_dodge(.9)
      )
    }
  } else {
    p <- ggplot() +
      geom_bar(aes(x = Sample, y = Value, fill = Grouping.var), data = .data, colour = "black", stat = "identity", position = position_bar)
  }

  p + theme_pubr(legend = "right") +
    labs(x = .labs[1], y = .labs[2], title = .title, subtitle = .subtitle, fill = .leg.title) +
    .colourblind_discrete(length(unique(.data$Grouping.var)))
}


#' Visualise results of the clonality analysis
#'
#' @concept clonality
#'
#' @name vis.immunr_clonal_prop
#'
#' @aliases vis.immunr_clonal_prop vis.immunr_homeo vis.immunr_top_prop vis.immunr_tail_prop
#'
#' @description An utility function to visualise the output from \code{\link{repClonality}}.
#'
#' @importFrom reshape2 melt
#' @importFrom scales percent
#'
#' @param .data Output from \code{\link{repClonality}}.
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as age, serostatus or hla.
#'
#' @param .errorbars A numeric vector of length two with quantiles for error bars
#' on sectors. Disabled if ".errorbars.off" is TRUE.
#'
#' @param .errorbars.off If TRUE then plot CI bars for distances between each group.
#' Disabled if no group passed to the ".by" argument.
#' @param .points A logical value defining whether points will be visualised or not.
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#' @param .signif.label.size An integer value defining the size of text for p-value.
#' @param ... Not used here.
#'
#' @details
#' If data is grouped, then statistical tests for comparing means of groups will be performed, unless \code{.test = FALSE} is supplied.
#' In case there are only two groups, the Wilcoxon rank sum test (https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) is performed
#' (R function \code{\link{wilcox.test}} with an argument \code{exact = FALSE}) for testing if there is a difference in mean rank values between two groups.
#' In case there more than two groups, the Kruskal-Wallis test (https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) is performed (R function \code{\link{kruskal.test}}), that is equivalent to ANOVA for ranks and it tests whether samples from different groups originated from the same distribution.
#' A significant Kruskal-Wallis test indicates that at least one sample stochastically dominates one other sample.
#' Adjusted for multiple comparisons P-values are plotted on the top of groups.
#' P-value adjusting is done using the Holm method (https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method) (also known as Holm-Bonferroni correction).
#' You can execute the command \code{?p.adjust} in the R console to see more.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \link{repClonality} \link{vis}
#'
#' @examples
#' data(immdata)
#' clp <- repClonality(immdata$data, "clonal.prop")
#' vis(clp)
#'
#' hom <- repClonality(immdata$data, "homeo")
#' # Remove p values and points from the plot
#' vis(hom, .by = "Status", .meta = immdata$meta, .test = FALSE, .points = FALSE)
#' @export
vis.immunr_clonal_prop <- function(.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .points = TRUE, .test = TRUE, .signif.label.size = 3.5, ...) {
  # ToDo: this and other repClonality and repDiversity functions doesn't work on a single repertoire. Fix it
  perc_value <- round(.data[1, 2][1])
  .data <- data.frame(Sample = row.names(.data), Value = .data[, 1])

  p <- vis_bar(
    .data = .data, .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
    .points = .points, .test = .test, .signif.label.size = .signif.label.size,
    .defgroupby = "Sample", .grouping.var = "Group",
    .labs = c(NA, "Clonal proportion"),
    .title = "Clonal proportion", .subtitle = paste0("Number of clonotypes occupying the ", perc_value, "% of repertoires"),
    .legend = NA, .leg.title = NA
  )

  p
}

#' @export
vis.immunr_homeo <- function(.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .stack = NA, .test = TRUE, .points = TRUE, ...) {
  melted <- reshape2::melt(.data)
  colnames(melted) <- c("Sample", "Clone.group", "Value")

  if (is.na(.by[1]) && is.na(.stack)) {
    .stack <- TRUE
  }

  p <- vis_bar(melted,
    .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = .stack,
    .grouping.var = "Clone.group", .points = .points,
    .labs = c(NA, "Relative abundance, perc"),
    .title = "Relative abundance", .subtitle = "Summary proportion of clonotypes with specific frequencies",
    .legend = NA, .leg.title = "Clonotype group", .test = .test
  ) +
    scale_y_continuous(labels = scales::percent)

  p
}

#' @export
vis.immunr_top_prop <- function(.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .stack = NA, .points = TRUE, .test = TRUE, ...) {
  tmp <- .data
  if (is.null(dim(tmp))) {
    tmp <- t(as.matrix(tmp))
    .data <- t(as.matrix(.data))
  }
  for (i in 2:ncol(.data)) {
    tmp[, i] <- .data[, i] - .data[, i - 1]
  }
  res <- tmp
  .head <- as.numeric(colnames(.data))
  colnames(res) <- paste0("[", c(1, .head[-length(.head)] + 1), ":", .head, ")")
  res <- as.data.frame(res)
  res$Sample <- row.names(res)
  res <- reshape2::melt(res)
  colnames(res) <- c("Sample", "Clone.index", "Value")

  if (is.na(.by[1]) && is.na(.stack)) {
    .stack <- TRUE
  }

  p <- vis_bar(res,
    .by = .by, .meta = .meta, .points = .points, .test = .test,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = .stack,
    .grouping.var = "Clone.index",
    .labs = c(NA, "Occupied repertoire space"),
    .title = "Top clonal proportion", .subtitle = "Summary proportion of clonotypes with specific indices",
    .legend = NA, .leg.title = "Clonotype indices"
  ) +
    scale_y_continuous(labels = scales::percent)

  p
}

#' @export
vis.immunr_rare_prop <- function(.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .stack = NA, .points = TRUE, .test = TRUE, ...) {
  tmp <- .data
  if (is.null(dim(tmp))) {
    tmp <- t(as.matrix(tmp))
    .data <- t(as.matrix(.data))
  }

  colnames_new <- c()
  prev_value <- "1"
  start_index <- 1
  if (colnames(.data)[1] == "1") {
    colnames_new <- "1"
    prev_value <- "2"
    start_index <- 2
  }
  for (i in start_index:ncol(.data)) {
    colnames_new <- c(colnames_new, paste0(prev_value, " - ", colnames(.data)[i]))
    if (colnames(.data)[i] != "MAX") {
      prev_value <- as.character(as.integer(colnames(.data)[i]) + 1)
    }
  }
  colnames(tmp) <- colnames_new

  for (i in 2:ncol(.data)) {
    tmp[, i] <- .data[, i] - .data[, i - 1]
  }
  res <- tmp

  res <- as.data.frame(res)
  res$Sample <- row.names(res)
  res <- reshape2::melt(res)
  colnames(res) <- c("Sample", "Counts", "Value")

  if (is.na(.by[1]) && is.na(.stack)) {
    .stack <- TRUE
  }

  p <- vis_bar(res,
    .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = .stack,
    .grouping.var = "Counts", .points = .points, .test = .test,
    .labs = c(NA, "Occupied repertoire space"),
    .title = "Rare clonal proportion", .subtitle = "Summary proportion of clonotypes with specific counts",
    .legend = NA, .leg.title = "Clonotype counts"
  ) +
    scale_y_continuous(labels = scales::percent)

  p
}


##### Diversity estimation & dodged bar plots #####


#' Bar plots
#'
#' @concept vis
#'
#' @importFrom ggpubr compare_means geom_signif stat_compare_means theme_pubr rotate_x_text
#'
#' @name vis_bar
#'
#' @param .data Data to visualise.
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as age, serostatus or hla.
#' @param .errorbars A numeric vector of length two with quantiles for error bars
#' on sectors. Disabled if ".errorbars.off" is TRUE.
#' @param .errorbars.off If TRUE then plot CI bars for distances between each group.
#' Disabled if no group passed to the ".by" argument.
#' @param .errorbar.width Numeric. Width for error bars.
#' @param .points A logical value defining whether points will be visualised or not.
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#' @param .signif.label.size An integer value defining the size of text for p-value.
#' @param .stack If TRUE and .errorbars.off is TRUE then plot stacked bar plots for each Group or Sample
#' @param .defgroupby A name for the column with sample names.
#' @param .grouping.var A name for the column to group by.
#' @param .subtitle The text for the plot's subtitle.
#' @param .leg.title The text for the plots's legend. Provide NULL to remove the legend's title completely.
#' @param .legend If TRUE then displays a legend, otherwise removes legend from the plot.
#' @param .labs A character vector of length two specifying names for x-axis and y-axis.
#' @param .title The text for the plot's title.
#' @param .legend.pos Positions of the legend: either "top", "bottom", "left" or "right".
#' @param .rotate_x How much the x tick text should be rotated? In angles.
#'
#' @return
#' A ggplot2 object.
#'
#' @examples
#' vis_bar(data.frame(Sample = c("A", "B", "C"), Value = c(1, 2, 3)))
#' @export
vis_bar <- function(.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .stack = FALSE,
                    .points = TRUE, .test = TRUE, .signif.label.size = 3.5, .errorbar.width = 0.2, .defgroupby = "Sample", .grouping.var = "Group",
                    .labs = c("X", "Y"), .title = "Barplot (.title argument)",
                    .subtitle = "Subtitle (.subtitle argument)",
                    .legend = NA, .leg.title = "Legend (.leg.title argument)", .legend.pos = "right",
                    .rotate_x = 90) {
  group_res <- process_metadata_arguments(.data, .by, .meta, .defgroupby)
  group_column <- group_res$name
  .data$Group <- group_res$group_column

  if (!(.grouping.var %in% names(.data))) {
    stop("Error: no column '", .grouping.var, "' in the input data frame.")
  }

  .data$Grouping.var <- .data[[.grouping.var]]

  .data_proc <- .data %>%
    dplyr::group_by(Grouping.var, Group) %>%
    dplyr::summarise(
      Value.mean = mean(Value, na.rm = TRUE),
      Value.min = quantile(Value, .errorbars[1], na.rm = TRUE),
      Value.max = quantile(Value, .errorbars[2], na.rm = TRUE)
    )

  if (is.na(.labs[1])) {
    if (group_res$is_grouped) {
      .labs[1] <- "Group"
    } else {
      .labs[1] <- "Sample"
    }
  }

  if (is.na(.leg.title)) {
    .leg.title <- group_column
    # if (group_res$is_grouped) {
    #   .leg.title = "Group"
    # } else {
    #   .leg.title = "Sample"
    # }
  } else if (group_res$is_grouped) {
    .leg.title <- group_column
  }

  position_name <- "dodge"
  if (group_res$is_grouped) {
    if (is.na(.stack)) {
      .stack <- FALSE
    }
  } else {
    if (is.na(.stack)) {
      .stack <- FALSE
    } else if (.stack) {
      .errorbars.off <- TRUE
      .points <- FALSE
      .test <- FALSE
    }
  }

  if (group_res$is_grouped) {
    if (.stack) {
      p <- ggplot() +
        geom_col(aes(x = Group, y = Value.mean, fill = Grouping.var),
          position = "stack", data = .data_proc, col = "black"
        )
    } else {
      p <- ggplot(aes(x = Grouping.var, y = Value, color = Group, fill = Group, group = Group), data = .data) +
        geom_col(aes(x = Grouping.var, y = Value.mean),
          position = position_name, data = .data_proc, col = "black"
        )
    }

    if (is.na(.legend)) {
      if (.grouping.var == "Group") {
        # p = p + guides(fill=FALSE)
        .legend.pos <- "none"
      }
    } else if (.legend == FALSE) {
      # p = p + guides(fill=FALSE)
      .legend.pos <- "none"
    }

    if (!.errorbars.off) {
      p <- p +
        geom_errorbar(aes(x = Grouping.var, y = Value.mean, ymin = Value.min, ymax = Value.max, color = Group),
          color = "black", data = .data_proc, width = .errorbar.width, position = position_dodge(.9)
        )
    }

    if (.points) {
      p <- p +
        geom_point(color = "black", position = position_jitterdodge(0.05), size = 1)
    }

    if (.test) {
      if (.grouping.var == "Group") {
        comparisons <- list()
        for (i in 1:(nrow(.data_proc) - 1)) {
          for (j in (i + 1):nrow(.data_proc)) {
            # if ((.data$Grouping.var[i] == .data$Grouping.var[i]) && (.data$Group[i] != .data$Group[j])) {
            comparisons <- c(comparisons, list(c(i, j)))
            # }
          }
        }

        p_df <- compare_means(Value ~ Group, .data, comparisons = comparisons, p.adjust.method = "holm")

        y_max <- max(.data$Value)
        p.value.y.coord <- rep(y_max, nrow(p_df))
        step.increase <- (1:nrow(p_df)) * (y_max / 10)
        p.value.y.coord <- p.value.y.coord + step.increase

        p_df <- p_df %>%
          mutate(
            y.coord = p.value.y.coord,
            p.adj = format.pval(p.adj, digits = 1)
          )

        p <- p + geom_signif(
          data = p_df,
          aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord),
          manual = TRUE, tip_length = 0.03, size = .5, inherit.aes = FALSE
        )
      } else {
        # Seems fine...
        # p_df = compare_means(Value ~ Group, group.by = "Grouping.var", method = "kruskal.test", .data, p.adjust.method = "holm")
        # print(p_df)

        p <- p +
          stat_compare_means(aes(label = ..p.adj..),
            bracket.size = .5, size = .signif.label.size,
            label.y = max(.data$Value, na.rm = TRUE) * 1.07
          )
      }
    }

    p <- p + .tweak_fill(length(unique(.data$Group)))
  } else {
    if (.stack) {
      p <- ggplot() +
        geom_bar(aes(x = Sample, y = Value, fill = Grouping.var), data = .data, colour = "black", stat = "identity", position = "stack")

      if (!is.na(.legend)) {
        if (!.legend) {
          # p = p + guides(fill=FALSE)
          .legend.pos <- "none"
        }
      }

      p <- p + .colourblind_discrete(length(unique(.data$Grouping.var)))
    } else {
      p <- ggplot() +
        geom_bar(aes(x = Grouping.var, y = Value, fill = Sample), position = position_name, stat = "identity", data = .data, col = "black")

      if (is.na(.legend)) {
        # p = p + guides(fill=FALSE)
        .legend.pos <- "none"
      }

      p <- p + .colourblind_discrete(length(unique(.data$Sample)))
    }
  }

  p <- p + theme_pubr(legend = .legend.pos) +
    labs(x = .labs[1], y = .labs[2], title = .title, subtitle = .subtitle, fill = .leg.title) +
    rotate_x_text() + theme_cleveland2()

  p
}


#' Visualise diversity.
#'
#' @concept diversity
#'
#' @aliases vis.immunr_chao1 vis.immunr_dxx vis.immunr_rarefaction vis.immunr_div vis.immunr_ginisimp vis.immunr_invsimp vis.immunr_hill
#' @description An utility function to visualise the output from \code{\link{repDiversity}}.
#'
#' @importFrom reshape2 melt
#'
#' @param .data Output from \code{\link{repDiversity}}.
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as age, serostatus or hla.
#' @param .errorbars A numeric vector of length two with quantiles for error bars
#' on sectors. Disabled if ".errorbars.off" is TRUE.
#' @param .errorbars.off If TRUE then plot CI bars for distances between each group.
#' Disabled if no group passed to the ".by" argument.
#' @param .points A logical value defining whether points will be visualised or not.
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#' @param .signif.label.size An integer value defining the size of text for p-value.
#' @param ... Not used here.
#'
#' @details
#' If data is grouped, then statistical tests for comparing means of groups will be performed, unless \code{.test = FALSE} is supplied.
#' In case there are only two groups, the Wilcoxon rank sum test (https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) is performed
#' (R function \code{\link{wilcox.test}} with an argument \code{exact = FALSE}) for testing if there is a difference in mean rank values between two groups.
#' In case there more than two groups, the Kruskal-Wallis test (https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) is performed (R function \code{\link{kruskal.test}}), that is equivalent to ANOVA for ranks and it tests whether samples from different groups originated from the same distribution.
#' A significant Kruskal-Wallis test indicates that at least one sample stochastically dominates one other sample.
#' Adjusted for multiple comparisons P-values are plotted on the top of groups.
#' P-value adjusting is done using the Holm method (https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method) (also known as Holm-Bonferroni correction).
#' You can execute the command \code{?p.adjust} in the R console to see more.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \link{repDiversity} \link{vis}
#'
#' @examples
#' data(immdata)
#' dv <- repDiversity(immdata$data, "chao1")
#' vis(dv)
#' @export
vis.immunr_chao1 <- function(.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = FALSE, .points = TRUE, .test = TRUE, .signif.label.size = 3.5, ...) {
  .data <- data.frame(Sample = row.names(.data), Value = .data[, 1])

  vis_bar(
    .data = .data, .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
    .points = .points, .test = .test, .signif.label.size = .signif.label.size,
    .defgroupby = "Sample", .grouping.var = "Group",
    .labs = c(NA, "Chao1"),
    .title = "Chao1", .subtitle = "Sample diversity estimation using Chao1",
    .legend = NA, .leg.title = NA
  )
}

#' @export
vis.immunr_hill <- function(.data, .by = NA, .meta = NA, .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                            .points = TRUE, .add.points = TRUE, .point.size = 1.5, .add.point.size = 1, .line.size = 0.75, .leg.title = NA, ...) {
  group_res <- process_metadata_arguments(.data, .by, .meta)
  .data$Group <- group_res$group_column

  if (is.na(.leg.title)) {
    if (group_res$is_grouped) {
      .leg.title <- "Group"
    } else {
      .leg.title <- "Sample"
    }
  }

  grouped_data <- .data %>%
    dplyr::group_by(Group, Q) %>%
    dplyr::summarise(
      Value.mean = mean(Value, na.rm = TRUE),
      Value.min = quantile(Value, .errorbars[1], na.rm = TRUE),
      Value.max = quantile(Value, .errorbars[2], na.rm = TRUE)
    )

  position_value <- "identity"
  if (group_res$is_grouped) {
    position_value <- position_jitter(width = 0.1)
  }

  p <- ggplot() +
    geom_line(aes(x = Q, y = Value.mean, col = Group), size = .line.size, data = grouped_data, stat = "identity")

  if (.points) {
    p <- p + geom_point(aes(x = Q, y = Value, col = Group),
      data = .data,
      position = position_value, alpha = .5, size = .point.size
    )
  }

  if (group_res$is_grouped) {
    if (.add.points) {
      p <- p +
        geom_point(aes(x = Q, y = Value.mean, col = Group), data = grouped_data, size = .add.point.size)
    }

    if (!.errorbars.off) {
      p <- p +
        geom_errorbar(aes(x = Q, ymin = Value.min, ymax = Value.max, group = Group, col = Group),
          data = grouped_data, width = .2, size = .line.size, position = position_dodge(0.1)
        )
    }
  }

  p <- p +
    theme_pubr(legend = "right") + theme_cleveland2() +
    labs(x = "Q values", y = "Diversity estimation", title = "Hill numbers", subtitle = "Sample diversity estimation using Hill numbers", color = .leg.title)

  p
}

#' @export
vis.immunr_div <- function(.data, .by = NA, .meta = NA,
                           .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                           .points = TRUE, .test = TRUE, .signif.label.size = 3.5,
                           .legend = NA, ...) {
  vis_bar(
    .data = .data, .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
    .points = .points, .test = .test, .signif.label.size = .signif.label.size,
    .defgroupby = "Sample", .grouping.var = "Group",
    .labs = c(NA, "Effective number of clonotypes"),
    .title = "True diversity", .subtitle = "Sample diversity estimation using the true diversity index",
    .legend = .legend, .leg.title = NA
  )
}

#' @export
vis.immunr_ginisimp <- function(.data, .by = NA, .meta = NA,
                                .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                                .points = TRUE, .test = TRUE, .signif.label.size = 3.5, ...) {
  vis_bar(
    .data = .data, .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
    .points = .points, .test = .test, .signif.label.size = .signif.label.size,
    .defgroupby = "Sample", .grouping.var = "Group",
    .labs = c(NA, "Gini-Simpson index"),
    .title = "Gini-Simpson index", .subtitle = "Sample diversity estimation using the Gini-Simpson index",
    .legend = NA, .leg.title = NA
  )
}

#' @export
vis.immunr_invsimp <- function(.data, .by = NA, .meta = NA,
                               .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                               .points = TRUE, .test = TRUE, .signif.label.size = 3.5,
                               .legend = NA, ...) {
  vis_bar(
    .data = .data, .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
    .points = .points, .test = .test, .signif.label.size = .signif.label.size,
    .defgroupby = "Sample", .grouping.var = "Group",
    .labs = c(NA, "Inverse Simpson index"),
    .title = "Inverse Simpson index", .subtitle = "Sample diversity estimation using the inverse Simpson index",
    .legend = .legend, .leg.title = NA
  )
}

#' @export
vis.immunr_dxx <- function(.data, .by = NA, .meta = NA,
                           .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                           .points = TRUE, .test = TRUE, .signif.label.size = 3.5,
                           .legend = NA, ...) {
  perc_value <- round(.data[1, 2][1])
  .data <- data.frame(Sample = row.names(.data), Value = .data[, 1])

  vis_bar(
    .data = .data, .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
    .points = .points, .test = .test, .signif.label.size = .signif.label.size,
    .defgroupby = "Sample", .grouping.var = "Group",
    .labs = c(NA, paste0("D", perc_value)),
    .title = paste0("D", perc_value, " diversity index"), .subtitle = paste0("Number of clonotypes occupying the ", perc_value, "% of repertoires"),
    .legend = .legend, .leg.title = NA
  )
}

#' @export
vis.immunr_rarefaction <- function(.data, .by = NA, .meta = NA,
                                   .mean = TRUE, .errors = TRUE, .log = FALSE,
                                   .labels = TRUE, ...) {
  .muc.res <- .data

  group_res <- process_metadata_arguments(.data, .by, .meta)
  data_groups <- group_res$groups
  group_column <- group_res$name
  names(data_groups) <- group_res$group_names
  is_grouped <- group_res$is_grouped
  .muc.res$Group <- group_res$group_column

  if ("Type" %in% colnames(.muc.res)) {
    .muc.res$Type <- factor(.muc.res$Type, levels = c("interpolation", "extrapolation"), ordered = TRUE)
  }

  p <- ggplot() +
    xlab("Sample size (clones)") +
    ylab("Estimated diversity (unique clonotypes)") +
    ggtitle("Rarefaction analysis") +
    theme_pubr(legend = "right") +
    theme_cleveland2()

  if (.mean && is_grouped) {
    data_proc <- .muc.res %>%
      group_by(Group, Size) %>%
      summarise(
        MinVal = quantile(Mean, .025, na.rm = TRUE),
        MaxVal = quantile(Mean, .975, na.rm = TRUE),
        MeanVal = mean(Mean, na.rm = TRUE)
      )

    p <- p + geom_line(aes(x = Size, y = MeanVal, color = Group), data = data_proc)

    if (.errors) {
      p <- p + geom_ribbon(aes(x = Size, group = Group, ymin = MinVal, ymax = MaxVal, fill = Group), data = data_proc, alpha = 0.15)
    }
  } else {
    colnames(.muc.res)[2] <- "Q1"
    colnames(.muc.res)[4] <- "Q2"
    if (.errors) {
      p <- p + geom_ribbon(aes(x = Size, group = Sample, ymin = Q1, ymax = Q2), col = "grey80", data = .muc.res, alpha = 0.2)
    }

    if ("Type" %in% colnames(.muc.res)) {
      p <- p +
        geom_line(aes(x = Size, y = Mean, colour = Group, linetype = Type), size = 1, data = .muc.res)
    } else {
      p <- p +
        geom_line(aes(x = Size, y = Mean, colour = Group), size = 1, data = .muc.res)
    }

    if (.labels) {
      tmp <- .muc.res[tapply(1:nrow(.muc.res), .muc.res$Sample, tail, 1), ]
      tmp <- tmp[order(tmp$Group), ]
      p <- p + ggrepel::geom_label_repel(aes(x = Size, y = Mean, label = Sample, fill = Group), data = tmp)
    }
  }

  if (.log) {
    p <- p + scale_x_log10()
  }

  p + .tweak_col(length(unique(.muc.res$Group))) + .tweak_fill(length(unique(.muc.res$Group)))
}


##### Exploratory analysis #####


#' Visualise results of the exploratory analysis
#'
#' @concept explore
#'
#' @aliases vis.immunr_exp_vol vis.immunr_exp_count vis.immunr_exp_len vis.immunr_exp_clones
#' @description An utility function to visualise the output from \code{\link{repExplore}}.
#'
#' @importFrom reshape2 melt
#'
#' @param .data Output from \code{\link{repExplore}}.
#' @param .by Pass NA if you want to plot samples without grouping.
#'
#' You can pass a character vector with one or several column names from ".meta"
#' to group your data before plotting. In this case you should provide ".meta".
#'
#' You can pass a character vector that exactly matches the number of samples in
#' your data, each value should correspond to a sample's property. It will be used
#' to group data based on the values provided. Note that in this case you should
#' pass NA to ".meta".
#'
#' @param .meta A metadata object. An R dataframe with sample names and their properties,
#' such as age, serostatus or hla.
#'
#' @param .errorbars A numeric vector of length two with quantiles for error bars
#' on sectors. Disabled if ".errorbars.off" is TRUE.
#'
#' @param .errorbars.off If TRUE then plot CI bars for distances between each group.
#' Disabled if no group passed to the ".by" argument.
#' @param .points A logical value defining whether points will be visualised or not.
#' @param .test A logical vector whether statistical tests should be applied. See "Details" for more information.
#' @param .signif.label.size An integer value defining the size of text for p-value.
#' @param ... Not used here.
#'
#' @details
#' If data is grouped, then statistical tests for comparing means of groups will be performed, unless \code{.test = FALSE} is supplied.
#' In case there are only two groups, the Wilcoxon rank sum test (https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) is performed
#' (R function \code{\link{wilcox.test}} with an argument \code{exact = FALSE}) for testing if there is a difference in mean rank values between two groups.
#' In case there more than two groups, the Kruskal-Wallis test (https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) is performed (R function \code{\link{kruskal.test}}), that is equivalent to ANOVA for ranks and it tests whether samples from different groups originated from the same distribution.
#' A significant Kruskal-Wallis test indicates that at least one sample stochastically dominates one other sample.
#' Adjusted for multiple comparisons P-values are plotted on the top of groups.
#' P-value adjusting is done using the Holm method (https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method) (also known as Holm-Bonferroni correction).
#' You can execute the command \code{?p.adjust} in the R console to see more.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \link{repExplore} \link{vis}
#'
#' @examples
#' data(immdata)
#' repExplore(immdata$data, "volume") %>% vis()
#' repExplore(immdata$data, "count") %>% vis()
#' repExplore(immdata$data, "len") %>% vis()
#' repExplore(immdata$data, "clones") %>% vis()
#' @export
vis.immunr_exp_vol <- function(.data, .by = NA, .meta = NA,
                               .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                               .points = TRUE, .test = TRUE,
                               .signif.label.size = 3.5, ...) {
  .data <- rename_column(.data, "Volume", "Value")

  vis_bar(
    .data = .data, .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
    .points = .points, .test = .test, .signif.label.size = .signif.label.size,
    .defgroupby = "Sample", .grouping.var = "Group",
    .labs = c(NA, "Clonotypes"),
    .title = "Number of clonotypes", .subtitle = "Number of unique clonotypes in the input data",
    .legend = NA, .leg.title = NA
  )
}

#' @export
vis.immunr_exp_count <- function(.data, .by = NA, .meta = NA, .logx = TRUE, .logy = TRUE, .labs = c("Abundance", "Number of clonotypes"),
                                 .title = "Distribution of clonotype abundances", .subtitle = NULL,
                                 .legend = TRUE, .leg.title = NA, ...) {
  group_res <- process_metadata_arguments(.data, .by, .meta)
  .data$Group <- group_res$group_column

  if (is.na(.leg.title)) {
    if (group_res$is_grouped) {
      .leg.title <- "Group"
    } else {
      .leg.title <- "Sample"
    }
  }

  p <- ggplot() +
    geom_line(aes(x = Clone.num, y = Clonotypes, col = Group, group = Sample), data = .data, stat = "identity")

  if (.logx) {
    p <- p + scale_x_log10()
  }
  if (.logy) {
    p <- p + scale_y_log10()
  }

  p <- p + labs(x = .labs[1], y = .labs[2], title = .title, subtitle = .subtitle, color = .leg.title) +
    theme_pubr(legend = "right") + rotate_x_text() + theme_cleveland2()

  p
}

#' @export
vis.immunr_exp_len <- function(.data, .by = NA, .meta = NA,
                               .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                               .points = TRUE, .test = TRUE,
                               .signif.label.size = 3.5, ...) {
  .data <- rename_column(.data, "Count", "Value")

  vis_bar(
    .data = .data, .by = .by, .meta = .meta,
    .errorbars = c(0.025, 0.975), .errorbars.off = .errorbars.off, .stack = FALSE,
    .points = .points, .test = .test, .signif.label.size = .signif.label.size,
    .defgroupby = "Sample", .grouping.var = "Length",
    .labs = c("CDR3 length", "Clonotypes"),
    .title = "Distribution of CDR3 lengths", .subtitle = NULL,
    .legend = TRUE, .leg.title = NA
  )
}

#' @export
vis.immunr_exp_clones <- function(.data, .by = NA, .meta = NA,
                                  .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                                  .points = TRUE, .test = TRUE,
                                  .signif.label.size = 3.5, ...) {
  .data <- rename_column(.data, "Clones", "Value")

  vis_bar(
    .data = .data, .by = .by, .meta = .meta,
    .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
    .points = .points, .test = .test, .signif.label.size = .signif.label.size,
    .defgroupby = "Sample", .grouping.var = "Group",
    .labs = c(NA, "Clones"),
    .title = "Number of clones", .subtitle = "Number of clones in the input data",
    .legend = NA, .leg.title = NA
  )
}


##### Kmer analysis & sequence logo plot #####


#' Most frequent kmers visualisation.
#'
#' @concept kmers
#'
#' @name vis.immunr_kmer_table
#'
#' @description
#' Plot a distribution (bar plot) of the most frequent kmers in a data.
#'
#' @param .data Data frame with two columns "Kmers" and "Count" or a list with such data frames. See Examples.
#' @param .head Number of the most frequent kmers to choose for plotting from each data frame.
#' @param .position Character vector of length 1. Position of bars for each kmers. Value for the \code{ggplot2} argument \code{position}.
#' @param .log Logical. If TRUE then plot log-scaled plots.
#' @param ... Not used here.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \code{get.kmers}
#'
#' @examples
#' # Load necessary data and package.
#' data(immdata)
#' # Get 5-mers.
#' imm.km <- getKmers(immdata$data[[1]], 5)
#' # Plots for kmer proportions in each data frame in immdata.
#' p1 <- vis(imm.km, .position = "stack")
#' p2 <- vis(imm.km, .position = "fill")
#' p1 + p2
#' @export
vis.immunr_kmer_table <- function(.data, .head = 100, .position = c("stack", "dodge", "fill"), .log = FALSE, ...) {
  .position <- switch(substr(.position[1], 1, 1),
    s = "stack",
    d = "dodge",
    f = "fill"
  )
  # .data[is.na(.data)] <- 0

  max_counts <- apply(.data[, -1], 1, max, na.rm = TRUE)
  max_indices <- head(order(max_counts, decreasing = TRUE), .head)
  n_samples <- ncol(.data) - 1
  .data <- .data[max_indices, ]

  melted <- reshape2::melt(.data, id.vars = "Kmer")
  colnames(melted) <- c("Kmer", "Sample", "Count")
  p <- ggplot() +
    geom_bar(aes(x = Kmer, y = Count, fill = Sample),
      col = "black",
      data = melted, stat = "identity", position = .position
    ) +
    theme_pubr(legend = "right") +
    rotate_x_text() +
    theme_cleveland2()
  if (.position[1] == "stack" || .position[1] == "dodge") {
    p <- p + ylab("Count")
  } else {
    p <- p + ylab("Proportions")
  }

  if (.log) {
    p <- p + scale_y_log10(expand = c(0, 0))
  } else {
    p <- p + scale_y_continuous(expand = c(0, 0))
  }

  p + scale_x_discrete(expand = c(0, 0)) + .tweak_fill(n_samples) + ggtitle("Kmers distribution")
}

#' Sequence logo plots for amino acid profiles.
#'
#' @concept kmers
#'
#' @importFrom ggseqlogo geom_logo theme_logo
#'
#' @aliases vis_seqlogo vis_textlogo
#'
#' @name vis_textlogo
#'
#' @usage
#' vis_textlogo(.data, .replace.zero.with.na = TRUE, .width = 0.1, ...)
#'
#' vis_seqlogo(.data, .scheme = "chemistry", ...)
#'
#' @description
#' Plot sequence logo plots for visualising of amino acid motif sequences / profiles.
#'
#' `vis_textlogo` plots sequences in a text format - each letter has the same height. Useful when there
#' are no big differences between occurences of amino acids in the motif.
#'
#' `vis_seqlogo` is a traditional sequence logo plots. Useful when there are one or two amino acids
#' with clear differences in their occurrences.
#'
#' @param .data Output from the \code{kmer.profile} function.
#' @param .replace.zero.with.na if TRUE then replace all zeros with NAs, therefore letters with
#' zero frequency wont appear at the plot.
#' @param .scheme Character. An argumentt passed to \link{geom_logo} specifying how to colour symbols.
#' @param .width Width for jitter, i.e., how much points will scatter around the verical line. Pass 0 (zero)
#' to plot points on the straight vertical line for each position.
#' @param ... Not used here.
#'
#' @return
#' A ggplot2 object.
#'
#' @seealso \link{getKmers}, \link{kmer_profile}
#'
#' @examples
#' data(immdata)
#' kmers <- getKmers(immdata$data[[1]], 5)
#' ppm <- kmer_profile(kmers, "prob")
#' vis(ppm, .plot = "text")
#' vis(ppm, .plot = "seq")
#'
#' d <- kmer_profile(c("CASLL", "CASSQ", "CASGL"))
#' vis_textlogo(d)
#' vis_seqlogo(d)
#' @export
vis_textlogo <- function(.data, .replace.zero.with.na = TRUE, .width = 0.1, ...) {
  # ToDo: make different color schemas, for type of aminoacids (polarity, etc), etc

  .data <- reshape2::melt(.data)
  if (.replace.zero.with.na) {
    .data$value[.data$value == 0] <- NA
  }
  colnames(.data) <- c("AA", "Position", "Freq")
  ggplot(aes(x = Position, y = Freq, colour = AA), data = .data) +
    geom_jitter(colour = "black", width = .width) +
    ggrepel::geom_label_repel(aes(label = AA), size = 5) +
    xlab("Position") +
    ylab("Proportion") +
    theme_pubr(legend = "none")
}

#' @export
vis_seqlogo <- function(.data, .scheme = "chemistry", ...) {
  ggplot() +
    ggseqlogo::geom_logo(.data, method = "custom", col_scheme = .scheme) +
    ggseqlogo::theme_logo()
}


#' Visualise kmer profiles
#'
#' @concept kmers
#'
#' @param .data Kmer data, an output from \link{kmer_profile}.
#' @param .plot String specifying the plot type:
#'
#' - "seqlogo" for traditional sequence logo plots using \link{vis_seqlogo};
#'
#' - "textlogo" for modified approach to sequence logo plots via text labels using \link{vis_textlogo};
#' @param ... Other arguments passed to \link{vis_textlogo} or \link{vis_seqlogo}, depending
#' on the ".plot" argument.
#'
#' @return
#' A ggplot2 object.
#'
#' @examples
#' data(immdata)
#' getKmers(immdata$data[[1]], 5) %>%
#'   kmer_profile() %>%
#'   vis("seqlogo")
#' @export
vis_immunr_kmer_profile_main <- function(.data, .plot, ...) {
  if (.plot[1] == "text") {
    vis_textlogo(.data, ...)
  } else if (.plot[1] == "seq") {
    vis_seqlogo(.data, ...)
  } else if (.plot[1] == "textlogo") {
    vis_textlogo(.data, ...)
  } else if (.plot[1] == "seqlogo") {
    vis_seqlogo(.data, ...)
  } else {
    stop('Error: unknown plot type. Please provide "seq" for seqlogo plots (vis_seqlogo) or "text" for textlogo plots (vis_textlogo) ')
  }
}

#' @export
vis.immunr_kmer_profile_pfm <- function(.data, .plot = c("textlogo", "seqlogo"), ...) {
  p <- vis_immunr_kmer_profile_main(.data, .plot, ...)
  p + ggtitle("Position frequency matrix")
}

#' @export
vis.immunr_kmer_profile_ppm <- function(.data, .plot = c("textlogo", "seqlogo"), ...) {
  p <- vis_immunr_kmer_profile_main(.data, .plot, ...)
  p + ggtitle("Position probability matrix")
}

#' @export
vis.immunr_kmer_profile_pwm <- function(.data, .plot = c("textlogo", "seqlogo"), ...) {
  p <- vis_immunr_kmer_profile_main(.data, .plot, ...)
  p + ggtitle("Position weight matrix")
}

#' @export
vis.immunr_kmer_profile_self <- function(.data, .plot = c("textlogo", "seqlogo"), ...) {
  p <- vis_immunr_kmer_profile_main(.data, .plot, ...)
  p + ggtitle("Position self-information matrix")
}


##### Other & WIP visualisations #####


#' Visualise clonotype dynamics
#'
#' @concept dynamics
#'
#' @importFrom data.table setnames melt.data.table
#' @importFrom ggalluvial geom_flow geom_stratum
#'
#' @name vis.immunr_dynamics
#'
#' @param .data Output from the \link{trackClonotypes} function.
#' @param .plot Character. Either "smooth", "area" or "line". Each specifies a type of plot for visualisation of clonotype dynamics.
#' @param .order Numeric or character vector. Specifies the order to samples, e.g., it used for ordering samples
#' by timepoints. Either See "Examples" below for more details.
#' @param .log Logical. If TRUE then use log-scale for the frequency axis.
#' @param ... Not used here.
#'
#' @return
#' A ggplot2 object.
#'
#' @examples
#' # Load an example data that comes with immunarch
#' data(immdata)
#'
#' # Make the data smaller in order to speed up the examples
#' immdata$data <- immdata$data[c(1, 2, 3, 7, 8, 9)]
#' immdata$meta <- immdata$meta[c(1, 2, 3, 7, 8, 9), ]
#'
#' # Option 1
#' # Choose the first 10 amino acid clonotype sequences
#' # from the first repertoire to track
#' tc <- trackClonotypes(immdata$data, list(1, 10), .col = "aa")
#' # Choose the first 20 nucleotide clonotype sequences
#' # and their V genes from the "MS1" repertoire to track
#' tc <- trackClonotypes(immdata$data, list("MS1", 20), .col = "nt+v")
#'
#' # Option 2
#' # Choose clonotypes with amino acid sequences "CASRGLITDTQYF" or "CSASRGSPNEQYF"
#' tc <- trackClonotypes(immdata$data, c("CASRGLITDTQYF", "CSASRGSPNEQYF"), .col = "aa")
#'
#' # Option 3
#' # Choose the first 10 clonotypes from the first repertoire
#' # with amino acid sequences and V segments
#' target <- immdata$data[[1]] %>%
#'   select(CDR3.aa, V.name) %>%
#'   head(10)
#' tc <- trackClonotypes(immdata$data, target)
#'
#' # Visualise the output regardless of the chosen option
#' # Therea are three way to visualise it, regulated by the .plot argument
#' vis(tc, .plot = "smooth")
#' vis(tc, .plot = "area")
#' vis(tc, .plot = "line")
#'
#' # Visualising timepoints
#' # First, we create an additional column in the metadata with randomly choosen timepoints:
#' immdata$meta$Timepoint <- sample(1:length(immdata$data))
#' immdata$meta
#' # Next, we create a vector with samples in the right order,
#' # according to the "Timepoint" column (from smallest to greatest):
#' sample_order <- order(immdata$meta$Timepoint)
#' # Sanity check: timepoints are following the right order:
#' immdata$meta$Timepoint[sample_order]
#' # Samples, sorted by the timepoints:
#' immdata$meta$Sample[sample_order]
#' # And finally, we visualise the data:
#' vis(tc, .order = sample_order)
#' @export
vis.immunr_dynamics <- function(.data, .plot = c("smooth", "area", "line"), .order = NA, .log = FALSE, ...) {
  .plot <- .plot[1]
  if (!(.plot %in% c("smooth", "area", "line"))) {
    stop("Error: unknown plot identifier \"", .plot, "\". Please provide one of the following: \"smooth\", \"area\" or \"line\".")
  }

  y_lab_title <- "Count"
  melted <- melt(.data) %>%
    lazy_dt() %>%
    rename(Count = value, Sample = variable) %>%
    collect()
  setDT(melted)
  if (max(melted[["Count"]]) <= 1) {
    y_lab_title <- "Proportion"
  }

  last_obj_column_i <- match("Sample", names(melted)) - 1
  column_names <- names(melted)[1:last_obj_column_i]
  melted[, Clonotype := do.call(paste, .SD), .SDcols = column_names]

  if (!is.na(.order[1])) {
    if (is.character(.order)) {
      ordered_sample_names <- .order
    } else {
      sample_names <- unique(melted$Sample)
      ordered_sample_names <- sample_names[.order]
    }
    melted <- melted[melted$Sample %in% ordered_sample_names, ]
    melted$Sample <- factor(melted$Sample, levels = ordered_sample_names, ordered = TRUE)
  }

  melted$Count <- melted$Count + min(melted$Count[melted$Count != 0]) / 1000
  p <- ggplot(melted, aes(
    x = Sample, fill = Clonotype, stratum = Clonotype,
    alluvium = Clonotype, y = Count, label = Clonotype
  ))

  if (.plot == "smooth") {
    p <- p +
      geom_flow() +
      geom_stratum()
  } else if (.plot == "area") {
    p <- p +
      geom_area(aes(group = Clonotype), color = "black")
  } else {
    p <- p +
      geom_line(aes(color = Clonotype, group = Clonotype))
  }

  if (.log) {
    p <- p + scale_y_log10()
  }

  p +
    ylab(y_lab_title) +
    ggtitle("Clonotype tracking") +
    theme_pubr(legend = "right") + rotate_x_text(90) + theme_cleveland2()
}



# vis.immunr_mutation_network <- function (.data) {
#   stop(IMMUNR_ERROR_NOT_IMPL)
# }


# vis.immunr_cdr_prop <- function (.data, .by = NA, .meta = NA, .plot = c("box", "hist")) {
#   stop(IMMUNR_ERROR_NOT_IMPL)
# }
