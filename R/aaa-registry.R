IMMUNARCH_METHOD_REGISTRY <- new.env(parent = emptyenv())

IMMUNARCH_VIS_REGISTRY <- new.env(parent = emptyenv())

IMMUNARCH_CLASS_PREFIX <- "immunarch_res"

IMMUNARCH_AUTOJOIN_OPTION <- "immunarch.autojoin"

IMMUNARCH_AUTOJOIN_DEFAULT <- TRUE


# ---------------------------------------------------------------------------- #
# --- Common arguments
# ---------------------------------------------------------------------------- #


#' Common arguments for immunarch helpers
#' @keywords internal
#' @param autojoin Logical. If TRUE, join repertoire metadata by the schema repertoire id.
#'  Change the default behaviour by calling `options(immunarch.autojoin = FALSE)`.
#' @param format String. One of `"long"` ("long" tibble with `imd_repertoire_id`, facet
#'   columns, and `value`; useful for visualizations) or `"wide"` (wide/unmelted table of features,
#'   with each row corresponding to a specific repertoire / pair of repertoires; useful for Machine Learning).
method_args_airr <- function(
  autojoin = getOption(IMMUNARCH_AUTOJOIN_OPTION, IMMUNARCH_AUTOJOIN_DEFAULT),
  format = c("long", "wide")
) {} # nocov


# ---------------------------------------------------------------------------- #
# --- Immunarch results attributes
# ---------------------------------------------------------------------------- #


im_norm <- function(x) {
  x <- tolower(x)
  gsub("[^a-z0-9]+", "_", x)
}


im_result_class <- function(family, name = NULL) {
  fam <- im_norm(family)
  if (is.null(name)) {
    paste0(IMMUNARCH_CLASS_PREFIX, "_", fam)
  } else {
    nm <- im_norm(name)
    paste0(IMMUNARCH_CLASS_PREFIX, "_", fam, "_", nm)
  }
}


im_as_result <- function(x, family, name) {
  # Wrap any object as an Immunarch result, preserving original classes
  cls_full <- im_result_class(family, name)
  cls_fam <- im_result_class(family, NULL)
  # TODO: maybe I need the "airr" or "receptor" instead of IMMUNARCH_CLASS_PREFIX?
  class(x) <- c(cls_full, cls_fam, IMMUNARCH_CLASS_PREFIX, class(x))
  x
}


# ---------------------------------------------------------------------------- #
# --- Immunarch methods
# ---------------------------------------------------------------------------- #


im_method <- function(
  core,
  family,
  name,
  required_cols = NULL,
  need_repertoires = TRUE,
  airr_args = FALSE,
  result_class = TRUE,
  vdj_args = FALSE
) {
  checkmate::assert_function(core, args = c("idata"))
  checkmate::assert_string(family)
  checkmate::assert_string(name)
  checkmate::assert_logical(need_repertoires)

  core_fmls <- formals(core)
  if (!"idata" %in% names(core_fmls)) {
    cli::cli_abort("Core method must declare an {.code idata} argument.")
  }
  if (any(c("autojoin", "format") %in% names(core_fmls))) {
    cli::cli_abort("Core method must not declare {.code autojoin} or {.code format}.")
  }

  # required_cols = columns expected in idata$annotations
  if (!is.null(required_cols)) {
    checkmate::assert_character(required_cols, any.missing = FALSE)
  }

  wrapper <- function() {}
  formals(wrapper) <- c(
    core_fmls,
    if (airr_args) formals(method_args_airr),
    if (vdj_args) alist(output = c("immundata", "data.frame"))
  )
  environment(wrapper) <- environment()

  body(wrapper) <- substitute(
    {
      if (airr_args) {
        format <- match.arg(format)
      }
      if (vdj_args) {
        output <- match.arg(output)
      }
      checkmate::assert_r6(idata, "ImmunData")

      if (need_repertoires) {
        if (is.null(idata$schema_repertoire)) {
          cli::cli_abort("Repertoire aggregation is needed for this function. Run {.code ?agg_repertoires} for more info.")
        }
      }

      # For each argument name in `required_cols`, fetch its runtime value
      # and ensure the referenced columns exist in idata$annotations.
      if (length(required_cols)) {
        ann_cols <- colnames(idata$annotations)
        for (.arg in required_cols) {
          .val <- get(.arg, inherits = TRUE)

          # Just in case
          if (is.null(.val) || (length(.val) == 1 && is.na(.val))) next

          # Allow symbols or character vectors
          if (rlang::is_symbol(.val)) {
            .val <- rlang::as_string(.val)
          }

          if (!is.character(.val)) {
            cli::cli_abort("Argument {.code {.arg}} must be a character (column name) or character vector; got a {.code {class(.val)[1]}}.")
          }

          missing <- setdiff(.val, ann_cols)
          if (length(missing)) {
            suggest <- function(x, pool, n = 3) {
              if (!length(pool)) {
                return(character())
              }
              d <- utils::adist(x, pool)
              pool[order(d)][seq_len(min(n, length(pool)))]
            }
            hints <- unique(unlist(lapply(missing, suggest, pool = ann_cols, n = 3)))
            cli::cli_abort("Passed column name(s) [{.code {missing}}] is not in the input ImmunData. Did you mean [{.code {hints}}]?")
          }
        }
      }

      # Call core with its own formals
      .core_args <- mget(names(core_fmls), inherits = TRUE)
      out <- do.call(core, .core_args)

      if (vdj_args && output == "data.frame") {
        out <- dplyr::collect(out)
      }

      # Autojoin: join repertoire metadata if requested and applicable (!)
      if (airr_args &&
        isTRUE(get("autojoin", inherits = FALSE)) &&
        !is.null(idata$repertoires)) {
        rep_col <- immundata::imd_schema("repertoire")
        if (!is.null(rep_col) && rep_col %in% names(out)) {
          out <- dplyr::left_join(out, idata$repertoires |> select(all_of(c(rep_col, idata$schema_repertoire))), by = rep_col)
        }
      }

      # Wrap the output to assign correct classes
      if (result_class) {
        out <- im_as_result(out, family, name)
      }

      out
    },
    list(
      core = core,
      core_fmls = core_fmls,
      required_cols = required_cols,
      airr_args = airr_args,
      result_class = result_class,
      vdj_args = vdj_args
    )
  )

  wrapper
}


#' Register an Immunarch method (developer)
#'
#' `r lifecycle::badge("experimental")`
#'
#' Wrap a core implementation into a user-facing function and (optionally)
#' register it in the in-memory method registry. Domain-specific registrars add
#' the appropriate arguments and post-processing while sharing the same safety checks.
#'
#' ## What your core must look like
#' * Signature: `function(idata, ...)`
#' * **Must not** declare `autojoin` or `format` - these are added to AIRR methods.
#'
#' ## What the wrapper adds
#' * AIRR args from `method_args_airr()`: `autojoin` and `format`
#'   (with `autojoin` default controlled by `getOption("immunarch.autojoin", TRUE)`).
#' * Validates `idata` is an [immundata::ImmunData] object.
#' * Ensures all columns in `required_cols` exist in `idata$annotations`.
#' * If `autojoin = TRUE` and the result is a data frame containing the repertoire id
#'   column (`immundata::imd_schema("repertoire")`), joins repertoire metadata from
#'   `idata$repertoires`.
#'
#' @param core A function with signature `function(idata, ...)`. This is your core
#'   implementation; it must accept an `ImmunData` as the first argument and **must not**
#'   declare `autojoin` or `format`.
#' @param family String. Method family name used for dispatch (e.g., `"airr_stats"`).
#' @param name String. Method name within the family (e.g., `"lengths"`).
#' @param register_family Logical (default `TRUE`). If `TRUE`, attempts to create/ensure
#'   the family environment by calling `register_method_family()` when available.
#' @param required_cols Character vector of column names that **must** be present in
#'   `idata$annotations`. Use this to declare the minimal input schema your core needs.
#' @param need_repertoires Logical. Use this to declare the necessity of having aggregated
#'  repertoires.
#' @param vdj_args Logical scalar. If `TRUE`, add the VDJ `output` argument and
#'   collect the result when `output = "data.frame"`. The default is `FALSE`.
#'
#' @return A **function** - the user-facing wrapper around `core`. Typical usage is to
#' assign it to the exported symbol of the method, e.g.:
#' `airr_stats_lengths <- register_method_airr(...)`.
#'
#' @examples
#' \dontrun{
#' # Minimal core implementation (must accept `idata`)
#' airr_stats_lengths_impl <- function(idata, seq_col = "cdr3_aa") {
#'   repertoire_col <- immundata::imd_schema("repertoire")
#'   idata$annotations |>
#'     dplyr::distinct(!!rlang::sym(repertoire_col), !!rlang::sym(seq_col)) |>
#'     dplyr::mutate(seq_len = nchar(!!rlang::sym(seq_col))) |>
#'     dplyr::count(!!rlang::sym(repertoire_col), seq_len, name = "n")
#' }
#'
#' # Register and expose a user-facing function
#' airr_stats_lengths <- register_method_airr(
#'   core = airr_stats_lengths_impl,
#'   family = "airr_stats",
#'   name = "lengths",
#'   required_cols = c("cdr3_aa", immundata::imd_schema("repertoire"))
#' )
#'
#' # Optional: call via dispatcher
#' # make_airr_dispatcher("airr_stats")(idata = immdata, method = "lengths")
#' }
#'
#' @keywords internal
register_method_impl <- function(
  core,
  family,
  name,
  register_family = TRUE,
  required_cols = NULL,
  need_repertoires = TRUE,
  airr_args = FALSE,
  result_class = TRUE,
  vdj_args = FALSE
) {
  fn <- im_method(
    core,
    family,
    name,
    required_cols = required_cols,
    need_repertoires = need_repertoires,
    airr_args = airr_args,
    result_class = result_class,
    vdj_args = vdj_args
  )

  if (isTRUE(register_family) && exists("register_method_family", mode = "function", inherits = TRUE)) {
    try(register_method_family(family), silent = TRUE)
  }

  if (exists("register_method_record", mode = "function", inherits = TRUE)) {
    try(register_method_record(
      family_name = family,
      method_name = name,
      fn = fn
    ), silent = TRUE)
  }

  # Link visualisation to a method if visualisation was already created
  im_ensure_vis_s3_for(family, name)

  fn
}


register_method_airr <- function(...) {
  register_method_impl(..., airr_args = TRUE)
}


register_method_annotate <- function(...) {
  register_method_impl(..., result_class = FALSE)
}


register_method_vdj <- function(...) {
  register_method_impl(..., vdj_args = TRUE, result_class = FALSE)
}


register_method_repsim <- function(...) {
  register_method_impl(...)
}


register_method_dist <- function(...) {
  register_method_impl(...)
}


register_method_clust <- function(...) {
  register_method_impl(...)
}


# ---------------------------------------------------------------------------- #
# --- Immunarch visualisations
# ---------------------------------------------------------------------------- #


IMMUNARCH_VIS_REGISTRY <- new.env(parent = emptyenv())

.im_ns <- function() asNamespace("immunarch")

im_vis_s3_exists <- function(class) {
  !is.null(utils::getS3method("vis", class, optional = TRUE))
}

im_ensure_vis_s3_for <- function(family, name) {
  cls <- im_result_class(family, name)
  fn <- IMMUNARCH_VIS_REGISTRY[[cls]]
  if (!is.function(fn)) {
    return(invisible(FALSE))
  }
  if (im_vis_s3_exists(cls)) {
    return(invisible(FALSE))
  }

  method <- function(.data, ...) {
    f <- IMMUNARCH_VIS_REGISTRY[[cls]]
    if (!is.function(f)) cli::cli_abort("Visualization for {.code {cls}} not found.")
    f(.data, ...)
  }

  base::registerS3method("vis", cls, method, envir = .im_ns())
  invisible(TRUE)
}

register_immunarch_visualisation <- function(fn, family, name) {
  checkmate::assert_function(fn, args = c(".data"))
  checkmate::assert_string(family)
  checkmate::assert_string(name, null.ok = TRUE)

  cls <- im_result_class(family, name)
  assign(cls, fn, envir = IMMUNARCH_VIS_REGISTRY)

  # immediate S3 registration (errors if vis generic not yet defined)
  if (!exists("vis", envir = .im_ns(), inherits = FALSE)) {
    stop("vis() generic must be defined before registering visualisations.")
  }
  im_ensure_vis_s3_for(family, name)
  invisible(cls)
}
