IMMUNARCH_METHOD_REGISTRY <- new.env(parent = emptyenv())


#' Common arguments for immundata helpers
#' @keywords internal
#' @param autojoin If TRUE, join repertoire metadata by the schema repertoire id.
#'   For `format="analysis"`, metadata is joined to the long table; for
#'   `format="ml"`, it’s joined after pivoting to wide. Defaults to
#'   `getOption("immundata.autojoin", FALSE)`.
#' @param format One of `"analysis"` (long tibble with `repertoire_id`, facet
#'   columns, and `value`) or `"ml"` (wide/unmelted table of features).
#' @param features Character vector of **feature keys** to keep when
#'   `format="ml"`. If `NULL`, features are derived from the data. A feature key
#'   looks like `family.method|facet1=...;facet2=...` (e.g.,
#'   `airr_stats.genes|v_call=TRBV7-2`).
im_common_args <- function(
    autojoin = getOption("immundata.autojoin", FALSE),
    format   = c("analysis", "ml"),
    features = NULL) {} # nocov


im_method <- function(core, family, name) {
  checkmate::assert_function(core, args = c("idata"))
  checkmate::assert_string(family)
  checkmate::assert_string(name)

  # Merge core formals with wrapper defaults to enable autocompletion
  core_fmls <- formals(core)
  if (!"idata" %in% names(core_fmls)) {
    cli::cli_abort("Core method must declare an {.code idata} argument.")
  }
  if (any(c("autojoin", "format", "features") %in% names(core_fmls))) {
    cli::cli_abort("Core method must not declare {.code autojoin}, {.code format}, or {.code features}.")
  }

  wrapper <- function() { }
  formals(wrapper) <- c(
    core_fmls,
    formals(im_common_args)
  )
  environment(wrapper) <- environment()

  body(wrapper) <- substitute(
    {
      format <- match.arg(format)

      # Pre: validate idata + schema
      checkmate::assert_r6(idata, "ImmunData")

      # Build argument list for core from our own formals (now visible to the user)
      .core_args <- mget(names(core_fmls), inherits = TRUE)
      out <- do.call(core, .core_args)

      out
    },
    list(core = core, core_fmls = core_fmls)
  )

  wrapper
}


register_immunarch_method <- function(core, family, name, register_family = TRUE) {
  fn <- im_method(core, family, name)

  if (isTRUE(register_family) && exists("register_airr_family", mode = "function", inherits = TRUE)) {
    try(register_airr_family(family), silent = TRUE)
  }

  if (exists("register_airr_method", mode = "function", inherits = TRUE)) {
    try(register_airr_method(
      family_name = family,
      method_name = name,
      fn = fn
    ), silent = TRUE)
  }

  fn
}
