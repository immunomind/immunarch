make_immunarch_method_record <- function(fn, method_name, ...) {
  list(
    fn = fn,
    method_name = method_name
  )
}

register_airr_family <- function(family_name) {
  fam_env <- IMMUNARCH_METHOD_REGISTRY[[family_name]]

  if (is.null(fam_env)) {
    fam_env <- new.env(parent = emptyenv())
    IMMUNARCH_METHOD_REGISTRY[[family_name]] <- fam_env
  }
}

register_airr_method <- function(family_name, method_name, fn) {
  # checkmate checks

  fam_env <- IMMUNARCH_METHOD_REGISTRY[[family_name]]

  assign(method_name,
    make_immunarch_method_record(
      fn = fn,
      method_name = method_name
    ),
    envir = fam_env
  )
}

get_airr_method <- function(family_name, method_name, verbose = TRUE) {
  # TODO: checkmate checks

  fam_env <- IMMUNARCH_METHOD_REGISTRY[[family_name]]

  if (is.null(method_name)) {
    if (verbose) {
      cli::cli_text("Available methods in {.code {family_name}}:")
      cli::cli_ol(ls(fam_env))
      return(invisible(ls(fam_env)))
    } else {
      return(ls(fam_env))
    }
  } else {
    # check and throw an error if no such object
    record <- get0(method_name, envir = fam_env, ifnotfound = NULL)

    if (is.null(record)) {
      cli::cli_abort("No such method: {.code {method_name}}. Available methods: {ls(fam_env)}")
    }

    record
  }
}

make_airr_dispatcher <- function(family_name) {
  function(idata = NULL, method = NULL, ...) {
    checkmate::assert_r6(idata, "ImmunData", null.ok = TRUE)
    checkmate::assert_character(method, null.ok = TRUE)

    if (is.null(idata)) {
      # list available methods for this family
      get_airr_method(family_name, NULL)
    } else if (is.null(method)) {
      cli::cli_abort("{.code idata} provided, but {.code method} is null, aborting the execution; please provide either both {.code idata} and {.code method}, or leave them as nulls to show the list of available methods")
    } else {
      checkmate::assert_data_frame(idata$repertoires, null.ok = FALSE)

      record <- get_airr_method(family_name, method)
      fn <- record$fn
      res <- fn(idata, ...)
      res
    }
  }
}
