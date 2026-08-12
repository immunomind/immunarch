make_immunarch_features <- function(ires, method_name = NULL, feature_col = NULL, value_col = NULL) {
  checkmate::assert_data_frame(ires)
  checkmate::assert_character(method_name, null.ok = TRUE)
  checkmate::assert_character(feature_col, null.ok = TRUE)
  checkmate::assert_character(value_col, null.ok = TRUE)

  to_rename <- c("feature" = feature_col, "value" = value_col)

  if (is.null(method_name)) {
    ires <- ires |>
      mutate(method = method_name)
  }

  ires |>
    rename(to_rename)
}

immunarch_methods <- function(family_name = NULL) {
  checkmate::assert_string(family_name, null.ok = TRUE)

  if (is.null(family_name)) {
    ls(IMMUNARCH_METHOD_REGISTRY)
  } else if (family_name %in% IMMUNARCH_METHOD_REGISTRY) {
    ls(IMMUNARCH_METHOD_REGISTRY[[family_name]])
  } else {
    cli::cli_abort("No such family name: {family_name}")
  }
}
