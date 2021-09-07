# meta_conditions must contain a string with dplyr conditions for metadata
# example usage: immdata <- filter_by_meta(immdata, "Age > 16 & Lane == \"C\"")
filter_by_meta <- function(.data, .meta_conditions) {
  filtered_meta <- filter(.data$meta, eval(rlang::parse_expr(.meta_conditions)))
  filtered_data <- .data$data[names(.data$data) %in% filtered_meta$Sample]

  # return filtered dataset
  list(data = filtered_data, meta = filtered_meta)
}
