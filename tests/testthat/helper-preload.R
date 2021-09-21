apply_DF_DT <- function(df_data, dt_data, .fun, ...) {
  res1 <- .fun(df_data, ...)
  res2 <- .fun(dt_data, ...)
  list(df = res1, dt = res2)
}

vis_results <- function(res_list, ...) {
  p_df <- vis(res_list$df, ...)
  p_dt <- vis(res_list$dt, ...)
  list(df = p_df, dt = p_dt)
}

check_for_mutation <- function(.frame, .table) {
  expect_equal(lapply(.frame, as.data.table), .table)
}

# add pretfix
ap <- function(.name, .prefix) {
  paste0(.prefix, .name)
}

add_mock_sample <- function(.immdata, .sample_name, .meta = NA) {
  # copy dataframe of 1st sample to the new sample
  .immdata$data[[.sample_name]] <- .immdata$data[[1]]
  .immdata$meta %<>% add_row(Sample = .sample_name)
  if (!is.na(.meta)) {
    # .meta must be a list of lists with 2 elements: meta_column and meta_value
    for (i in seq_along(.meta)) {
      meta_column <- .meta[[i]][["meta_column"]]
      meta_value <- .meta[[i]][["meta_value"]]
      # set meta value in the last (just added) row
      .immdata$meta[nrow(.immdata$meta), meta_column] <- meta_value
    }
  }
  return(.immdata)
}

data(immdata)
frame_data <- immdata$data
table_data <- lapply(frame_data, as.data.table)
