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

data(immdata)
frame_data <- immdata$data
table_data <- lapply(frame_data, as.data.table)
