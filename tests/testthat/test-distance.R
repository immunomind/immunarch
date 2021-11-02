data(immdata)
short_immdata <- immdata$data
short_immdata <- map(short_immdata, ~ .x %>% head(3)) # less sample size saves time in computations

f <- function(x, y) {
  res <- matrix(nrow = length(x), ncol = length(y))
  for (i in seq_along(x)) {
    res[i, ] <- abs(nchar(x[i]) - nchar(y))
  }
  dimnames(res) <- list(x, y)
  return(as.dist(res))
}

# Arrange
positive_test_cases <- list(
  "Changing method" = list(
    args = list(
      .data = short_immdata[1],
      .method = "lv"
    ),
    result = c(18, 20, 14)
  ),
  "Changing column" = list(
    args = list(
      .data = short_immdata[1],
      .col = "CDR3.aa",
      .method = "lv"
    ),
    result = c(10, 12, 6)
  ),
  "Custom_func" = list(args = list(.data = short_immdata[1], .method = f), result = c(6, 3, 9)),
  "Custom_func_col" = list(args = list(.data = short_immdata[1], .col = "CDR3.aa", .method = f), result = c(2, 1, 3))
)

negative_test_cases <- list(
  "Immdata test" = list(
    args = list(
      .data = immdata
    )
  ),
  "Wrong col" = list(
    args = list(
      .data = short_immdata,
      .col = "aa"
    )
  ),
  "Wrong method" = list(
    args = list(
      .data = short_immdata,
      .method = "ddddd"
    )
  )
)

# Act
args <- map(positive_test_cases, "args")
results <- map(positive_test_cases, "result")
positive_act_result <- map(args, ~ do.call(seqDist, .x)[[1]] %>% as.numeric())
positive_test_values <- list(names(positive_test_cases), positive_act_result, results)
negative_args <- map(negative_test_cases, "args")
# Assert
pmap(positive_test_values, ~ test_that(..1, expect_equal(..2, ..3)))
## for negative tests act can be done only with assert
map2(names(negative_test_cases), negative_args, ~ test_that(.x, expect_error(do.call(seqDist, .y))))
