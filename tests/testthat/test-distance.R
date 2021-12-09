data(immdata)
library(purrr)
short_immdata <- map(immdata$data, ~ .x %>% head(1000)) # smaller sample size saves time in computations

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
    result = c(0, 6, 6, 0)
  ),
  "Changing column" = list(
    args = list(
      .data = short_immdata[1],
      .col = "CDR3.aa",
      .method = "lv"
    ),
    result = c(0, 6, 6, 0)
  ),
  "Custom_func" = list(args = list(.data = short_immdata[1], .method = f, .group_by_seqLength = FALSE), result = c(0, 3, 3, 0)),
  "Group_by changing" = list(args = list(.data = short_immdata[1], .group_by = "V.name"), result = c(0, 21, 19, 18, 13, 21, 0, 18, 23, 19, 19, 18, 0, 16, 15, 18, 23, 16, 0, 15, 13, 19, 15, 15, 0))
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
    ),
    "Wrong group_by" = list(
      args = list(
        .data = short_immdata,
        .group_by = "ddddd"
      )
    )
  )
)

# Act
args <- map(positive_test_cases, "args")
results <- map(positive_test_cases, "result")
positive_act_result <- map(args, ~ do.call(seqDist, .x)[[1]][10]) %>%
  map(1) %>%
  map(., ~ as.matrix(.x) %>% as.numeric())
positive_test_values <- list(names(positive_test_cases), positive_act_result, results)
negative_args <- map(negative_test_cases, "args")
# Assert
pmap(positive_test_values, ~ test_that(..1, expect_equal(..2, ..3)))
## for negative tests act can be done only with assert
map2(names(negative_test_cases), negative_args, ~ test_that(.x, expect_error(do.call(seqDist, .y))))
