data(immdata)

dist_test <- seqDist(immdata$data[1])


# Arrange
negative_test_cases <- list(
  "Wrong data" = list(
    args = list(
      .data = immdata[1:2],
      .dist = dist_test
    )
  ),
  "No matching col" = list(
    args = list(
      .data = immdata$data$`A2-i131` %>% mutate("CDR3.nt" = NULL),
      .dist = dist_test
    )
  ),
  "Wrong samples" = list(
    args = list(
      .data = immdata$data[3],
      .dist = dist_test
    )
  ),
  "Multiple thresholds" = list(
    args = list(
      .data = immdata$data[1],
      .dist = dist_test,
      .fixed_threshold = 10,
      .nt_similarity = 3
    )
  ),
  "No thresholds" = list(
    args = list(
      .data = immdata$data[1],
      .dist = dist_test,
      .fixed_threshold = NULL
    )
  )
)

# Act

negative_args <- map(negative_test_cases, "args")

# Assert

map2(names(negative_test_cases), negative_args, ~ test_that(.x, expect_error(do.call(seqCluster, .y))))
