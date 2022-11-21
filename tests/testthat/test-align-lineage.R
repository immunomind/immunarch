data(bcrdata)
test_bcr_data <- bcrdata$data %>% top(1000)
test_input <- test_bcr_data %>%
  seqCluster(seqDist(test_bcr_data), .fixed_threshold = 3) %>%
  repGermline(.threads = 1) %>%
  suppressWarnings()

positive_test_cases <- list(
  "Not empty result" = list(
    args = list(
      .data = test_input,
      .min_lineage_sequences = 2,
      .prepare_threads = 1,
      .align_threads = 1
    ),
    assert_function = function(result) {
      expect_equal(immunarch:::has_no_data(result), FALSE)
      expect_equal(nrow(result[["full_clones"]]) > 0, TRUE)
    }
  ),
  "Multiple threads" = list(
    args = list(
      .data = test_input,
      .min_lineage_sequences = 2
    ),
    assert_function = function(result) {
      expect_equal(immunarch:::has_no_data(result), FALSE)
      expect_equal(nrow(result[["full_clones"]]) > 0, TRUE)
    }
  ),
  "Dataframe only" = list(
    args = list(
      .data = test_input[["full_clones"]],
      .min_lineage_sequences = 2,
      .prepare_threads = 1,
      .align_threads = 1
    ),
    assert_function = function(result) {
      expect_equal(immunarch:::has_no_data(result), FALSE)
      expect_equal(nrow(result) > 0, TRUE)
    }
  )
)

for (i in seq_along(positive_test_cases)) {
  # Arrange
  test_name <- names(positive_test_cases)[i]
  args <- positive_test_cases[[i]][["args"]]
  assert_function <- positive_test_cases[[i]][["assert_function"]]

  # Act
  result <- suppressWarnings(do.call(repAlignLineage, args))

  # Assert
  test_that(
    test_name,
    assert_function(result)
  )
}

negative_test_cases <- list(
  "List of lists" = list(
    args = list(
      .data = bcrdata
    )
  ),
  "Missing columns" = list(
    args = list(
      .data = test_bcr_data[["full_clones"]],
      .prepare_threads = 1,
      .align_threads = 1
    )
  ),
  "Missing Cluster column" = list(
    args = list(
      .data = subset(test_input[["full_clones"]], select = -c(get("Cluster"))),
      .prepare_threads = 1,
      .align_threads = 1
    )
  ),
  "Missing Germline.sequence column" = list(
    args = list(
      .data = subset(test_input[["full_clones"]], select = -c(get("Germline.sequence"))),
      .prepare_threads = 1,
      .align_threads = 1
    )
  )
)

for (i in seq_along(negative_test_cases)) {
  # Arrange
  test_name <- names(negative_test_cases)[i]
  args <- negative_test_cases[[i]][["args"]]

  # Act, Assert
  test_that(
    test_name,
    expect_error(suppressWarnings(do.call(repAlignLineage, args)))
  )
}
