data(bcrdata)
test_bcr_data <- bcrdata$data %>% top(1000)

positive_test_cases <- list(
  "Not empty result" = list(
    args = list(
      .data = test_bcr_data,
      .threads = 1
    ),
    assert_function = function(result) {
      expect_equal(immunarch:::has_no_data(result), FALSE)
      expect_equal(nrow(result[["full_clones"]]) > 0, TRUE)
    }
  ),
  "Multiple threads" = list(
    args = list(
      .data = test_bcr_data
    ),
    assert_function = function(result) {
      expect_equal(immunarch:::has_no_data(result), FALSE)
      expect_equal(nrow(result[["full_clones"]]) > 0, TRUE)
    }
  ),
  "Dataframe only" = list(
    args = list(
      .data = test_bcr_data[["full_clones"]],
      .threads = 1
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
  result <- suppressWarnings(do.call(repGermline, args))

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
  "Missing column" = list(
    args = list(
      .data = subset(test_bcr_data[["full_clones"]], select = -c(FR1.nt)),
      .threads = 1
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
    expect_error(suppressWarnings(do.call(repGermline, args)))
  )
}
