data(bcrdata)
test_bcr_data <- bcrdata$data %>% top(1000)
test_clusters <- seqCluster(test_bcr_data, seqDist(test_bcr_data), .fixed_threshold = 3)

positive_test_cases <- list(

)

for (i in seq_along(positive_test_cases)) {
  # Arrange

  # Act

  # Assert

}

negative_test_cases <- list(
  "List of lists" = list(
    args = list(
      .data = bcrdata
    )
  ),
  "Missing column" = list(
    args = list(
      .data =  subset(test_clusters[["full_clones"]], select = -c(FR1.nt))
    )
  )
)

for (i in seq_along(negative_test_cases)) {
  # Arrange
  test_name <- names(negative_test_cases)[i]
  args <- negative_test_cases[[i]][["args"]]

  # Act, Assert
  test_that(test_name, {
    expect_error(suppressWarnings(do.call(repGermline, args)))
  })
}
