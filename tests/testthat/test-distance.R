context("Distance computing base tests")
data(immdata)
short_immdata <- immdata$data
short_immdata <- map(short_immdata, ~ .x %>% head(3)) # less sample size saves time in computations

test_that("Input response is sanitized", {
  test_cases <- list(
    list(),
    list(data.frame()),
    list(data.table()),
    c("a", "b", "c")
  )
  for (i in test_cases) {
    i %>%
      seqDist() %>%
      expect_error()
  }
})

test_that("Argument values are sanitized", {
  args_cases <- list(
    list(.col = "aaa"),
    list(.method = "aaa")
  )

  for (i in args_cases) {
    i[".data"] <- short_immdata[1]
    do.call(seqDist, i) %>% expect_error()
  }
})

test_that("Return values are correct", {
  args_case <- list(.method = "lv", .data = short_immdata[1]) # it is easier to test without Inf values with lv method
  result <- do.call(seqDist, args_case)
  result %>% expect_type("list")
  result[[1]] %>% expect_type("double")
  result[[1]] %>%
    as.numeric() %>%
    expect_equal(c(18, 20, 14))
})

test_that("Custom function returning value is correct", {
  f <- function(x, y) {
    res <- matrix(nrow = length(x), ncol = length(y))
    for (i in 1:length(x)) {
      res[i, ] <- abs(nchar(x[i]) - nchar(y))
    }
    dimnames(res) <- list(x, y)
    return(as.dist(res))
  }
  result <- seqDist(short_immdata[1], .method = f)
  result %>% expect_type("list")
  result[[1]] %>% expect_type("integer")
  result[[1]] %>%
    as.numeric() %>%
    expect_equal(c(6, 3, 9))
})
