# TODO: implement these functions in filters.R and use them in this test
include <- function(...) {}
exclude <- function(...) {}
lessthan <- function(...) {}
morethan <- function(...) {}
interval <- function(...) {}

# TODO: implement filters, check number of samples after filtering the test data and fill expected_samples
test_cases <- list(
  list(
    method = "by.meta", query = list(Status = include("C")),
    expected_samples = 0
  ),
  list(
    method = "by.meta", method = "by.meta", query = list(Lane = exclude("A")),
    expected_samples = 0
  ),
  list(
    method = "by.meta", query = list(Age = lessthan(23)),
    expected_samples = 0
  ),
  list(
    method = "by.meta", query = list(Age = interval(15, 23)),
    expected_samples = 0
  ),
  list(
    method = "by.meta", query = list(Lane = include("B", "C")),
    expected_samples = 0
  ),
  list(
    method = "by.meta", query = list(Lane = exclude("A", "C")),
    expected_samples = 0
  ),
  list(
    method = "by.meta", query = list(Lane = include("A", "B"), Age = morethan(15)),
    expected_samples = 0
  ),
  list(
    method = "by.repertoire", query = list(n_clonotypes = morethan(1000)),
    expected_samples = 0
  ),
  list(
    method = "by.repertoire", query = list(n_clones = lessthan(100)),
    expected_samples = 0
  ),
  list(
    method = "by.rep", query = list(n_clonotypes = morethan(1000)),
    expected_samples = 0
  ),
  list(
    method = "by.rep", query = list(n_clones = lessthan(100)),
    expected_samples = 0
  ),
  list(
    method = "by.clonotype", query = list(CDR3.aa = exclude("partial", "out_of_frame")),
    expected_samples = 0
  ),
  list(
    method = "by.clonotype", query = list(Clones = interval(5, 50)),
    expected_samples = 0
  ),
  list(
    method = "by.cl", query = list(CDR3.aa = exclude("partial", "out_of_frame")),
    expected_samples = 0
  ),
  list(
    method = "by.cl", query = list(Clones = interval(5, 50)),
    expected_samples = 0
  ),
  list(
    method = "by.gene", query = list(V = exclude("TRBV1", "TRGV11")),
    expected_samples = 0
  ),
  list(
    method = "by.gene", query = list(V = include("TRAV9"), J = include("TRAJ11")),
    expected_samples = 0
  )
)

for (i in seq_along(test_cases)) {
  method <- test_cases[[i]][["method"]]
  query <- test_cases[[i]][["query"]]
  expected_samples <- test_cases[[i]][["expected_samples"]]
  test_name <- paste0("method:", method, ".case:", i)

  # TODO: implement all filters and remove skip()
  skip("Not implemented")

  compute_res <- apply_DF_DT(frame_with_meta, table_with_meta,
    repFilter,
    .method = method,
    .query = query
  )

  test_that(ap(test_name, "compute_"), {
    expect_equal(compute_res[[1]]$data %>% length(), expected_samples)
    expect_equal(compute_res[[1]]$meta %>% nrow(), expected_samples)
    expect_equal(lapply(compute_res[[1]]$data, as.data.table), compute_res[[2]]$data)
  })
}
