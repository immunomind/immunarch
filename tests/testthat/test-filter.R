# TODO: implement these functions in filters.R and use them in this test
include <- function(...) {}
exclude <- function(...) {}
lessthan <- function(...) {}
morethan <- function(...) {}
interval <- function(...) {}

data(original_data)
original_samples_count <- nrow(original_data$data)

test_cases <- list()

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1", list(
      list(meta_column = "Status", meta_value = "N")
    ))
  },
  method = "by.meta", query = list(Status = include("N")),
  expected_samples = 1
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1", list(
      list(meta_column = "Lane", meta_value = "D")
    ))
  },
  method = "by.meta", query = list(Lane = exclude("D")),
  expected_samples = original_samples_count
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1", list(
      list(meta_column = "Age", meta_value = "1")
    ))
    test_data %<>% add_mock_sample("S2", list(
      list(meta_column = "Age", meta_value = "2")
    ))
  },
  method = "by.meta", query = list(Age = lessthan(5)),
  expected_samples = 2
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1", list(
      list(meta_column = "Age", meta_value = "95")
    ))
    test_data %<>% add_mock_sample("S2", list(
      list(meta_column = "Age", meta_value = "99")
    ))
  },
  method = "by.meta", query = list(Age = interval(95, 100)),
  expected_samples = 2
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1", list(
      list(meta_column = "Age", meta_value = "100")
    ))
  },
  method = "by.meta", query = list(Age = interval(95, 100)),
  expected_samples = 0
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1", list(
      list(meta_column = "Lane", meta_value = "D")
    ))
    test_data %<>% add_mock_sample("S2", list(
      list(meta_column = "Lane", meta_value = "E")
    ))
  },
  method = "by.meta", query = list(Lane = include("D", "E")),
  expected_samples = 2
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1", list(
      list(meta_column = "Lane", meta_value = "D")
    ))
    test_data %<>% add_mock_sample("S2", list(
      list(meta_column = "Lane", meta_value = "E")
    ))
  },
  method = "by.meta", query = list(Lane = exclude("D", "E")),
  expected_samples = original_samples_count
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1", list(
      list(meta_column = "Lane", meta_value = "D"),
      list(meta_column = "Age", meta_value = 95)
    ))
    test_data %<>% add_mock_sample("S2", list(
      list(meta_column = "Lane", meta_value = "E"),
      list(meta_column = "Age", meta_value = 96)
    ))
  },
  method = "by.meta", query = list(Lane = include("D", "E"), Age = morethan(95)),
  expected_samples = 1
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1")
    test_data[["S1"]] %<>% rbind(test_data[["S1"]][rep(1, 10000), ])
    return(test_data)
  },
  method = "by.repertoire", query = list(n_clonotypes = morethan(10000)),
  expected_samples = 1
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    data(test_data)
    test_data %<>% add_mock_sample("S1")
    test_data$data[["S1"]] <- test_data$data[["S1"]][1, ]
    test_data$data[["S1"]][["Clones"]] <- 50
    test_data %<>% add_mock_sample("S2")
    test_data$data[["S2"]] <- test_data$data[["S2"]][1:2, ]
    test_data$data[["S2"]][1, ][["Clones"]] <- 50
    test_data$data[["S2"]][2, ][["Clones"]] <- 50
    return(test_data)
  },
  method = "by.repertoire", query = list(n_clones = lessthan(100)),
  expected_samples = 1
)

# repeat the last 2 test cases, but abbreviate method as "by.rep"
for (i in 1:2) {
  test_cases[[length(test_cases) + 1]] <- test_cases[[length(test_cases) - 1]]
  test_cases[[length(test_cases)]][["method"]] <- "by.rep"
}

# test_cases[[length(test_cases) + 1]] <- list(
#   data_factory = function() {
#     data(test_data)
#     return(test_data)
#   },
#   method = "by.clonotype", query = list(CDR3.aa = exclude("partial", "out_of_frame")),
#   expected_samples = 12
# )

#   list(
#     data_factory = function() {
#       data(test_data)
#       return(test_data)
#     },
#     method = "by.clonotype", query = list(Clones = interval(150, 200)),
#     expected_samples = 4
#   ),
#   list(
#     data_factory = function() {
#       data(test_data)
#       return(test_data)
#     },
#     method = "by.cl", query = list(CDR3.aa = exclude("partial", "out_of_frame")),
#     expected_samples = 12
#   ),
#   list(
#     data_factory = function() {
#       data(test_data)
#       return(test_data)
#     },
#     method = "by.cl", query = list(Clones = interval(150, 200)),
#     expected_samples = 4
#   ),
#   list(
#     data_factory = function() {
#       data(test_data)
#       return(test_data)
#     },
#     method = "by.gene", query = list(V = exclude("TRBV1", "TRGV11")),
#     expected_samples = 12
#   ),
#   list(
#     data_factory = function() {
#       data(test_data)
#       return(test_data)
#     },
#     method = "by.gene", query = list(V = include("TRAV9"), J = include("TRAJ11")),
#     expected_samples = 0
#   )
# )

for (i in seq_along(test_cases)) {
  # Arrange
  data_factory <- test_cases[[i]][["data_factory"]]
  method <- test_cases[[i]][["method"]]
  query <- test_cases[[i]][["query"]]
  expected_samples <- test_cases[[i]][["expected_samples"]]
  test_name <- paste0("method:", method, ".case:", i)
  immdata <- data_factory()
  frame_with_meta <- immdata
  table_with_meta <- list(data = lapply(immdata$data, as.data.table), meta = immdata$meta)

  # TODO: implement all filters and remove skip()
  skip("Not implemented")

  # Act
  compute_res <- apply_DF_DT(frame_with_meta, table_with_meta,
    repFilter,
    .method = method,
    .query = query
  )

  # Assert
  test_that(ap(test_name, "compute_"), {
    expect_equal(compute_res[[1]]$data %>% length(), expected_samples)
    expect_equal(compute_res[[1]]$meta %>% nrow(), expected_samples)
    expect_equal(lapply(compute_res[[1]]$data, as.data.table), compute_res[[2]]$data)
  })
}
