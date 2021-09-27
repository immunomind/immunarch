test_cases <- list()

original_samples_count <- nrow(immdata$data)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", list(Status = "N"))
  },
  method = "by.meta", query = list(Status = include("N")),
  expected_samples = 1
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", list(Lane = "D"))
  },
  method = "by.meta", query = list(Lane = exclude("D")),
  expected_samples = original_samples_count
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", list(Age = 1))
    immdata %<>% add_mock_sample("S2", list(Age = 2))
  },
  method = "by.meta", query = list(Age = lessthan(5)),
  expected_samples = 2
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", list(Age = 95))
    immdata %<>% add_mock_sample("S2", list(Age = 99))
  },
  method = "by.meta", query = list(Age = interval(95, 100)),
  expected_samples = 2
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", list(Age = 100))
  },
  method = "by.meta", query = list(Age = interval(95, 100)),
  expected_samples = 0
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", list(Lane = "D"))
    immdata %<>% add_mock_sample("S2", list(Lane = "E"))
  },
  method = "by.meta", query = list(Lane = include("D", "E")),
  expected_samples = 2
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", list(Lane = "D"))
    immdata %<>% add_mock_sample("S2", list(Lane = "E"))
  },
  method = "by.meta", query = list(Lane = exclude("D", "E")),
  expected_samples = original_samples_count
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", list(Lane = "D", Age = 95))
    immdata %<>% add_mock_sample("S2", list(Lane = "E", Age = 96))
  },
  method = "by.meta", query = list(Lane = include("D", "E"), Age = morethan(95)),
  expected_samples = 1
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1")
    immdata[["S1"]] %<>% rbind(immdata[["S1"]][rep(1, 10000), ])
    return(immdata)
  },
  method = "by.repertoire", query = list(n_clonotypes = morethan(10000)),
  expected_samples = 1
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1")
    immdata$data[["S1"]] <- immdata$data[["S1"]][1, ]
    immdata$data[["S1"]][["Clones"]] <- 50
    immdata %<>% add_mock_sample("S2")
    immdata$data[["S2"]] <- immdata$data[["S2"]][1:2, ]
    immdata$data[["S2"]][1, ][["Clones"]] <- 50
    immdata$data[["S2"]][2, ][["Clones"]] <- 50
    return(immdata)
  },
  method = "by.repertoire", query = list(n_clones = lessthan(100)),
  expected_samples = 1
)

# repeat the last 2 test cases, but abbreviate method as "by.rep"
for (i in 1:2) {
  test_cases[[length(test_cases) + 1]] <- test_cases[[length(test_cases) - 1]]
  test_cases[[length(test_cases)]][["method"]] <- "by.rep"
}

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", .empty = TRUE)
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(CDR3.aa = "partial")))
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(CDR3.aa = "out_of_frame")))
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(CDR3.aa = "other")))
    return(immdata)
  },
  method = "by.clonotype", query = list(CDR3.aa = exclude("partial", "out_of_frame")),
  expected_samples = original_samples_count + 1,
  expected_sample_rows = list(S1 = 1)
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", .empty = TRUE)
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(Clones = 1000)))
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(Clones = 1500)))
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(Clones = 2000)))
    return(immdata)
  },
  method = "by.clonotype", query = list(Clones = interval(1000, 2000)),
  expected_samples = 1,
  expected_sample_rows = list(S1 = 2)
)

mock_genes <- function() {
  immdata %<>% add_mock_sample("S1", .empty = TRUE)
  # delete all other samples
  immdata$data <- immdata$data[names(immdata$data) == "S1"]
  immdata$meta %<>% filter(Sample == "S1")

  immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(V.name = "TRBV1")))
  immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(V.name = "TRAV1")))
  immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(V.name = "TRAV2")))
  immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(V.name = "TRAV11")))
  return(immdata)
}

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = mock_genes,
  method = "by.clonotype", query = list(V.name = exclude("TRBV1", "TRAV1")),
  expected_samples = 1,
  expected_sample_rows = list(S1 = 2)
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = mock_genes,
  method = "by.clonotype", query = list(V.name = exclude("TRBV1", "TRAV1")), match = "exact",
  expected_samples = 1,
  expected_sample_rows = list(S1 = 2)
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = mock_genes,
  method = "by.clonotype", query = list(V.name = exclude("TRBV1", "TRAV1")), match = "startswith",
  expected_samples = 1,
  expected_sample_rows = list(S1 = 1)
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = mock_genes,
  method = "by.clonotype", query = list(V.name = include("TRBV1", "TRAV1")), match = "startswith",
  expected_samples = 1,
  expected_sample_rows = list(S1 = 3)
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = mock_genes,
  method = "by.clonotype", query = list(V.name = exclude("V1")), match = "substring",
  expected_samples = 1,
  expected_sample_rows = list(S1 = 1)
)

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = mock_genes,
  method = "by.clonotype", query = list(V.name = include("V1")), match = "substring",
  expected_samples = 1,
  expected_sample_rows = list(S1 = 3)
)

# repeat the last 8 test cases, but abbreviate method as "by.cl"
for (i in 1:8) {
  test_cases[[length(test_cases) + 1]] <- test_cases[[length(test_cases) - 7]]
  test_cases[[length(test_cases)]][["method"]] <- "by.cl"
}

test_cases[[length(test_cases) + 1]] <- list(
  data_factory = function() {
    immdata %<>% add_mock_sample("S1", .empty = TRUE)
    immdata$data <- immdata$data[names(immdata$data) == "S1"]
    immdata$meta %<>% filter(Sample == "S1")
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(V.name = "TRBV1", J.name = "TRAJ1")))
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(V.name = "TRAV1", J.name = "TRAJ11")))
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(V.name = "TRAV2", J.name = "TRAJ1")))
    immdata$data[["S1"]] %<>% bind_rows(as.data.frame(list(V.name = "TRAV11", J.name = "TRBJ2")))
    return(immdata)
  },
  method = "by.clonotype",
  query = list(V.name = include("AV1"), J.name = include("AJ")),
  match = "substring",
  expected_samples = 1,
  expected_sample_rows = list(S1 = 1)
)

for (i in seq_along(test_cases)) {
  # Arrange
  data_factory <- test_cases[[i]][["data_factory"]]
  method <- test_cases[[i]][["method"]]
  query <- test_cases[[i]][["query"]]
  match <- test_cases[[i]][["match"]]
  expected_samples <- test_cases[[i]][["expected_samples"]]
  expected_sample_rows <- test_cases[[i]][["expected_sample_rows"]]
  # expected_sample_rows is named list, contains sample names and expected rows;
  # if not specified, don't check sample rows
  if (is.null(expected_sample_rows)) {
    expected_sample_rows <- list()
  }

  test_name <- paste0("method:", method, ".case:", i)
  immdata <- data_factory()
  frame_with_meta <- immdata
  table_with_meta <- list(data = lapply(immdata$data, as.data.table), meta = immdata$meta)

  # TODO: implement all filters and remove skip()
  skip("Not implemented")

  # Act
  if (is.null(match)) {
    compute_res <- apply_DF_DT(frame_with_meta, table_with_meta,
      repFilter,
      .method = method, .query = query
    )
  } else {
    compute_res <- apply_DF_DT(frame_with_meta, table_with_meta,
      repFilter,
      .method = method, .query = query, .match = match
    )
  }

  # Assert
  test_that(ap(test_name, "compute_"), {
    expect_equal(compute_res[[1]]$data %>% length(), expected_samples)
    expect_equal(compute_res[[1]]$meta %>% nrow(), expected_samples)
    for (j in seq_along(expected_sample_rows)) {
      sample_name <- names(expected_sample_rows)[[j]]
      expected_rows <- expected_sample_rows[[j]]
      expect_equal(compute_res[[1]]$data[[sample_name]] %>% nrow(), expected_rows)
    }
    expect_equal(lapply(compute_res[[1]]$data, as.data.table), compute_res[[2]]$data)
  })
}
