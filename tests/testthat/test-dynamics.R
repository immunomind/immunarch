#'
#'
#'
#' trackClonotypes <- function (.data, .which = list(1, 15), .col = "aa", .norm = T) {
#'
#'
#'   tc = trackClonotypes(immdata$data, list(1, 10), .col = "aa")
#'   tc = trackClonotypes(immdata$data, list("MS1", 20), .col = "nt+v")
#'
#'   tc = trackClonotypes(immdata$data, c("CASRGLITDTQYF", "CSASRGSPNEQYF"), .col = "aa")
#'
#'   target = immdata$data[[1]] %>% select(CDR3.aa, V.name) %>% head(10)
#'   tc = trackClonotypes(immdata$data, target)


for (norm_val in c(TRUE, FALSE)) {
  for (method in list(list(1, 10), list(names(frame_data)[1], 20))) {
    for (column in c("nt", "aa", "nt+v", "aa+v", "nt+j", "aa+j", "aa+v+j")) {
      test_name <- paste0("method:1.norm:", norm_val, ".column:", column)

      compute_res <- apply_DF_DT(frame_data, table_data,
        trackClonotypes,
        .which = method, .col = column, .norm = norm_val
      )
      test_that(ap(test_name, "compute_"), {
        expect_equal(compute_res[[1]], compute_res[[2]])
      })
    }
  }

  target <- c("CASSLEETQYF", "CASSDSSGGANEQFF", "CASSLQETQYF", "CASSLDRETQYF", "CASSPGGGNQPQHF")
  column <- "aa"
  test_name <- paste0("methid:2.norm:", norm_val, ".column:", column)

  compute_res <- apply_DF_DT(frame_data, table_data,
    trackClonotypes,
    .which = target, .col = column, .norm = norm_val
  )
  test_that(ap(test_name, "compute_"), {
    expect_equal(compute_res[[1]], compute_res[[2]])
  })

  method_list <- list()
  method_list[["nt+v"]] <- immdata$data[[1]] %>%
    select(CDR3.nt, V.name) %>%
    head(10)
  method_list[["aa+v"]] <- immdata$data[[1]] %>%
    select(CDR3.aa, V.name) %>%
    head(15)
  method_list[["nt+v+j"]] <- immdata$data[[1]] %>%
    select(CDR3.nt, V.name, J.name) %>%
    head(10)
  method_list[["aa"]] <- immdata$data[[1]] %>%
    select(CDR3.aa) %>%
    head(10)
  method_list[["nt"]] <- immdata$data[[1]] %>%
    select(CDR3.nt) %>%
    head(10)

  for (method_name in names(method_list)) {
    test_name <- paste0("methid:3.norm:", norm_val, ".which:", method_name)

    compute_res <- apply_DF_DT(frame_data, table_data,
      trackClonotypes,
      .which = method_list[[method_name]], .col = column, .norm = norm_val
    )
    test_that(ap(test_name, "compute_"), {
      expect_equal(compute_res[[1]], compute_res[[2]])
    })
  }
}
