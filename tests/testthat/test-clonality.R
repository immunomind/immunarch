for (method in c("clonal.prop", "homeo", "top", "tail")) {
  test_name <- paste0("method:", method)

  compute_res <- apply_DF_DT(frame_data, table_data,
    repClonality,
    .method = method
  )
  test_that(ap(test_name, "compute_"), {
    expect_equal(compute_res[[1]], compute_res[[2]])
  })

  # vis_res = vis_results(compute_res)
  # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
}
