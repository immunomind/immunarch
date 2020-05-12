for (column in c("nt", "aa", "nt+v", "aa+v", "nt+j", "aa+j", "aa+v+j")) {
  for (quant in c("count", "prop")) {
    for (coding_val in c(TRUE, FALSE)) {
      for (min_samples in c(1, 3)) {
        for (max_samples in c(5, NA)) {
          test_name <- paste0("column:", column, ".quant:", quant, ".coding:", coding_val, ".min:", min_samples, ".max:", max_samples)

          compute_res <- apply_DF_DT(frame_data, table_data,
            pubRep,
            .col = column, .quant = quant, .coding = coding_val, .min.samples = min_samples, .max.samples = max_samples, .verbose = F
          )
          test_that(ap(test_name, "compute_"), {
            expect_equal(compute_res[[1]], compute_res[[2]])
          })

          # vis_res = vis_results(compute_res)
          # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
        }
      }
    }
  }
}
