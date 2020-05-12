for (coding_val in c(TRUE, FALSE)) {
  for (k_size in c(1, 3, 5, 10)) {
    for (column in c("nt", "aa")) {
      test_name <- paste0("kmer:", k_size, ":", column, ".coding:", coding_val)

      compute_res <- apply_DF_DT(frame_data[[1]], table_data[[1]],
        getKmers,
        .k = k_size, .col = column, .coding = coding_val
      )
      test_that(ap(test_name, "compute_"), {
        expect_equal(compute_res[[1]], compute_res[[2]])
      })

      # vis_res = vis_results(compute_res)
      # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })

      if (column == "aa") {
        for (method in c("freq", "prob", "wei", "self")) {
          for (remove_stop in c(TRUE, FALSE)) {
            test_name <- paste0("kmer:", k_size, ".method:", method, ".stop:", remove_stop, ".coding:", coding_val)

            profile_res1 <- kmer_profile(compute_res[[1]], .method = method, .remove.stop = remove_stop)
            profile_res2 <- kmer_profile(compute_res[[2]], .method = method, .remove.stop = remove_stop)

            test_that(ap(test_name, "compute_"), {
              expect_equal(profile_res1, profile_res2)
            })

            # vis_res = vis_results(compute_res)
            # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
          }
        }
      }
    }
  }
}
