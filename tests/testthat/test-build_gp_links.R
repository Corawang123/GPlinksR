test_that("build_gp_links returns valid outputs for KPMP example peaks and genes", {
  skip_on_cran()  # skip heavy online test for CRAN checks

  # ---- Prepare input (real KPMP-style examples) ----
  file_path <- get_peregrine_file(19)
  test_peaks <- c("chr1:819770-822338", "chr1:983871-984475")
  test_genes <- c("TTLL10", "PERM1")

  # ---- Run ----
  result <- build_gp_links(
    pk = test_peaks,
    gn = test_genes,
    enh_file = file_path
  )

  # ---- Validate structure ----
  expect_s3_class(result, "data.frame")
  expect_true(all(c("Peak", "Gene", "Src") %in% colnames(result)))

  # ---- Validate contents ----
  expect_gt(nrow(result), 0)
  expect_false(any(is.na(result$Gene)))
  expect_true(all(unique(result$Src) %in% c("enh", "prom", "clo")))

  # ---- Check integrity ----
  expect_true(all(result$Peak %in% test_peaks))

  # ---- Optional: check at least one enhancer link found ----
  expect_true(any(result$Src == "enh"))
})
