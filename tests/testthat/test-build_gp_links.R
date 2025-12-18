context("build_gp_links: valid usage and user-facing error messages")

# 1. Valid KPMP-style example (success case)

test_that("build_gp_links returns valid outputs for KPMP example peaks and genes", {

  skip_on_cran()  # skip heavy online test for CRAN checks

  # ---- Prepare input (real KPMP-style examples) ----
  file_path <- get_peregrine_file(19)

  test_peaks <- c(
    "chr1:819770-822338",
    "chr1:983871-984475"
  )

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



# 2. Missing gene input (gn)


test_that("build_gp_links errors when gn is missing", {

  expect_error(
    build_gp_links(pk = "chr1:1000-2000"),
    "'gn' must be provided"
  )
})



# 3. gn is not a character vector


test_that("build_gp_links errors when gn is not a character vector", {

  expect_error(
    build_gp_links(
      pk = "chr1:1000-2000",
      gn = 123
    ),
    "'gn' must be a character vector"
  )
})



# 4. gn contains NA or empty strings


test_that("build_gp_links errors when gn contains NA or empty strings", {

  expect_error(
    build_gp_links(
      pk = "chr1:1000-2000",
      gn = c("TP53", NA)
    ),
    "contains missing or empty entries"
  )

  expect_error(
    build_gp_links(
      pk = "chr1:1000-2000",
      gn = c("TP53", "")
    ),
    "contains missing or empty entries"
  )
})


# 5. Missing peak input (pk)


test_that("build_gp_links errors when pk is missing", {

  expect_error(
    build_gp_links(gn = "TP53"),
    "'pk' must be provided"
  )
})



# 6. pk is wrong type


test_that("build_gp_links errors when pk is not character or data.frame", {

  expect_error(
    build_gp_links(
      pk = 123,
      gn = "TP53"
    ),
    "'pk' must be a character vector or a data.frame"
  )
})


# 7. Invalid peak coordinate format


test_that("build_gp_links errors for malformed peak strings", {

  expect_error(
    build_gp_links(
      pk = "1:1000-2000",
      gn = "TP53"
    ),
    "Invalid peak format"
  )
})


# 8. pk data.frame missing required columns

test_that("build_gp_links errors when pk data.frame lacks required columns", {

  bad_pk <- data.frame(a = 1, b = 2)

  expect_error(
    build_gp_links(
      pk = bad_pk,
      gn = "TP53"
    ),
    "Peak data.frame must contain"
  )
})


# 9. Invalid PEREGRINE enhancer version


test_that("build_gp_links errors for invalid PEREGRINE version", {

  expect_error(
    build_gp_links(
      pk = "chr1:1000-2000",
      gn = "TP53",
      ver = 16
    ),
    "'ver' must be one of 17, 18, or 19"
  )
})


# 10. Nonexistent enhancer file

test_that("build_gp_links errors when enh_file does not exist", {

  expect_error(
    build_gp_links(
      pk = "chr1:1000-2000",
      gn = "TP53",
      enh_file = "not_a_real_file.tsv"
    ),
    "does not exist"
  )
})


# 11. Gene symbols cannot be mapped to HGNC


test_that("build_gp_links errors when gene symbols cannot be mapped to HGNC", {

  skip_on_cran()  # biomaRt lookup

  expect_error(
    build_gp_links(
      pk = "chr1:1000-2000",
      gn = "THISISNOTAREALGENE"
    ),
    "No HGNC IDs found"
  )
})
