# tests/testthat/test-build_gp_links.R

.local_peregrine_files <- function() {
    enh_file <- tempfile(fileext = ".tsv")
    coords_file <- tempfile()

    write.table(
        data.frame(
            enhancer = "EH1",
            gene = "gene;HGNC=11998",
            stringsAsFactors = FALSE
        ),
        file = enh_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        data.frame(
            chr = "chr1",
            start = 1050L,
            end = 1150L,
            enhancer = "EH1"
        ),
        file = coords_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    list(enh_file = enh_file, coords_file = coords_file)
}


# 1. Valid deterministic example (success case)

test_that(
    "build_gp_links returns valid outputs for local example peaks and genes",
    {
        local_files <- .local_peregrine_files()
        testthat::local_mocked_bindings(
            .map_genes_to_hgnc = function(gn) {
                data.frame(
                    hgnc_symbol = "TP53",
                    hgnc_id = "HGNC:11998",
                    stringsAsFactors = FALSE
                )
            },
            .get_peregrine_enhancer_coords_file = function() {
                local_files$coords_file
            },
            .build_promoter_links = function(pk_gr, gn, edb) {
                data.frame(
                    Peak = character(0),
                    Gene = character(0),
                    Src = character(0),
                    stringsAsFactors = FALSE
                )
            },
            .build_closest_links = function(pk_gr, gn, edb) {
                data.frame(
                    Peak = character(0),
                    Gene = character(0),
                    Src = character(0),
                    stringsAsFactors = FALSE
                )
            },
            .package = "GPlinksR"
        )

        test_peaks <- c(
            "chr1:1000-1100",
            "chr1:5000-5100"
        )

        test_genes <- "TP53"

        result <- build_gp_links(
            pk = test_peaks,
            gn = test_genes,
            enh_file = local_files$enh_file
        )

        expect_s3_class(result, "data.frame")
        expect_true(all(c("Peak", "Gene", "Src") %in% colnames(result)))
        expect_gt(nrow(result), 0)
        expect_false(any(is.na(result$Gene)))
        expect_true(all(unique(result$Src) %in% c("enh", "prom", "clo")))
        expect_true(all(result$Peak %in% test_peaks))
        expect_true(any(result$Src == "enh" & result$Gene == "TP53"))
    }
)


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
    local_files <- .local_peregrine_files()
    testthat::local_mocked_bindings(
        .map_genes_to_hgnc = function(gn) {
            stop("No HGNC IDs found for the provided gene symbols.")
        },
        .get_peregrine_enhancer_coords_file = function() {
            local_files$coords_file
        },
        .package = "GPlinksR"
    )

    expect_error(
        build_gp_links(
            pk = "chr1:1000-2000",
            gn = "THISISNOTAREALGENE",
            enh_file = local_files$enh_file
        ),
        "No HGNC IDs found"
    )
})
