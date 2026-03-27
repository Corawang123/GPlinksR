#' Build Gene-Peak Links from MAE or SCE Input
#'
#' Wrapper around [build_gp_links()] that extracts peak coordinates and gene
#' symbols from a `MultiAssayExperiment` or `SingleCellExperiment`, then runs
#' the standard GPlinksR workflow.
#'
#' For `MultiAssayExperiment` input, supply the experiment names containing peak
#' and gene features.
#'
#' For `SingleCellExperiment` input, the wrapper can use the main experiment and
#' one `altExp()` to represent the gene and peak feature spaces. By default, the
#' main experiment is treated as genes and `peak_experiment` is taken from an
#' `altExp()` name.
#'
#' Peak features are extracted in this order:
#' 1. `rowRanges()`
#' 2. `rowData()` column `PeakRegion`
#' 3. `rowData()` columns `chr`, `start`, `end`
#' 4. row names in `"chr:start-end"` format
#'
#' Gene symbols are extracted in this order:
#' 1. `gene_col`, if supplied
#' 2. `rowData()` columns `symbol`, `gene_name`, `hgnc_symbol`,
#'    `external_gene_name`, or `Gene`
#' 3. row names
#'
#' @param x A `MultiAssayExperiment`, `SingleCellExperiment`,
#'   `SummarizedExperiment`, or `RangedSummarizedExperiment`.
#' @param peak_experiment For `MultiAssayExperiment`, the experiment name
#'   holding peak features. For `SingleCellExperiment`, the `altExp()` name
#'   holding peak features.
#' @param gene_experiment For `MultiAssayExperiment`, the experiment name
#'   holding gene features. Ignored for `SingleCellExperiment`, where the main
#'   experiment is used for genes.
#' @param peak_col Optional `rowData()` column containing peak coordinates in
#'   `"chr:start-end"` format.
#' @param gene_col Optional `rowData()` column containing gene symbols.
#' @param ver Integer (17, 18, or 19). PEREGRINE enhancer-gene link version.
#' @param enh_file Optional local path to enhancer-gene link file (`.tsv`).
#'
#' @return A data.frame with columns returned by [build_gp_links()].
#' @export
build_gp_links_wrapper <- function(
    x,
    peak_experiment = NULL,
    gene_experiment = NULL,
    peak_col = NULL,
    gene_col = NULL,
    ver = 19,
    enh_file = NULL
) {
    inputs <- .extract_wrapper_inputs(
        x = x,
        peak_experiment = peak_experiment,
        gene_experiment = gene_experiment,
        peak_col = peak_col,
        gene_col = gene_col
    )

    build_gp_links(
        pk = inputs$pk,
        gn = inputs$gn,
        ver = ver,
        enh_file = enh_file
    )
}


.extract_wrapper_inputs <- function(x,
                                    peak_experiment = NULL,
                                    gene_experiment = NULL,
                                    peak_col = NULL,
                                    gene_col = NULL) {
    if (methods::is(x, "MultiAssayExperiment")) {
        return(.extract_mae_inputs(
            x = x,
            peak_experiment = peak_experiment,
            gene_experiment = gene_experiment,
            peak_col = peak_col,
            gene_col = gene_col
        ))
    }

    if (methods::is(x, "SingleCellExperiment")) {
        return(.extract_sce_inputs(
            x = x,
            peak_experiment = peak_experiment,
            peak_col = peak_col,
            gene_col = gene_col
        ))
    }

    if (methods::is(x, "SummarizedExperiment") ||
        methods::is(x, "RangedSummarizedExperiment")) {
        return(.extract_se_inputs(
            x = x,
            peak_col = peak_col,
            gene_col = gene_col
        ))
    }

    stop(.wrapper_class_message())
}


.extract_mae_inputs <- function(x,
                                peak_experiment,
                                gene_experiment,
                                peak_col,
                                gene_col) {
    exp_list <- MultiAssayExperiment::experiments(x)
    exp_names <- names(exp_list)
    available <- paste(exp_names, collapse = ", ")

    if (is.null(peak_experiment) || !peak_experiment %in% exp_names) {
        msg <- paste0(
            "For MultiAssayExperiment input, 'peak_experiment' must ",
            "match one of: ",
            available
        )
        stop(msg)
    }

    if (is.null(gene_experiment) || !gene_experiment %in% exp_names) {
        msg <- paste0(
            "For MultiAssayExperiment input, 'gene_experiment' must ",
            "match one of: ",
            available
        )
        stop(msg)
    }

    list(
        pk = .extract_peak_regions_wrapper(
            exp_list[[peak_experiment]],
            peak_col = peak_col
        ),
        gn = .extract_gene_symbols_wrapper(
            exp_list[[gene_experiment]],
            gene_col = gene_col
        )
    )
}


.extract_sce_inputs <- function(x, peak_experiment, peak_col, gene_col) {
    alt_names <- SingleCellExperiment::altExpNames(x)
    available_alt <- paste(alt_names, collapse = ", ")

    if (is.null(peak_experiment)) {
        if (length(alt_names) == 1) {
            peak_experiment <- alt_names[[1]]
        } else {
            msg <- paste0(
                "For SingleCellExperiment input, provide 'peak_experiment' as ",
                "an altExp() name. Available altExp names: ", available_alt
            )
            stop(msg)
        }
    }

    if (!peak_experiment %in% alt_names) {
        msg <- paste0(
            "'peak_experiment' was not found among altExp() names: ",
            available_alt
        )
        stop(msg)
    }

    list(
        pk = .extract_peak_regions_wrapper(
            SingleCellExperiment::altExp(x, peak_experiment),
            peak_col = peak_col
        ),
        gn = .extract_gene_symbols_wrapper(x, gene_col = gene_col)
    )
}


.extract_se_inputs <- function(x, peak_col, gene_col) {
    if (is.null(peak_col)) {
        msg <- paste0(
            "For SummarizedExperiment input, supply 'peak_col' if the object ",
            "stores peaks, or use a SingleCellExperiment/MultiAssayExperiment ",
            "so genes and peaks can both be extracted."
        )
        stop(msg)
    }

    list(
        pk = .extract_peak_regions_wrapper(x, peak_col = peak_col),
        gn = .extract_gene_symbols_wrapper(x, gene_col = gene_col)
    )
}


.wrapper_class_message <- function() {
    paste0(
        "'x' must be a MultiAssayExperiment, SingleCellExperiment, ",
        "SummarizedExperiment, or RangedSummarizedExperiment."
    )
}


.extract_peak_regions_wrapper <- function(x, peak_col = NULL) {
    pk <- .peak_regions_from_ranges(x)
    if (length(pk) == 0) {
        pk <- .peak_regions_from_metadata(x = x, peak_col = peak_col)
    }

    pk <- unique(stats::na.omit(as.character(pk)))
    pk <- pk[nzchar(pk)]

    if (length(pk) == 0) {
        stop("Could not extract peak coordinates from the supplied object.")
    }

    .validate_wrapper_peak_regions(pk)
    pk
}


.peak_regions_from_ranges <- function(x) {
    rr <- tryCatch(
        SummarizedExperiment::rowRanges(x),
        error = function(e) NULL
    )

    if (is.null(rr) || length(rr) == 0) {
        return(character(0))
    }

    paste0(
        as.character(GenomicRanges::seqnames(rr)),
        ":",
        GenomicRanges::start(rr),
        "-",
        GenomicRanges::end(rr)
    )
}


.peak_regions_from_metadata <- function(x, peak_col = NULL) {
    rd <- tryCatch(
        SummarizedExperiment::rowData(x),
        error = function(e) NULL
    )
    rd_names <- if (is.null(rd)) character(0) else colnames(rd)

    if (!is.null(peak_col)) {
        if (!peak_col %in% rd_names) {
            available_cols <- paste(rd_names, collapse = ", ")
            msg <- paste0(
                "'peak_col' was not found in rowData(). Available columns: ",
                available_cols
            )
            stop(msg)
        }
        return(as.character(rd[[peak_col]]))
    }

    if ("PeakRegion" %in% rd_names) {
        return(as.character(rd[["PeakRegion"]]))
    }

    if (all(c("chr", "start", "end") %in% rd_names)) {
        return(paste0(rd[["chr"]], ":", rd[["start"]], "-", rd[["end"]]))
    }

    rownames(x)
}


.validate_wrapper_peak_regions <- function(pk) {
    is_valid <- grepl("^chr[0-9XYMT]+:[0-9]+-[0-9]+$", pk)
    if (all(is_valid)) {
        return(invisible(NULL))
    }

    bad <- pk[!is_valid][1]
    msg <- paste0(
        "Extracted peak coordinates are not in 'chr:start-end' format. ",
        "First bad value: ", bad
    )
    stop(msg)
}


.extract_gene_symbols_wrapper <- function(x, gene_col = NULL) {
    rd <- tryCatch(
        SummarizedExperiment::rowData(x),
        error = function(e) NULL
    )

    rd_names <- if (is.null(rd)) character(0) else colnames(rd)
    candidate_cols <- c(
        "symbol",
        "gene_name",
        "hgnc_symbol",
        "external_gene_name",
        "Gene"
    )

    if (!is.null(gene_col)) {
        if (!gene_col %in% rd_names) {
            available_cols <- paste(rd_names, collapse = ", ")
            msg <- paste0(
                "'gene_col' was not found in rowData(). Available columns: ",
                available_cols
            )
            stop(msg)
        }

        gn <- as.character(rd[[gene_col]])
    } else {
        matched_col <- candidate_cols[candidate_cols %in% rd_names][1]

        if (!is.na(matched_col)) {
            gn <- as.character(rd[[matched_col]])
        } else {
            gn <- rownames(x)
        }
    }

    gn <- unique(stats::na.omit(as.character(gn)))
    gn <- gn[nzchar(gn)]

    if (length(gn) == 0) {
        stop("Could not extract gene symbols from the supplied object.")
    }

    gn
}
