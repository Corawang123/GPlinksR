#' Download Enhancer-Gene Link File from PANTHER Peregrine Database
#'
#' @param version Integer (17, 18, or 19). Which PANTHER enhancer-gene link
#'   version to download.
#' @return The full path to the downloaded `.tsv` file stored in the
#'   `BiocFileCache`.
#'
#' @export
#'
#' @examples
#' if (interactive()) {
#'     f <- get_peregrine_file(19)
#'     f
#'     readLines(f, n = 3)
#' }
get_peregrine_file <- function(version = 19) {
    if (!requireNamespace("BiocFileCache", quietly = TRUE)) {
        stop("Package 'BiocFileCache' must be installed.")
    }

    if (!version %in% c(17, 18, 19)) {
        stop("Version must be 17, 18, or 19.")
    }

    url <- paste0(
        "https://data.pantherdb.org/ftp/peregrine_data/",
        "enhancer_gene_link_", version, ".tsv"
    )

    bfc <- BiocFileCache::BiocFileCache(ask = FALSE)

    rname <- paste0("peregrine_v", version)

    query <- BiocFileCache::bfcquery(
        bfc,
        rname,
        field = "rname"
    )

    if (nrow(query) == 0) {
        message("Downloading PEREGRINE enhancer file...")

        rid <- names(
            BiocFileCache::bfcadd(
                bfc,
                rname = rname,
                fpath = url
            )
        )
    } else {
        message("Using cached PEREGRINE file.")
        rid <- query$rid[1]
    }

    BiocFileCache::bfcpath(bfc, rid)
}


.get_peregrine_enhancer_coords_file <- function() {
    if (!requireNamespace("BiocFileCache", quietly = TRUE)) {
        stop("Package 'BiocFileCache' must be installed.")
    }

    url <- paste0(
        "https://data.pantherdb.org/ftp/peregrine_data/",
        "PEREGRINEenhancershg38"
    )
    bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
    rname <- "peregrine_enhancers_hg38"

    query <- BiocFileCache::bfcquery(
        bfc,
        rname,
        field = "rname"
    )

    if (nrow(query) == 0) {
        message("Downloading PEREGRINE enhancer coordinates...")
        rid <- names(
            BiocFileCache::bfcadd(
                bfc,
                rname = rname,
                fpath = url
            )
        )
    } else {
        message("Using cached PEREGRINE enhancer coordinates.")
        rid <- query$rid[1]
    }

    BiocFileCache::bfcpath(bfc, rid)
}
