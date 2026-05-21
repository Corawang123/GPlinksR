utils::globalVariables(c("Gene"))

#' Build Gene-Peak Network (Enhancer, Promoter, Closest)
#'
#' @param pk Either:
#'   (1) A character vector of genomic coordinates ("chr1:1000-1200"), or
#'   (2) A data.frame containing peak coordinates with columns
#'       "chr","start","end" or a "PeakRegion" column
#'       ("chr1:1000-1200" style).
#' @param gn Character vector of gene symbols
#'   (HGNC official symbols, e.g. "TP53", "FMNL2").
#' @param ver Integer (17, 18, or 19) - PANTHER enhancer-gene link version.
#' @param enh_file Optional local path to enhancer-gene link file (.tsv).
#'
#' @return A data.frame with columns: Peak, Gene, Src.
#'   Closest-gene mappings are based on distance from each peak to gene TSS.
#' @export
#'
#'
#' @examples
#' data("gp_example_inputs", package = "GPlinksR")
#' data("gp_example_links", package = "GPlinksR")
#' pk <- gp_example_inputs$pk[seq_len(3)]
#' gn <- gp_example_inputs$gn[seq_len(4)]
#'
#' head(gp_example_links)
#'
#' if (interactive()) {
#'     gp <- build_gp_links(pk, gn)
#'     head(gp)
#' }
build_gp_links <- function(pk, gn, ver = 19, enh_file = NULL) {
    message("Checking inputs...")
    .validate_build_gp_inputs(pk = pk, gn = gn, ver = ver, enh_file = enh_file)
    msg <- paste0(
        "Input validation passed: ", length(pk),
        " peaks and ", length(gn), " genes provided."
    )
    message(msg)
    message("Preparing peaks...")
    pk_gr <- .prepare_peak_granges(pk)
    message("Building enhancer-based links...")
    if (is.null(enh_file)) {
        enh_file <- get_peregrine_file(ver)
    }
    df_enh <- .build_enhancer_links(pk_gr = pk_gr, gn = gn, enh_file = enh_file)
    msg <- paste0("Enhancer-gene links found: ", nrow(df_enh))
    message(msg)
    message("Building promoter-based links...")
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    df_prom <- .build_promoter_links(pk_gr = pk_gr, gn = gn, edb = edb)
    msg <- paste0("Promoter-gene links found: ", nrow(df_prom))
    message(msg)
    message("Building closest-gene links...")
    df_clo <- .build_closest_links(pk_gr = pk_gr, gn = gn, edb = edb)
    msg <- paste0("Closest-gene links found: ", nrow(df_clo))
    message(msg)
    df_all <- .combine_link_tables(df_enh, df_prom, df_clo)
    msg <- paste0("Total combined links: ", nrow(df_all))
    message(msg)
    df_all
}


.validate_build_gp_inputs <- function(pk, gn, ver, enh_file) {
    .validate_gene_input(gn)
    .validate_peak_input(pk)
    if (!ver %in% c(17, 18, 19)) {
        stop("'ver' must be one of 17, 18, or 19 (PEREGRINE version numbers).")
    }

    if (!is.null(enh_file) && !file.exists(enh_file)) {
        msg <- paste0("The provided 'enh_file' path does not exist: ", enh_file)
        stop(msg)
    }
}


.validate_gene_input <- function(gn) {
    if (missing(gn) || is.null(gn) || length(gn) == 0) {
        msg <- paste0(
            "'gn' must be provided - a non-empty character vector ",
            "of gene symbols."
        )
        stop(msg)
    }

    if (!is.character(gn)) {
        stop("'gn' must be a character vector (e.g., c('TP53', 'FMNL2')).")
    }

    if (any(is.na(gn) | gn == "")) {
        msg <- paste0(
            "'gn' contains missing or empty entries. ",
            "Please remove or replace them."
        )
        stop(msg)
    }
}


.validate_peak_input <- function(pk) {
    if (missing(pk) || is.null(pk) || length(pk) == 0) {
        msg <- paste0(
            "'pk' must be provided - either a character vector ",
            "or a data.frame."
        )
        stop(msg)
    }

    if (is.character(pk)) {
        .validate_peak_strings(pk)
        return(invisible(NULL))
    }

    if (is.data.frame(pk)) {
        .validate_peak_dataframe(pk)
        return(invisible(NULL))
    }

    pk_class <- paste(class(pk), collapse = ", ")
    msg <- paste0(
        "'pk' must be a character vector or a data.frame, not ", pk_class
    )
    stop(msg)
}


.validate_peak_strings <- function(pk) {
    valid_peak <- grepl("^chr[0-9XYMT]+:[0-9]+-[0-9]+$", pk)
    if (all(valid_peak)) {
        return(invisible(NULL))
    }

    bad <- pk[!valid_peak][1]
    msg <- paste0(
        "Invalid peak format: '", bad,
        "'. Expected 'chr#:start-end' (e.g., 'chr1:1000-1200')."
    )
    stop(msg)
}


.validate_peak_dataframe <- function(pk) {
    cols <- colnames(pk)
    has_required <- "PeakRegion" %in% cols ||
        all(c("chr", "start", "end") %in% cols)

    if (has_required) {
        return(invisible(NULL))
    }

    msg <- paste0(
        "Peak data.frame must contain either 'PeakRegion' or columns ",
        "'chr','start','end'."
    )
    stop(msg)
}


.prepare_peak_granges <- function(pk) {
    if (is.data.frame(pk)) {
        pk_dt <- .peak_dataframe_to_dt(pk)
    } else if (is.character(pk)) {
        pk_dt <- .peak_vector_to_dt(pk)
    } else {
        msg <- paste0(
            "Unsupported peak input type. Provide a character vector ",
            "or data.frame."
        )
        stop(msg)
    }

    pk_dt$start <- as.integer(pk_dt$start)
    pk_dt$end <- as.integer(pk_dt$end)
    pk_dt$id <- paste0(pk_dt$chr, ":", pk_dt$start, "-", pk_dt$end)
    pk_dt$chr <- ifelse(
        !grepl("^chr", pk_dt$chr),
        paste0("chr", pk_dt$chr),
        pk_dt$chr
    )

    valid_chr <- paste0("chr", c(seq_len(22), "X", "Y", "M", "MT"))
    pk_dt <- pk_dt[pk_dt$chr %in% valid_chr, ]

    GenomicRanges::makeGRangesFromDataFrame(
        pk_dt,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = TRUE
    )
}


.peak_dataframe_to_dt <- function(pk) {
    if ("PeakRegion" %in% colnames(pk)) {
        pk_split <- data.table::tstrsplit(pk$PeakRegion, "[:-]")
        return(data.table::data.table(
            chr = pk_split[[1]],
            start = pk_split[[2]],
            end = pk_split[[3]]
        ))
    }

    if (all(c("chr", "start", "end") %in% colnames(pk))) {
        return(data.table::data.table(
            chr = pk$chr,
            start = pk$start,
            end = pk$end
        ))
    }

    msg <- paste0(
        "Peak input must contain either 'PeakRegion' or columns ",
        "'chr','start','end'."
    )
    stop(msg)
}


.peak_vector_to_dt <- function(pk) {
    pk_split <- data.table::tstrsplit(pk, "[:-]")
    data.table::data.table(
        chr = pk_split[[1]],
        start = pk_split[[2]],
        end = pk_split[[3]]
    )
}


.build_enhancer_links <- function(pk_gr, gn, enh_file) {
    enh_raw <- data.table::fread(
        enh_file,
        sep = "\t",
        header = TRUE,
        fill = TRUE,
        quote = ""
    )
    enh_raw$HGNC <- sub(".*HGNC=([0-9]+).*", "\\1", enh_raw$gene)
    enh_raw$enhancer <- as.character(enh_raw$enhancer)

    gn_map <- .map_genes_to_hgnc(gn)
    gn_ids <- sub("HGNC:", "", gn_map$hgnc_id)
    enh_sub <- enh_raw[enh_raw$HGNC %in% gn_ids, ]

    message("Loading enhancer coordinates...")
    enh_coords_file <- .get_peregrine_enhancer_coords_file()
    enh_reg <- data.table::fread(
        enh_coords_file,
        sep = "\t",
        header = FALSE,
        col.names = c("chr", "start", "end", "enhancer")
    )

    enh_full <- merge(enh_sub, enh_reg, by = "enhancer", all.x = TRUE)
    enh_full$sym <- gn_map$hgnc_symbol[match(enh_full$HGNC, gn_ids)]

    enh_gr <- GenomicRanges::makeGRangesFromDataFrame(
        enh_full,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = TRUE
    )

    enh_gr <- .align_seqlevels(enh_gr, pk_gr)

    .overlap_df(
        query_gr = pk_gr,
        subject_gr = enh_gr,
        subject_gene = S4Vectors::mcols(enh_gr)$sym,
        src = "enh"
    )
}


.map_genes_to_hgnc <- function(gn) {
    mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gn_map <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "hgnc_id"),
        filters = "hgnc_symbol",
        values = gn,
        mart = mart
    )

    if (nrow(gn_map) == 0) {
        stop("No HGNC IDs found for the provided gene symbols.")
    }

    gn_map
}


.build_promoter_links <- function(pk_gr, gn, edb) {
    tx_gr <- ensembldb::transcripts(
        edb,
        columns = c("tx_id", "gene_id", "gene_name", "seq_name")
    )

    prom_gr <- .promoter_ranges(tx_gr)
    prom_gr$symbol <- tx_gr$gene_name
    prom_gr <- prom_gr[prom_gr$symbol %in% gn]
    prom_gr <- .align_seqlevels(prom_gr, pk_gr)

    .overlap_df(
        query_gr = pk_gr,
        subject_gr = prom_gr,
        subject_gene = prom_gr$symbol,
        src = "prom",
        ignore_strand = TRUE
    )
}


.promoter_ranges <- function(tx_gr) {
    prom_gr <- withCallingHandlers(
        GenomicRanges::promoters(tx_gr, upstream = 2000, downstream = 200),
        warning = function(w) {
            if (grepl("out-of-bound", conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        }
    )
    GenomicRanges::trim(prom_gr)
}


.build_closest_links <- function(pk_gr, gn, edb) {
    tx_gene <- ensembldb::genes(
        edb,
        columns = c("gene_id", "gene_name", "seq_name")
    )
    tx_gene$symbol <- tx_gene$gene_name
    tx_gene <- tx_gene[tx_gene$symbol %in% gn]

    tx_gene <- .align_seqlevels(tx_gene, pk_gr)
    pk_gr4 <- .subset_common_seqlevels(pk_gr, tx_gene)

    pk_gr4 <- GenomicRanges::trim(pk_gr4)
    tx_gene <- GenomicRanges::trim(tx_gene)

    if (length(tx_gene) == 0 || length(pk_gr4) == 0) {
        return(.empty_link_df())
    }

    tss_gr <- GenomicRanges::resize(tx_gene, width = 1, fix = "start")
    tss_gr$symbol <- tx_gene$symbol

    nearest_hits <- GenomicRanges::nearest(pk_gr4, tss_gr, ignore.strand = TRUE)

    data.frame(
        Peak = pk_gr4$id,
        Gene = tss_gr$symbol[nearest_hits],
        Src = "clo",
        stringsAsFactors = FALSE
    ) |>
        dplyr::filter(!is.na(Gene)) |>
        dplyr::distinct()
}


.align_seqlevels <- function(subject_gr, query_gr) {
    subject_gr <- .set_seqlevels_style(subject_gr, query_gr)

    common <- intersect(
        GenomeInfoDb::seqlevels(query_gr),
        GenomeInfoDb::seqlevels(subject_gr)
    )
    GenomeInfoDb::keepSeqlevels(subject_gr, common, pruning.mode = "coarse")
}


.set_seqlevels_style <- function(subject_gr, query_gr) {
    withCallingHandlers(
        {
            target_style <- GenomeInfoDb::seqlevelsStyle(query_gr)
            GenomeInfoDb::seqlevelsStyle(subject_gr) <- target_style
        },
        warning = function(w) {
            if (grepl("cannot switch some .* seqlevels", conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        }
    )

    subject_gr
}


.subset_common_seqlevels <- function(query_gr, subject_gr) {
    common <- intersect(
        GenomeInfoDb::seqlevels(query_gr),
        GenomeInfoDb::seqlevels(subject_gr)
    )
    GenomeInfoDb::keepSeqlevels(query_gr, common, pruning.mode = "coarse")
}


.overlap_df <- function(query_gr,
                        subject_gr,
                        subject_gene,
                        src,
                        ignore_strand = FALSE) {
    query_gr2 <- .subset_common_seqlevels(query_gr, subject_gr)
    hits <- GenomicRanges::findOverlaps(
        query_gr2,
        subject_gr,
        ignore.strand = ignore_strand
    )

    if (length(hits) == 0) {
        return(.empty_link_df())
    }

    data.frame(
        Peak = query_gr2$id[S4Vectors::queryHits(hits)],
        Gene = subject_gene[S4Vectors::subjectHits(hits)],
        Src = src,
        stringsAsFactors = FALSE
    ) |>
        dplyr::distinct()
}


.empty_link_df <- function() {
    data.frame(
        Peak = character(0),
        Gene = character(0),
        Src = character(0),
        stringsAsFactors = FALSE
    )
}


.combine_link_tables <- function(df_enh, df_prom, df_clo) {
    dplyr::bind_rows(df_enh, df_prom, df_clo) |>
        dplyr::distinct()
}
