#' Build Gene–Peak Network (Enhancer, Promoter, Closest)
#'
#' @param pk Either:
#'   (1) A character vector of genomic coordinates ("chr1:1000-1200"), or
#'   (2) A data.frame containing columns "chr","start","end" or "PeakRegion".
#' @param gn Character vector of HGNC gene symbols (e.g., "TP53", "FMNL2").
#' @param ver Integer (17, 18, or 19) — PANTHER enhancer–gene link version.
#' @param enh_file Optional path to enhancer–gene link file (.tsv).
#'
#' @return Data frame with columns: Peak, Gene, Src (enh/prom/clo)
#' @export
build_gp_links <- function(pk, gn, ver = 19, enh_file = NULL) {
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(EnsDb.Hsapiens.v86)
    library(ensembldb)
    library(biomaRt)
    library(dplyr)
    library(data.table)
  })

  # ensure clean thread safety
  data.table::setDTthreads(1)

  # ---- 0. Enhancer–gene link file ----
  if (is.null(enh_file)) enh_file <- get_peregrine_file(ver)

  # ---- 1. Prepare peaks ----
  message("Preparing peaks...")
  if (is.data.frame(pk)) {
    if ("PeakRegion" %in% names(pk)) {
      pk_dt <- data.table::as.data.table(data.table::tstrsplit(pk$PeakRegion, "[:-]"))
    } else if (all(c("chr","start","end") %in% names(pk))) {
      pk_dt <- data.table::as.data.table(pk[, c("chr","start","end")])
    } else stop("Peak input must have PeakRegion or chr/start/end.")
  } else if (is.character(pk)) {
    pk_dt <- data.table::as.data.table(data.table::tstrsplit(pk, "[:-]"))
  } else stop("Unsupported peak input type — provide vector or data.frame.")

  data.table::setnames(pk_dt, c("chr","start","end"))
  pk_dt[, start := as.integer(start)]
  pk_dt[, end := as.integer(end)]
  pk_dt[, id := paste0(chr, ":", start, "-", end)]

  # normalize chr names and remove non-standard contigs
  pk_dt[, chr := ifelse(!grepl("^chr", chr), paste0("chr", chr), chr)]
  valid_chr <- paste0("chr", c(1:22, "X", "Y", "M", "MT"))
  pk_dt <- pk_dt[chr %in% valid_chr]

  pk_gr <- GenomicRanges::makeGRangesFromDataFrame(
    pk_dt, seqnames.field = "chr", start.field = "start",
    end.field = "end", keep.extra.columns = TRUE
  )

  # ---- 2. Enhancer links ----
  message("Loading enhancer links...")
  enh_raw <- data.table::fread(enh_file, sep = "\t", header = TRUE, fill = TRUE)
  enh_raw[, HGNC := sub(".*HGNC=([0-9]+).*", "\\1", gene)]
  enh_raw[, enhancer := as.character(enhancer)]

  # Map gene symbols to HGNC IDs
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gn_map <- biomaRt::getBM(
    attributes = c("hgnc_symbol", "hgnc_id"),
    filters = "hgnc_symbol",
    values = gn,
    mart = mart
  )
  gn_ids <- sub("HGNC:", "", gn_map$hgnc_id)
  enh_sub <- enh_raw[HGNC %in% gn_ids]

  # Load enhancer coordinates
  message("Loading enhancer coordinates...")
  enh_reg <- data.table::fread(
    "https://data.pantherdb.org/ftp/peregrine_data/PEREGRINEenhancershg38",
    sep = "\t", header = FALSE,
    col.names = c("chr", "start", "end", "enhancer")
  )
  enh_full <- merge(enh_sub, enh_reg, by = "enhancer", all.x = TRUE)
  enh_full[, sym := gn_map$hgnc_symbol[match(HGNC, gn_ids)]]

  enh_gr <- GenomicRanges::makeGRangesFromDataFrame(
    enh_full,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )

  # Overlaps
  common <- intersect(GenomicRanges::seqlevels(pk_gr), GenomicRanges::seqlevels(enh_gr))
  pk_gr  <- GenomicRanges::keepSeqlevels(pk_gr, common, pruning.mode = "coarse")
  enh_gr <- GenomicRanges::keepSeqlevels(enh_gr, common, pruning.mode = "coarse")

  hit_enh <- GenomicRanges::findOverlaps(pk_gr, enh_gr)
  df_enh <- data.table::data.table(
    Peak = pk_gr$id[queryHits(hit_enh)],
    Gene = mcols(enh_gr)$sym[subjectHits(hit_enh)],
    Src  = "enh"
  ) %>% dplyr::distinct()

  # ---- 3. Promoters ----
  message("Mapping promoters...")
  edb <- EnsDb.Hsapiens.v86
  tx <- ensembldb::transcripts(edb, columns = c("tx_id", "gene_name", "seq_name", "strand"))
  prom <- GenomicRanges::promoters(tx, upstream = 2000, downstream = 200)
  prom$gn <- tx$gene_name
  prom <- prom[prom$gn %in% gn]
  suppressWarnings(GenomicRanges::seqlevelsStyle(prom) <- GenomicRanges::seqlevelsStyle(pk_gr))
  hit_prom <- GenomicRanges::findOverlaps(pk_gr, prom, ignore.strand = TRUE)
  df_prom <- data.table::data.table(
    Peak = pk_gr$id[queryHits(hit_prom)],
    Gene = prom$gn[subjectHits(hit_prom)],
    Src  = "prom"
  ) %>% dplyr::distinct()

  # ---- 4. Closest genes ----
  message("Finding closest genes...")
  gns <- ensembldb::genes(edb, columns = c("gene_id", "gene_name", "seq_name"))
  gns <- gns[gns$gene_name %in% gn]
  suppressWarnings(GenomicRanges::seqlevelsStyle(gns) <- GenomicRanges::seqlevelsStyle(pk_gr))
  nearest_hit <- GenomicRanges::nearest(pk_gr, gns, ignore.strand = TRUE)
  df_clo <- data.table::data.table(
    Peak = pk_gr$id,
    Gene = gns$gene_name[nearest_hit],
    Src  = "clo"
  ) %>%
    dplyr::filter(!is.na(Gene)) %>%
    dplyr::distinct()

  # ---- 5. Merge all results ----
  df_all <- dplyr::bind_rows(df_enh, df_prom, df_clo) %>%
    dplyr::distinct(Peak, Gene, Src)
  df_all$Src <- factor(df_all$Src, levels = c("enh", "prom", "clo"))

  message("Final combined links: ", nrow(df_all))
  return(as.data.frame(df_all))
}

