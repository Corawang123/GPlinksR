#' Build Gene–Peak Network (Enhancer, Promoter, Closest)
#'
#' @param pk Either:
#'   (1) A character vector of genomic coordinates ("chr1:1000-1200"), or
#'   (2) A data.frame containing peak coordinates with columns "chr","start","end"
#'       or a "PeakRegion" column ("chr1:1000-1200" style).
#' @param gn Character vector of gene symbols (HGNC official symbols, e.g. "TP53", "FMNL2").
#' @param ver Integer (17, 18, or 19) — PANTHER enhancer–gene link version.
#' @param enh_file Optional local path to enhancer–gene link file (.tsv).
#'
#' @return A data.frame with columns: Peak, Gene, Source.
#' @export
#'
#' @examples
#' gp <- build_gp_links(c("chr1:819770-822338", "chr1:983871-984475"), c("TTLL10", "PERM1"))
#' head(gp)
build_gp_links <- function(pk, gn, ver = 19, enh_file = NULL) {

  suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(EnsDb.Hsapiens.v86)
    library(ensembldb)
    library(biomaRt)
    library(dplyr)
  })


  # 0. Input validation

  message("Checking inputs...")

  # --- Check gene input ---
  if (missing(gn) || is.null(gn) || length(gn) == 0)
    stop("'gn' must be provided — a non-empty character vector of gene symbols.")

  if (!is.character(gn))
    stop("'gn' must be a character vector (e.g., c('TP53', 'FMNL2')).")

  if (any(is.na(gn) | gn == ""))
    stop("'gn' contains missing or empty entries. Please remove or replace them.")

  # --- Check peak input ---
  if (missing(pk) || is.null(pk) || length(pk) == 0)
    stop("'pk' must be provided — either a character vector or a data.frame.")

  if (is.character(pk)) {
    if (!all(grepl("^chr[0-9XYMT]+:[0-9]+-[0-9]+$", pk))) {
      bad <- pk[!grepl("^chr[0-9XYMT]+:[0-9]+-[0-9]+$", pk)][1]
      stop(paste0("Invalid peak format: '", bad, "'. Expected 'chr#:start-end' (e.g., 'chr1:1000-1200')."))
    }
  } else if (is.data.frame(pk)) {
    cols <- colnames(pk)
    if (!("PeakRegion" %in% cols || all(c("chr","start","end") %in% cols))) {
      stop("Peak data.frame must contain either 'PeakRegion' or columns 'chr','start','end'.")
    }
  } else {
    stop("'pk' must be a character vector or a data.frame, not ", class(pk))
  }

  # --- Check enhancer version ---
  if (!ver %in% c(17, 18, 19))
    stop("'ver' must be one of 17, 18, or 19 (PEREGRINE version numbers).")

  # --- Check enhancer file existence ---
  if (!is.null(enh_file) && !file.exists(enh_file))
    stop("The provided 'enh_file' path does not exist: ", enh_file)

  message(paste0(" Input validation passed: ", length(pk), " peaks and ", length(gn), " genes provided."))


  # 1. Prepare Peak GRanges

  message("Preparing peaks...")

  if (is.data.frame(pk)) {
    if ("PeakRegion" %in% colnames(pk)) {
      pk_split <- tstrsplit(pk$PeakRegion, "[:-]")
      pk_dt <- data.table(chr = pk_split[[1]], start = pk_split[[2]], end = pk_split[[3]])
    } else if (all(c("chr", "start", "end") %in% colnames(pk))) {
      pk_dt <- data.table(chr = pk$chr, start = pk$start, end = pk$end)
    } else {
      stop("Peak input must contain either 'PeakRegion' or columns 'chr','start','end'.")
    }
  } else if (is.character(pk)) {
    pk_split <- tstrsplit(pk, "[:-]")
    pk_dt <- data.table(chr = pk_split[[1]], start = pk_split[[2]], end = pk_split[[3]])
  } else {
    stop("Unsupported peak input type. Provide a character vector or data.frame.")
  }

  pk_dt$start <- as.integer(pk_dt$start)
  pk_dt$end   <- as.integer(pk_dt$end)
  pk_dt$id    <- paste0(pk_dt$chr, ":", pk_dt$start, "-", pk_dt$end)
  pk_dt$chr   <- ifelse(!grepl("^chr", pk_dt$chr), paste0("chr", pk_dt$chr), pk_dt$chr)

  valid_chr <- paste0("chr", c(1:22, "X", "Y", "M", "MT"))
  pk_dt <- pk_dt[pk_dt$chr %in% valid_chr, ]

  pk_gr <- makeGRangesFromDataFrame(pk_dt,
                                    seqnames.field = "chr",
                                    start.field = "start",
                                    end.field = "end",
                                    keep.extra.columns = TRUE)


  # 2. Enhancer-based links (PEREGRINE)

  message("Building enhancer-based links...")
  if (is.null(enh_file)) enh_file <- get_peregrine_file(ver)

  enh_raw <- fread(enh_file, sep = "\t", header = TRUE, fill = TRUE, quote = "")
  enh_raw$HGNC <- sub(".*HGNC=([0-9]+).*", "\\1", enh_raw$gene)
  enh_raw$enhancer <- as.character(enh_raw$enhancer)

  # Map user genes to HGNC IDs
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gn_map <- getBM(attributes = c("hgnc_symbol", "hgnc_id"),
                  filters = "hgnc_symbol", values = gn, mart = mart)
  if (nrow(gn_map) == 0)
    stop("No HGNC IDs found for the provided gene symbols.")
  gn_ids <- sub("HGNC:", "", gn_map$hgnc_id)

  enh_sub <- enh_raw[enh_raw$HGNC %in% gn_ids, ]

  message("Loading enhancer coordinates...")
  enh_url <- "https://data.pantherdb.org/ftp/peregrine_data/PEREGRINEenhancershg38"
  enh_reg <- fread(enh_url, sep = "\t", header = FALSE,
                   col.names = c("chr", "start", "end", "enhancer"))
  enh_full <- merge(enh_sub, enh_reg, by = "enhancer", all.x = TRUE)
  enh_full$sym <- gn_map$hgnc_symbol[match(enh_full$HGNC, gn_ids)]

  enh_gr <- makeGRangesFromDataFrame(enh_full,
                                     seqnames.field = "chr",
                                     start.field = "start",
                                     end.field = "end",
                                     keep.extra.columns = TRUE)

  suppressWarnings(seqlevelsStyle(enh_gr) <- seqlevelsStyle(pk_gr))
  common <- intersect(seqlevels(pk_gr), seqlevels(enh_gr))
  pk_gr2 <- keepSeqlevels(pk_gr, common, pruning.mode = "coarse")
  enh_gr <- keepSeqlevels(enh_gr, common, pruning.mode = "coarse")

  hits_enh <- findOverlaps(pk_gr2, enh_gr)
  df_enh <- data.frame(
    Peak = pk_gr2$id[queryHits(hits_enh)],
    Gene = mcols(enh_gr)$sym[subjectHits(hits_enh)],
    Src  = "enh",
    stringsAsFactors = FALSE
  ) %>% distinct()

  message("Enhancer–gene links found: ", nrow(df_enh))


  # 3. Promoter-based links (EnsDb)

  message("Building promoter-based links...")

  edb <- EnsDb.Hsapiens.v86
  tx_gr <- transcripts(edb, columns = c("tx_id","gene_id","gene_name","seq_name","strand"))
  prom_gr <- promoters(tx_gr, upstream = 2000, downstream = 200)
  prom_gr$symbol <- tx_gr$gene_name
  prom_gr <- prom_gr[prom_gr$symbol %in% gn]
  prom_gr$PromoterRegion <- paste0(seqnames(prom_gr), ":", start(prom_gr), "-", end(prom_gr))

  suppressWarnings(seqlevelsStyle(prom_gr) <- seqlevelsStyle(pk_gr))
  common <- intersect(seqlevels(pk_gr), seqlevels(prom_gr))
  pk_gr3 <- keepSeqlevels(pk_gr, common, pruning.mode = "coarse")
  prom_gr <- keepSeqlevels(prom_gr, common, pruning.mode = "coarse")

  hits_prom <- findOverlaps(pk_gr3, prom_gr, ignore.strand = TRUE)
  df_prom <- data.frame(
    Peak = pk_gr3$id[queryHits(hits_prom)],
    Gene = prom_gr$symbol[subjectHits(hits_prom)],
    Src  = "prom",
    stringsAsFactors = FALSE
  ) %>% distinct()

  message("Promoter–gene links found: ", nrow(df_prom))


  # 4. Closest-gene links (EnsDb)

  message("Building closest-gene links...")

  tx_gene <- genes(edb, columns = c("gene_id","gene_name","seq_name"))
  tx_gene$symbol <- tx_gene$gene_name
  tx_gene <- tx_gene[tx_gene$symbol %in% gn]

  suppressWarnings(seqlevelsStyle(tx_gene) <- seqlevelsStyle(pk_gr))
  common <- intersect(seqlevels(pk_gr), seqlevels(tx_gene))
  pk_gr4 <- keepSeqlevels(pk_gr, common, pruning.mode = "coarse")
  tx_gene <- keepSeqlevels(tx_gene, common, pruning.mode = "coarse")

  # trim out-of-bound GRanges
  pk_gr4 <- trim(pk_gr4)
  tx_gene <- trim(tx_gene)

  nearest_hits <- nearest(pk_gr4, tx_gene, ignore.strand = TRUE)

  df_clo <- data.frame(
    Peak = pk_gr4$id,
    Gene = tx_gene$symbol[nearest_hits],
    Src  = "clo",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(Gene)) %>%
    dplyr::distinct()

  message("Closest–gene links found: ", nrow(df_clo))


  # 5. Combine and return

  df_all <- bind_rows(df_enh, df_prom, df_clo) %>% distinct()
  message("Total combined links: ", nrow(df_all))

  return(df_all)
}

