#' Example Link Table for the Packaged GPlinksR Demo Subset
#'
#' A precomputed example gene-peak link table corresponding to the curated
#' leading entries in `gp_example_inputs`. The table includes enhancer-,
#' promoter-, and closest-based links and is used in the vignette to display
#' representative output without rerunning external downloads.
#'
#' @format A data.frame with 7 rows and 3 columns:
#' \describe{
#'   \item{Peak}{Peak coordinate in `"chr:start-end"` format.}
#'   \item{Gene}{Gene symbol linked to the peak.}
#'   \item{Src}{Link source, one of `"enh"`, `"prom"`, or `"clo"`.}
#' }
"gp_example_links"
