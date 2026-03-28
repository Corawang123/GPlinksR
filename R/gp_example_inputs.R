#' Example Gene and Peak Inputs for GPlinksR
#'
#' A small example dataset for demonstrating `GPlinksR` workflows. The object
#' contains 300 peak coordinates and 100 gene symbols suitable for examples and
#' vignettes. The first few entries are curated so that the vignette examples
#' produce enhancer-, promoter-, and closest-based links from the packaged data.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{pk}{A character vector of 300 peak coordinates in
#'   `"chr:start-end"` format.}
#'   \item{gn}{A character vector of 100 gene symbols.}
#' }
"gp_example_inputs"
