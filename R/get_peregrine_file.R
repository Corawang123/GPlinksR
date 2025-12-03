#' Download Enhancer–Gene Link File from PANTHER Peregrine Database
#'
#' @param version Integer (17, 18, or 19). Which PANTHER enhancer–gene link version to download.
#' @param destdir Directory to save the file. Defaults to a temporary directory.
#' @return The full path to the downloaded .tsv file.
#' @export
#'
#' @examples
#' file_path <- get_peregrine_file(19)
#' head(readLines(file_path), 3)
get_peregrine_file <- function(version = 19, destdir = tempdir()) {
  if (!version %in% c(17, 18, 19))
    stop("Version must be one of 17, 18, or 19.", call. = FALSE)

  base_url <- "https://data.pantherdb.org/ftp/peregrine_data/"
  fname <- paste0("enhancer_gene_link_", version, ".tsv")
  file_url <- paste0(base_url, fname)
  dest_file <- file.path(destdir, fname)

  if (!file.exists(dest_file)) {
    message("Downloading ", fname, " from PANTHER FTP...")

    # extend timeout temporarily
    old_timeout <- getOption("timeout")
    options(timeout = max(600, old_timeout))  # extend to 10 minutes

    tryCatch({
      utils::download.file(file_url, destfile = dest_file, mode = "wb", quiet = FALSE)
      message("Download complete: ", dest_file)
    }, error = function(e) {
      stop("Failed to download enhancer file from ", file_url, call. = FALSE)
    }, finally = {
      options(timeout = old_timeout)
    })
  } else {
    message("Using cached file: ", dest_file)
  }

  return(dest_file)
}


