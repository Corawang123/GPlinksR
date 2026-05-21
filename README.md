# GPlinksR

**GPlinksR** builds **gene-peak networks** for multi-omics integration by combining enhancer-, promoter-, and proximity-based (closest-gene) mapping.  
It uses public enhancer-gene link data from **PANTHER / PEREGRINE (AnnoQ)**.

---

## Installation

Install **GPlinksR** with `BiocManager`:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("GPlinksR")
```

Until the package is available from Bioconductor, install the development
version from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("Corawang123/GPlinksR")
```
