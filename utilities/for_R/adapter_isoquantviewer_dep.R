library(reticulate)
library(rhdf5)

## allow character or NULL for optional paths
setClassUnion("characterORNULL", c("character", "NULL"))
setClassUnion("listORNULL",      c("list", "NULL"))

setClass(
  "isoquantViewer",
  slots = list(
    output_directory = "character",         # required
    prefix           = "character",         # required
    reference_gtf    = "characterORNULL",   # optional
    gene_db          = "characterORNULL",   # optional
    mudata_path      = "characterORNULL"    # optional
  ),
  prototype = list(
    output_directory = NULL,
    prefix       = NULL,
    reference_gtf   = NULL,
    gene_db          = NULL,
    mudata_path         = NULL
  )
)

Make_isoquantViewer <- function(output_directory,
                  prefix,
                  reference_gtf = NULL,
                  gene_db       = NULL,
                  mudata_path   = NULL) {

  new("isoquantViewer",
      output_directory = output_directory,
      prefix           = prefix,
      reference_gtf    = reference_gtf,
      gene_db          = gene_db,
      mudata_path      = mudata_path)
}


setMethod("show", "Make_isoquantViewer", function(object) {
  cat("An object of class 'isoquantViewer'\n")
  cat("  output_directory:", object@output_directory, "\n")
  cat("  prefix          :", object@prefix, "\n")
  cat("  reference_gtf   :", if (is.null(object@reference_gtf)) "<NULL>" else object@reference_gtf, "\n")
  cat("  gene_db         :", if (is.null(object@gene_db)) "<NULL>" else object@gene_db, "\n")
  cat("  mudata_path     :", if (is.null(object@mudata_path)) "<NULL>" else object@mudata_path, "\n")
 
})

