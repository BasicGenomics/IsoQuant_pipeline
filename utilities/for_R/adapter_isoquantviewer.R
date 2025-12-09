library(logger)
library(glue)
## Import python class
p_isoquantViewer <- import_from_path(
  "isoquantViewer",
  path = "~/Documents/git_hub/isoquant/IsoQuant_pipeline/utilities/"
)
p_plot_func <- import_from_path(
  "plot_func",
  path = "~/Documents/git_hub/isoquant/IsoQuant_pipeline/utilities/"
)


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
    gene_db_model    = 'characterORNULL',   # optional
    mudata_path      = "characterORNULL",   # optional
    plot_dir         = "characterORNULL",   # optional
    python_viewer    = "ANY"
  ),
  prototype = list(
    output_directory = NULL,
    prefix       = NULL,
    reference_gtf   = NULL,
    gene_db = NULL,
    gene_db_model = NULL,
    mudata_path = NULL,
    plot_dir = NULL,
    python_viewer = NULL
    )
)

isoquantViewer <- function(output_directory,
                  prefix,
                  reference_gtf = NULL,
                  gene_db = NULL,
                  gene_db_model = NULL,
                  mudata_path = NULL,
                  plot_dir = NULL
                  ) {

  new("isoquantViewer",
      output_directory = output_directory,
      prefix           = prefix,
      reference_gtf    = reference_gtf,
      gene_db = gene_db,
      gene_db_model = gene_db_model,
      mudata_path = mudata_path,
      plot_dir = plot_dir
      )
}

setMethod("show", "isoquantViewer", function(object) {
  cat("An object of class 'isoquantViewer'\n")
  cat("  output_directory:", object@output_directory, "\n")
  cat("  prefix          :", object@prefix, "\n")
  cat("  reference_gtf   :", if (is.null(object@reference_gtf)) "<NULL>" else object@reference_gtf, "\n")
  cat("  gene_db   :", if (is.null(object@gene_db)) "<NULL>" else object@gene_db, "\n")
  cat("  gene_db_model   :", if (is.null(object@gene_db_model)) "<NULL>" else object@gene_db_model, "\n")
  cat("  mudata_path   :", if (is.null(object@mudata_path)) "<NULL>" else object@mudata_path, "\n")
  cat("  plot_dir   :", if (is.null(object@plot_dir)) "<NULL>" else object@plot_dir, "\n")
}
)

setMethod(
  "initialize",
  "isoquantViewer",
  function(.Object, ...) {

    ## Fill the S4 slots from new(...)
    .Object <- callNextMethod()

    ## Only create Python viewer if we have the required info
    if (.Object@output_directory != "" && .Object@prefix != "") {

      .Object@python_viewer <- p_isoquantViewer$isoquantViewer(
        output_directory = .Object@output_directory,
        prefix           = .Object@prefix,
        # NOTE: Python arg name looks like 'referene_gtf' in your code
        referene_gtf     = .Object@reference_gtf,
        gene_db          = .Object@gene_db,
        gene_db_model    = .Object@gene_db_model,
        mudata_path      = .Object@mudata_path,
        plot_dir         = .Object@plot_dir
      )
    } else {
      .Object@python_viewer <- NULL
    }

    logger::log_info(glue("Plots are saved in {.Object@python_viewer$plot_dir}"))

    .Object
  }
)



### function that work on class isoquantViewer


setGeneric(
  "plot_transcript_map",
  function(object,...) {
    standardGeneric("plot_transcript_map")
  }
)


setMethod(
  "plot_transcript_map",
  signature(object = "isoquantViewer"),
  function(object,
          use_transcript_model=FALSE,
          Ensembl_ID=NULL,
          gene_names=NULL,
          figsize=c(8,3.5),
          savefig=TRUE,
          save_format='png') {
            p_plot_func$plot_transcript_map(obj=object@python_viewer,
            use_transcript_model=use_transcript_model,
            Ensembl_ID=Ensembl_ID,
            gene_names=gene_names,
            figsize=figsize,
            savefig=savefig,
            save_format=save_format)
  }
  )



setGeneric(
  "plot_genomic_region",
  function(object,...) {
    standardGeneric("plot_genomic_region")
  }
)


setMethod(
  "plot_genomic_region",
  signature(object = "isoquantViewer"),
  function(object,
          region = NULL,
          Ensembl_ID = NULL,
          gene_name = NULL,
          plot_tracks = c('ref','model','reads'),
          tracks_params = c(),
          tracks_file = NULL,
          padding = 1000,
          figsize=c(10,18),
          label_fontsize = 5,
          savefig=TRUE,
          save_format='png'
          ) {
            p_plot_func$plot_genomic_region(obj=object@python_viewer,
            region = region,
            Ensembl_ID = Ensembl_ID,
            gene_name = gene_name,
            plot_tracks = plot_tracks,
            tracks_params = tracks_params,
            tracks_file = tracks_file,
            padding = padding,
            figsize = figsize,
            label_fontsize = label_fontsize,
            savefig = savefig,
            save_format = save_format)
            }
  )



setGeneric(
  "plot_pie_assignment",
  function(object,...) {
    standardGeneric("plot_pie_assignment")
  }
)


setMethod(
  "plot_pie_assignment",
  signature(object = "isoquantViewer"),
  function(object,
          feature_to_plot="ANY",
          figsize = c(5,4),
          savefig = TRUE,
          save_format = 'png') {
            p_plot_func$plot_pie_assignment(obj=object@python_viewer,
            feature_to_plot = feature_to_plot,
            figsize = figsize,
            savefig = savefig,
            save_format = save_format)
            }
  )


setGeneric(
  "plot_count_bar",
  function(object,...) {
    standardGeneric("plot_count_bar")
  }
)

setMethod(
  "plot_count_bar",
  signature(object = "isoquantViewer"),
  function(object,
          layer="count",
          sample_id = NULL,
          gene_list = NULL,
          isoform_list = NULL,
          use_transcript_model = FALSE, 
          savefig = TRUE,
          save_format= 'png') {
          p_plot_func$plot_count_bar(obj=object@python_viewer,
          layer = layer,
          sample_id = sample_id,
          gene_list = gene_list,
          isoform_list = isoform_list,
          use_transcript_model = use_transcript_model,
          savefig = savefig,
          save_format = save_format)
            }
  )