library(polars)
library(ggplot2)
library(dplyr)
library(scales)
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
    gene_db          = "characterORNULL",            # optional
    mudata_path      = "characterORNULL",
    python_viewer    = "ANY",
    reads_assignment = "ANY"

  ),
  prototype = list(
    output_directory = NULL,
    prefix       = NULL,
    reference_gtf   = NULL,
    gene_db = NULL,
    mudata_path = NULL,
    python_viewer = NULL
    )
)

isoquantViewer <- function(output_directory,
                  prefix,
                  reference_gtf = NULL,
                  gene_db = NULL,
                  mudata_path = NULL
                  ) {

  new("isoquantViewer",
      output_directory = output_directory,
      prefix           = prefix,
      reference_gtf    = reference_gtf,
      gene_db = gene_db,
      mudata_path = mudata_path
      )

}

setMethod("show", "isoquantViewer", function(object) {
  cat("An object of class 'isoquantViewer'\n")
  cat("  output_directory:", object@output_directory, "\n")
  cat("  prefix          :", object@prefix, "\n")
  cat("  reference_gtf   :", if (is.null(object@reference_gtf)) "<NULL>" else object@reference_gtf, "\n")
  cat("  gene_db   :", if (is.null(object@gene_db)) "<NULL>" else object@gene_db, "\n")
  cat("  mudata_path   :", if (is.null(object@gene_db)) "<NULL>" else object@gene_db, "\n")
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
        mudata_path      = .Object@mudata_path
      )
    } else {
      .Object@python_viewer <- NULL
    }

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


# setMethod(
#   "plot_pie_assignment",
#   signature(object = "isoquantViewer"),
#   function(object,
#           feature_to_plot="ANY",
#           figsize = c(5,4),
#           savefig = TRUE,
#           save_format = 'png', 
#           return_data = FALSE
#           ) {
#             p_plot_func$plot_pie_assignment(obj=object@python_viewer,
#             feature_to_plot = feature_to_plot,
#             figsize = figsize,
#             savefig = savefig,
#             save_format = save_format,
#             return_data = return_data)
#             }
#   )

setMethod(
  "plot_pie_assignment",
  signature(object = "isoquantViewer"),
  function(object,
          feature_to_plot="ANY",
          figsize = c(5,4),
          savefig = TRUE,
          save_format = 'png', 
          return_data = FALSE
          ) {
          
          read_assign_fname <- glue("{object@output_directory}{object@prefix}/{object@prefix}.read_assignments.tsv.gz")
          scan <-pl$scan_csv(
          path.expand(read_assign_fname),
          skip_rows=2,
          separator="\t",
          has_header=TRUE,
          infer_schema_length=0)

        if (feature_to_plot == 'assignment_type'){
            scan <- scan$with_columns(
            pl$col(feature_to_plot)$
            alias('feature'))
        }
        else {
          
          pattern <- sprintf("(?:^|;)%s=([^;]*)", feature_to_plot)
          scan <- scan$with_columns(
            pl$col("additional_info")$
            str$extract(pattern, 1L)$
            alias('feature'))
        }

          df<-scan$select(
            pl$col('feature')$
            value_counts())$
            unnest('feature')$
            collect(engine = "streaming")

          df<- as.data.frame(df) %>%
          mutate(
            frac = count / sum(count),
            label = percent(frac, accuracy = 0.1),
            feature_label = glue("{feature} (n={count})"),
            ymax = cumsum(frac),
            ymin = c(0, head(ymax, n = -1)),
            mid = (ymin + ymax) / 2
          )

          total = sum(df$count)

          p <- ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0,
                      fill = feature_label)) +
          geom_rect() +
          geom_text(
            aes(x = 0.5, y = mid, label = label),
            size = 4
          ) +
          coord_polar(theta = "y") +
          theme_void() +
          ggtitle(paste0("Assignment\n total no. reads: ", total)) +
          theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "right",
            legend.title = element_blank()
          )
        print (p)

        if (savefig) {
        plot_output <- "./plot_output"
        if (!dir.exists(plot_output)) dir.create(plot_output, recursive = TRUE)
        
        ofname <- file.path(plot_output,
                            paste0("Pie_", feature_to_plot, ".", save_format))
        ggsave(ofname, p, width = figsize[1], height = figsize[2])
      }
      
      invisible(p)
    }

            
  )


# setMethod(
#   "plotfunc",
#   signature(object = "isoquantViewer"),
#   function(
#           object, 
#           plot_type,
#           use_transcript_model,
#           Ensembl_ID,
#           gene_names,
#           figsize,
#           savefig=,

#           ) {

#     if (plot_type == 'pie_assignment'){

#     }

#     if (plot_type == 'count_bar'){

#     }

#      if (plot_type == 'transcript_map'){

#       p_plot_func$plot_transcript_map(  obj=,
#             use_transcript_model=,
#             Ensembl_ID=,
#             gene_names=,
#             figsize=,
#             savefig=,
#             save_format=)

#     }

#     if (plot_type == 'genomic_region'){

#     }

   
#   }
# )