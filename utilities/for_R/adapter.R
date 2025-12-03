
library(reticulate)
library(glue)
library(rhdf5)
library(Rgff)


read_mudata_layer <- function(
    isoquantViewer,
    use_transcript_model,
    layer
    ) {

    if(use_transcript_model) {
        mod <-'isoform'
    } else {  
    mod <-'reference_isoform'
    }

    toread <-glue('mod/reference_isoform/{mod}/{layer}')

    read_layer <- h5read(isoquantViewer@mudata_path,toread)

    }


read_mudata_var <- function(
    isoquantViewer,
    use_transcript_model
    ) {

    if(use_transcript_model) {
        mod <-'isoform'
    } else {  
    mod <-'reference_isoform'
    }

    toread <-glue('mod/reference_isoform/{mod}/var/')
    read_layer <- h5read(isoquantViewer@mudata_path,toread)

    }
    

# plot_pie_assignment <-function
# obj,
#     feature_to_plot: Optional[
#         Literal["assignment_type", "gene_assignment", "Classification"]
#     ],
#     figsize:tuple[float,float]=(5,4),
#     savefig:bool=True

plot_transcript_map <- function(use_transcript_model = FALSE,
                   Ensembl_ID = NULL,
                   gene_names = NULL,
                   figsize = c(8, 3.5),
                   savefig = TRUE) {
        
    plot_func <- reticulate::import("../plot_func")
    plot_func$plot_transcript_map(
    use_transcript_model = use_transcript_model,
    Ensembl_ID =Ensembl_ID,
    gene_names = gene_names,
    figsize = figsize,
    savefig = savefig)
                   }


plot_genomic_region <- function(
    mudata_path = "",
    region       = NULL,
    Ensembl_ID   = NULL,
    gene_name    = NULL,
    plot_tracks  = c("ref", "model", "reads"),
    tracks_params = NULL,
    tracks_file   = list(),
    padding       = 1000L,
    figsize       = c(10, 18),
    label_fontsize = 5,
    savefig        = TRUE,
    save_format    = "png") {

        plot_func <- reticulate::import("../plot_func")
        plot_func$plot_genomic_region(
            obj = mudata_path,
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
