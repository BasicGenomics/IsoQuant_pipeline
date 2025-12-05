import matplotlib.pyplot as plt
from typing import List,Optional,Dict,Literal
from collections import Counter
import numpy as np
import matplotlib.ticker as ticker
import pandas as pd
from pygenometracks.tracks.BigWigTrack import BigWigTrack
from pygenometracks.tracks.GtfTrack import GtfTrack
from pygenometracks.tracks.ScaleBarTrack import ScaleBarTrack
from pygenometracks.tracks.BedTrack import BedTrack
import mudata

from pybedtools import BedTool
import math
import os
import logging
from default_tracks import bed_props, gtf_props, ref_props, bam_props
from pathlib import Path
import gffutils
from gffutils.exceptions import FeatureNotFoundError
from isoquantViewer import isoquantViewer,create_db_genedict_from_geneid

def plot_pie_assignment(
    obj,
    feature_to_plot: Optional[
        Literal["assignment_type", "gene_assignment", "Classification"]
    ],
    figsize:tuple[float,float]=(5,4),
    savefig:bool=True,
    save_format:str='png',
    return_data:bool=False
):
    

    if hasattr(obj, "reads_assignment"):
        df = obj.reads_assignment
        count = df[feature_to_plot].value_counts()
        

    else:
        import polars as pl
        import gzip

        read_assign_fname =  f'{obj.output_directory}/{obj.prefix}/{obj.prefix}.read_assignments.tsv.gz'
        with gzip.open(read_assign_fname, "rt") as f:
            comment_lines = []
            for line in f:
                if line.startswith("#"):
                    comment_lines.append(line)
                else:
                    break
    
        header = comment_lines[-1].lstrip("#").strip().split()
        n_comment_lines = len(comment_lines)
        scan = pl.scan_csv(
        read_assign_fname,
        separator="\t",
        has_header=False,
        new_columns=header,
        comment_prefix="#",      
        infer_schema_length=0,
    )
        
        if feature_to_plot != 'assignment_type':    
            scan = scan.with_columns(
                pl.col("additional_info").str.extract(rf"(?:^|;){feature_to_plot}=([^;]*)", 1).alias(f"{feature_to_plot}")
            )

        count = scan.select(pl.col(f"{feature_to_plot}").value_counts()).unnest(f"{feature_to_plot}").collect(engine = "streaming")
    
    labels = count[feature_to_plot].to_list()
    sizes = count['count'].to_list()
   

    total = sum(sizes)
    percentages = [(s / total) * 100 for s in sizes]

    # Legend format: Label – 54.3% (12345)
    legend_labels = [
        f"{l} ({c})" 
        for l, c in zip(labels, sizes)
    ]

    if return_data:
        data = pd.DataFrame(data={'labels':labels,
                                  'sizes':sizes})
        return data

    fig, ax = plt.subplots(figsize=figsize)

    # Pie: only percent on wedges, no label text
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=None,
        autopct="%1.1f%%",
        pctdistance=0.8,
        textprops={'fontsize': 8}
    )

    # Add legend with full info
    ax.legend(
        wedges,
        legend_labels,
        loc="center left",
        bbox_to_anchor=(1, 0.5),frameon=False,
        fontsize=8
    )

    ax.set_title(f"Assignment \n total no. reads:{total}",fontsize=10)

    if savefig:
        plot_output = './plot_output'
        if os.path.exists(plot_output) == False:
            os.makedirs(plot_output,exist_ok=True)
        ofname = os.path.join(plot_output,f'Pie_{feature_to_plot}.{save_format}')
        plt.savefig(ofname, format=save_format, dpi=144, bbox_inches='tight')



def plot_count_bar(obj,
            layer:str='count',
            sample_id:List[str]=None,
            gene_list:List[str]=None,
            isoform_list:List[str]=None,
            use_transcript_model:bool=False,
            savefig:bool=True,
            save_format:str='png',
            return_data=False
):

    """
    Plot bar chart of counts for either specified genes or isoforms.
    If providing gene_list, plots stacked bar of isoforms per gene.
    If providing isoform_list, plots bar chart of specified isoforms.
    
    Parameters:
    ----------
    obj: IsoQuant object containing mdata and gene dictionaries
    layer: data layer to use for counts (e.g., 'count', 'tpm').
    sample_id:  list of sample IDs to include in the plot,  if None, use all samples
    gene_list: list of gene names to plot
    isoform_list: list of isoform IDs to plot
    use_transcript_model: bool, whether to use transcript model from IsoQuant
    """

# --- loading gene dict
    print(gene_list,sample_id)
    if use_transcript_model:
        # if not hasattr(obj, "gene_dict_model"):
        #     obj.parse_input_gtf(use_ref=False)
        # gene_dict = obj.gene_dict_model

        gene_var = obj.mdata['gene'].var
    else:  
        # if not hasattr(obj, "gene_dict_ref"):
        #     obj.parse_input_gtf(use_ref=True)
        # gene_dict = obj.gene_dict_ref

        gene_var = obj.mdata['reference_gene'].var
        
    def is_nonempty(x):
        return x is not None and len(x) > 0

    # Validate inputs: require at least one non-empty
    if not (is_nonempty(gene_list) or is_nonempty(isoform_list)):
        raise ValueError("Provide non-empty 'gene_list' or 'isoform_list' for plotting.")

    # Only unique() when non-empty; otherwise keep as None
    gene_list = np.unique(gene_list) if is_nonempty(gene_list) else None
    isoform_list = np.unique(isoform_list) if is_nonempty(isoform_list) else None

    if sample_id is None:
        sample_id = list(obj.mdata.obs_names)

    if is_nonempty(gene_list):

        for g in gene_list:

            # isoforms = np.array(list(gene_dict[g]['transcripts'].keys()))

            isoforms = np.array(gene_var.var[g,:]['transcripts'].split(','))

            if use_transcript_model:
                X = obj.mdata['isoform'][sample_id,isoforms].to_df(layer)
            else:
                X = obj.mdata['reference_isoform'][sample_id,isoforms].to_df(layer)
                
            ylabel = 'Transcripts per million' if layer=='tpm' else 'Counts'

            X_sorted = X.reindex(X.sum().sort_values(ascending=False).index, axis=1)

            if return_data:
                return X_sorted
            
            xticks = X_sorted.index

            #Adjusting the figure width based on numbers of samples
            n_vars = len(xticks)
            base_width = 0.6 
            width = max(4, n_vars * base_width) 
            height = width*2/3
            
            fig, ax = plt.subplots(figsize=(width, height))

            X_sorted.plot.bar(ax=ax,stacked=True)
            ax.legend(loc='center left', 
            bbox_to_anchor=(1, 0.5),
            title=f'Transripts of {g}',
            frameon=False)

            ax.set_ylabel(ylabel)    
            ax.set_xticklabels(xticks, rotation=40, ha='right')
            ax.set_xlabel('Sample')

            if savefig:
                plot_output = './plot_output'
                if os.path.exists(plot_output) == False:
                    os.makedirs(plot_output,exist_ok=True)
                ofname = os.path.join(plot_output,f'Bar_{g}_TranscriptModel_{use_transcript_model}.{save_format}')
                plt.savefig(ofname, format=save_format, dpi=144, bbox_inches='tight')
                plt.show()
                plt.close()
            
    if is_nonempty(isoform_list):

        if use_transcript_model:
            X = obj.mdata['isoform'][sample_id,isoform_list].to_df(layer)
        else: 
            X = obj.mdata['reference_isoform'][sample_id,isoform_list].to_df(layer)
        
        X_sorted = X.T.reindex(X.T.sum().sort_values(ascending=False).index, axis=1)

        xticks = X_sorted.index
        ylabel = 'Counts' if layer=='count' else 'Transcripts per million'
            
        #Adjusting the figure width based on numbers of isoforms
        n_vars = len(xticks)
        base_width = 0.6 
        width = max(4, n_vars * base_width) 
        height = width*2/3

        fig, ax = plt.subplots(figsize=(width, height))

        X_sorted.plot.bar(ax=ax)
        ax.legend(loc='center left', 
                bbox_to_anchor=(1, 0.5),
                frameon=False,
                title='Sample')
       
        ax.set_ylabel(ylabel)    
        ax.set_xticklabels(xticks, rotation=40, ha='right')
        ax.set_xlabel('Sample')

        if savefig:
            plot_output = './plot_output'
            if os.path.exists(plot_output) == False:
                os.makedirs(plot_output,exist_ok=True)
            isoform_list_str = '_'.join(isoform_list)
            ofname = os.path.join(plot_output,f'Bar_isoforms_{isoform_list_str}_TranscriptModel_{use_transcript_model}.{save_format}')
            plt.show()
            plt.savefig(ofname, format=save_format, dpi=144, bbox_inches='tight')
            plt.close()


# --- helper function for plot_transcript_map, plot one gene per ax ---
def _draw_gene_on_ax(ax,
    gene_id:str,
    gene_data: Dict,
    show_xlabel: bool = False,
    strand_markers: bool = True,
):

    transcripts = gene_data["transcripts"]
   
    num_transcripts = len(transcripts)
    ax.set_title(
        f"Transcripts of Gene: {gene_data['name']} ({gene_id}) on Chromosome {gene_data['chromosome']}"
    )

    ax.set_ylabel("Transcripts")
    ax.set_yticks(range(num_transcripts))
    ax.set_yticklabels(
        [
            f"{transcript_id}"
            for transcript_id in gene_data["transcripts"].keys()
        ]
    )

    ax.xaxis.set_major_locator(
        ticker.MaxNLocator(integer=True)
    )  # Ensure genomic positions are integers
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, pos: f"{int(x)}")
    )  # Format x-axis ticks as integers

    # Plot each transcript
    for i, (transcript_id, transcript_info) in enumerate(
        gene_data["transcripts"].items()
    ):
        # Determine the direction based on the gene's strand information
        direction_marker = ">" if gene_data["strand"] == "+" else "<"
        marker_pos = (
            transcript_info["end"] + 100
            if gene_data["strand"] == "+"
            else transcript_info["start"] - 100
        )
        ax.plot(
            marker_pos,
            i,
            marker=direction_marker,
            markersize=5,
            color="blue",
        )

        # Draw the line for the whole transcript
        ax.plot(
            [transcript_info["start"], transcript_info["end"]],
            [i, i],
            color="grey",
            linewidth=2,
        )

        # Exon blocks
        for exon in transcript_info["exons"]:
            exon_length = exon["end"] - exon["start"]
            ax.add_patch(
                plt.Rectangle(
                    (exon["start"], i - 0.4),
                    exon_length,
                    0.8,
                    color="skyblue",
                )
            )

    ax.set_xlim(gene_data["start"], gene_data["end"])
    ax.tick_params(axis='both', which='major', labelsize=8)

    ax.invert_yaxis()

    if show_xlabel:
        ax.set_xlabel("Chromosomal position")
    else:
        ax.set_xlabel("")

def _plot_transcript_map_helper(obj:isoquantViewer|str,
            use_transcript_model:bool=False,
            Ensembl_ID:List[str]=None,
            gene_names:List[str]=None):
        
        if use_transcript_model:
            if not os.path.exists(obj.genedb_filename_model):
                obj.create_db(use_ref=True)
            db = gffutils.FeatureDB(obj.genedb_filename_model)
            dname = 'transcript model databse'
        else:
            if not os.path.exists(obj.genedb_filename):
                obj.create_db(use_ref=False)
            db = gffutils.FeatureDB(obj.genedb_filename)
            dname = 'reference transcripts database'
        
        if isinstance(obj,str):
            if os.path.exists(obj):
                mdata = mudata.read_h5ad(obj)

        elif hasattr(obj, "mdata"):
            mdata = obj.mdata 

        else:
            raise ValueError(f"Cant locate mdata.")
        
        if (Ensembl_ID is None) and (gene_names is None):
            raise ValueError(f"Need to provide either  Ensembl_ID or gene_names.")

        if gene_names is not None:
            if use_transcript_model:
                b00l = np.isin(mdata['gene'].var['name'],gene_names)
                tmp_dict = mdata['gene'].var.loc[b00l,:]['name'].to_dict()
                
            else: 
                b00l = np.isin(mdata['reference_gene'].var['name'],gene_names)
                tmp_dict = mdata['reference_gene'].var.loc[b00l,:]['name'].to_dict()
            
            Ensembl_ID = list(tmp_dict.keys())

            notin = np.isin(np.array(list(tmp_dict.values())),np.array(gene_names),invert=True)
            notin = np.array(list(tmp_dict.values()))[notin]
            if len(notin)>0:
                logging.info(f'Skipping {notin}, not found in {dname} .var.')


        # genes = [g for g in (Ensembl_ID or []) if g in gene_dict]
        # if not genes:
        #     if gene_names is not None:
        #         raise ValueError(f"None of the provided gene_names were found in ({dname}).")
        #     raise ValueError(f"None of the provided Ensembl_ID were found in ({dname}).")

        # notin = np.isin(np.array(Ensembl_ID),np.array(genes),invert=True)
        # notin = np.array(Ensembl_ID)[notin]
        # if len(notin)>0:
        #     if gene_names is not None:
        #         notin = [tmp_dict[i] for i in notin]
        #     logging.info(f'Skipping {notin}, as they are not in the ({dname}).')

        return Ensembl_ID,db


def plot_transcript_map(
            obj:isoquantViewer|str,
            use_transcript_model:bool=False,
            Ensembl_ID:List[str]=None,
            gene_names:List[str]=None,
            figsize:tuple[float,float]=(8,3.5),
            savefig:bool=True,
            save_format:str='png'):
        """
        Plot transcript structures for specified genes.
        Parameters:
        ----------
        obj: IsoQuant object containing mdata and gene dictionaries
        use_transcript_model: bool, whether to use transcript model from IsoQuant
        Ensembl_ID: list of Ensembl ID to plot (Provide either Ensembl_ID or gene_names)
        gene_names: list of gene names to plot (Provide either Ensembl_ID or gene_names)
        """
        # Adapted from IsoQuant PlotOutputs.py  

        Ensembl_ID_,db = _plot_transcript_map_helper(obj,use_transcript_model,Ensembl_ID,gene_names)
        

        for g in Ensembl_ID_:
            gene_dict = create_db_genedict_from_geneid(g,db)
            n_tx = len(gene_dict[g]["transcripts"])
            total_height = max(3.0, n_tx * 0.3)


            fig_width = figsize[0]
            fig_height = figsize[1] if figsize[1] > total_height else total_height
            fig,ax = plt.subplots(figsize=(fig_width, fig_height))

            _draw_gene_on_ax(
                ax=ax,
                gene_id = g,
                gene_data=gene_dict[g],
                show_xlabel=True,
                strand_markers=True,
            )

            ax.set_frame_on(False)
            ax.set_xlim(gene_dict[g]["start"], gene_dict[g]["end"])


            if savefig:
                plot_output = './plot_output'
                if os.path.exists(plot_output) == False:
                    os.makedirs(plot_output,exist_ok=True)
                ofname = os.path.join(plot_output,f'TranscriptMap_{g}_TranscriptModel_{use_transcript_model}.{save_format}')
                plt.savefig(ofname, format=save_format, dpi=144, bbox_inches='tight')
                plt.show()
                plt.close()


# --- helpfer functions for genomic region plotting ---
def _parse_region(s):
    """Parse a genomic region string into chromosome, start, and end."""
    chrom, positions = s.split(":")
    start, end = positions.split("-")
    return chrom, int(start), int(end)

def genomic_formatter(x, pos):
    if x >= 1_000_000:
        return f"{x/1_000_000:.4f} Mb"
    elif x >= 1_000:
        return f"{x/1_000:.3f} kb"
    return f"{int(x):,}"

def plot_genomic_region(
        obj:isoquantViewer|str,
        region: Optional[str] = None,
        Ensembl_ID:Optional[str] = None,
        gene_name: Optional[str] = None,
        plot_tracks: List[str]=['ref','model','reads'],
        tracks_params:Optional[Dict] = None,
        tracks_file:Optional[Dict] = dict(),
        padding:int=1000,
        figsize:tuple[float,float]=(10,18),
        label_fontsize: int=5,
        savefig:bool=True,
        save_format:str='png'
    ):
    """
    Plot a genomic region using pyGenomeTracks classes

    Parameters:
    ----------
    - obj: IsoQuant object
    - region: str, genomic region in the format "chr:start-end" e.g. '1:6625000-6635500'
    - plot_tracks: list of str, order of tracks to plot
        - 'ref' the reference GTF
        - 'model' the IsoQuant transcript model GTF
        - 'reads' IsoQuant corrected alignment BED
    - tracks_params: dict, optional, to override default track properties by setting key-value pairs for any supported track property.
        e.g. tracks_params = {'reads': {
                            'arrowhead_fraction':0.5
                            }}
    - tracks_file: dict, optional, If provided, the track will be read from this dictionary instead.
        e.g. tracks_file = {'ref': 'ref.gff3'}
    - padding: int, padding to add if Ensembl_ID or gene_name is used
    - figsize : tuple, figure size (width, height)
    - label_fontsize: int for x-axis and y-axis labels
    
    """


    if (region is None) and (Ensembl_ID is None) and (gene_name is None):
        raise ValueError(f"Need to provide either region,Ensembl_ID or gene_names")
    
    if isinstance(obj,str):
        if os.path.exists(obj):
            mdata = mudata.read_h5ad(obj)   
    elif isinstance(obj,isoquantViewer) and hasattr(obj, "mdata"):
        mdata = obj.mdata
    else:
        raise ValueError(f"Cant locate mdata.")


    if region is None:
        if gene_name is not None: ### check in reference first, if doesnt exist go to the transcript model
            b00l = np.isin(mdata['reference_gene'].var['name'],[gene_name])
            if np.sum(b00l) == 0: 
                b00l_1 = np.isin(mdata['gene'].var['name'],[gene_name])
                if np.sum(b00l_1) == 0: 
                    raise ValueError(f"Can not find {gene_name} in either reference_gene.var or gene.var")
                else:
                    Ensembl_ID = mdata['gene'].var.loc[b00l_1,:].index.values[0]
                    if np.sum(b00l_1)>1:
                        logging.info(f'More than one entry identified for {gene_name}, proceed with the entry with Ensembl ID: {Ensembl_ID}.')
                    # entry = obj.gene_dict_model.get(Ensembl_ID)
                    entry = mdata['gene'].var.loc[Ensembl_ID,:]
            else:
                Ensembl_ID = mdata['reference_gene'].var.loc[b00l,:].index.values[0]
                if np.sum(b00l)>1:
                        logging.info(f'More than one entry identified for {gene_name}, proceed with the entry with Ensembl ID: {Ensembl_ID}.')
                # entry = obj.gene_dict_ref.get(Ensembl_ID)
                entry = mdata['reference_gene'].var.loc[Ensembl_ID,:]

        elif Ensembl_ID is not None:
            
            # entry = obj.gene_dict_ref.get(Ensembl_ID)
            entry = mdata['reference_gene'].var.loc[Ensembl_ID,:]
            if entry is None:
                # entry = obj.gene_dict_model.get(Ensembl_ID)
                entry = mdata['gene'].var.loc[Ensembl_ID,:]
                if entry is None:
                    #  raise ValueError(f"Can not find {Ensembl_ID} in either reference (gene_dict_ref) or transcript model (gene_dict_model).")
                    raise ValueError(f"Can not find {Ensembl_ID} in either reference (reference_gene.var) or transcript model (gene.var)")


        chrom, start, end = entry['chromosome'],entry['start'],entry['end']
        region_str = f'{chrom}:{int(start)}-{int(end)}'
        start-=padding
        end+=padding
    
    else:
        region_str = region
        chrom, start, end = _parse_region(region)

    plot_ref = False
    plot_model = False
    plot_reads = False

    gtf_ref = obj.referene_gtf
    print('gtf_ref:',gtf_ref)
    gtf_model = obj.transcript_model
    print('gtf_model:',gtf_model)
    bedfile = f'{obj.output_directory}/{obj.prefix}/{obj.prefix}.corrected_reads.bed.gz'

    if 'ref' in plot_tracks: 
        plot_ref = True
        tmp = gtf_ref if os.path.exists(gtf_ref) else tracks_file.get('ref')
        if tmp is None or not os.path.exists(tmp):
            raise ValueError(f"Neither {gtf_ref} nor {tmp} exists.")

    if 'model' in plot_tracks:
        plot_model = True
        tmp = gtf_model if os.path.exists(gtf_model) else tracks_file.get('model')
        if tmp is None or not os.path.exists(tmp):
            raise ValueError(f"Neither {gtf_model} nor {tmp} exists.")

    if 'reads' in plot_tracks:
        plot_reads = True
        tmp = bedfile if os.path.exists(bedfile) else tracks_file.get('reads')
        if tmp is None or not os.path.exists(tmp):
            raise ValueError(f"Neither {bedfile} nor {tmp} exists.")


# --- Extract region slices for faster plotting ---
    tmp_gtf_ref = 'ref.slice.gtf'
    tmp_gtf_model = 'model.slice.gtf'
    tmp_bed = 'reads.slice.bed'
   
    str_ = f'{chrom} {int(start)} {int(end)}'
    region = BedTool(str_, from_string=True)

    print('str: ',str_)
    if plot_ref is True:
        BedTool(gtf_ref).intersect(region).saveas(tmp_gtf_ref)
    if plot_model is True:
        BedTool(gtf_model).intersect(region).saveas(tmp_gtf_model)
    if plot_reads is True:
        BedTool(bedfile).intersect(region,wa=True).saveas(tmp_bed)

    if tracks_params is not None:
        bed_props.update(tracks_params.get('reads',{}))
        gtf_props.update(tracks_params.get('model',{}))
        ref_props.update(tracks_params.get('ref',{}))

# --- Define track properties ---
    if plot_reads is True:
        bed_props['file'] =  tmp_bed
    
    if plot_model is True:
        gtf_props['file'] = tmp_gtf_model

    if plot_ref is True:
        ref_props['file'] = tmp_gtf_ref

    bed_track = BedTrack(bed_props) if plot_reads is True else None
    gtf_track = GtfTrack(gtf_props) if plot_model is True else None
    ref_track = GtfTrack(ref_props) if plot_ref is True else None

# --- Ordering tracks and count itemd for each track---
    tracks_dict = {
        'ref': ref_track,
        'model': gtf_track,
        'reads': bed_track,
        'scale': None}

    tracks = []
    for order in plot_tracks:
        track = tracks_dict[order]
        if track is not None:
            tracks.append((order,track))
    tracks.append(("scale", None))   # always last

    def safe_count(track):
        if track is None:
            return 0
        _, n = track.get_bed_handler()
        return n
    
    n_ref = safe_count(tracks_dict['ref'])
    n_gtf = safe_count(tracks_dict['model'])
    n_bed = safe_count(tracks_dict['reads'])

    if bed_track is not None:
        if n_bed > 50:
            bed_track.properties["labels"] = False      # turn labels off
            bed_track.properties["max_labels"] = 0      # (optional extra safety)
        else:
            bed_track.properties["labels"] = True
            # you can also limit how many labels at most:
            bed_track.properties["max_labels"] = 50

# --- Adjust the height ratio of each track based on the number of transcripts ---
    def adjust(ns):
        length = min(3.0, min(ns * 0.3, 12.0))
        return length

    def adjust_bed(ns):
        length = 3.0 +   math.log2(ns + 1)
        return length

    scale_units = 0.3

    height_units = []
    for ttype, track in tracks:
        if ttype == "ref":
            height_units.append(adjust(n_ref))
        elif ttype == "model":
            height_units.append(adjust(n_gtf))
        elif ttype == "reads":
            height_units.append(adjust_bed(n_bed))
        elif ttype == "scale":
            height_units.append(scale_units)

    total_units = sum(height_units)
    height_ratios = [u / total_units for u in height_units]


# --- Create figure and axes ---
    # total_height = sum(height_per_ax) 
    fig_height_in = 5.0 + 0.15 * total_units  
    fig_height_in = min(fig_height_in, 18.0)


    fig_width = figsize[0]
    fig_height = figsize[1] if figsize[1] > fig_height_in else fig_height_in

    fig = plt.figure(figsize=(fig_width, fig_height),dpi=300)

    gs = fig.add_gridspec(
        nrows=len(tracks), ncols=1,
        height_ratios=height_ratios,
        hspace=0.2
    )
    axes = [fig.add_subplot(gs[i, 0]) for i in range(len(tracks))]

# --- Plot actual tracks ---
    for (ttype, track), ax in zip(tracks, axes):
        if ttype in ("ref", "model", "reads"):
            track.plot(ax, chrom, start, end)
            ax.set_ylabel(track.properties.get('title',''),
                        labelpad=5, 
                        va='center',
                        fontsize=label_fontsize)

    for a in axes[1:]:
        a.sharex(axes[0])

    for ax in axes:
        ax.set_xlim(start, end)
    
# ---  scalebar  ---
    ax_scale = axes[-1]
    ax_scale.xaxis.set_major_formatter(ticker.FuncFormatter(genomic_formatter))
    ax_scale.set_xlabel(f'chr {chrom}',fontsize=label_fontsize)
    ax_scale.tick_params(axis='x', labelsize=5,width=0.1)

  
# ---  fix ticks and others  ---
    for ax in axes[:-1]:  # all but last
        ax.tick_params(labelbottom=False,length=0)
        ax.set_xlabel("") 

    for ax in axes:  
        ax.set_yticks([])
        ax.set_frame_on(False)

    if savefig:
        plot_output = './plot_output'
        if os.path.exists(plot_output) == False:
            os.makedirs(plot_output,exist_ok=True)
        ofname = os.path.join(plot_output,f'GenomeTrack_{region_str}.{save_format}')
        plt.savefig(ofname, format=save_format, dpi=144, bbox_inches='tight')
        plt.show()
        plt.close()

    for f in [tmp_gtf_ref,tmp_gtf_model,tmp_bed]:
        Path(f).unlink(missing_ok=True)