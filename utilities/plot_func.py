import matplotlib.pyplot as plt
from typing import List,Optional,Dict,Literal
from collections import Counter
import numpy as np
import matplotlib.ticker as mticker

from pygenometracks.tracks.BigWigTrack import BigWigTrack
from pygenometracks.tracks.GtfTrack import GtfTrack
from pygenometracks.tracks.ScaleBarTrack import ScaleBarTrack
from pygenometracks.tracks.BedTrack import BedTrack

from pybedtools import BedTool
import math
from adjustText import adjust_text

def plot_pie_assignment(
    df,
    feature_to_plot: Optional[
        Literal["assignment_type", "gene_assignment", "Classification"]
    ],
):
    count = df[feature_to_plot].value_counts()
    labels = count.index
    sizes = count.values

    total = sum(sizes)
    percentages = [(s / total) * 100 for s in sizes]

    # Legend format: Label – 54.3% (12345)
    legend_labels = [
        f"{l} ({c})" 
        for l, c in zip(labels, sizes)
    ]

    fig, ax = plt.subplots(figsize=(5, 4))

    # Pie: only percent on wedges, no label text
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=None,
        autopct="%1.1f%%",
        pctdistance=0.8,
    )

    # Add legend with full info
    ax.legend(
        wedges,
        legend_labels,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
    )

    ax.set_title(f"Assignment \n total no. reads:{total}",)

def plot_count_bar(obj,
            layer:str,
            sample_id:List[str]=None,
            gene_list:List[str]=None,
            isoform_list:List[str]=None,
            use_transcript_model:bool=False,
):

    """
    Plot bar chart of counts for either specified genes or isoforms.
    If providing gene_list, plots stacked bar of isoforms per gene.
    If providing isoform_list, plots bar chart of specified isoforms.
    
    Parameters:
    ----------
    obj: IsoQuant object containing mdata and gene dictionaries
    layer: data layer to use for counts (e.g., 'count', 'tpm'). If using transcript model, only 'count' is used.
    sample_id:  list of sample IDs to include in the plot,  if None, use all samples
    gene_list: list of gene names to plot
    isoform_list: list of isoform IDs to plot
    use_transcript_model: bool, whether to use transcript model from IsoQuant
    """

# --- loading gene dict
    if use_transcript_model:
        if not hasattr(obj, "gene_dict_model"):
            obj.parse_input_gtf(use_ref=False)
            
        gene_dict = obj.gene_dict_model
    else:  
        if not hasattr(obj, "gene_dict_ref"):
            obj.parse_input_gtf(use_ref=True)

        gene_dict = obj.gene_dict_ref
        
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

            isoforms = np.array(list(gene_dict[g]['transcripts'].keys()))

            if use_transcript_model:
                if layer == 'count':
                    layer = ''
                X = obj.mdata['isoform'][sampld_id,isoforms].to_df(layer)
            else:
                X = obj.mdata['reference_isoform'][sample_id,isoforms].to_df(layer)
                
            xticks = X.index
            ylabel = 'Transcripts per million'
            
            
            fig, ax = plt.subplots(figsize=(3, 4))
            
            X.plot.bar(ax=ax,stacked=True)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=f'Transripts of {g}')
       
            

    elif is_nonempty(isoform_list):

        if use_transcript_model:
            if layer == 'count':
                layer = ''
            X = obj.mdata['isoform'][sample_id,isoform_list].to_df(layer)
        else: 
            X = obj.mdata['reference_isoform'][sample_id,isoform_list].to_df(layer)
        xticks = X.T.index
        ylabel = 'Counts'
            
        #Adjusting the figure width based on numbers of isoforms
        n_vars = len(xticks)
        base_width = 0.6 
        width = max(6, n_vars * base_width) 
        height = width*2/3
    
        fig, ax = plt.subplots(figsize=(width, height))

        X.T.plot.bar(ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
       
    ax.set_ylabel(ylabel)    
    ax.set_xticklabels(xticks, rotation=40, ha='right')
    ax.set_xlabel('Sample')

# --- helper function for plot_transcript_map, plot one gene per ax ---
def _draw_gene_on_ax(ax,
    obj,
    gene_data: Dict,
    filter_transcripts: Optional[float] = None,
    show_xlabel: bool = False,
    strand_markers: bool = True,
):

    transcripts = gene_data["transcripts"]
    if filter_transcripts is not None:
        pass

    num_transcripts = len(transcripts)
    ax.set_title(
        f"Transcripts of Gene: {gene_data['name']} on Chromosome {gene_data['chromosome']}"
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
    ax.invert_yaxis()

    if show_xlabel:
        ax.set_xlabel("Chromosomal position")
    else:
        ax.set_xlabel("")


def plot_transcript_map(obj,
            use_transcript_model:bool=False,
            gene_names:List[str]=None,
            filter_transcripts:Optional[float]=None):
        """
        Plot transcript structures for specified genes.
        Parameters:
        ----------
        obj: IsoQuant object containing mdata and gene dictionaries
        use_transcript_model: bool, whether to use transcript model from IsoQuant
        gene_names: list of gene names to plot
        filter_transcripts: float, optional threshold to filter transcripts by expression
        """
        # Adapted from IsoQuant PlotOutputs.py  

        # --- loading gene dict ---
        if use_transcript_model:
            if not hasattr(obj, "gene_dict_model"):
                obj.parse_input_gtf(use_ref=False)
                
            gene_dict = obj.gene_dict_model
        else:  
            if not hasattr(obj, "gene_dict_ref"):
                obj.parse_input_gtf(use_ref=True)

            gene_dict = obj.gene_dict_ref
        genes = [g for g in (gene_names or []) if g in gene_dict]
        if not genes:
            raise ValueError("None of the provided gene_names were found in the gene dictionary.")

        # height per axis proportional to number of transcripts
        height_per_ax = []
        for g in genes:
            n_tx = len(gene_dict[g]["transcripts"])
            height_per_ax.append(max(3.0, n_tx * 0.3))

        total_height = sum(height_per_ax) + 0.5 * (len(genes) - 1)
        fig = plt.figure(figsize=(12, total_height))
        gs = fig.add_gridspec(
            nrows=len(genes),
            ncols=1,
            height_ratios=height_per_ax,
            hspace=0.4,
        )

        axes = []
        for i, g in enumerate(genes):
            ax = fig.add_subplot(gs[i, 0])
            axes.append(ax)

            _draw_gene_on_ax(
                ax=ax,
                obj=obj,
                gene_data=gene_dict[g],
                filter_transcripts=filter_transcripts,
                show_xlabel=(i == len(genes) - 1),
                strand_markers=True,
            )
            ax.set_frame_on(False)
            ax.set_xlim(gene_dict[g]["start"], gene_dict[g]["end"])

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
        obj,
        region: str,
        plot_model: bool = True,
        plot_reads: bool = False,
        plot_BCreads: bool = False,
        labelsize: int = 8,
        title_fontsize: int= 10,
        ylab_fontsize: int=10,
        fig_width: int=10
    ):
    """
    Plot a genomic region using pyGenomeTracks classes

    Parameters:
    ----------
    - obj: IsoQuant object
    - region: str, genomic region in the format "chr:start-end"
    - plot_model: bool, whether to plot the IsoQuant transcript model GTF
    - plot_reads: bool, whether to plot the IsoQuant corrected alignment BED
    - labelsize: int 
    - title_fontsize: int
    - ylab_fontsize: int
    - fig_width: int
    """
    
    chrom, start, end = _parse_region(region)

    gtf_ref = obj.referene_gtf
    gtf_model = obj.transcript_model
    bedfile = f'{obj.output_directory}/{obj.prefix}/{obj.prefix}.corrected_reads.bed.gz'
    # dir_ =  Path(obj.output_directory).parent
    # basecode_stitched_bam = dir_/f'{obj.prefix}.stitched.molecules.sorted.bam'



# --- Extract region slices for faster plotting ---
    tmp_gtf_ref = 'ref.slice.gtf'
    tmp_gtf_model = 'model.slice.gtf'
    tmp_bed = 'reads.slice.bed'
    tmp_bam =  'reads.sorted.bam'
   
    str_ = f'{chrom} {start} {end}'
    region = BedTool(str_, from_string=True)
    if plot_model is True:
        BedTool(gtf_model).intersect(region).saveas(tmp_gtf_model)
    if plot_reads is True:
        BedTool(bedfile).intersect(region,wa=True).saveas(tmp_bed)
    if plot_BCreads is True:
        BedTool(basecode_stitched_bam).intersect(region,wa=True).saveas(tmp_bam)
    

    BedTool(gtf_ref).intersect(region).saveas(tmp_gtf_ref)

# --- Define track properties ---
    bed_props = {
        "file": tmp_bed,
        "title": "Isoquant Corrected Alignment",
        "file_type": "bed",
        'height':7,
        'arrow_interval':1000,
        'color_arrow':'red',
        'style':'UCSC',
        "merge_transcripts": True,
        "merge_overlapping_exons": True,
        'fontsize':labelsize
    } if plot_reads is True else None

    gtf_props = {
        "file": tmp_gtf_model,
        "title": "IsoQuant Transcript Model",
        "height": 0.3,
        "prefered_name": "transcript_id",
        "style": "UCSC",
        "labels": True,
        "display": "stacked",
        "fontsize":labelsize
    } if plot_model is True else None

    ref_props = {
        "file": tmp_gtf_ref,
        "title": "Reference",
        "height": 0.5,
        "style": "UCSC",
        "prefered_name": "transcript_id",
        "labels": True,
        "display": "stacked",
        'fontsize':labelsize
    }

    bam_props = {
        "file": tmp_bam,
        "title": "BaseCode Alignment",
        "file_type": "bam",
        'height':7,
        'arrow_interval':1000,
        'color_arrow':'red',
        'style':'UCSC',
        "merge_transcripts": True,
        "merge_overlapping_exons": True,
        'fontsize':labelsize
    } if plot_BCreads is True else None


    bed_track = BedTrack(bed_props) if plot_reads is True else None
    gtf_track = GtfTrack(gtf_props) if plot_model is True else None
    # bam_tract = BamTrack
    ref_track = GtfTrack(ref_props)  

# --- Ordering tracks and count itemd for each track---
    tracks = []
    tracks.append(("ref", ref_track))
    if gtf_track is not None:
        tracks.append(("gtf", gtf_track))
    if bed_track is not None:
        tracks.append(("bed", bed_track))
    tracks.append(("scale", None))   # always last

    def safe_count(track):
        if track is None:
            return 0
        _, n = track.get_bed_handler()
        return n
    
    n_ref = safe_count(ref_track)
    n_gtf = safe_count(gtf_track)
    n_bed = safe_count(bed_track)

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
        length = max(3.0, min(ns * 0.3, 12.0))
        return length

    def adjust_bed(ns):
        length = 3.0 +   math.log2(ns + 1)
        return length

    scale_units = 0.3

    height_units = []
    for ttype, track in tracks:
        if ttype == "ref":
            height_units.append(adjust(n_ref))
        elif ttype == "gtf":
            height_units.append(adjust(n_gtf))
        elif ttype == "bed":
            height_units.append(adjust_bed(n_bed))
        elif ttype == "scale":
            height_units.append(scale_units)

    total_units = sum(height_units)
    height_ratios = [u / total_units for u in height_units]


# --- Create figure and axes ---
    # total_height = sum(height_per_ax) 
    fig_height_in = 5.0 + 0.15 * total_units  
    fig_height_in = min(fig_height_in, 18.0)

    fig = plt.figure(figsize=(fig_width, fig_height_in),dpi=300)

    gs = fig.add_gridspec(
        nrows=len(tracks), ncols=1,
        height_ratios=height_ratios,
        hspace=0.2
    )
    axes = [fig.add_subplot(gs[i, 0]) for i in range(len(tracks))]

# --- Plot actual tracks ---
    for (ttype, track), ax in zip(tracks, axes):
        if ttype in ("ref", "gtf", "bed"):
            track.plot(ax, chrom, start, end)
            ax.set_ylabel(track.properties.get('title',''),
                        labelpad=5, 
                        va='center',
                        fontsize=ylab_fontsize)

    for a in axes[1:]:
        a.sharex(axes[0])
    
# ---  scalebar  ---
    ax_scale = axes[-1]
    ax_scale.xaxis.set_major_formatter(mticker.FuncFormatter(genomic_formatter))
    ax_scale.set_xlabel(f'chr {chrom}',fontsize=ylab_fontsize)
    ax_scale.tick_params(axis='x', labelsize=5,width=0.1)

  
# ---  fix ticks and others  ---
    for ax in axes[:-1]:  # all but last
        ax.tick_params(labelbottom=False,length=0)
        ax.set_xlabel("") 

    for ax in axes:  
        ax.set_yticks([])
        ax.set_frame_on(False)


    # fig.suptitle(f"chr{chrom}:{start:,}-{end:,}",fontsize=title_fontsize)
    plt.show()

