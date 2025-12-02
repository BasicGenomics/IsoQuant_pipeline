# --- Define track properties ---
bed_props = {
    "file": None,
    "title": "Isoquant Corrected Alignment",
    "file_type": "bed",
    'height':3,
    'arrow_interval':1000,
    'color_arrow':'black',
    'style':'UCSC',
    "merge_transcripts": True,
    "merge_overlapping_exons": True,
    'fontsize':5,
    'line_width':0.1,
    'color_arrow':'black',
    'color':'#1f78b4'

} 

gtf_props = {
    "file": None,
    "title": "IsoQuant Transcript Model",
    "height": 0.3,
    "prefered_name": "transcript_id",
    "style": "UCSC",
    "labels": True,
    "display": "stacked",
    "fontsize":5
} 

ref_props = {
    "file": None,
    "title": "Reference",
    "height": 0.5,
    "style": "UCSC",
    "prefered_name": "transcript_id",
    "labels": True,
    "display": "stacked",
    'fontsize':5
}

bam_props = {
    "file":  None,
    "title": "BaseCode Alignment",
    "file_type": "bam",
    'height':7,
    'arrow_interval':1000,
    'color_arrow':'red',
    'style':'UCSC',
    "merge_transcripts": True,
    "merge_overlapping_exons": True,
    'fontsize':5
} 
