from typing import List,Optional
from mudata import MuData
import mudata
import subprocess
import sys
import os
import gzip
import pandas as pd
import pickle
from anndata import AnnData
import numpy as np
from collections import Counter
from argparse import Namespace
import tempfile
import gffutils
from pathlib import Path
import polars as pl
import logging

class SafeUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module.startswith("src"):
            # replace with a harmless placeholder
            import types
            return types.SimpleNamespace
        return super().find_class(module, name)


def is_in_set(a, b):
    """Check if elements of array a are in set b, return boolean array."""
    set_b = set(b)
    return np.array([x in set_b for x in a])

class isoquantViewer:

    def __init__(self,
                 output_directory:Path|str,
                 prefix:str,
                 referene_gtf: Optional[Path] = None,
                 gene_db: Optional[Path] = None,
                 mudata_path: Optional[Path] = None):

        """Class to handle IsoQuant output files and data parsing.
        Parameters:
        -----------
        output_directory: Path or str
            Directory where IsoQuant output files are located.
        prefix: str
            Prefix used for IsoQuant output files.
        referene_gtf: Optional[Path]
            Path to the reference GTF file used by IsoQuant. If not provided, it will be loaded from the .params file used during IsoQuant run.
        gene_db: Optional[Path]
            Path to the gene database file used by IsoQuant. If not provided, it will be loaded from the .params file used during IsoQuant run.
        mudata_path: Optional[Path]
            Path to the MuData file containing counts and metadata. If not provided, it will be inferred from the output directory and prefix.
        """
        
        self.output_directory = Path(output_directory).expanduser()
        self.prefix = prefix
        self.referene_gtf = referene_gtf
        self.gene_db = gene_db
        self.plot_output = self.output_directory/'plot_output'
        self.plot_output = Path(self.plot_output).expanduser()
        transcript_model_reads_fname = self.output_directory/f'{prefix}/{prefix}.transcript_model_reads.tsv.gz'
        transcript_model_reads_fname = Path(transcript_model_reads_fname).expanduser()

        if mudata_path is None:
            mudata_path = Path(self.output_directory)/f"mudata/{prefix}_counts.h5mu"
            mudata_path = mudata_path.expanduser()
        self.mdata= mudata.read_h5mu(mudata_path)


        with gzip.open(transcript_model_reads_fname, "rt") as f: 
            comment_lines = []
            for line in f:
                if line.startswith("#"):
                    comment_lines.append(line)
                else:
                    break
        header_line = comment_lines[-1].lstrip("#").strip().split()

        self.transcript_model_reads = pd.read_csv(transcript_model_reads_fname, sep="\t", comment="#", header=None, compression="gzip")
        self.transcript_model_reads.columns = header_line

        self.transcript_model = f'{self.output_directory}/{self.prefix}/{self.prefix}.transcript_models.gtf'
        self.transcript_model = Path(self.transcript_model).expanduser()

        self.genedb_filename_model = f'{self.output_directory}/{self.prefix}/{self.prefix}.transcript_models.db'
        self.genedb_filename_model = Path(self.genedb_filename_model).expanduser()

        self._load_params_file()

        self.create_db(use_ref=True)
        self.create_db(use_ref=False)

        # self.get_assignment_df()


    def _load_params_file(self):
        """Load the .params file for necessary configuration and commands. From IsoQuant code."""
        params_path = self.output_directory/".params"
        params_path = Path(params_path).expanduser()
        assert os.path.exists(params_path), f"Params file not found: {params_path}"
        try:
            with open(params_path, "rb") as f:
                params = SafeUnpickler(f).load()
                if isinstance(params, Namespace):
                    self._process_params(vars(params))
                else:
                    print("Unexpected params format.")

        except Exception as e:
            raise ValueError(f"An error occurred while loading params: {e}")

    def _process_params(self, params):
        """Process parameters loaded from the .params file. From IsoQuant code."""

        self.genedb_filename = self.gene_db or params.get("genedb_filename")
        if os.path.exists(self.genedb_filename) == False:
            self.genedb_filename = self.output_directory/'geneannotations.db'
            self.genedb_filename = Path(self.genedb_filename).expanduser()
            
        self.referene_gtf = self.referene_gtf or params.get("genedb")
        self.referene_gtf = Path(self.referene_gtf).expanduser()

    def create_db(self,use_ref:bool=False):
        """Creating database based on GTF file using gffutils.
        Parameters:
        -----------
        use_ref: bool
            If True, parse the reference GTF; if False, parse the transcript model GTF.
        """

        if use_ref:
            if not os.path.exists(self.genedb_filename):
                logging.info("Creating reference GTF database")
                # convert GTF to DB if we use previous IsoQuant runs
                input_gtf_path = str(self.referene_gtf)
                gffutils.create_db(
                    input_gtf_path,
                    dbfn=str(self.genedb_filename),
                    force=True,
                    keep_order=True,
                    merge_strategy="merge",
                    sort_attribute_values=True,
                    disable_infer_genes=True,
                    disable_infer_transcripts=True,
                )
    
            # try:
            #     # Create a database without using a context manager
            #     db = gffutils.FeatureDB(self.genedb_filename, keep_order=True)  # read-only
    
            #     gene_dict = self._create_db_genedict(db)
    
            #     self.gene_dict_ref = gene_dict
    
            # except Exception as e:
            #     raise Exception(f"Error parsing GTF file: {str(e)}")

        else:
            if not os.path.exists(self.genedb_filename_model):
                logging.info("Creating transcript model GTF database")
                input_gtf_path = str(self.transcript_model)
                gffutils.create_db(
                                input_gtf_path,
                                dbfn=str(self.genedb_filename_model),
                                force=True,
                                keep_order=True,
                                merge_strategy="merge",
                                sort_attribute_values=True,
                                disable_infer_genes=True,
                                disable_infer_transcripts=True,
                            )


            # gene_dict = self._create_db_genedict(db)
    
            # self.gene_dict_model = gene_dict
    
        # if return_gene_dict:
        #     return gene_dict
        

    def get_assignment_df(self,
                          return_df:bool=False):
        """Load the read assignment TSV file into a DataFrame.
        Parameters:
        -----------
        return_df: bool
            If True, return the DataFrame containing read assignments.
        Returns:
        --------
        df: pl.DataFrame
            DataFrame containing read assignments (if return_df is True).
        """

        logging.info(f"Loading: read assignments")
        read_assign_fname =  f'{self.output_directory}/{self.prefix}/{self.prefix}.read_assignments.tsv.gz'


        # read comment header lines once
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

        scan = scan.with_columns([
            pl.col("additional_info").str.extract(r"(?:^|;)gene_assignment=([^;]*)", 1).alias("gene_assignment"),
            pl.col("additional_info").str.extract(r"(?:^|;)PolyA=([^;]*)", 1).alias("PolyA"),
            pl.col("additional_info").str.extract(r"(?:^|;)Classification=([^;]*)", 1).alias("Classification"),
        ])
    
        # df = scan.collect().to_pandas(use_pyarrow_extension_array=True)

        
        self.reads_assignment = scan.collect()
        if return_df:
            return scan

    def write_inut_for_drimseq(self):
        """Generate a transcript-level count file formatted for DRIMSeq analysis."""
    
        transcript_dict ={}

        db = gffutils.FeatureDB(self.genedb_filename_model)
        gene_dict = self._create_db_genedict(db)
        for _, vals in gene_dict.items():
            for t, val_t in vals['transcripts'].items():
                transcript_dict[t] = val_t
    
        count_fname = f'{self.output_directory}/{self.prefix}/{self.prefix}.transcript_model_grouped_tpm.tsv'
        count = pd.read_csv(count_fname,sep= '\t')
    
        geneid = [transcript_dict[i]['gene_id'] for i in count['#feature_id'].values]
        count['gene_id'] = geneid
    
        # organise the columns
        cols = count.columns
        b00l = np.isin(cols,['#feature_id','gene_id'],invert=True)
        cols_reorder = ['#feature_id','gene_id'] + cols[b00l].to_list()
        
        count = count.loc[:,cols_reorder]
        count = count.rename(columns={"#feature_id": 'feature_id'})
        count.to_csv(f'{self.output_directory}/{self.prefix}/{self.prefix}.transcript_model_grouped_tpm_drimseq_formatted.tsv',index=False,sep='\t')

    def extract_read_bam(self,
                        isoform_id:List[str],
                        use_transcript_model:bool=False,
                        output:str='read_extract.bam',
                       ):
        """Extract reads from the basecode.stitched.molecules.sorted.bam based on isoform IDs.
        Parameters:     
        ----------- 
        isoform_id: List[str]
            List of isoform IDs to extract reads for.
        use_transcript_model: bool
            If True, use isoform IDs based on the transcript model; if False, use isoform IDs based on the reference.
        output: str
            Path to the output BAM file where extracted reads will be saved.
        """

        if not hasattr(self, "reads_assignment") or self.reads_assignment is None:
            self.get_assignment_df()

        def find_read_ids(use_transcript_model:bool,isoform_id:list=None):
            """Helper function to find read IDs based on isoform IDs."""
            if use_transcript_model:   ### use the isoform ID based on the transcript model
                b00l = is_in_set(self.transcript_model_reads['transcript_id'],isoform_id)
                read_id = self.transcript_model_reads.loc[b00l,:]['read_id'].values
            else:                     ### use the isoform ID based on the reference
                # b00l = is_in_set(self.reads_assignment['isoform_id'],isoform_id)
                # read_id = self.reads_assignment.loc[b00l,'read_id'].values
                read_id = self.reads_assignment.filter(
                    pl.col("isoform_id").is_in(isoform_id)
                )["read_id"].to_list()

            return read_id

        read_id = find_read_ids(use_transcript_model,isoform_id)
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as f:
            for r in read_id:
                f.write(r)
                f.write('\n')
            readID_fname = f.name

        dir_ =  Path(self.output_directory).parent
        basecode_stitched_bam = dir_/f'{self.prefix}.stitched.molecules.sorted.bam'
        
        cmd_view = [
            "samtools", "view", "-b", "-N", readID_fname,
            "-o", output, basecode_stitched_bam
        ]
        cmd_index = ["samtools", "index", output]

        # Run commands in bash
        try:
            subprocess.run(cmd_view, check=True)
            subprocess.run(cmd_index, check=True)
        finally:
            if os.path.exists(readID_fname):
                os.remove(readID_fname)

        logging.info(f" Extracted reads written to: {output}")
        
        
    def _create_db_genedict(self,db):
        """Create a nested dictionary from a gffutils FeatureDB.
        Parameters:         
        -----------
        db: gffutils.FeatureDB
            The gffutils FeatureDB object to parse.     
        Returns:    
        --------
        gene_dict: dict
            A nested dictionary containing gene, transcript, and exon information.
        """
        gene_dict = {}
    
        # --- PASS 1: genes
        for g in db.features_of_type("gene"):
            attrs = g.attributes
            gene_dict[g.id] = {
                "chromosome": g.seqid,
                "start": g.start,
                "end": g.end,
                "strand": g.strand,
                "name": (attrs.get("gene_name") or [""])[0],
                "biotype": (attrs.get("gene_biotype") or [""])[0],
                "transcripts": {},
            }
    
        # --- PASS 2: transcripts
        t2g = {}
        for t in db.features_of_type("transcript"):
            attrs = t.attributes
            gene_id = (attrs.get("gene_id") or [None])[0]
            if gene_id is None:
                parents = list(db.parents(t, featuretype="gene"))
                gene_id = parents[0].id if parents else None
            if gene_id is None or gene_id not in gene_dict:
                continue  
    
            t2g[t.id] = gene_id
            gene_dict[gene_id]["transcripts"][t.id] = {
                "start": t.start,
                "end": t.end,
                "name": (attrs.get("transcript_name") or [""])[0],
                "biotype": (attrs.get("transcript_biotype") or [""])[0],
                "exons": [],
                "tags": ((attrs.get("tag") or [""])[0]).split(","),
                'gene_id':gene_id
            }
    
        # --- PASS 3: exons
        for e in db.features_of_type("exon"):
            attrs = e.attributes
            tid = (attrs.get("transcript_id") or [None])[0]
            if tid is None:
                continue
            gid = t2g.get(tid)
            if gid is None:
                continue 
            gene_dict[gid]["transcripts"][tid]["exons"].append({
                "exon_id": e.id,
                "start": e.start,
                "end": e.end,
                "number": (attrs.get("exon_number") or [""])[0],
            })
        return gene_dict
    

def create_db_genedict_from_geneid(geneid,
                                   db):
        """Create a nested dictionary from a gffutils FeatureDB.
        Parameters:         
        -----------
        db: gffutils.FeatureDB
            The gffutils FeatureDB object to parse.     
        Returns:    
        --------
        gene_dict: dict
            A nested dictionary containing gene, transcript, and exon information.
        """
        
        gene = db[geneid]
       
        gene_dict = {}
    
        # --- PASS 1: genes
        attrs = gene.attributes
        gene_dict[gene.id] = {
            "chromosome": gene.seqid,
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "name": (attrs.get("gene_name") or [""])[0],
            "biotype": (attrs.get("gene_biotype") or [""])[0],
            "transcripts": {},
        }
    
        # --- PASS 2: transcripts
        t2g = {}
        transcript = db.children(geneid, featuretype='transcript')

        for t in transcript:
            attrs = t.attributes
            gene_id = (attrs.get("gene_id") or [None])[0]
            if gene_id is None:
                parents = list(db.parents(t, featuretype="gene"))
                gene_id = parents[0].id if parents else None
            if gene_id is None or gene_id not in gene_dict:
                continue  
    
            t2g[t.id] = gene_id
            gene_dict[gene_id]["transcripts"][t.id] = {
                "start": t.start,
                "end": t.end,
                "name": (attrs.get("transcript_name") or [""])[0],
                "biotype": (attrs.get("transcript_biotype") or [""])[0],
                "exons": [],
                "tags": ((attrs.get("tag") or [""])[0]).split(","),
                'gene_id':gene_id
            }
    
        # --- PASS 3: exons
        exon = db.children(geneid, featuretype='exon')
        for e in exon:
            attrs = e.attributes
            tid = (attrs.get("transcript_id") or [None])[0]
            if tid is None:
                continue
            gid = t2g.get(tid)
            if gid is None:
                continue 
            gene_dict[gid]["transcripts"][tid]["exons"].append({
                "exon_id": e.id,
                "start": e.start,
                "end": e.end,
                "number": (attrs.get("exon_number") or [""])[0],
            })
        return gene_dict