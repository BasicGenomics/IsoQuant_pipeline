import argparse
import gffutils
import tempfile
import os
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import sparse
from mudata import MuData
import anndata
anndata.settings.allow_write_nullable_strings = True

version = "1.1"

def parse_gtf(gtffile,genedb=None):
    if genedb is None:
        tmp_file = tempfile.NamedTemporaryFile(suffix=".db")
        db = gffutils.create_db(
                        gtffile,
                        dbfn=tmp_file.name,
                        force=True,
                        keep_order=True,
                        merge_strategy="merge",
                        sort_attribute_values=True,
                        disable_infer_genes=True,
                        disable_infer_transcripts=True,
                    )
    else:
        try:
            db = gffutils.FeatureDB(genedb, 
                                    sort_attribute_values=True,
                                    keep_order=True) 
        except Exception as e:
            return None,None
    gene_dict = {}
    transcript_dict = {}
    for gene in db.features_of_type("gene"):
        gene_id = gene.id
        gene_dict[gene_id] = {
            "chromosome": gene.seqid,
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "name": gene.attributes.get("gene_name", [""])[0],
            "biotype": gene.attributes.get("gene_biotype", [""])[0],
            "transcripts": [],
        }
        for transcript in db.children(gene, featuretype="transcript"):
            transcript_id = transcript.id
            gene_dict[gene_id]['transcripts'].append(transcript_id)
            transcript_dict[transcript_id] = {
                "chromosome": transcript.seqid,
                "start": transcript.start,
                "end": transcript.end,
                "gene_id": gene_id,
                "name": transcript.attributes.get("transcript_name", [""])[0],
                "biotype": transcript.attributes.get("transcript_biotype", [""])[0]
            }
        gene_dict[gene_id]['n_transcripts'] = len(gene_dict[gene_id]['transcripts'])
        gene_dict[gene_id]['transcripts'] = ','.join(gene_dict[gene_id]['transcripts'])
    return gene_dict, transcript_dict

def generate_anndata(X, obs_df, var_df):
    ad = anndata.AnnData(X, obs_df, var_df)
    return ad

def make_count_adata(count_tsv, tpm_tsv):
    count = pd.read_csv(count_tsv,sep='\t',index_col=0, comment=None)
    tpm = pd.read_csv(tpm_tsv,sep='\t',index_col=0, comment=None)
    count.index.name, tpm.index.name = 'feature_id','feature_id'
    tpm = tpm.loc[count.index,:]
    tpm = tpm.loc[:,count.columns]
    assert np.array_equal(count.index,tpm.index)
    assert np.array_equal(count.columns,tpm.columns)
    obs_names = count.columns
    feature_name = count.index
    adata = anndata.AnnData(X=count.values.T,
            layers={'count':count.T,
                    'tpm':tpm.T})
    adata.obs_names = obs_names
    adata.var_names = feature_name
    return adata

def load_matrix_modality(count_tsv, tpm_tsv, obs_names, var_dict):
    adata = make_count_adata(count_tsv, tpm_tsv)
    adata = adata[obs_names, :]
    assert np.array_equal(adata.obs_names, obs_names)
    if var_dict is not None:
        var_df = pd.DataFrame.from_dict(var_dict, orient='index').reindex(adata.var_names)
        adata.var = var_df
    return adata

def main():
    parser = argparse.ArgumentParser(description='Generate MuData object', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--disc-transcript-counts', required=True, type=str,
                        help='Discovered transcript grouped counts (.linear.tsv)')
    parser.add_argument('--disc-transcript-tpm', required=True, type=str,
                        help='Discovered transcript grouped TPM matrix (.tsv)')
    parser.add_argument('--disc-gene-counts', required=True, type=str,
                        help='Discovered gene grouped counts matrix (.tsv)')
    parser.add_argument('--disc-gene-tpm', required=True, type=str,
                        help='Discovered gene grouped TPM matrix (.tsv)')
    parser.add_argument('--ref-transcript-counts', type=str, default=None,
                        help='Reference transcript grouped counts matrix (.tsv)')
    parser.add_argument('--ref-transcript-tpm', type=str, default=None,
                        help='Reference transcript grouped TPM matrix (.tsv)')
    parser.add_argument('--ref-gene-counts', type=str, default=None,
                        help='Reference gene grouped counts matrix (.tsv)')
    parser.add_argument('--ref-gene-tpm', type=str, default=None,
                        help='Reference gene grouped TPM matrix (.tsv)')
    parser.add_argument('-m','--models', required=True, type=str,
                        help='Discovered transcript models (.gtf)')
    parser.add_argument('-g','--genedb', required=True, type=str,
                        help='IsoQuant reference gene database (.db)')
    parser.add_argument('-o','--output', required=True, type=str, help='Output .h5mu file')

    args = parser.parse_args()

    outfile = args.output

    _, disc_transcript_dict = parse_gtf(args.models)
    ref_gene_dict, ref_transcript_dict = parse_gtf(args.models, args.genedb)

    transcript_var_df = pd.DataFrame.from_dict(disc_transcript_dict, orient='index')

    linear_transcript_count_df = pd.read_csv(args.disc_transcript_counts, sep='\t')
    linear_transcript_count_sum_df = linear_transcript_count_df.groupby('feature_id').sum()
    transcript_subset = linear_transcript_count_sum_df.index[(linear_transcript_count_sum_df['count'] > 0.0).values.reshape(-1)]
    linear_transcript_count_df = linear_transcript_count_df[linear_transcript_count_df.apply(lambda row: row['feature_id'] in transcript_subset, axis=1)]
    linear_transcript_count_df = linear_transcript_count_df[linear_transcript_count_df['count'] > 0]
    linear_transcript_count_df = linear_transcript_count_df[~linear_transcript_count_df.apply(lambda row: 'Unassigned' in row['group_id'], axis=1)]
    transcript_var_df = transcript_var_df.reindex(linear_transcript_count_df['feature_id'].unique())
    transcript_to_col_dict = {transcript: j for (j, transcript) in enumerate(transcript_var_df.index)}
    cell_to_row_dict = {cell: k for (k, cell) in enumerate(linear_transcript_count_df['group_id'].unique())}
    linear_transcript_count_df['row_ind'] = linear_transcript_count_df.apply(lambda row: cell_to_row_dict[row['group_id']], axis=1)
    linear_transcript_count_df['col_ind'] = linear_transcript_count_df.apply(lambda row: transcript_to_col_dict[row['feature_id']], axis=1)
    obs_df = pd.DataFrame(index=linear_transcript_count_df['group_id'].unique())
    transcript_X = sparse.csr_matrix((linear_transcript_count_df['count'].values,(linear_transcript_count_df['row_ind'].values, linear_transcript_count_df['col_ind'].values)))

    transcript_adata = generate_anndata(transcript_X, obs_df, transcript_var_df)
    mdata = MuData({'isoform': transcript_adata})
 
    tpm = pd.read_csv(args.disc_transcript_tpm, sep='\t', index_col=0, comment=None)
    tpm = tpm.loc[mdata['isoform'].var_names, :]
    tpm = tpm.loc[:, mdata.obs_names]
    assert np.array_equal(tpm.index, mdata['isoform'].var_names)
    assert np.array_equal(tpm.columns, mdata['isoform'].obs_names)
    mdata['isoform'].layers['count'] = mdata['isoform'].X
    mdata['isoform'].layers['tpm'] = tpm.T

    obs_names = mdata.obs_names

    gene_adata = load_matrix_modality(args.disc_gene_counts, args.disc_gene_tpm, obs_names, ref_gene_dict)

    reference_isoform = None
    if args.ref_transcript_counts and args.ref_transcript_tpm:
        reference_isoform = load_matrix_modality(args.ref_transcript_counts, args.ref_transcript_tpm,
                                                 obs_names, ref_transcript_dict)

    reference_gene = None
    if args.ref_gene_counts and args.ref_gene_tpm:
        reference_gene = load_matrix_modality(args.ref_gene_counts, args.ref_gene_tpm,
                                              obs_names, ref_gene_dict)

    mods = {
        'gene': gene_adata,
        'isoform': mdata.mod['isoform'].copy(),
    }
    if reference_gene is not None:
        mods['reference_gene'] = reference_gene
    if reference_isoform is not None:
        mods['reference_isoform'] = reference_isoform

    mdata = MuData(mods)

    Path('/'.join(outfile.split('/')[:-1])).mkdir(parents=True, exist_ok=True)
    mdata.write(outfile)


if __name__ == '__main__':
    main()
    