import argparse
import gffutils
import tempfile
import pandas as pd
from scipy import sparse
import anndata
from pathlib import Path
from mudata import MuData

def parse_gtf(gtffile):
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
                "biotype": transcript.attributes.get(
                    "transcript_biotype", [""]
                )[0]
            }
        gene_dict[gene_id]['n_transcripts'] = len(gene_dict[gene_id]['transcripts'])
        gene_dict[gene_id]['transcripts'] = ','.join(gene_dict[gene_id]['transcripts'])
    return gene_dict, transcript_dict

def generate_anndata(X, obs_df, var_df):
    ad = anndata.AnnData(X, obs_df, var_df)
    return ad

def main():
    parser = argparse.ArgumentParser(description='Generate MuData object', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c','--counts',metavar='counts', type=str, help='Input counts file')
    parser.add_argument('-m','--models',metavar='models', type=str, help='Input gtf transcript models')
    parser.add_argument('-o','--output',metavar='output', type=str, help='Output .h5mu file')

    args = parser.parse_args()

    transcript_count_file = args.counts
    gtffile = args.models
    outfile = args.output

    gene_dict, transcript_dict = parse_gtf(gtffile)

    gene_var_df = pd.DataFrame.from_dict(gene_dict, orient='index')
    transcript_var_df = pd.DataFrame.from_dict(transcript_dict, orient='index')

    linear_transcript_count_df = pd.read_csv(transcript_count_file, sep='\t')

    linear_transcript_count_sum_df = linear_transcript_count_df.groupby('#feature_id').sum()
    transcript_subset = linear_transcript_count_sum_df.index[(linear_transcript_count_sum_df['count'] > 0.0).values.reshape(-1)]
    linear_transcript_count_df = linear_transcript_count_df[linear_transcript_count_df.apply(lambda row: row['#feature_id'] in transcript_subset, axis=1)]

    linear_transcript_count_df['gene_id'] = linear_transcript_count_df.apply(lambda row: transcript_dict[row['#feature_id']]['gene_id'], axis=1)
    linear_gene_count_df = linear_transcript_count_df.groupby('gene_id').apply(lambda gene_df: gene_df.groupby('group_id').sum().reset_index())

    linear_transcript_count_df = linear_transcript_count_df[linear_transcript_count_df['count'] > 0]
    linear_gene_count_df = linear_gene_count_df[linear_gene_count_df['count'] > 0]

    linear_transcript_count_df = linear_transcript_count_df[~linear_transcript_count_df.apply(lambda row: 'Unassigned' in row['group_id'], axis=1)]
    linear_gene_count_df = linear_gene_count_df[~linear_gene_count_df.apply(lambda row: 'Unassigned' in row['group_id'], axis=1)]

    transcript_var_df = transcript_var_df.reindex(linear_transcript_count_df['#feature_id'].unique())
    transcript_to_col_dict = {transcript: j for (j, transcript) in enumerate(transcript_var_df.index)}
    gene_var_df = gene_var_df.reindex(linear_gene_count_df['gene_id'].unique())
    gene_to_col_dict = {gene: i for (i, gene) in enumerate(gene_var_df.index)}

    cell_to_row_dict = {cell: k for (k, cell) in enumerate(linear_transcript_count_df['group_id_modified'].unique())}

    linear_transcript_count_df['row_ind'] = linear_transcript_count_df.apply(lambda row: cell_to_row_dict[row['group_id_modified']], axis=1)
    linear_transcript_count_df['col_ind'] = linear_transcript_count_df.apply(lambda row: transcript_to_col_dict[row['#feature_id']], axis=1)

    linear_gene_count_df['row_ind'] = linear_gene_count_df.apply(lambda row: cell_to_row_dict[row['group_id_modified']], axis=1)
    linear_gene_count_df['col_ind'] = linear_gene_count_df.apply(lambda row: gene_to_col_dict[row['gene_id']], axis=1)

    obs_df = pd.DataFrame(index=linear_transcript_count_df['group_id_modified'].unique())

    transcript_X = sparse.csr_matrix((linear_transcript_count_df['count'].values,(linear_transcript_count_df['row_ind'].values, linear_transcript_count_df['col_ind'].values)))
    gene_X = sparse.csr_matrix((linear_gene_count_df['count'].values,(linear_gene_count_df['row_ind'].values, linear_gene_count_df['col_ind'].values)))

    transcript_adata = generate_anndata(transcript_X, obs_df, transcript_var_df)
    gene_adata = generate_anndata(gene_X, obs_df, gene_var_df)
    mdata = MuData({'gene': gene_adata, 'isoform': transcript_adata})

    Path('/'.join(outfile.split('/')[:-1])).mkdir(parents=True, exist_ok=True)

    mdata.write(outfile)

if __name__ == '__main__':
    main()
    