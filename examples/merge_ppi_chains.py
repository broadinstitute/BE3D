"""
File: merge_ppi_chains.py
Description:
    Post-processing utility to merge results from two or more be3d_local.py PPI-mode
    runs that each targeted a different chain of the same complex (e.g. KBTBD4_chain_B.yaml
    and HDAC1_chain_C.yaml, both run against the same user_pdb).

    Produces, per score type requested:
      - a combined TSV stacking each run's own-chain rows
      - a combined PDB with B-factors set per-chain from each run's own score table
        (built off one run's processed PDB, since it already contains every chain's
        coordinates from the shared input complex structure)

Usage:
    python merge_ppi_chains.py merge_config.yaml
"""

import os
import sys
import yaml
import pandas as pd
from biopandas.pdb import PandasPdb


def load_config(config_yaml):
    with open(config_yaml, "r") as file:
        return yaml.safe_load(file)


def signed_score(df, neg_col, pos_col):
    neg = pd.to_numeric(df[neg_col], errors='coerce')
    pos = pd.to_numeric(df[pos_col], errors='coerce')
    return neg.where(neg.notna(), pos).fillna(0.0)


def merge_tables(runs, output_tsv):
    tagged = []
    for run in runs:
        df = pd.read_csv(run['table'], sep='\t')
        df = df[df['chain'] == run['chain']].copy()
        df['merged_score'] = signed_score(df, run['neg_col'], run['pos_col'])
        df['gene'] = run['gene']
        tagged.append(df)

    merged = pd.concat(tagged, ignore_index=True, sort=False)
    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)
    merged.to_csv(output_tsv, sep='\t', index=False)
    return merged


def write_bfactor_pdb(base_pdb_path, out_pdb_path, merged_df):
    value_by_chain_pos = {
        (row['chain'], int(row['unipos'])): float(row['merged_score'])
        for _, row in merged_df.iterrows()
    }

    ppdb = PandasPdb().read_pdb(base_pdb_path)
    atom_df = ppdb.df['ATOM'].copy()
    atom_df['b_factor'] = atom_df.apply(
        lambda row: value_by_chain_pos.get((row['chain_id'], int(row['residue_number'])), 0.0),
        axis=1,
    )
    ppdb.df['ATOM'] = atom_df

    os.makedirs(os.path.dirname(out_pdb_path), exist_ok=True)
    ppdb.to_pdb(out_pdb_path)


def main(config):
    output_dir = config['output_dir']
    base_pdb = config['base_pdb']
    os.makedirs(output_dir, exist_ok=True)

    for score_cfg in config['score_types']:
        score_type = score_cfg['name']
        merged = merge_tables(score_cfg['runs'], os.path.join(output_dir, f'merged_{score_type}.tsv'))
        write_bfactor_pdb(base_pdb, os.path.join(output_dir, f'merged_{score_type}.pdb'), merged)
        print(f'[{score_type}] wrote merged_{score_type}.tsv and merged_{score_type}.pdb '
              f'({len(merged)} residues across {len(score_cfg["runs"])} chains)')


if __name__ == '__main__':
    config_yaml = sys.argv[1]
    config = load_config(config_yaml)
    main(config)
