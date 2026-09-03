"""
File: calculate_lfc3d.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Calculates LFC3D scores by aggregating local neighborhood mutation effects.
             Translated from Notebook 3.3
"""

import pandas as pd
import numpy as np
from pathlib import Path
import os
import warnings
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# THERE ARE 2 MEAN FUNCTIONS, #
# MEAN FOR CALCULATING LFC3D WHICH IS TUNABLE, #
# AND MEAN FOR AVG RANDOMIZATIONS WHICH IS NOT TUNABLE #
def calculate_lfc3d(
    df_struc, 
    df_edits_list, 
    df_rand_list, 
    workdir, 
    input_gene, 
    screen_names, 
    nRandom=1000, 
    muttype='Missense', 
    function_type_lfc='mean', 
    function_type_lfc3d='mean',
    LFC_only=False, 
    conserved_only=False,
    gene_type='Human',
    target_gene_chain = 'A',
    ppi_chain_gene_dict = {}, # {'GENE1':'B','GENE2':'C'}
    ppi_gene_edits_dict = {}, # {'GENE1': edits_dict, 'GENE2': edits_dict}
    func_map={'mean':np.mean, 'median':np.median, 'sum':np.sum, 'min':np.min, 'max':np.max},
    min_neighbors=0,
):
    """
    Calculates LFC3D scores using structural data. 

    Parameters
    ----------
    df_struc : pd.DataFrame
        DataFrame containing structural conservation data for residues. 
        Must include columns 'unipos', 'unires', 'chain', 'Naa_pos', 'Naa_chain'.

    df_edits_list : list of pd.DataFrame
        List of mutation DataFrames for each screen. 

    df_rand_list : list of pd.DataFrame
        List of randomized mutation DataFrames for each screen.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in df_edits_list and df_rand_list.

    nRandom : int, optional (default=1000)
        Number of randomizations per screen for calculating randomized LFC and LFC3D scores.

    muttype : str, optional (default='Missense')
        Type of mutation to focus on (e.g., 'Missense', 'Nonsense', etc.).

    function_type_lfc : str, optional (default='mean')
        String label for the type of aggregation function used to compute LFC3D scores.

    function_type_lfc3d : str, optional (default='mean')
        String label for the type of aggregation function used to compute LFC3D scores.        

    LFC_only : bool, optional (default=False)
        If True, skips the LFC3D computation.

    conserved_only : bool, optional (default=False)
        If True, calculates LFC3D only for residues marked as 'conserved' in the conservation data.
        Non-conserved residues will be skipped (set to NaN or '-').

    min_neighbors : int, optional (default=0)
        Minimum number of structural *neighbours* (non-self sources) that must actually
        contribute a value for a residue's LFC3D to be emitted.

        The LFC3D score for a residue is an aggregation (default: mean) of the residue's OWN
        1-D LFC together with the LFCs of its structural neighbours. When a residue has NO
        contributing neighbour (e.g. positions outside the supplied structure, or in a
        fragment/disordered region with no resolved contacts), the aggregation collapses to a
        single value -- the residue's own LFC -- so ``LFC3D == LFC`` exactly. Such a "self-leak"
        residue is really a 1-D measurement masquerading as a 3-D (spatial) one, and it can be
        spuriously flagged as a 3-D hotspot. Stress-testing on structure-restricted / large
        proteins showed this produces many false hotspots at out-of-structure and poorly-packed
        residues.

        ``min_neighbors=0`` (the default) preserves the original behaviour exactly: every residue
        with at least one contributing source (self is always a source) gets an LFC3D value, and
        existing LFC3D values are unchanged. Set ``min_neighbors=1`` (or higher) for
        structure-restricted or large-protein runs so that residues with fewer than
        ``min_neighbors`` contributing neighbours have their LFC3D set to the missing sentinel
        ``'-'`` (and their randomised LFC3D likewise), and are therefore NOT emitted as spatial
        hotspots. The aggregation math for residues that DO have neighbours is never changed.

        Regardless of ``min_neighbors``, an additive per-residue, per-screen column
        ``{screen}_LFC3D_n_neighbors`` is written recording how many neighbours contributed, so
        downstream filtering (e.g. a validation-shortlist utility) can always exclude self-only
        calls even on data generated with the default settings.

    Returns
    -------
    df_struct_3d : pd.DataFrame
        DataFrame containing the structural data, LFC, LFC3D, and randomized scores.
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'LFC3D'):
        os.mkdir(working_filedir / 'LFC3D')

    # CHECK INPUTS ARE SELF CONSISTENT #
    for str_cons_df, str_cons_rand_df in zip(df_edits_list, df_rand_list): 
        assert len(df_struc) == len(str_cons_df) == len(str_cons_rand_df)
    assert 'unipos' in df_struc.columns and 'unires' in df_struc.columns and 'chain' in df_struc.columns
    # COLUMNS FOR SMOOTHING ACROSS RESIDUES AND#
    assert 'Naa_pos' in df_struc.columns
    assert 'Naa_chain' in df_struc.columns
    structure_columns = ['Naa_pos', 'Naa_chain']
    core_columns = ['unipos', 'unires', 'chain']

    assert len(df_edits_list) == len(df_rand_list) == len(screen_names)

    df_struct_3d = df_struc[core_columns + structure_columns].copy()
    df_struc = df_struc.fillna('-')
    
    naa_pos_chain_dict = dict()
    for idx, row in df_struc.iterrows():
        if row['Naa_pos'] == '-':
            naa_pos_chain_dict[f"{row['chain']}_{idx}"] = np.nan
        else:
            naa_chain_list = row['Naa_chain'].split(';')
            naa_pos_list = row['Naa_pos'].split(';')
            naa_chain_pos_list = list()
            for naa_chain, naa_pos in zip(naa_chain_list,naa_pos_list):
                naa_chain_pos_list.append(f'{naa_chain}_{naa_pos}')

            naa_pos_chain_dict[f"{row['chain']}_{idx}"] = ';'.join(naa_chain_pos_list)

    # MAP AGGREGATION FUNCTION #
    function_aggr_lfc3d = func_map[function_type_lfc3d]
    
    assert function_type_lfc3d in func_map.keys()
    # FOR EVERY SCREEN #
    for screen_name, df_edits, df_rand in zip(screen_names, df_edits_list, df_rand_list):
        ppi_edits_dict2 = dict()
        taa_conserv_dict = df_edits['conservation'].to_dict() ###
        # ADD LFC COLUMNS FROM DF #
        lfc_colname = f'{function_type_lfc}_{muttype}_LFC'
        df_struct_3d = pd.concat([df_struct_3d, 
                                  df_edits[[lfc_colname]].rename(columns={lfc_colname: f"{screen_name}_LFC"}), 
                                  df_edits[[f'{lfc_colname}_Z']].rename(columns={f'{lfc_colname}_Z': f"{screen_name}_LFC_Z"})], axis=1)
        if ppi_gene_edits_dict:
            for gene_identifier, be3d_dir in ppi_gene_edits_dict.items():
                _gene = gene_identifier.split('_')[0]
                ppi_edits_dict2[gene_identifier] = pd.read_csv(os.path.join(be3d_dir,'screendata_sequence',f'{_gene}_{screen_name}_protein_edits.tsv'),sep='\t')[lfc_colname].to_dict()

        # PRECOMPUTE NEIGHBOR SOURCES ONCE PER SCREEN, SINCE NEITHER THE NEIGHBOR TOPOLOGY NOR #
        # THE conserved_only/CHAIN-MATCHING GATING DEPENDS ON r; ONLY THE LOOKED-UP VALUES DO #
        if not LFC_only:
            aa_eligible = [True] * len(df_edits)
            aa_sources = [None] * len(df_edits)
            for aa in range(len(df_edits)):
                if conserved_only and taa_conserv_dict[aa] != 'conserved':  ###
                    aa_eligible[aa] = False
                    continue
                aa_sources[aa] = _resolve_neighbor_sources(
                    target_gene_chain, aa, taa_conserv_dict,
                    naa_pos_chain_dict[f'{target_gene_chain}_{aa}'],
                    conserved_only, ppi_chain_gene_dict,
                )

        # CALCULATE LFC3D, IF LFC_only SKIP OVER #
        if not LFC_only:
            aggr_vals = []
            n_neighbors_vals = []
            # WHICH RESIDUES FALL BELOW min_neighbors ON THE REAL DATA -- REUSED TO GATE THE #
            # RANDOMIZED PASSES SO THE NULL STAYS CONSISTENT WITH THE OBSERVED SCORE #
            aa_below_min = [False] * len(df_edits)
            taa_LFC_dict = df_edits[lfc_colname].to_dict() ###

            for aa in range(len(df_edits)): # FOR EVERY RESIDUE # ###
                # RESIDUE NEEDS TO BE CONSERVED #
                if not aa_eligible[aa]:  ###
                    aggr_vals.append('-')
                    n_neighbors_vals.append(0)
                    continue
                # COUNT CONTRIBUTING (NON-SELF) NEIGHBORS; SELF IS aa_sources[aa][0] #
                n_neigh = _count_contributing_neighbors(aa_sources[aa], taa_LFC_dict, ppi_edits_dict2)
                n_neighbors_vals.append(n_neigh)
                below_min = (min_neighbors >= 1) and (n_neigh < min_neighbors)
                aa_below_min[aa] = below_min
                # CALCULATE LFC3D #
                taa_naa_LFC_vals = _gather_values(aa_sources[aa], taa_LFC_dict, ppi_edits_dict2)
                if len(taa_naa_LFC_vals) == 0 or below_min:
                    # NO CONTRIBUTING SOURCE, OR TOO FEW NEIGHBORS: NOT A 3-D MEASUREMENT #
                    aggr_vals.append('-')
                else:
                    aggr_vals.append(str(function_aggr_lfc3d(taa_naa_LFC_vals)))

            df_struct_3d = pd.concat([df_struct_3d, pd.DataFrame({
                f"{screen_name}_LFC3D": aggr_vals,
                f"{screen_name}_LFC3D_n_neighbors": n_neighbors_vals,
            })], axis=1)
            del taa_LFC_dict, aggr_vals, n_neighbors_vals

        # REPEAT LFC LFC3D CALCULATIONS FOR RANDOMIZED DATA #

        dict_temp = {}
        for r in range(nRandom):
            # ADD LFC RANDOMIZATION COLUMNS FROM DF #
            dict_temp[f"{screen_name}_LFCr{str(r+1)}"] = df_rand[f'{lfc_colname}r{str(r+1)}']

            # CALCULATE LFC3D RANDOMIZATION, IF LFC_only SKIP OVER #
            if not LFC_only:
                aggr_vals = []
                taa_LFC_rand_dict = df_rand[f'{lfc_colname}r{str(r+1)}'].to_dict() ###

                for aa in range(len(df_rand)): # FOR EVERY RESIDUE # ###
                    # RESIDUE NEEDS TO BE CONSERVED #
                    if not aa_eligible[aa]: ###
                        aggr_vals.append('-')
                        continue
                    # RESIDUE BELOW min_neighbors ON REAL DATA IS EXCLUDED FROM THE NULL TOO #
                    if aa_below_min[aa]:
                        aggr_vals.append('-')
                        continue
                    # CALCULATE LFC3D RANDOMIZATION #
                    taa_naa_LFC_vals = _gather_values(aa_sources[aa], taa_LFC_rand_dict, ppi_edits_dict2)
                    if len(taa_naa_LFC_vals) == 0:
                        aggr_vals.append('-')
                    else:
                        aggr_vals.append(function_aggr_lfc3d(taa_naa_LFC_vals))

                dict_temp[f"{screen_name}_LFC3Dr{str(r+1)}"] = aggr_vals
                del aggr_vals

        df_struct_3d = pd.concat((df_struct_3d, pd.DataFrame(dict_temp)), axis=1)
        # CONVERT '-' TO NAN FOR EASIER CALCULATIONS #
        df_struct_3d = df_struct_3d.replace('-', np.nan).infer_objects(copy=False)
        df_struct_3d = df_struct_3d.apply(lambda col: pd.to_numeric(col, errors='coerce'))
        del dict_temp

        # AVG OVER LFC RANDOMIZATION COLUMNS FROM DF #
        LFC_colnames = [f"{screen_name}_LFCr{str(r+1)}" for r in range(nRandom)]
        df_struct_3d[f"{screen_name}_AVG_LFCr"]     = df_struct_3d[LFC_colnames].mean(axis=1) # AVG ALL
        df_struct_3d[f"{screen_name}_AVG_LFCr_neg"] = (df_struct_3d[LFC_colnames]
                                                       .apply(lambda col: col.map(lambda x: x if x < 0 else np.nan))
                                                       .mean(axis=1)) # AVG NEG
        df_struct_3d[f"{screen_name}_AVG_LFCr_pos"] = (df_struct_3d[LFC_colnames]
                                                       .apply(lambda col: col.map(lambda x: x if x > 0 else np.nan))
                                                       .mean(axis=1)) # AVG POS
        
        # AVG OVER LFC3D RANDOMIZATION COLUMNS FROM DF #
        if not LFC_only: 
            LFC3D_colnames = [f"{screen_name}_LFC3Dr{str(r+1)}" for r in range(nRandom)]
            df_struct_3d[f"{screen_name}_AVG_LFC3Dr"]     = df_struct_3d[LFC3D_colnames].mean(axis=1) # AVG ALL
            df_struct_3d[f"{screen_name}_AVG_LFC3Dr_neg"] = (df_struct_3d[LFC3D_colnames]
                                                            .apply(lambda col: col.map(lambda x: x if x < 0 else np.nan))
                                                            .mean(axis=1)) # AVG NEG
            df_struct_3d[f"{screen_name}_AVG_LFC3Dr_pos"] = (df_struct_3d[LFC3D_colnames]
                                                            .apply(lambda col: col.map(lambda x: x if x > 0 else np.nan))
                                                            .mean(axis=1)) # AVG POS
        
        # CONVERT NAN TO '-' FOR REPRESENTATION #
        df_struct_3d = df_struct_3d.fillna('-')
        print('Calculated LFC3D for', screen_name)

    df_struct_3d[core_columns + structure_columns] = df_struc[core_columns + structure_columns]
    out_filename = working_filedir / f"LFC3D/{gene_type}_{input_gene}_LFC_LFC3D_LFC3Dr.tsv.gz"
    df_struct_3d.to_csv(out_filename, sep = '\t', index=False, compression="gzip")

    return df_struct_3d

def _resolve_neighbor_sources(
    target_gene_chain, # should be main target gene
    aa, # should be main target gene
    df_struc_edits_dict, # only for main chain or target gene
    naa_chain_pos_str,  # only for main target gene
    conserved_only, # only for main target gene
    ppi_chain_gene_dict,
):
    """
    Resolves, for one residue, the fixed list of value sources (self + eligible neighbors)
    that feed into its LFC3D aggregation. Everything here is independent of which randomization
    column is being processed, so it's computed once per screen and reused for the real data
    pass and every one of the nRandom randomized passes (see _gather_values).

    Each source is either ('local', idx) meaning "look up taa_LFC_dict[idx]" (covers both the
    residue itself and same-chain neighbors), or ('cross', gene_identifier, idx) meaning
    "look up ppi_edits_dict[gene_identifier][idx]" (cross-chain PPI neighbor).
    """
    # naa IS NEIGHBORING AMINO ACIDS #
    # taa IS THIS AMINO ACID #
    # VALUE FOR THIS RESIDUE: caller already skips aa entirely when conserved_only excludes it #
    is_ppi_mode = isinstance(ppi_chain_gene_dict, dict)
    sources = [('local', aa)]

    # CHECK NEIGHBORING RESIDUES #
    if isinstance(naa_chain_pos_str, str):  ###
        naa_chain_pos_list = naa_chain_pos_str.split(';') ###
        for naa_chain_pos in naa_chain_pos_list:  ###
            naa_chain, naa_pos = naa_chain_pos.split('_')
            naa_idx = int(naa_pos) - 1

            if is_ppi_mode: # For PPIs
                if naa_chain == target_gene_chain:
                    if not conserved_only or df_struc_edits_dict[naa_idx] == 'conserved': ###
                        sources.append(('local', naa_idx))
                elif naa_chain in ppi_chain_gene_dict:
                    # CROSS-CHAIN PPI NEIGHBORS ARE NEVER conserved_only-GATED #
                    sources.append(('cross', ppi_chain_gene_dict[naa_chain], naa_idx))
                # ELSE: naa_chain IS A REAL CHAIN IN THE PDB BUT NOT LISTED IN ppi_chain_gene_dict --
                # IGNORE IT ENTIRELY RATHER THAN ERROR, SO CALLERS CAN OPT INTO A SUBSET OF CHAINS #
            else: # For Monomer
                if not conserved_only or df_struc_edits_dict[naa_idx] == 'conserved': ###
                    if naa_chain == target_gene_chain:
                        sources.append(('local', naa_idx))

    return sources

def _count_contributing_neighbors(sources, taa_LFC_dict, ppi_edits_dict):
    """
    Counts how many NON-SELF neighbor sources actually contribute a value for one residue.

    ``sources`` is the list produced by ``_resolve_neighbor_sources``; by construction its first
    element is always the residue itself (``('local', aa)``) and every remaining element is a
    structural neighbor. A neighbor "contributes" only if its looked-up LFC value is not the
    missing sentinel ``'-'`` (mirroring the filtering done in ``_gather_values``), so a residue
    surrounded only by missing/undefined neighbors correctly counts as zero.

    A count of 0 means the residue's LFC3D collapses to its own 1-D LFC (a self-leak), i.e. it is
    not really a 3-D measurement. This is what the ``{screen}_LFC3D_n_neighbors`` column records
    and what the ``min_neighbors`` gate thresholds on.
    """
    count = 0
    for source in sources[1:]:  # sources[0] IS SELF #
        if source[0] == 'cross':
            _, gene_identifier, idx = source
            val = ppi_edits_dict[gene_identifier][idx]
        else:
            _, idx = source
            val = taa_LFC_dict[idx]
        if val != '-':
            count += 1
    return count

def _gather_values(sources, taa_LFC_dict, ppi_edits_dict):
    """
    Fetches the actual LFC values for a precomputed source list. This is the only part of the
    per-residue computation that depends on r (which randomization column's values are used),
    including which entries happen to be '-' after permutation.
    """
    taa_naa_LFC_vals = []
    for source in sources:
        if source[0] == 'cross':
            _, gene_identifier, idx = source
            val = ppi_edits_dict[gene_identifier][idx] ###
        else:
            _, idx = source
            val = taa_LFC_dict[idx] ###
        if val != '-':
            taa_naa_LFC_vals.append(float(val))

    return taa_naa_LFC_vals
