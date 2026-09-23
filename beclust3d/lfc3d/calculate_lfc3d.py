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
    skip_no_coords=True,
    gene_type='Human',
    target_gene_chain = 'A',
    ppi_chain_gene_dict = {}, # {'GENE1':'B','GENE2':'C'}
    ppi_gene_edits_dict = {}, # {'GENE1': edits_dict, 'GENE2': edits_dict}
    func_map={'mean':np.mean, 'median':np.median, 'sum':np.sum, 'min':np.min, 'max':np.max},
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

    skip_no_coords : bool, optional (default=True)
        If True, blanks out both the 1D LFC and the LFC3D arms at residues that have no
        resolved xyz coordinates in df_struc (x_coord == '-'), setting LFC, LFC_Z, every LFCr,
        LFC3D and every LFC3Dr to '-' instead. Requires an 'x_coord' column.

        Such residues have no 3D neighborhood at all, so the neighbor search returns only the
        residue itself and LFC3D silently collapses to an exact copy of that residue's own 1D
        LFC -- a sequence-level score carrying a 3D label. Left ungated it propagates into the
        z-scores, p-values, meta-aggregation and union hit calls.

        The LFC arm is gated as well so that a residue absent from the model cannot be called
        a hit by ANY arm of a structure-based run: with both arms '-' their z/p/psig become
        '-', find_union returns '-', and they drop out of the per-screen and Meta union hit
        lists instead of sneaking back in through LFC alone.

        Gating LFC does not perturb any surviving residue's scores. A residue with no xyz
        never appears in another residue's Naa_pos, so it is never a neighbor and never feeds
        another residue's LFC3D; and the LFC null (mu, sigma) comes from the No_Mutation
        controls, not from these columns. Residues that do have coordinates are bit-identical
        with the flag on and off.

        Defaults to True: a residue that is not in the model should not be scored or called
        a hit by any arm, and the old default silently produced both. Because the gate needs
        to know which residues have coordinates, df_struc must carry an 'x_coord' column --
        pass the full *_coord_struc_features.tsv table. Set skip_no_coords=False to reproduce
        a run made before this became the default, or when df_struc genuinely has no
        coordinate columns.

        On an AlphaFold model every residue is coordinated, so the flag is a no-op there and
        only changes results for experimental PDBs with unresolved regions.

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

    # RESIDUES WITH NO RESOLVED xyz HAVE NO 3D NEIGHBORHOOD: THEIR Naa_pos IS '-', SO #
    # _resolve_neighbor_sources RETURNS ONLY THE SELF-ENTRY AND LFC3D DEGENERATES TO AN EXACT #
    # COPY OF THE RESIDUE'S OWN 1D LFC. WE BLANK BOTH ARMS -- LFC/LFC_Z/LFCr AS WELL AS #
    # LFC3D/LFC3Dr -- SO SUCH A RESIDUE CANNOT BE CALLED A HIT BY EITHER ARM, AND IN #
    # PARTICULAR CANNOT RE-ENTER THE union / Meta-union HIT LISTS THROUGH ITS LFC ARM. #
    # THIS PROPAGATES CLEANLY INTO THE z/p, NonAggr, META AND CLUSTERING STAGES. #
    no_coord_idx = set()
    if skip_no_coords:
        if 'x_coord' not in df_struc.columns:
            raise ValueError(
                "skip_no_coords=True requires an 'x_coord' column in df_struc; got columns "
                f"{list(df_struc.columns)}. Pass the full *_coord_struc_features.tsv table, "
                "or set skip_no_coords=False (note: True is the default as of 2026-09-23)."
            )
        for _idx, _val in enumerate(df_struc['x_coord']):
            if str(_val).strip() in ('-', '', 'nan', 'NaN', 'None'):
                no_coord_idx.add(_idx)
    # POSITIONAL MASK FOR THE 1D LFC COLUMNS: df_edits/df_rand ARE ROW-ALIGNED WITH df_struc #
    # (ASSERTED ABOVE), SO INDEX-FREE POSITIONAL MASKING IS SAFE. EMPTY UNLESS skip_no_coords. #
    no_coord_mask = np.zeros(len(df_struc), dtype=bool)
    if no_coord_idx:
        no_coord_mask[list(no_coord_idx)] = True

    
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
                                  _blank_no_coords(df_edits[lfc_colname], no_coord_mask).rename(f"{screen_name}_LFC"), 
                                  _blank_no_coords(df_edits[f'{lfc_colname}_Z'], no_coord_mask).rename(f"{screen_name}_LFC_Z")], axis=1)
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
                # NO COORDINATES -> NO NEIGHBORHOOD -> NO MEANINGFUL LFC3D #
                if aa in no_coord_idx:
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
            taa_LFC_dict = df_edits[lfc_colname].to_dict() ###

            for aa in range(len(df_edits)): # FOR EVERY RESIDUE # ###
                # RESIDUE NEEDS TO BE CONSERVED #
                if not aa_eligible[aa]:  ###
                    aggr_vals.append('-')
                    continue
                # CALCULATE LFC3D #
                taa_naa_LFC_vals = _gather_values(aa_sources[aa], taa_LFC_dict, ppi_edits_dict2)
                if len(taa_naa_LFC_vals) == 0:
                    aggr_vals.append('-')
                else:
                    aggr_vals.append(str(function_aggr_lfc3d(taa_naa_LFC_vals)))

            df_struct_3d = pd.concat([df_struct_3d, pd.DataFrame({f"{screen_name}_LFC3D": aggr_vals})], axis=1)
            del taa_LFC_dict, aggr_vals

        # REPEAT LFC LFC3D CALCULATIONS FOR RANDOMIZED DATA #

        dict_temp = {}
        for r in range(nRandom):
            # ADD LFC RANDOMIZATION COLUMNS FROM DF #
            dict_temp[f"{screen_name}_LFCr{str(r+1)}"] = _blank_no_coords(
                df_rand[f'{lfc_colname}r{str(r+1)}'], no_coord_mask)

            # CALCULATE LFC3D RANDOMIZATION, IF LFC_only SKIP OVER #
            if not LFC_only:
                aggr_vals = []
                taa_LFC_rand_dict = df_rand[f'{lfc_colname}r{str(r+1)}'].to_dict() ###

                for aa in range(len(df_rand)): # FOR EVERY RESIDUE # ###
                    # RESIDUE NEEDS TO BE CONSERVED #
                    if not aa_eligible[aa]: ###
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

def _blank_no_coords(series, no_coord_mask):
    """
    Blanks a 1D LFC column at the residues with no resolved xyz coordinates (see
    calculate_lfc3d's skip_no_coords). Returns the series untouched when nothing is masked,
    so the skip_no_coords=False path, and any fully-coordinated model, is a no-op.

    NaN rather than '-' is written because calculate_lfc3d normalizes the whole table with
    replace('-', nan) -> to_numeric -> fillna('-') later in the same loop iteration; the value
    reaches the output file as '-' either way.
    """
    if not no_coord_mask.any():
        return series
    return series.where(~no_coord_mask, np.nan)

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
