"""
File: preprocess_data.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
    Parses raw base editing screen data into per-mutation-type DataFrames for each screen.
"""

import os
import warnings
import pandas as pd
from pathlib import Path

from .preprocess_data_helpers import *

def parse_be_data(
    workdir, 
    input_dfs, 
    input_gene, 
    screen_names, 
    mut_col='Mutation category', 
    val_col='logFC', 
    gene_col='Target Gene Symbol', 
    edits_col='Amino Acid Edits', 
    mut_categories=["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"],
    mut_delimiter=',',
    conserv_dfs=[],
    conserv_col='mouse_res_pos',
    v_score_threshold=3, ### conserv_col
    gene_list=False,
    mutation_priority=None,
    control_category='No Mutation',
):
    """
    Parses raw base editing screen data into per-mutation-type DataFrames for each screen.

    Optionally filters mutations based on conservation scores from conservation DataFrames produced by conservation().

    Parameters
    ----------
    
    workdir : str
        Path to the working directory where output files and results will be saved.

    input_dfs : list of pd.DataFrame
        List of input dataframes, one for each screen, each containing mutation category, gene, and value columns.
    
    input_gene : str
        Name of the gene being processed (e.g., 'DNMT3A', 'MEN1'). 
        
    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in input_dfs, used in plot labels and output filenames.

    mut_col : str, optional (default='Mutation category')
        Column name in input_dfs specifying the mutation category (e.g., 'Missense', 'Nonsense').

    val_col : str, optional (default='logFC')
        Column name in input_dfs specifying the value measurement (e.g., log fold-change).

    gene_col : str, optional (default='Target Gene Symbol')
        Column name specifying the target gene name in input_dfs.

    edits_col : str, optional (default='Amino Acid Edits')
        Column name specifying the amino acid edits or mutation information in input_dfs.

    mut_categories : list of str, optional
        List of mutation categories to extract separately. Each category is
        written to its own ``{gene}_{screen}_{category}.tsv`` file, so this list
        also defines the vocabulary of edit types the parser will keep.

        The default ``["Nonsense", "Splice Site", "Missense", "No Mutation",
        "Silent"]`` matches classic base-editing dropout screens. Prime editing
        and other emerging screens install edit types outside this set --
        insertions, deletions, frameshifts and arbitrary stop-gain
        substitutions (Anzalone 2019; Erwood 2022; Cirincione 2024). To keep
        those, pass an extended vocabulary, e.g.
        ``[..., "Frameshift", "Insertion", "Deletion", "Stop-gain"]``.

        Any category token that appears in the input data but is NOT listed
        here is dropped -- but the parser now emits a ``warnings.warn`` naming
        the dropped categories instead of discarding them silently, so
        unexpected prime-editing categories are visible.

    control_category : str, optional (default='No Mutation')
        Name of the mutation category that represents the neutral / no-effect
        control group for this screen. In classic base-editing dropout screens
        this is guides with no coding consequence (``'No Mutation'``), but
        differently tokened screens may label their controls otherwise
        (e.g. ``'UTR'`` or ``'Intron'``). Downstream steps read the control
        distribution from ``{gene}_{screen}_{control_category}.tsv``; setting
        this to match your screen's token lets those steps find the control
        file instead of failing with FileNotFoundError. This category must also
        appear in ``mut_categories`` so that its file is written; a warning is
        raised if it does not. When left at the default the behavior is
        byte-identical to previous releases.

    mut_delimiter : str, optional (default=',')
        Delimiter used to separate multiple mutations within the edits_col field.

    conserv_dfs : list of pd.DataFrame, optional (default=[])
        List of conservation DataFrames, one per screen, used to optionally filter mutations based on conserved residues.

    conserv_col : str, optional (default='mouse_res_pos')
        Column name in conserv_dfs containing residue positions to filter on.
        
    v_score_threshold : int, optional (default=3)
        Conservation score for filtering. Scores are: -1 (not conserved), 1 (weakly similar), 2 (similar), 3 (conserved).
        
    gene_list : bool, optional (default=False)
        If True, processes a list of genes rather than a single gene.

    mutation_priority : list of str or None, optional (default=None)
        Priority order (most to least deleterious) used to collapse a guide's
        mut_col value into a single category when it is a delimiter-joined
        list of per-edit categories (e.g. 'Silent;Missense;'). If None, mut_col
        values are used as-is, assuming they are already single categories.

    Returns
    -------
    mut_dfs : dict
        Nested dictionary where:
          - Keys are screen names (from screen_names)
          - Values are dictionaries mapping mutation types (e.g., 'Missense') to processed DataFrames
            containing parsed mutation information and LFC values.
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'screendata'):
        os.mkdir(working_filedir / 'screendata')
    if not os.path.exists(working_filedir / 'screendata/plots'):
        os.mkdir(working_filedir / 'screendata/plots')

    # CHECK INPUTS ARE SELF CONSISTENT #
    for df in input_dfs: 
        assert mut_col in df.columns, 'Check [mut_col] input'
        assert val_col in df.columns, 'Check [val_col] input'
        assert gene_col in df.columns, 'Check [gene_col] input'
        assert edits_col in df.columns, 'Check [edits_col] input'
    # for df in conserv_dfs: 
    #     if df is not None: 
    #         assert conserv_col in df.columns, 'Check [conserv_col] input'
    #         assert conserv_score_col in df.columns, 'Check [conserv_col] input'

    assert len(input_dfs) == len(screen_names) == len(conserv_dfs), 'Lengths of [input_dfs] and [screen_names] and [conservation_dfs] must match'

    # THE NEUTRAL CONTROL MUST BE PART OF THE PARSED VOCABULARY SO ITS FILE IS WRITTEN #
    if control_category not in mut_categories:
        warnings.warn(
            f"control_category '{control_category}' is not in mut_categories "
            f"{mut_categories}; the neutral control file "
            f"'{{gene}}_{{screen}}_{control_category.replace(' ', '_')}.tsv' "
            f"will not be written and downstream prioritization may fail with "
            f"FileNotFoundError. Add '{control_category}' to mut_categories."
        )

    mut_dfs = {}
    # OUTPUT TSV BY INDIVIDUAL SCREENS #
    for input_gene, screen_df, screen_name, conserv_df in zip(gene_list, input_dfs, screen_names, conserv_dfs): 
        print('Processing', screen_name)
        # IF WE LOOK AT CONSERVATION #
        if conserv_df is not None:
            conserv_df['v_score'] = conserv_df['v_score'].astype(int)
            conserv_list = [str(x) for x in conserv_df[conserv_df['v_score']>=v_score_threshold][conserv_col].tolist()]
        # NARROW DOWN TO INPUT_GENE #
        df_gene = screen_df.loc[screen_df[gene_col] == input_gene, ]
        if mutation_priority:
            df_gene = df_gene.copy()
            df_gene[mut_col] = df_gene[mut_col].apply(
                lambda x: reduce_mutation_type(x, mut_delimiter, mutation_priority))

        # WARN (DON'T SILENTLY DROP) CATEGORIES PRESENT IN DATA BUT ABSENT FROM mut_categories #
        present_categories = set(df_gene[mut_col].dropna().unique())
        dropped_categories = present_categories - set(mut_categories)
        if dropped_categories:
            warnings.warn(
                f"parse_be_data: categories present in screen '{screen_name}' "
                f"but not in mut_categories were dropped: "
                f"{sorted(dropped_categories)}. Add them to mut_categories to "
                f"retain them (e.g. prime-editing Frameshift/Insertion/Deletion/"
                f"Stop-gain categories)."
            )

        mut_dfs[screen_name] = {}

        # NARROW DOWN TO EACH MUTATION TYPE #
        gene_mut_df = {}
        for mut in mut_categories: 

            # MAKE SURE MUT CATEGORY APPEARS IN DF #
            if not mut in df_gene[mut_col].unique(): 
                warnings.warn(f'{mut} not in Dataframe')
                continue

            # IF USER WANTS TO CATEGORIZE BY ONE SINGLE MUTATION PER GUIDE OR MULTIPLE MUTATIONS PER GUIDE #
            df_mut = df_gene.loc[df_gene[mut_col] == mut, ]
            df_mut = df_mut.reset_index(drop=True)
            gene_mut_df[mut] = len(df_mut)

            # ASSIGN position refAA altAA #
            df_mut[edits_col] = df_mut[edits_col].str.strip(mut_delimiter) # CLEAN
            df_mut[edits_col] = df_mut[edits_col].str.split(mut_delimiter) # STR to LIST
            df_mut[edits_col] = df_mut[edits_col].apply(lambda xs: identify_mutations(xs)) # FILTER FOR MUTATIONS ONLY #

            df_exploded = df_mut.explode(edits_col) # EACH ROW IS A MUTATION #
            df_exploded['edit_pos'] = df_exploded[edits_col].str.extract('(\d+)')
            df_exploded['refAA'] = df_exploded[edits_col].str.extract('([A-Za-z*]+)')
            df_exploded['altAA'] = df_exploded[edits_col].str.extract('[A-Za-z]+\d+([A-Za-z*]+)$')
            # IF 3 LETTER CODES ARE USED, TRANSLATE TO 1 LETTER CODE #
            df_exploded['refAA'] = df_exploded['refAA'].str.upper().apply(lambda x: aa_map.get(x, x))
            df_exploded['altAA'] = df_exploded['altAA'].str.upper().apply(lambda x: aa_map.get(x, x))

            # FILTER OUT SCORES WHERE POS DOES NOT APPEAR IN CONSERVED #
            if conserv_df is not None: 
                df_exploded = df_exploded[df_exploded['edit_pos'].isin(conserv_list) | df_exploded['edit_pos'].isna() | (df_exploded['edit_pos'] == "")]

            df_subset = df_exploded[[edits_col, 'edit_pos', 'refAA', 'altAA', val_col]]
            df_subset = df_subset.rename(columns={edits_col: 'this_edit', val_col: 'LFC'})

            # FOR PARTICULAR MUTATIONS, NEED TO SUBSET FURTHER #
            if mut == 'Missense': 
                df_subset = df_subset[(df_subset['refAA'] != df_subset['altAA']) & (df_subset['altAA'] != '*')]
            elif mut == 'Silent': # SILENT BEFORE NONSENSE (ie *248* MUTATION IS SILENT NOT NONSENSE)
                df_subset = df_subset[df_subset['refAA'] == df_subset['altAA']]
            elif mut == 'Nonsense': 
                df_subset = df_subset[df_subset['altAA'] == '*']
            else: 
                df_subset = df_subset[df_subset['LFC'] != df_subset['LFC'].shift()]
                df_subset = df_subset['LFC']

            # WRITE LIST OF MUT AND THEIR LFC VALUES #
            screen_name_nospace, mut_nospace = screen_name.replace(' ','_'), mut.replace(' ','_')
            edits_filename = f"screendata/{input_gene}_{screen_name_nospace}_{mut_nospace}.tsv"
            df_subset.to_csv(working_filedir / edits_filename, sep='\t')
            
            mut_dfs[screen_name][mut] = df_subset
        
        print(gene_mut_df) # OUTPUT COUNTS PER GENE FOR EA MUTATION #

    return mut_dfs

def sanitary_check(df_struc, df_missense_list, mute=True, screen_names=None, warn_fraction=0.05):
    """
        Check how the number of missense edits mapped to the target protein.

        Parameters
        ----------
        df_struc : pd.DataFrame
            Dataframe for target structure-sequence information.

        df_missense_list : list of pd.DataFrame
            List of missense dataframes, one for each screen. Screens numbered on another
            species' sequence (cross-species) should be left out, since they are only mapped
            onto df_struc later, in prioritize_by_sequence.

        mute : bool, optional (default=True)
            If False, prints the mapped / unmapped counts for every screen.

        screen_names : list of str or None, optional
            Names for df_missense_list, used in the printed report and warnings.

        warn_fraction : float, optional (default=0.05)
            Warn when more than this fraction of a screen's distinct missense edits have a
            reference residue/position that is not in df_struc -- usually a sign that the
            reference sequence (structure/UniProt/user_fasta) is not the one the screen
            library was designed on.

        Returns
        -------
        """    
    struc_refAA_pos_set = set((df_struc['unires']+df_struc['unipos'].astype(str)).to_list())
    if screen_names is None: 
        screen_names = [f'screen {i+1}' for i in range(len(df_missense_list))]

    for screen_name, each_df_missense in zip(screen_names, df_missense_list):
        missense_refAA_pos_set = set(each_df_missense['this_edit'].str[:-1].to_list())
        not_mapped = missense_refAA_pos_set.difference(struc_refAA_pos_set)
        if not mute: 
            print('-----[SANITARY CHECK]-----')
            print(f'{screen_name}: #of missense edits:{len(missense_refAA_pos_set)},\
                  #of mapped missense edits:{len(missense_refAA_pos_set) - len(not_mapped)},\
                  #of not mapped missense edits:{len(not_mapped)},\
                  list of not mapped missense edits: {list(not_mapped)}')
        if missense_refAA_pos_set and len(not_mapped) / len(missense_refAA_pos_set) > warn_fraction: 
            # PRINTED, NOT warnings.warn: SEVERAL PIPELINE MODULES TURN ALL WARNINGS OFF AT IMPORT #
            print(f'WARNING: {screen_name}: {len(not_mapped)}/{len(missense_refAA_pos_set)} missense edits do not match '
                  f'the reference sequence (e.g. {sorted(not_mapped)[:5]}); check that the reference '
                  'sequence (see sequence_source / user_fasta) matches the screen library numbering')

def check_screen_residues(
    workdir, 
    input_gene, 
    screen_names, 
    pdb_processed_file, 
    target_chainid, 
    df_struc, 
    gene_list=None, 
    conserv_dfs=None, 
    mut_categories=('Missense', 'Silent', 'Nonsense'), 
    on_mismatch='warn', 
): 
    """
    Description
        Guardrail comparing every reference residue/position a screen edits (e.g. the 'I257'
        of I257V, from parse_be_data's screendata/ tables) against the residue the PDB has at
        that position on target_chainid, and against the reference sequence in df_struc.
        Any disagreement means the screen library, the reference sequence and the structure
        are not numbered on the same sequence, and those edits' scores land on the wrong
        residues (or on none).

        Writes sequence_check/{input_gene}_screen_vs_structure.tsv (every edited position
        that is not a clean match) and sequence_check/{input_gene}_screen_vs_structure_summary.tsv
        (per-screen counts), and prints a WARNING per screen with any discrepancy.

    Parameters
    ----------
    workdir : str
        Output directory parse_be_data wrote screendata/ into.

    input_gene : str
        Gene whose own-species screens are checked.

    screen_names : list of str
        Screen identifiers, as passed to parse_be_data.

    pdb_processed_file : str
        Processed PDB from sequence_structural_features(_lite) (sequence_structure/{structureid}_processed.pdb).

    target_chainid : str
        Chain ID of input_gene in the PDB structure.

    df_struc : pd.DataFrame
        Residue table with 'unipos' and 'unires' (the reference sequence).

    gene_list : list of str or None, optional
        Per-screen gene symbol, as passed to parse_be_data. Screens whose entry is not
        input_gene (cross-species screens, still in the other species' numbering here) are
        listed in the summary as skipped. None checks every screen.

    conserv_dfs : list of (pd.DataFrame or None) or None, optional (default=None)
        Per-screen residue map, as passed to parse_be_data. A screen that has one is
        still in the alternative sequence's numbering at this point, so it is listed
        as skipped. This covers the case gene_list cannot: an alternative carrying the
        same gene symbol as input_gene -- another isoform of the same protein, say --
        where the two names are equal and the gene_list test never fires. None checks
        every screen.

    mut_categories : list of str, optional (default=('Missense', 'Silent', 'Nonsense'))
        parse_be_data categories whose tables carry per-edit refAA/edit_pos columns.

    on_mismatch : str, optional (default='warn')
        'warn' prints a warning and continues; 'error' raises ValueError after writing the report.

    Returns
    -------
    df_summary : pd.DataFrame
        One row per screen with match / discrepancy counts.

    Status values in the detailed table
    -----------------------------------
    pdb_mismatch       : the PDB has a different residue at this position
    reference_mismatch : the position is not resolved in the PDB, and the reference sequence
                         has a different residue there
    outside_reference  : the position is beyond the end of the reference sequence
    not_in_structure   : the reference sequence agrees but the PDB has no residue there
                         (reported, not counted as a discrepancy)
    """
    from .structure_helpers import extract_sequence_from_pdb

    assert on_mismatch in ('warn', 'error'), f"on_mismatch must be 'warn' or 'error', got '{on_mismatch}'"
    working_filedir = Path(workdir)
    os.makedirs(working_filedir / 'sequence_check', exist_ok=True)
    if gene_list is None: 
        gene_list = [input_gene] * len(screen_names)
    if conserv_dfs is None: 
        conserv_dfs = [None] * len(screen_names)
    assert len(conserv_dfs) == len(screen_names), '[conserv_dfs] must match [screen_names]'
    assert len(gene_list) == len(screen_names), '[gene_list] must match [screen_names]'

    resnums, pdb_seq, _ = extract_sequence_from_pdb(pdb_processed_file, target_chainid)
    pdb_res = dict(zip(resnums, pdb_seq))
    if 'chain' in df_struc.columns: # A COMPLEX'S TABLE ALSO CARRIES THE OTHER CHAINS' RESIDUES #
        df_struc = df_struc[df_struc['chain'].astype(str) == str(target_chainid)]
    ref_res = dict(zip(df_struc['unipos'].astype(int), df_struc['unires'].astype(str)))
    ref_len = max(ref_res) if ref_res else 0

    detail_rows, summary_rows = [], []
    for screen_name, screen_gene, conserv_df in zip(screen_names, gene_list, conserv_dfs): 
        if screen_gene != input_gene: 
            summary_rows.append({'gene': input_gene, 'chain': target_chainid, 'screen': screen_name,
                                 'checked': False, 'note': f'skipped: numbered on {screen_gene}'})
            continue
        if conserv_df is not None: 
            summary_rows.append({'gene': input_gene, 'chain': target_chainid, 'screen': screen_name,
                                 'checked': False,
                                 'note': 'skipped: numbered on the alternative sequence'})
            continue

        # DISTINCT (refAA, position) PAIRS THIS SCREEN EDITS, OVER EVERY PER-EDIT CATEGORY TABLE #
        edits = {}
        for mut in mut_categories: 
            path = working_filedir / f"screendata/{screen_gene}_{screen_name.replace(' ', '_')}_{mut.replace(' ', '_')}.tsv"
            if not os.path.exists(path): 
                continue
            df = pd.read_csv(path, sep='\t')
            if not {'this_edit', 'edit_pos', 'refAA'}.issubset(df.columns): 
                continue
            df = df.dropna(subset=['edit_pos', 'refAA'])
            for edit, pos, ref in zip(df['this_edit'], df['edit_pos'], df['refAA']): 
                edits.setdefault((str(ref), int(pos)), set()).add(str(edit))

        counts = dict.fromkeys(['match', 'pdb_mismatch', 'reference_mismatch', 'outside_reference', 'not_in_structure', 'stop_codon'], 0)
        for (ref, pos), edit_set in sorted(edits.items(), key=lambda x: x[0][1]): 
            if ref == '*': # STOP-LOSS EDITS (e.g. *1033R) SIT ON THE STOP CODON, WHICH HAS NO RESIDUE #
                counts['stop_codon'] += 1
                continue
            if pos in pdb_res: 
                status = 'match' if pdb_res[pos] == ref else 'pdb_mismatch'
            elif pos > ref_len: 
                status = 'outside_reference'
            elif ref_res.get(pos) != ref: 
                status = 'reference_mismatch'
            else: 
                status = 'not_in_structure'
            counts[status] += 1
            if status != 'match': 
                detail_rows.append({'gene': input_gene, 'chain': target_chainid, 'screen': screen_name,
                                    'position': pos, 'screen_res': ref, 'pdb_res': pdb_res.get(pos, '-'),
                                    'reference_res': ref_res.get(pos, '-'), 'status': status,
                                    'edits': ';'.join(sorted(edit_set))})

        n_discrepant = counts['pdb_mismatch'] + counts['reference_mismatch'] + counts['outside_reference']
        summary_rows.append({'gene': input_gene, 'chain': target_chainid, 'screen': screen_name, 'checked': True,
                             'n_positions': len(edits), 'n_discrepant': n_discrepant, **counts, 'note': ''})
        if n_discrepant: 
            examples = [f"{r['screen_res']}{r['position']} (PDB {r['pdb_res']}, reference {r['reference_res']})"
                        for r in detail_rows if r['screen'] == screen_name and r['status'] != 'not_in_structure'][:5]
            # PRINTED, NOT warnings.warn: SEVERAL PIPELINE MODULES TURN ALL WARNINGS OFF AT IMPORT #
            print(f'WARNING: {input_gene} chain {target_chainid}, {screen_name}: {n_discrepant}/{len(edits)} edited '
                  f'positions disagree with the structure/reference residue '
                  f'(pdb_mismatch={counts["pdb_mismatch"]}, reference_mismatch={counts["reference_mismatch"]}, '
                  f'outside_reference={counts["outside_reference"]}); e.g. {", ".join(examples)}')

    detail_cols = ['gene', 'chain', 'screen', 'position', 'screen_res', 'pdb_res', 'reference_res', 'status', 'edits']
    detail_file = working_filedir / f'sequence_check/{input_gene}_screen_vs_structure.tsv'
    summary_file = working_filedir / f'sequence_check/{input_gene}_screen_vs_structure_summary.tsv'
    pd.DataFrame(detail_rows, columns=detail_cols).to_csv(detail_file, sep='\t', index=False)
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(summary_file, sep='\t', index=False)

    n_total = int(df_summary['n_discrepant'].sum()) if 'n_discrepant' in df_summary else 0
    n_checked = int(df_summary['checked'].sum()) if 'checked' in df_summary else 0
    if n_total: 
        print(f'WARNING: {input_gene}: see {detail_file} for every discrepant position')
        if on_mismatch == 'error': 
            raise ValueError(f'{input_gene} chain {target_chainid}: {n_total} edited position(s) disagree with the '
                             f'structure/reference residue (on_residue_mismatch: error); see {detail_file}')
    elif n_checked: 
        print(f'{input_gene} chain {target_chainid}: every edited residue matches the structure/reference sequence')
    else: 
        # NOTHING WAS COMPARED, SO SAY SO RATHER THAN REPORT A CLEAN BILL OF HEALTH #
        print(f'{input_gene} chain {target_chainid}: no screen checked '
              f'(all {len(screen_names)} numbered on another sequence); see {summary_file}')
    return df_summary
