"""
File: af_structural_features.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
"""

import os
from pathlib import Path
import shutil

from .structure_helpers import *

def sequence_structural_features(
        workdir, 
        input_gene, input_uniprot, structureid, chains=['A'], radius=6.0, 
        user_fasta=None, user_pdb=None, user_dssp=None, 
        domains_dict=None, 
        # CHAINS REFERS TO CHAIN OF THE TARGET PROTEIN #
): 
    """
    Queries Uniprot, AlphaFold, DSSP, and domain features 
    Generate a combined sequence-structure feature table.

    Parameters
    ----------
    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    input_uniprot : str
        Uniprot of the gene being processed. 

    structureid : str
        Identifier for the protein structure (used for naming processed files).

    chains : list of str, optional (default=['A'])
        Only residues from these chains are included in downstream structural analyses.

    radius : float, optional (default=6.0)
        Radius in Angstroms used when counting neighboring amino acids.

    user_fasta : str or None, optional
        Path to a user-supplied UniProt FASTA file. If provided, bypasses querying UniProt online.

    user_pdb : str or None, optional
        Path to a user-supplied AlphaFold PDB file. If provided, bypasses querying AlphaFold online.

    user_dssp : str or None, optional
        Path to a user-supplied DSSP secondary structure file. If provided, bypasses running DSSP locally.

    domains_dict : dict or None, optional
        Dictionary specifying domain annotations (e.g., Pfam, InterPro) for the protein.
        If None, domain information will be queried automatically from UniProt or another source.

    Returns
    -------
    df_coord_dssp : pd.DataFrame
        DataFrame containing structural features like coordinates, secondary structure, surface exposure, etc.
    """

    # NAME VARIABLES, PATHS, CREATE DIRECTORIES #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'sequence_structure'):
        os.mkdir(working_filedir / 'sequence_structure')
    
    # UNIPROT #
    out_fasta = working_filedir / f"sequence_structure/{input_gene}_{input_uniprot}.tsv"
    if user_fasta is not None: # USER INPUT FOR UNIPROT #
        assert os.path.isfile(user_fasta), f'{user_fasta} does not exist'
        uFasta_file = user_fasta
    else: # QUERY DATABASE #
        uFasta_file = query_uniprot(working_filedir, input_uniprot)
    parse_uniprot(uFasta_file, out_fasta)

    # DOMAINS #
    domains_filename = f"sequence_structure/{input_gene}_{input_uniprot}_domains.tsv"
    if domains_dict is not None: 
        parse_domains(working_filedir, out_fasta, domains_filename, domains_dict)
    else: 
        query_domains(working_filedir, input_uniprot, domains_filename)

    # STRUCTURE #
    af_filename = f"sequence_structure/AF_{input_uniprot}.pdb"
    af_processed_filename = f"sequence_structure/{structureid}_processed.pdb"
    if user_pdb is not None: # USER INPUT FOR ALPHAFOLD #
        assert os.path.isfile(working_filedir / user_pdb), f'{user_pdb} does not exist'
        af_filename = user_pdb
    else: # QUERY DATABASE #
        query_af(working_filedir, af_filename, structureid)
    parse_af(working_filedir, af_filename, af_processed_filename)

    coord_filename = f"sequence_structure/{structureid}_coord.tsv"
    parse_coord(working_filedir, af_processed_filename, out_fasta, coord_filename, chains)

    # SECONDAY STRUCTURE DSSP #
    dssp_filename = f"sequence_structure/{structureid}_processed.dssp"
    if user_dssp is not None: # USER INPUT FOR DSSP #
        assert os.path.isfile(working_filedir / user_dssp), f'{user_dssp} does not exist'
        if str(user_dssp) != str(dssp_filename): 
            shutil.copy2(working_filedir / user_dssp, working_filedir / dssp_filename)
    else: # QUERY DATABASE #
        run_dssp(working_filedir, af_filename, dssp_filename)

    dssp_parsed_filename = f"sequence_structure/{structureid}_dssp_parsed.tsv"
    parse_dssp(working_filedir, dssp_filename, out_fasta, dssp_parsed_filename, chains)

    # UNIPROT AND PDB SHOULD MATCH AND THERE ARE CHECKS #
    # DSSP IS BASED ON PDB SO UNIPROT AND DSSP SHOULD ALSO MATCH #

    df_dssp = pd.read_csv(working_filedir / dssp_parsed_filename, sep = '\t')
    coord_radius_filename = f"sequence_structure/{structureid}_coord_radius.tsv"
    df_coord = count_aa_within_radius(working_filedir, coord_filename, coord_radius_filename, radius=radius)

    coord_dssp_filename = f"sequence_structure/{structureid}_coord_struc_features.tsv"
    df_coord_dssp = degree_of_burial(df_dssp, df_coord, working_filedir, coord_dssp_filename)
    return df_coord_dssp
