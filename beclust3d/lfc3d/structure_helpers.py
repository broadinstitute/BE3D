"""
File: structure_features.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
"""

import os
import math
import wget
import warnings
import requests
import pandas as pd
import shutil
import subprocess
import csv
import tempfile
import numpy as np
from pathlib import Path

from biopandas.pdb import PandasPdb

from Bio.PDB.DSSP import DSSP
from DSSPparser import parseDSSP
from Bio.SeqUtils import seq1

aamap = {
    'A': {'max_asa': 129.0, 'aa3cap': 'ALA'}, 
    'R': {'max_asa': 247.0, 'aa3cap': 'ARG'}, 
    'N': {'max_asa': 195.0, 'aa3cap': 'ASN'}, 
    'D': {'max_asa': 193.0, 'aa3cap': 'ASP'}, 
    'C': {'max_asa': 167.0, 'aa3cap': 'CYS'}, 
    'E': {'max_asa': 223.0, 'aa3cap': 'GLU'}, 
    'Q': {'max_asa': 225.0, 'aa3cap': 'GLN'}, 
    'G': {'max_asa': 104.0, 'aa3cap': 'GLY'}, 
    'H': {'max_asa': 224.0, 'aa3cap': 'HIS'}, 
    'I': {'max_asa': 197.0, 'aa3cap': 'ILE'}, 
    'L': {'max_asa': 201.0, 'aa3cap': 'LEU'}, 
    'K': {'max_asa': 236.0, 'aa3cap': 'LYS'}, 
    'M': {'max_asa': 224.0, 'aa3cap': 'MET'}, 
    'F': {'max_asa': 240.0, 'aa3cap': 'PHE'}, 
    'P': {'max_asa': 159.0, 'aa3cap': 'PRO'}, 
    'S': {'max_asa': 155.0, 'aa3cap': 'SER'}, 
    'T': {'max_asa': 172.0, 'aa3cap': 'THR'}, 
    'W': {'max_asa': 285.0, 'aa3cap': 'TRP'}, 
    'Y': {'max_asa': 263.0, 'aa3cap': 'TYR'}, 
    'V': {'max_asa': 174.0, 'aa3cap': 'VAL'}, 
}

dssp_dict = {'H':'H', 'G':'H', 'I':'H', 'P':'H',   # alpha-helix, 3-10 helix, pi-helix, polyproline helix
             'B':'B', 'E':'B', 'S':'C', 'T':'C', } # beta-bridge, strand, bend, turn/loop

# QUERY UNIPROT AND PARSE IT INTO A TSV FILE #

def query_uniprot(
    working_filedir, 
    input_uniprot, 
): 
    """
    Description
        A function to query UniProt for the protein sequence         
    """

    # QUERY FASTA FILE #
    ffile = input_uniprot + '.fasta'
    if not os.path.exists(os.path.join(working_filedir, 'sequence_structure', ffile)): 
        _ = wget.download(f'https://rest.uniprot.org/uniprotkb/{ffile}', 
                          out=str(working_filedir / 'sequence_structure'))
    else: print(f'sequence_structure/{ffile} exists')

    uFasta_file = os.path.join(working_filedir, f'sequence_structure/{ffile}')
    return uFasta_file

def parse_uniprot(
    uFasta_file, 
    out_fasta, 
): 
    """
    Description
        A function to process UniProt .fasta into 
        a list of positions and amino acids
    """
    
    # OPEN INPUT AND OUTPUT FILES #
    uFasta_list = open(out_fasta, "w")
    uFasta_list.write('unipos\tunires\n')

    uFasta = open(uFasta_file, "r")
    header = uFasta.readline() # skip header

    # READ FASTA SEQUENCE, AND WRITES POS, AMINO ACID #
    j = 0
    for fasta_line in uFasta:
        fasta_line = fasta_line.strip()
        for i in range(len(fasta_line)):
            uFasta_list.write("%d\t%s\n" % (j+1, fasta_line[i]))
            j += 1

    uFasta.close()
    uFasta_list.close()
    return None

# QUERY PDB FILE AND PARSE IT INTO A TSV FILE #

def query_af(
    working_filedir, 
    af_filename, 
    structureid, 
): 
    """
    Description
        A function to query AlphaFold for the protein structure
    """

    # QUERY ALPHAFOLD #
    affile = structureid + '.pdb'
    if not os.path.exists(working_filedir / af_filename): 
        _ = wget.download(f'https://alphafold.ebi.ac.uk/files/{affile}', out=str(working_filedir))
        os.rename(working_filedir / affile, working_filedir / af_filename)
    return None

def parse_af(
    working_filedir, 
    af_filename, 
    af_processed_filename, 
): 
    """
    Description
        Process AlphaFold structure for all atoms and their information
    """

    # PREPROCESS AF TO KEEP ATOMS #
    af_file = open(working_filedir / af_filename, "r")
    af_lines = af_file.readlines()
    af_file.close()

    af_processed_lines = [idx for idx in af_lines if idx[0:4] == "ATOM"]
    af_processed_lines = ["HEADER\n"] + ["CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"] + af_processed_lines ###
    af_processed_file = open(working_filedir / af_processed_filename, "w")
    af_processed_file.writelines(af_processed_lines)
    af_processed_file.close()

    return None

def parse_coord(
    working_filedir, 
    af_processed_filename, 
    fastalist_filename, 
    coord_filename, 
    target_chainid, 
):
    """
    Take in processed AlphaFold PDB and processed FASTA list, and parse
    X/Y/Z coordinates (CA), chain, and pLDDT (stored in B-factor for AF).
    Writes a tab-separated file with columns:
    unipos  unires  x_coord y_coord z_coord chain  bfactor_pLDDT
    """

    # --- Load inputs ---
    ppdb = PandasPdb().read_pdb(str(working_filedir / af_processed_filename))
    atom_df = ppdb.df["ATOM"]                        # only ATOM; HETATM is in ppdb.df["HETATM"]
    fasta_df = pd.read_csv(fastalist_filename, sep="\t")

    # expected fasta columns
    if not {"unipos", "unires"}.issubset(fasta_df.columns):
        raise ValueError("fasta_df must contain columns: 'unipos' and 'unires'")

    # ensure integer residue numbers for matching
    atom_df = atom_df.copy()
    atom_df["residue_number"] = atom_df["residue_number"].astype(int)
    fasta_df = fasta_df.copy()
    fasta_df["unipos"] = fasta_df["unipos"].astype(int)

    # Subset ATOM to the requested chains and CA atoms
    atom_ca = atom_df[(atom_df["atom_name"] == "CA") & (atom_df["chain_id"].isin([target_chainid]))]

    # Build quick lookup by residue_number → rows (permitted to be >1; we handle below)
    # (A merge would also work; we’ll keep the loop structure you used.)
    rows = []
    for i in range(len(fasta_df)):
        unipos = fasta_df.at[i, "unipos"]
        uniaa  = fasta_df.at[i, "unires"]

        # entries for that residue number across selected chains
        ca_entry = atom_ca.loc[atom_ca["residue_number"] == int(unipos)]

        # defaults if not found / mismatch
        x_coord = y_coord = z_coord = b_factor = "-"
        chain_id = target_chainid
        if len(ca_entry) == 1:
            row0 = ca_entry.iloc[0]
            aa_at_ca = row0["residue_name"]

            # You used aamap[str(uniaa)]['aa3cap'] — keep it if available; else trust AA
            try:
                expected_aa3 = aamap[str(uniaa)]["aa3cap"]
            except Exception:
                expected_aa3 = aa_at_ca  # fall back (or set to None and skip check)

            if aa_at_ca == expected_aa3:
                x_coord  = round(float(row0["x_coord"]), 3)
                y_coord  = round(float(row0["y_coord"]), 3)
                z_coord  = round(float(row0["z_coord"]), 3)
                chain_id = str(row0.get("chain_id", "-"))
                try:
                    b_factor = round(float(row0["b_factor"]), 3)
                except Exception:
                    b_factor = "-"
            else:
                warnings.warn(f"PDB and UNIPROT residue mismatch at pos {unipos}: {aa_at_ca} vs {expected_aa3}")

        elif len(ca_entry) > 1:
            warnings.warn(f"Multiple CA rows for residue {unipos} in chains {[target_chainid]}; skipping this residue.")

        # Collect one row per FASTA position
        rows.append({
            "unipos": unipos,
            "unires": uniaa,
            "x_coord": x_coord,
            "y_coord": y_coord,
            "z_coord": z_coord,
            "chain": target_chainid,
            "bfactor_pLDDT": b_factor,
        })

    df_loop = pd.DataFrame(
        rows,
        columns=["unipos","unires","x_coord","y_coord","z_coord","chain","bfactor_pLDDT"]
        )

    # --- Handle "others" (chains NOT in `chains`) at residue-level (CA only) ---
    other_ca = atom_df.loc[
        (atom_df["atom_name"] == "CA") & (~atom_df["chain_id"].isin([target_chainid])),
        ["residue_number","residue_name","x_coord","y_coord","z_coord","chain_id","b_factor"]
    ].copy()

    # Align to df_loop schema
    other_ca.rename(columns={
        "residue_number": "unipos",
        "residue_name":  "unires",
        "chain_id":      "chain",
        "b_factor":      "bfactor_pLDDT",
    }, inplace=True)

    # 3-letter -> 1-letter AA (unknowns -> 'X')
    def _to_aa1(x):
        try:
            return seq1(str(x)) if pd.notnull(x) else "X"
        except Exception:
            return "X"
    other_ca["unires"] = other_ca["unires"].apply(_to_aa1)

    # If multiple CA for same (chain,unipos) (e.g., altLoc), keep first
    other_ca.sort_values(["chain","unipos"], inplace=True)
    other_ca.drop_duplicates(subset=["chain","unipos"], keep="first", inplace=True)

    # Round numeric cols safely
    for c in ["x_coord","y_coord","z_coord","bfactor_pLDDT"]:
        other_ca[c] = pd.to_numeric(other_ca[c], errors="coerce").round(3)

    # Combine (both are residue-level now)
    final_df = pd.concat([df_loop, other_ca], ignore_index=True)

    # Ensure string-safe for TSV (optional)
    for c in ["unipos","unires","x_coord","y_coord","z_coord","chain","bfactor_pLDDT"]:
        final_df[c] = final_df[c].astype(str)

    # Write
    out_path = working_filedir / coord_filename
    final_df.to_csv(out_path, sep="\t", index=False)
    return None


# RUN DSSP AND PARSE IT INTO A TSV FILE #

def run_dssp(
    working_filedir, 
    af_filename, 
    dssp_filename, 
): 
    # os.environ["LIBCIFPP_DATA_DIR"] = "src/helpers/libcifpp_data"

    if not os.path.exists(working_filedir / dssp_filename): 
        pdb_file_path = str(working_filedir / af_filename)
        out_file_path = str(os.path.join(working_filedir / dssp_filename))
        # dic_file_path = "src/helpers/mmcif_pdbx_v50.dic"

        # RUN DSSP COMMAND #
        if shutil.which('dssp') is None: 
            subprocess.run(["mkdssp", pdb_file_path, out_file_path, 
                            "--output-format", "dssp"], 
                            check=True)
        else: 
            subprocess.run(["dssp", pdb_file_path, out_file_path, 
                            "--output-format", "dssp"], 
                            check=True)

def parse_dssp(
    working_filedir, 
    alphafold_dssp_filename, 
    fastalist_filename, 
    dssp_parsed_filename, 
    target_chainid, 
): 
    """
    Description
        A function to parse .dssp file for burial, phi, psi, etc
    """

    # PARSE DSSP #
    parser = parseDSSP(working_filedir / alphafold_dssp_filename)
    parser.parse()
    pddict = parser.dictTodataframe()
    pddict_ch = pddict.loc[pddict['chain'].isin([target_chainid])]
    pddict_ch = pddict_ch.fillna('-')
    pddict_ch = pddict_ch.replace(r'^\s*$', '-', regex=True)

    # READ FASTA AND DSSP, WRITE PROCESSED DSSP #
    fasta_df = pd.read_csv(fastalist_filename, sep = '\t')
    dssp_output_file = open(working_filedir / dssp_parsed_filename, 'w')
    dssp_output_file.write('\t'.join(['unipos', 'unires', 'SS9', 'SS3', 'ACC', 'RSA', 
                                      'exposure', 'PHI', 'normPHI', 'PSI', 'normPSI']) + '\n')
    
    output_data = []
    unipos_dict = fasta_df['unipos'].to_dict()
    uniaa_dict = fasta_df['unires'].to_dict()

    for i in range(len(fasta_df)): 
        unipos = unipos_dict[i]
        uniaa = uniaa_dict[i]
        pddict_ch_entry = pddict_ch.loc[pddict_ch['inscode'] == str(unipos), ] ###

        if len(pddict_ch_entry) == 0:
            dssp_SS9, dssp_ASA, dssp_Phi, dssp_Psi, dssp_SS3 = '-', '-', '-', '-', '-'
            norm_ASA, exposure, norm_Phi, norm_Psi = '-', '-', '-', '-'
        elif len(pddict_ch_entry) == 1:
            dssp_SS9 = pddict_ch_entry['struct'].iloc[0].strip() ###
            if dssp_SS9 == "-":               dssp_SS9 = "L"
            if dssp_SS9 in dssp_dict.keys():  dssp_SS3 = dssp_dict[dssp_SS9]
            else:                             dssp_SS3 = "C"

            dssp_ASA = pddict_ch_entry['acc'].iloc[0]
            Gly_X_Gly = aamap[str(uniaa)]['max_asa']
            norm_ASA = round(float(dssp_ASA) / float(Gly_X_Gly), 2)
            if           norm_ASA < 0.05:  exposure = "core"
            elif 0.05 <= norm_ASA < 0.25:  exposure = "buried"
            elif 0.25 <= norm_ASA < 0.50:  exposure = "medburied"
            elif 0.50 <= norm_ASA < 0.75:  exposure = "medexposed"
            else:                          exposure = "exposed"

            dssp_Phi = pddict_ch_entry['phi'].iloc[0]
            norm_Phi = round(float(dssp_Phi) / 180.0, 2)
            dssp_Psi = pddict_ch_entry['psi'].iloc[0]
            norm_Psi = round(float(dssp_Psi) / 180.0, 2)
        else:
            warnings.warn(pddict_ch_entry)

        out = '\t'.join([str(unipos), str(uniaa), str(dssp_SS9), str(dssp_SS3), str(dssp_ASA), str(norm_ASA), 
                         str(exposure), str(dssp_Phi), str(norm_Phi), str(dssp_Psi), str(norm_Psi)])
        output_data.append(out)

    output_data_all = "\n".join(output_data)
    dssp_output_file.writelines(output_data_all)
    dssp_output_file.close()
    del fasta_df, pddict, pddict_ch
    return None

# QUERY DOMAINS AND PARSE IT INTO A TSV FILE #

def query_domains(
    working_filedir, 
    uniprot_id, 
    output_file
):
    """
    Fetches domain annotations from UniProt given a Uniprot ID and saves as a TSV file
    """

    # FETCH DATA #
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Error fetching data: {response.status_code}")
        return

    # SEQUENCE #
    data = response.json()
    sequence = data.get("sequence", {}).get("value", "")
    if not sequence:
        print("No sequence found.")
        return
    
    # DOMAINS #
    features = data.get("features", [])
    domain_map = {}
    for feature in features:
        if feature.get("type") in {"Domain", "Repeat"}:
            domain_name = feature.get("description", "Unknown domain")
            begin = int(feature["location"]["start"]["value"])
            end = int(feature["location"]["end"]["value"])
            for pos in range(begin, end + 1):
                domain_map[pos] = domain_name

    # WRITE #
    with open(working_filedir / output_file, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["Position", "Residue", "Domain"])
        for i, residue in enumerate(sequence, start=1):
            domain = domain_map.get(i, "None")
            writer.writerow([i, residue, domain])
    return

def parse_domains(
    working_filedir, 
    out_fasta, 
    domains_filename, 
    domains_dict
): 
    df_sequence = pd.read_csv(out_fasta, sep='\t')
    df_sequence = df_sequence.rename(columns={'unipos':'Position', 'unires':'Residue'})

    df_sequence['Domain'] = df_sequence['Position'].apply(lambda x: get_domain(x, domains_dict))
    df_sequence.to_csv(working_filedir / domains_filename, sep='\t',index=False)
    return None

def get_domain(
    pos, 
    domains, 
):
    for name, (start, end) in domains.items():
        if start <= pos <= end:
            return name
    return 'None'

# OTHER PREPROCESS #

def count_aa_within_radius(
    working_filedir, 
    coord_filename, 
    coord_radius_filename, 
    target_chain,
    radius=6.0, 
): 
    """
    Description
        Count the number of residues within [radius] Angstroms of the focal residue
    """

    # COUNT AMINO ACIDS IN 6A DISTANCE AND TEIR IDENTITY #
    df_coord = pd.read_csv(working_filedir / coord_filename, sep = "\t")    
    # PRE EXTRACT VALUES TO AVOID DF CALLS #
    x_coords = df_coord["x_coord"].values
    y_coords = df_coord["y_coord"].values
    z_coords = df_coord["z_coord"].values
    unires_all = df_coord["unires"].values
    unipos_all = df_coord["unipos"].values
    chains_all = df_coord["chain"].values

    incomplete_structure = False
    taa_count, taa_naa, taa_naa_positions, taa_naa_chains = [], [], [], []

    for taa in range(len(df_coord)):
        t_chain = chains_all[taa]
        t_xcoord = x_coords[taa]
        t_ycoord = y_coords[taa]
        t_zcoord = z_coords[taa]
        if t_chain == target_chain:
            # print(t_chain,taa)
            dis_count, naas, naas_positions, chains = 0, [], [], []
            # IF STRUCTURE IS INCOMPLETE #
            if t_xcoord == '-' or t_ycoord == '-' or t_zcoord == '-': 
                taa_count.append(dis_count)
                taa_naa.append(';'.join(naas))
                taa_naa_positions.append(';'.join(naas_positions))
                taa_naa_chains.append(';'.join(chains))
                incomplete_structure = True
                continue

            t_xcoord = float(t_xcoord)
            t_ycoord = float(t_ycoord)
            t_zcoord = float(t_zcoord)

            for naa in range(len(df_coord)):
                if taa == naa:
                    continue
                xcoord = x_coords[naa]
                ycoord = y_coords[naa]
                zcoord = z_coords[naa]

                if xcoord == '-' or ycoord == '-' or zcoord == '-':
                    continue
                dx = float(xcoord) - t_xcoord
                dy = float(ycoord) - t_ycoord
                dz = float(zcoord) - t_zcoord
                pairwise_dist = math.sqrt(dx**2 + dy**2 + dz**2)

                # ADD TO LIST IF WITHIN RADIUS CUTOFF #
                if pairwise_dist <= radius:
                    dis_count += 1
                    naas.append(unires_all[naa])
                    naas_positions.append(str(unipos_all[naa]))
                    chains.append(chains_all[naa])
            
            taa_count.append(dis_count)
            taa_naa.append(';'.join(naas))
            taa_naa_positions.append(';'.join(naas_positions))
            taa_naa_chains.append(';'.join(chains))
        
    if incomplete_structure: warnings.warn(f"Incomplete Structure Detected")
    df_coord = df_coord[df_coord['chain']==target_chain]
    df_coord['Naa_count'] = taa_count
    df_coord['Naa'] = taa_naa
    df_coord['Naa_pos'] = taa_naa_positions
    df_coord['Naa_chain'] = taa_naa_chains

    df_coord.to_csv(working_filedir / coord_radius_filename, sep='\t')
    return df_coord

def degree_of_burial(
    df_dssp, 
    df_coord, 
    working_filedir, 
    coord_dssp_filename,
    target_chainid
): 
    """
    Description
        Calculate the degree of burial per residue with maxRSA metric
    """
    df_coord_dssp = pd.merge(df_coord, df_dssp, on=["unipos", "unires"])
    for colname in ['unipos', 'unires', 'SS9', 'SS3', 'ACC', 'RSA', 'exposure', 'PHI', 'normPHI', 'PSI', 'normPSI']: 
        assert colname in df_coord_dssp

    ### 'dBurial', 'normSumdBurial', 'pLDDT_dis'

    df_dssp_rsa = df_dssp["RSA"].tolist()
    maxRSA = float(max([x for x in df_dssp_rsa if x != '-']))

    # Create a numeric version of RSA for calculation
    rsa_num = pd.to_numeric(df_coord_dssp['RSA'], errors='coerce')

    # Compute dBurial only for numeric RSA values
    df_coord_dssp['dBurial'] = np.where(
        rsa_num.notna(),
        (maxRSA - rsa_num).round(3),
       0.0   # keep '0.0' where RSA was '-'
    )

    # CALCULATE DEGREE OF BURIAL PER RESIDUE normSumdBurial AND CATEGORY pLDDT_dis #
    aa_wise_cdBurial = []
    arr_pLDDT_discrete = []
    naa_list_dict = df_coord_dssp['Naa'].to_dict()
    naa_pos_list_dict = df_coord_dssp['Naa_pos'].to_dict()
    naa_chain_list_dict = df_coord_dssp['Naa_chain'].to_dict()
    taa_dBurial_dict = df_coord_dssp['dBurial'].to_dict()
    pLDDT_dict = df_coord_dssp['bfactor_pLDDT'].to_dict()

    for i in range(len(df_coord_dssp)): 
        taa_dBurial  = taa_dBurial_dict[i]
        if taa_dBurial == '-': taa_dBurial = 0.0
        naa_list     = naa_list_dict[i].split(';') # NEIGHBORING AAs
        naa_pos_list = naa_pos_list_dict[i].split(';')
        naa_chain_list = naa_chain_list_dict[i].split(';')

        # CALCULATE #
        sum_dBurial = 0
        for naa_chain, naa_pos in zip(naa_chain_list, naa_pos_list): 
            if naa_pos != '' and naa_chain in [target_chainid]: 
                dburial = taa_dBurial_dict[int(naa_pos)-1]
                if dburial != '-': 
                    sum_dBurial += round(dburial, 2)
                 ### 250414 only errored on men1 pdb ###
        norm_sum_dBurial = round(sum_dBurial / len(naa_list), 2)
        aa_wise_cdBurial.append(round(norm_sum_dBurial * taa_dBurial, 3))

        # CATEGORIZE #
        pLDDT = pLDDT_dict[i]
        if pLDDT == '-': pLDDT_discrete = '-'
        else: 
            pLDDT = float(pLDDT)
            if         pLDDT < 50:  pLDDT_discrete = 'very low'
            elif 50 <= pLDDT < 70:  pLDDT_discrete = 'low'
            elif 70 <= pLDDT < 90:  pLDDT_discrete = 'confident'
            else:                   pLDDT_discrete = 'high'
        arr_pLDDT_discrete.append(pLDDT_discrete)

    df_coord_dssp['normSumdBurial'] = aa_wise_cdBurial
    df_coord_dssp['pLDDT_dis'] = arr_pLDDT_discrete
    df_coord_dssp['dBurial'] = df_coord_dssp['Naa_count'].map(lambda x: '-' if x == 0 else x)    

    df_coord_dssp.to_csv(working_filedir / coord_dssp_filename, sep="\t", index=False)
    return df_coord_dssp

def infer_element_symbol(
    atom_name, 
):
    """
    Infers the chemical element symbol from the atom name.
    """
    atom_name = atom_name.strip()
    if not atom_name:
        return "  "
    # If the first character is a digit (e.g., '1H'), the element is the second character
    if atom_name[0].isdigit():
        return atom_name[1].upper().rjust(2)
    # If the first two characters are letters and the second is lowercase, it's a two-letter element
    if len(atom_name) >= 2 and atom_name[1].islower():
        return atom_name[:2].capitalize().rjust(2)
    # Otherwise, it's a one-letter element
    return atom_name[0].upper().rjust(2)

def update_pdb_element_symbols(
    input_pdb_path, 
    output_pdb_path=None, 
):
    if output_pdb_path is None or input_pdb_path == output_pdb_path:
        # In-place edit: use a temp file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_outfile:
            tmp_path = tmp_outfile.name
            with open(input_pdb_path, 'r') as infile:
                for line in infile:
                    if line.startswith(('ATOM  ', 'HETATM')):
                        atom_name = line[12:16]
                        element_symbol = infer_element_symbol(atom_name)
                        line = line.rstrip('\n')
                        if len(line) < 78:
                            line = line.ljust(78)
                        updated_line = line[:76] + element_symbol + line[78:] + '\n'
                        tmp_outfile.write(updated_line)
                    else:
                        tmp_outfile.write(line)
        shutil.move(tmp_path, input_pdb_path)  # Replace original file
    else:
        # Normal input → output
        with open(input_pdb_path, 'r') as infile, open(output_pdb_path, 'w') as outfile:
            for line in infile:
                if line.startswith(('ATOM  ', 'HETATM')):
                    atom_name = line[12:16]
                    element_symbol = infer_element_symbol(atom_name)
                    line = line.rstrip('\n')
                    if len(line) < 78:
                        line = line.ljust(78)
                    updated_line = line[:76] + element_symbol + line[78:] + '\n'
                    outfile.write(updated_line)
                else:
                    outfile.write(line)
                    
def count_residue_contacts_all_atoms_single(
    pdb_path,
    coord_filename,
    coord_radius_filename,
    target_chain=None,
    radius=6.0,
    include_hetero=False,
):
    """
    Read a PDB and, for each residue, count neighboring residues whose **minimum atom–atom
    distance** to the focal residue is <= radius (Å). Outputs a TSV and returns a pandas
    DataFrame with:
        ['unires','unipos','chain','x_coord','y_coord','z_coord',
         'Naa_count','Naa','Naa_pos','Naa_chain']

    Notes:
    - `unipos` is an **integer** (resseq). Insertion codes are ignored for `unipos`.
    - Cα coordinates are reported in x/y/z columns; if a residue lacks a Cα, values are NaN.
    - Neighbors can be from any chain; output rows can be restricted via `target_chain`.
    """

    # ---------- helpers (scoped locally) ----------
    def _aa3_to_aa1():
        return {
            'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
            'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
            'THR':'T','TRP':'W','TYR':'Y','VAL':'V','SEC':'U','PYL':'O'
        }

    def _one_letter(resname):
        return _aa3_to_aa1().get(resname.upper(), 'X')

    def _parse_pdb_minimal(pdb_path_, include_hetero_):
        """
        Minimal PDB parser (ATOM/HETATM) producing a list of residue dicts:
        {chain,resname,resseq(int),icode,atoms=[(x,y,z),...], ca=(x,y,z)|None}
        Skips altlocs other than '' or 'A'.
        """
        residues = []
        key_to_idx = {}
        with open(pdb_path_, 'r') as fh:
            for line in fh:
                rec = line[:6]
                if rec not in ('ATOM  ', 'HETATM'):
                    continue
                if rec == 'HETATM' and not include_hetero_:
                    continue
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    continue
                altloc = line[16].strip()
                if altloc not in ('', 'A'):
                    continue
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip() or '_'
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    continue
                icode = line[26].strip()
                key = (chain, resname, resseq, icode)
                if key not in key_to_idx:
                    key_to_idx[key] = len(residues)
                    residues.append({
                        'chain': chain, 'resname': resname,
                        'resseq': resseq, 'icode': icode,
                        'atoms': [], 'ca': None
                    })
                idx = key_to_idx[key]
                residues[idx]['atoms'].append((x, y, z))
                if atom_name == 'CA' and residues[idx]['ca'] is None:
                    residues[idx]['ca'] = (x, y, z)
        return residues

    def _safe_altloc(atom):
        try:
            altloc = atom.get_altloc()
        except Exception:
            altloc = getattr(atom, 'altloc', ' ')
        return altloc in (' ', '', 'A')

    # ---------- main logic ----------
    pdb_path = Path(pdb_path)
    coord_radius_filename = Path(coord_radius_filename)
    coord_radius_filename.parent.mkdir(parents=True, exist_ok=True)

    # Try Biopython fast path
    use_biopy = True
    try:
        from Bio.PDB import PDBParser, NeighborSearch
        from Bio.PDB.Polypeptide import is_aa
    except Exception:
        use_biopy = False

    incomplete_structure = False

    if use_biopy:
        # ------ Biopython branch ------
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("x", str(pdb_path))

        # collect residues and atoms (+ Cα)
        res_records = []  # (cid, res, resname, resseq_int, atoms_list, ca_xyz|None)
        for model in structure:
            for chain in model:
                cid = chain.id
                for res in chain:
                    hetflag, resseq, icode = res.id  # resseq int
                    if not include_hetero and hetflag != ' ':
                        continue
                    if not is_aa(res, standard=False) and not include_hetero:
                        continue
                    atoms = [a for a in res.get_atoms() if _safe_altloc(a)]
                    if not atoms:
                        incomplete_structure = True
                        # continue
                    ca_atom = next((a for a in atoms if a.get_id() == 'CA'), None)
                    ca_xyz = tuple(map(float, ca_atom.coord)) if ca_atom is not None else None
                    if ca_xyz is None:
                        incomplete_structure = True
                    res_records.append((
                        cid, res, res.get_resname().strip(), int(resseq), atoms, ca_xyz
                    ))

        if not res_records:
            raise ValueError("No residues with atoms found in the PDB after filters.")

        # build neighbor search over all atoms
        all_atoms = []
        for _, _, _, _, atoms, _ in res_records:
            all_atoms.extend(atoms)
        ns = NeighborSearch(all_atoms)

        rows = []
        for cid, res, resname, pos_i, atoms, ca_xyz in res_records:
            neighbor_residues = set()
            if atoms:
                for a in atoms:
                    for na in ns.search(a.coord, radius, level='A'):
                        r = na.get_parent()
                        if r is None or r is atoms[0].get_parent():
                            continue
                        neighbor_residues.add(r)

            neighbors = []
            for nr in neighbor_residues:
                n_cid = nr.get_parent().id
                n_resname = nr.get_resname().strip()
                n_resseq = int(nr.id[1])
                neighbors.append((n_cid, _one_letter(n_resname), n_resseq))
            neighbors.sort(key=lambda t: (t[0], t[2]))

            x, y, z = (np.nan, np.nan, np.nan) if ca_xyz is None else ca_xyz
            rows.append({
                'unires': _one_letter(resname),
                'unipos': int(pos_i),
                'chain': cid,
                'x_coord': float(x),
                'y_coord': float(y),
                'z_coord': float(z),
                'Naa_count': len(neighbors),
                'Naa': ';'.join([n[1] for n in neighbors]) if neighbors else '-',
                'Naa_pos': ';'.join(str(n[2]) for n in neighbors) if neighbors else '-',
                'Naa_chain': ';'.join([n[0] for n in neighbors]) if neighbors else '-',
            })

        df = pd.DataFrame(rows)

    else:
        # ------ Minimal parser branch ------
        residues = _parse_pdb_minimal(str(pdb_path), include_hetero_=include_hetero)
        if not residues:
            raise ValueError("No residues with atoms found in the PDB after filters.")

        # gather all atoms for brute-force neighbor test
        atom_coords = []
        atom_res_idx = []
        for i, r in enumerate(residues):
            if not r['atoms']:
                incomplete_structure = True
                # continue
            for (x, y, z) in r['atoms']:
                atom_coords.append((x, y, z))
                atom_res_idx.append(i)

        rows = []
        r2aa = [_one_letter(r['resname']) for r in residues]
        r2pos = [int(r['resseq']) for r in residues]  # integer positions
        r2chain = [r['chain'] for r in residues]

        r2 = radius * radius

        for i, r in enumerate(residues):
            if not r['atoms']:
                ca_xyz = None
                neigh = []
            else:
                ca_xyz = r['ca']  # (x,y,z) or None
                neighbor_idx = set()
                for ax, ay, az in r['atoms']:
                    for j, (bx, by, bz) in enumerate(atom_coords):
                        if atom_res_idx[j] == i:
                            continue
                        dx = bx - ax; dy = by - ay; dz = bz - az
                        if (dx*dx + dy*dy + dz*dz) <= r2:
                            neighbor_idx.add(atom_res_idx[j])
                neigh = [(r2chain[k], r2aa[k], r2pos[k]) for k in neighbor_idx]
                neigh.sort(key=lambda t: (t[0], t[2]))

            x, y, z = (np.nan, np.nan, np.nan) if ca_xyz is None else ca_xyz
            rows.append({
                'unires': r2aa[i],
                'unipos': r2pos[i],
                'chain': r2chain[i],
                'x_coord': float(x),
                'y_coord': float(y),
                'z_coord': float(z),
                'Naa_count': len(neigh),
                'Naa': ';'.join([n[1] for n in neigh]),
                'Naa_pos': ';'.join(str(n[2]) for n in neigh),
                'Naa_chain': ';'.join([n[0] for n in neigh]),
            })

        df = pd.DataFrame(rows)

    if incomplete_structure:
        warnings.warn("Incomplete Structure Detected (some residues missing atoms and/or Cα)")

    if target_chain is not None:
        df = df[df['chain'] == target_chain].copy()

    df['unique_id'] = df['chain'] + '_' + df['unipos'].astype(str)
    coord_pd = pd.read_csv(coord_filename,sep='\t')
    coord_pd['unique_id'] = coord_pd['chain'] + "_" + coord_pd['unipos'].astype(str)
    id_list = coord_pd['unique_id'].to_list()

    result_pd = pd.merge(coord_pd, df[['unique_id','Naa_count','Naa','Naa_pos','Naa_chain']], how='outer', on='unique_id', suffixes=('_A','_B'))
    result_pd['unipos'] = result_pd['unipos'].astype(int)
    result_pd['unique_id'] = pd.Categorical(result_pd['unique_id'], categories=id_list, ordered=True)
    result_pd = result_pd.sort_values('unique_id')
    result_pd = result_pd.drop('unique_id',axis=1)
    result_pd = result_pd[result_pd['chain']==target_chain]
    result_pd = result_pd.reset_index(drop=True)
    result_pd['Naa_count'] = result_pd['Naa'].map(
        lambda x: 0 if (pd.isna(x) or x in ['-', '']) else len(str(x).split(';'))
    )

    # Write TSV
    result_pd = result_pd.fillna('-')
    result_pd.to_csv(coord_radius_filename, sep='\t', index=False)

    return result_pd
