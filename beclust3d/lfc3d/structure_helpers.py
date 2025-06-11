"""
File: af_structural_features.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
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

from biopandas.pdb import PandasPdb

from Bio.PDB.DSSP import DSSP
from DSSPparser import parseDSSP

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
        working_filedir, input_uniprot
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
    uFasta_file, out_fasta
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
    working_filedir, af_filename, structureid
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
    af_filename, af_processed_filename, 
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
    fastalist_filename, coord_filename, 
    chains, 
): 
    """
    Description
        Take in processed AlphaFold and processed fasta and parse
        coordinates and values
    """

    # GET COORDS AND CONFID VALS FROM PROCESSED AF #
    ### chain?
    alphafold_pdb = PandasPdb()
    alphafold_pdb.read_pdb(str(working_filedir / af_processed_filename))
    atom_df = alphafold_pdb.df['ATOM']
    fasta_df = pd.read_csv(fastalist_filename, sep = '\t')
    coord_file = open(working_filedir / coord_filename, 'w')
    coord_file.write('\t'.join(['unipos', 'unires', 
                                'x_coord', 'y_coord', 'z_coord', 
                                'chain', 'bfactor_pLDDT']) + '\n')

    # PARSE OUT X Y Z B DATA FROM PROCESSED FASTA, PROCESSED AF #
    output_data = []
    unipos_dict = fasta_df['unipos'].to_dict()
    uniaa_dict = fasta_df['unires'].to_dict()
    
    for i in range(len(fasta_df)):
        unipos = unipos_dict[i]
        uniaa = uniaa_dict[i]
        entry = atom_df.loc[atom_df['residue_number'] == int(unipos), ] ###
        ca_entry = entry.loc[(entry['atom_name'] == "CA") & (entry['chain_id'].isin(chains)), ] ###

        x_coord, y_coord, z_coord, chain_id, b_factor = "-", "-", "-", "-", "-"
        
        if len(ca_entry) == 1: 
            row = ca_entry.iloc[0]
            aa_at_ca = row['residue_name']
            uni_res = aamap[str(uniaa)]['aa3cap']

            if aa_at_ca == uni_res: 
                x_coord = round(float(row['x_coord']), 3)
                y_coord = round(float(row['y_coord']), 3)
                z_coord = round(float(row['z_coord']), 3)

                try: chain_id = row['chain_id']
                except Exception: pass
                try: b_factor = round(float(row['b_factor']), 3)
                except Exception: pass
            else: 
                warnings.warn(f"PDB and UNIPROT residue mismatch {aa_at_ca}: {uni_res}")
        elif len(ca_entry) > 1: 
            warnings.warn("PROBLEM - CHECK PDB")

        output_data_entry = '\t'.join([str(unipos), str(uniaa), 
                                       str(x_coord), str(y_coord), str(z_coord), 
                                       str(chain_id), str(b_factor)])+'\n'
        output_data.append(output_data_entry)

    del alphafold_pdb, atom_df, fasta_df
    coord_file.writelines(output_data)
    coord_file.close()
    return None

# RUN DSSP AND PARSE IT INTO A TSV FILE #

def run_dssp(
    working_filedir, af_filename, dssp_filename, 
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
        alphafold_dssp_filename, fastalist_filename, 
        dssp_parsed_filename, chains, 
): 
    """
    Description
        A function to parse .dssp file for burial, phi, psi, etc
    """

    # PARSE DSSP #
    parser = parseDSSP(working_filedir / alphafold_dssp_filename)
    parser.parse()
    pddict = parser.dictTodataframe()
    pddict_ch = pddict.loc[pddict['chain'].isin(chains)]
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

def query_domains(working_filedir, uniprot_id, output_file):
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

def parse_domains(working_filedir, out_fasta, domains_filename, domains_dict): 
    df_sequence = pd.read_csv(out_fasta, sep='\t')
    df_sequence = df_sequence.rename(columns={'unipos':'Position', 'unires':'Residue'})

    df_sequence['Domain'] = df_sequence['Position'].apply(lambda x: get_domain(x, domains_dict))
    df_sequence.to_csv(working_filedir / domains_filename, sep='\t',index=False)
    return None

def get_domain(pos, domains):
    for name, (start, end) in domains.items():
        if start <= pos <= end:
            return name
    return 'None'

# OTHER PREPROCESS #

def count_aa_within_radius(
        working_filedir, coord_filename, coord_radius_filename, 
        radius=6.0, 
): 
    """
    Description
        Count the number of residues within [radius] Angstroms
        of the focal residue
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
        t_xcoord = x_coords[taa]
        t_ycoord = y_coords[taa]
        t_zcoord = z_coords[taa]

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

    df_coord['Naa_count'] = taa_count
    df_coord['Naa'] = taa_naa
    df_coord['Naa_pos'] = taa_naa_positions
    df_coord['Naa_chain'] = taa_naa_chains

    df_coord.to_csv(working_filedir / coord_radius_filename, sep='\t')
    return df_coord

def degree_of_burial(
        df_dssp, df_coord, 
        working_filedir, coord_dssp_filename, 
): 
    """
    Description
        Calculate the degree of burial per residue with maxRSA metric
    """
    # df_coord_dssp = df_coord.copy()
    # dssp_columns = ['SS9', 'SS3', 'ACC', 'RSA', 'exposure', 'PHI', 'normPHI', 'PSI', 'normPSI']
    # df_coord_dssp[dssp_columns] = df_dssp[dssp_columns]
    df_coord_dssp = pd.merge(df_coord, df_dssp, on=["unipos", "unires"])
    for colname in ['unipos', 'unires', 'SS9', 'SS3', 'ACC', 'RSA', 'exposure', 'PHI', 'normPHI', 'PSI', 'normPSI']: 
        assert colname in df_coord_dssp

    ### 'dBurial', 'normSumdBurial', 'pLDDT_dis'

    df_dssp_rsa = df_dssp["RSA"].tolist()
    maxRSA = float(max([x for x in df_dssp_rsa if x != '-']))
    df_coord_dssp['dBurial'] = [round(maxRSA-float(x), 3) if x != '-' else x for x in df_dssp_rsa] 

    # CALCULATE DEGREE OF BURIAL PER RESIDUE normSumdBurial AND CATEGORY pLDDT_dis #
    aa_wise_cdBurial = []
    arr_pLDDT_discrete = []
    naa_list_dict = df_coord_dssp['Naa'].to_dict()
    naa_pos_list_dict = df_coord_dssp['Naa_pos'].to_dict()
    taa_dBurial_dict = df_coord_dssp['dBurial'].to_dict()
    pLDDT_dict = df_coord_dssp['bfactor_pLDDT'].to_dict()

    for i in range(len(df_coord_dssp)): 
        taa_dBurial  = taa_dBurial_dict[i]
        if taa_dBurial == '-': taa_dBurial = 0.0
        naa_list     = naa_list_dict[i].split(';') # NEIGHBORING AAs
        naa_pos_list = naa_pos_list_dict[i].split(';')

        # CALCULATE #
        sum_dBurial = 0
        for naa_pos in naa_pos_list: 
            if naa_pos != '': 
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

    df_coord_dssp.to_csv(working_filedir / coord_dssp_filename, sep="\t", index=False)
    return df_coord_dssp

def infer_element_symbol(atom_name):
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

def update_pdb_element_symbols(input_pdb_path, output_pdb_path=None):
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