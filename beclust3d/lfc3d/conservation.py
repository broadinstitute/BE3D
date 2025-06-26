"""
File: conservation.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 2

"""

import pandas as pd
from pathlib import Path
import requests
import sys, os
import time
import warnings
import subprocess

from Bio import AlignIO

cons_dict = {
    '*': ('conserved', 3),
    ':': ('similar', 2),
    '.': ('weakly_similar', 1),
    ' ': ('not_conserved', -1),
}

def conservation(
        workdir, 
        input_gene, alt_input_gene, 
        input_uniprot, alt_input_uniprot, 
        alignment_filename=None, mode='run', 
        title=None, email=None, wait_time=30, 
): 
    """
    Generate dataframes of sequence conservation for each residue. 
    The default is to query the alignment from protein sequences, but the user can also upload an alignment file.
    The alternate gene and uniprot is for if the data corresponds to a different sequence, 
    whether it's a different gene or isoform of species, to be aligned to the sequence and structure. 

    Parameters
    ----------
    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    alt_input_gene : str
        Name of the alternate gene (e.g., mouse ortholog or isoform) to align against.

    input_uniprot : str
        Uniprot of the gene being processed. 

    alt_input_uniprot : str
        UniProt for the alternate gene.

    alignment_filename : str or Path or None, optional (default=None)
        If provided, path to a precomputed sequence alignment file (.afa or other format).
        If None, the function generates the alignment automatically.

    mode : {'run', 'query'}, optional (default='run')
        - 'run': Perform local multiple sequence alignment using MUSCLE.
        - 'query': Submit a remote MUSCLE API job (requires `email` and `title`).

    title : str or None, optional
        Title for the MUSCLE remote query (required if `mode='query'`).

    email : str or None, optional
        Email address required for using the remote MUSCLE API (required if `mode='query'`).

    wait_time : int, optional (default=30)
        Time in seconds to wait before re-polling the MUSCLE API when submitting a remote alignment request.

    Returns
    -------
    df_alignconserv : pd.DataFrame
        DataFrame containing per-residue conservation scores based on the sequence alignment.

    df_residuemap : pd.DataFrame
        DataFrame mapping residue indices between the primary and alternate sequences 
    """
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'conservation'):
        os.mkdir(working_filedir / 'conservation')

    assert mode=='run' or mode=='query'
    if mode=='query': 
        assert title is not None and email is not None

    # QUERY UNIPROT AND WRITE TO sequences.fasta #
    if alignment_filename is None: 
        request_filename_original, request_filename_alternate = f"{input_uniprot}.fasta", f"{alt_input_uniprot}.fasta"
        original_seq = query_protein_fasta(working_filedir, request_filename_original)
        alternate_seq = query_protein_fasta(working_filedir, request_filename_alternate)
        seqs_filename = f"conservation/sequences.fasta"
        with open(working_filedir / seqs_filename, "w") as text_file:
            text_file.write(original_seq)
            text_file.write(alternate_seq)

        # MUSCLE ALIGNMENT #
        muscle_output_filename = f"conservation/{input_gene}_{alt_input_gene}.afa"
        align_filename = f"conservation/{input_gene}_{alt_input_gene}.align"
        # RUN ALIGNMENT LOCALLY #
        if mode=='run': run_muscle(working_filedir, seqs_filename, muscle_output_filename, align_filename)
        # IF MUSCLE CANT BE RUN LOCALLY, QUERY API #
        if mode=='query': query_muscle(working_filedir, seqs_filename, align_filename, email, title, wait_time)
    else: 
        align_filename = alignment_filename
    
    # PARSE ALIGNMENT #
    alignconserv_filename = f"conservation/{input_gene}{input_uniprot}_{alt_input_gene}{alt_input_uniprot}_align_conservation.tsv"
    residuemap_filename =  f"conservation/{input_gene}{input_uniprot}_{alt_input_gene}{alt_input_uniprot}_residuemap_conservation.tsv"
    df_alignconserv, df_residuemap = parse_alignment(working_filedir, align_filename, alignconserv_filename, residuemap_filename)

    return df_alignconserv, df_residuemap

def query_protein_fasta(
    edits_filedir, request_filename
): 
    """
    Description
        Query a Uniprot ID for .fasta file
    """

    # QUERY UNIPROT #
    url = "https://rest.uniprot.org/uniprotkb/"
    requestURL = url+request_filename
    response = requests.get(requestURL)
    if not response.ok:
        response.raise_for_status()
        sys.exit()

    # RESPONSE TEXT #
    response_body = response.text
    with open(edits_filedir / 'conservation' / request_filename, "w") as text_file:
        text_file.write(response_body)
    return response_body

def query_muscle(
    edits_filedir, 
    seqs_filename, align_filename, 
    email, title, wait_time, 
): 
    """
    Description
        Query ebi.ac.uk for MUSCLE alignment
    """
    # POST MUSCLE #
    url = 'https://www.ebi.ac.uk/Tools/services/rest/muscle/run'
    files = {'sequence': open(edits_filedir / seqs_filename, 'rb')}
    data = {
        'email': email, 'title': title, 'format': 'clw', 'tree': 'none',
    }
    response = requests.post(url, data=data, files=files)
    if response.status_code != 200: 
        warnings.warn('Error with MUSCLE query')
    time.sleep(wait_time)

    # GET MUSCLE #
    job_url_code = response.text.split(' ')[-1]
    print(f'Job ID: {job_url_code}')
    url = f'https://www.ebi.ac.uk/Tools/services/rest/muscle/result/{job_url_code}/aln-clustalw'
    response = requests.get(url)

    with open(edits_filedir / align_filename, "wb") as f:
        f.write(response.content)
    return job_url_code

def run_muscle(
    edits_filedir, 
    seqs_filename, afa_filename, align_filename, 
): 
    """
    Description
        Run local MUSCLE for sequence alignment
    """
    # Run MUSCLE alignment
    subprocess.run([
        "muscle", 
        "-align", str(edits_filedir / seqs_filename), 
        "-output", str(edits_filedir / afa_filename), 
        "-threads", "1", 
        ], check=True)

    # Run Clustal Omega for formatting
    subprocess.run([
        "clustalo",
        "--in", edits_filedir / afa_filename,
        "--infmt=sta",
        "--outfmt=clustal",
        "-o", edits_filedir / align_filename,
        "--force"
    ], check=True)

    # # "chmod +x src/helpers/align/muscle-osx-arm64.v5.3"
    # muscle_exe = "src/helpers/align/muscle-osx-arm64.v5.3"
    # in_file = str(edits_filedir / seqs_filename)
    # inter_file = str(edits_filedir / afa_filename)
    # out_file = str(edits_filedir / align_filename)

    # subprocess.run([
    #     muscle_exe,
    #     "-align", in_file,
    #     "-output", inter_file
    # ], check=True)
    
def parse_alignment(
    edits_filedir, 
    align_filename, alignconserv_filename, residuemap_filename, 
): 

    align = AlignIO.read(edits_filedir / align_filename, "clustal")
    human_align_res = align[0].seq
    mouse_align_res = align[1].seq
    score = align.column_annotations['clustal_consensus']

    index  = [i+1 for i in range(len(human_align_res))]
    dis, v = [cons_dict[s][0] for s in score], [cons_dict[s][1] for s in score]
    colnames = ['alignment_pos', 'human_aligned_res', 'mouse_aligned_res', 'score', 'dis_score', 'v_score']
    colvals  = [index, 
                [c for c in human_align_res], 
                [c for c in mouse_align_res], 
                [s for s in score], 
                dis, v]

    df_alignconserv = pd.DataFrame()
    for name, col in zip(colnames, colvals): 
        df_alignconserv[name] = col
    df_alignconserv.to_csv(edits_filedir / alignconserv_filename, sep='\t', index=False)

    # COUNT CONSERVATION POSITIONS, ADD TO DATAFRAME #
    i, j = 0, 0
    human_res_pos, mouse_res_pos = [], []
    for s in human_align_res: 
        if s != '-': 
            i += 1
        human_res_pos.append(i)
    for s in mouse_align_res: 
        if s != '-': 
            j += 1
        mouse_res_pos.append(j)

    colnames = ['alignment_pos', 'human_res_pos', 'human_res', 'mouse_res_pos', 'mouse_res', 'conservation']
    colvals  = [index, human_res_pos, human_align_res, mouse_res_pos, mouse_align_res, dis]
    df_residuemap = pd.DataFrame()
    for name, col in zip(colnames, colvals): 
        df_residuemap[name] = col
    df_residuemap = df_residuemap[df_residuemap['human_res'] != '-']
    df_residuemap.to_csv(edits_filedir / residuemap_filename, sep='\t', index=False)

    del index, human_align_res, mouse_align_res, score, dis, v, human_res_pos, mouse_res_pos, colnames, colvals

    return df_alignconserv, df_residuemap
