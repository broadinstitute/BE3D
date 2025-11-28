import os
import sys
import pandas as pd
import numpy as np
import yaml
import glob
import shutil

def load_config(config_yaml):
    with open(config_yaml, "r") as file:
        config = yaml.safe_load(file)

    return config

def main(**kwargs):
    ## REQUIRED
    input_gene=kwargs['input_gene']
    input_uniprot=kwargs['input_uniprot']
    input_chain=kwargs['input_chain']
    screen_dir=kwargs['screen_dir']
    screens=kwargs['screens']
    mut_col=kwargs['mut_col']
    val_col=kwargs['val_col']
    gene_col=kwargs['gene_col']
    edits_col=kwargs['edits_col']
    output_dir=kwargs['output_dir']
    mut_categories=kwargs['mut_categories']
    mut_delimiter=kwargs['mut_delimiter']
    conservation_run=kwargs['conservation_run']
    alt_gene_name=kwargs['alt_gene_name']
    alt_uniprot_id=kwargs['alt_uniprot_id']
    alt_screen_start=kwargs['alt_screen_start']
    
    ## OPTIONAL
    user_fasta=kwargs['user_fasta']
    user_pdb=kwargs['user_pdb']
    user_dssp=kwargs['user_dssp']
    function_for_lfc=kwargs['function_for_lfc']
    function_for_lfc3d=kwargs['function_for_lfc3d']
    function_for_meta=kwargs['function_for_meta']
    nRandom=kwargs['nRandom']
    single_screen_pthr=kwargs['single_screen_pthr']
    multi_screen_pthr=kwargs['multi_screen_pthr']
    structure_radius=kwargs['structure_radius']    
    clustering_radius=kwargs['clustering_radius']
    qa_passed_only=kwargs['qa_passed_only']
    qa_only=kwargs['qa_only']
    qa_controls=kwargs['qa_controls']
    qa_cases=kwargs['qa_cases']
    priority_on_alternative=kwargs['priority_on_alternative']
    ppi_chain_gene_dict=kwargs['ppi_chain_gene_dict']
    ppi_gene_edits_dict=kwargs['ppi_gene_edits_dict']
    v_score_threshold=kwargs['v_score_threshold']
    atom_level_naa=kwargs['atom_level_naa']

    if user_pdb:
        structureid = f'PDB-{input_uniprot}'
    else:
        structureid = f'AF-{input_uniprot}-F1-model_v6'
        
    single_screen_pthr_str = str(single_screen_pthr)
    multi_screen_pthr_str = str(multi_screen_pthr)
    # pthr_str_short = str(pthr).split('.')[1]
    pdb_file = os.path.join(output_dir, "sequence_structure", f"{structureid}_processed.pdb")

    def find_union(input, pthr_str):
        if input[0] == f'p<{pthr_str}' or input[1] == f'p<{pthr_str}':
            return f'p<{pthr_str}'
        else:
            return f'p>={pthr_str}'

    ## ASSIGN VALUES FOR EMPTY VARIABLES
    # if output_dir == '' or not output_dir:
    #     output_dir = os.path.join(workdir, f'{gene}')
    os.makedirs(output_dir, exist_ok=True)
    shutil.copy2(config_yaml, os.path.join(output_dir,os.path.basename(config_yaml)))
    print('All results will be saved in the following directory:')
    print(output_dir)

    screens = [screen.strip() for screen in screens.split(',')]
    screen_names = [s.split('.')[0] for s in screens]
    input_dfs = [pd.read_csv(os.path.join(screen_dir,s), sep='\t') for s in screens]
                          
    df_residuemap = pd.DataFrame()
    conserv_dfs = list()
    gene_list = list()

    if conservation_run:
        _, df_residuemap = conservation(output_dir,
                    input_gene, alt_gene_name,
                    input_uniprot, alt_uniprot_id,
        )
    #     for screen_name in screen_names: # mCherry-high or mCherry-low
    #             conserv_dfs.append(df_residuemap)
    #             gene_list.append(alt_gene_name)

    # else:
    #     for screen_name in screen_names: # mCherry-high or mCherry-low
    #         conserv_dfs.append(None)
    #         gene_list.append(input_gene)

    # TODO: conservation mode : Human + Mouse to Human # SETDDB1
    # TODO: transfer mode : Mouse to Human
    for screen_name in screen_names: # 
        if alt_gene_name and screen_name.startswith(alt_screen_start):
            conserv_dfs.append(df_residuemap)
            gene_list.append(alt_gene_name)
        else:
            conserv_dfs.append(None)
            gene_list.append(input_gene)

    # Only for original Gene
    sequence_structural_features(
        output_dir,
        input_gene, input_uniprot, structureid,
        user_fasta=user_fasta,
        user_pdb=user_pdb,
        user_dssp=user_dssp,
        target_chainid=input_chain,
        radius=structure_radius,
        atom_level_naa=atom_level_naa
        )

    # For All gene
    hypothesis_test(
        output_dir,
        input_dfs, screen_names,
        cases=qa_cases,
        controls=qa_controls,
        comp_name=f"{'_'.join(qa_cases)}_vs_{'_'.join(qa_controls)}",
        mut_col=mut_col,
        val_col=val_col,
        gene_col=gene_col,
        save_type='svg'
        )

    if qa_only:
        sys.exit()
        
    if qa_passed_only:
        h2_ks_test_pd = pd.read_csv(f'{output_dir}/hypothesis_qc/KolmogorovSmirnov_hypothesis2.tsv',sep='\t')
        h2_ks_test_pd = h2_ks_test_pd.replace(-999,None)
        white_screen_list = h2_ks_test_pd[(h2_ks_test_pd[f"p_{'_'.join(qa_cases)}_vs_{'_'.join(qa_controls)}"]<0.05)&(h2_ks_test_pd['gene_name'].isin(gene_list))]['screenid'].to_list()
        print(f'original screen size: {len(screen_names)}, screen white list size: {len(white_screen_list)}, QA-passed screen size: {len(list(set(screen_names).intersection(white_screen_list)))}')
        screen_names = list(set(screen_names).intersection(white_screen_list))
        input_dfs = [pd.read_csv(os.path.join(screen_dir,f'{s}.tsv'), sep='\t') for s in screen_names]
        conserv_dfs = list()
        gene_list = list() # where we need only have genes passed QA, so OVERWRITE the previous raw gene list.
        for screen_name in screen_names:
            if screen_name.startswith(alt_screen_start):
                conserv_dfs.append(df_residuemap)
                gene_list.append(alt_gene_name)
            else:
                conserv_dfs.append(None)
                gene_list.append(input_gene)

    # sys.exit()
    ## WHERE WE CAN HAVE A BARRICADE FOR FILTERING QA_PASSED or ALL screens # Re-organise input_dfs, conservation_dfs

    # For All    
    parse_be_data(
        output_dir,
        input_dfs, input_gene, screen_names,
        mut_col=mut_col,
        val_col=val_col,
        gene_col=gene_col,
        edits_col=edits_col,
        mut_categories = mut_categories,
        mut_delimiter = mut_delimiter,
        conserv_dfs = conserv_dfs,
        conserv_col='alternative_res_pos',
        gene_list=gene_list,
        v_score_threshold=v_score_threshold
        )

    # For All
    plot_rawdata(
        output_dir,
        input_dfs,
        screen_names,
        mut_col=mut_col,
        val_col=val_col,
        gene_col=gene_col,
        mut_categories = mut_categories,
        save_type='svg'
        )

    # For All
    df_missense_list = [
        pd.read_csv(f'{output_dir}/screendata/{gene}_{screen_name}_Missense.tsv',
                    sep='\t') for gene, screen_name in zip(gene_list, screen_names)
    ]

    # For All
    for df_missense, screen_name, gene in zip(df_missense_list, screen_names, gene_list):
        randomize_data(
            df_missense,
            output_dir, gene,
            screen_name,
            nRandom=nRandom,
            seed=True,
            )

    ## PRIORITIZE
    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

    # Human and Non-human mixed
    for gene, screen_name, df_consrv in zip(gene_list, screen_names, conserv_dfs):
        df_control = pd.read_csv(f'{output_dir}/screendata/{gene}_{screen_name}_No_Mutation.tsv', sep='\t', index_col=0)
        df_dict = {}

        for mut in ['Missense', 'Silent', 'Nonsense']:
            filepath = f'{output_dir}/screendata/{gene}_{screen_name}_{mut}.tsv'
            if os.path.exists(filepath):
                df_dict[mut] = pd.read_csv(filepath, sep='\t', index_col=0)
        
        if df_consrv is not None:
            target_res_pos, target_res = 'original_res_pos', 'unires'
            alternate_res_pos, alternate_res = 'alternative_res_pos','alternative_res'        
            df_missense = prioritize_by_sequence(
                df_dict,
                df_struc, df_consrv, df_control,
                output_dir,
                gene, screen_name,
                target_res_pos=target_res_pos, 
                alt_res_pos=alternate_res_pos,
                alt_res=alternate_res
            )
        else:
            df_missense = prioritize_by_sequence(
                df_dict,
                df_struc, df_consrv, df_control,
                output_dir,                
                gene, screen_name,
            )

        df_rand = pd.read_csv(f'{output_dir}/screendata_rand/{gene}_{screen_name}_Missense_rand.tsv.gz', sep='\t')

        # For ALL
        randomize_sequence(
            df_missense, df_rand,
            output_dir,
            gene, screen_name,
            nRandom=nRandom, conservation=False,
            muttype='Missense',
            function_name=function_for_lfc,
            target_pos='unipos', target_res=None,
            )

        plot_screendata_sequence(
            df_missense,
            output_dir,
            gene, screen_name, function_name=function_for_lfc, muttype='Missense',
        )

    def run_Clust3D_per_species(gene, screen_names, gene_type, pdb_file):    
        ## CLUSTERING ON PRIORITIZED
        df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')
        ## LFC3D 
        df_edits_list = []
        df_rand_list = []
        
        for screen_name in screen_names:
            df_missense = pd.read_csv(f'{output_dir}/screendata_sequence/{gene}_{screen_name}_protein_edits.tsv', sep='\t')
            df_missense[f'{function_for_lfc}_Missense_LFC_plab_input'] = df_missense[f'{function_for_lfc}_Missense_LFC_p'].apply(lambda x: f'p<{single_screen_pthr_str}' if x < single_screen_pthr else f'p>={single_screen_pthr_str}')

            df_hits_clust, distances, yvalues = clustering(
                df_struc, df_missense,
                output_dir, gene,
                psig_columns=[f'{function_for_lfc}_Missense_LFC_plab_input'],
                pthr_cutoffs=[f'p<{single_screen_pthr_str}'], # TODO:
                screen_name=screen_name, score_type='LFC',
                max_distances=20, merge_cols=['unipos', 'chain'],
                atom_level= pdb_file if atom_level_naa == True else False
            )

            df_protein_edits = pd.read_csv(f'{output_dir}/screendata_sequence/{gene}_{screen_name}_protein_edits.tsv', sep='\t')
            df_edits_list.append(df_protein_edits)
            df_protein_edits_rand = pd.read_csv(f'{output_dir}/screendata_sequence_rand/{gene}_{screen_name}_Missense_protein_edits_rand.tsv.gz', sep='\t')
            df_rand_list.append(df_protein_edits_rand)

        df_LFC_LFC3D = calculate_lfc3d(
            df_struc, df_edits_list, df_rand_list,
            output_dir, gene, screen_names,
            nRandom=nRandom,  muttype='Missense',
            function_type_lfc=function_for_lfc,
            function_type_lfc3d=function_for_lfc3d,
            conserved_only=False, gene_type=gene_type,
            target_gene_chain=input_chain,
            ppi_chain_gene_dict=ppi_chain_gene_dict,
            ppi_gene_edits_dict=ppi_gene_edits_dict
        )
        
        # LFC #
        df_bidir = average_split_score(
            df_LFC_LFC3D,
            output_dir, gene, screen_names,
            score_type='LFC',gene_type=gene_type
        )
        df_dis, _, _ = bin_score(
            df_bidir,
            output_dir, gene, screen_names,
            score_type='LFC',gene_type=gene_type
        )
        
        znorm_score(
            df_bidir, output_dir, gene, screen_names,
            pthrs=[0.05, 0.01, 0.001], score_type='LFC', gene_type=gene_type
        )
        df_lfc = pd.read_csv(f'{output_dir}/LFC/{gene_type}_{gene}_NonAggr_LFC.tsv', sep='\t')

        for screen_name in screen_names:
            average_split_bin_plots(
                df_lfc,
                workdir = output_dir,
                input_gene = gene,
                screen_name=screen_name, # BLANK FOR META #
                func='', # BLANK FOR NON AGGR #
                pthr=single_screen_pthr,
                score_type='LFC',
                aggregate_dir='LFC',
                save_type='svg'
                )
        
        df_bidir = average_split_score(
            df_LFC_LFC3D,
            output_dir, gene, screen_names,
            score_type='LFC3D',gene_type=gene_type
        )
        df_dis, _, _ = bin_score(
            df_bidir,
            output_dir, gene, screen_names,
            score_type='LFC3D',gene_type=gene_type
        )

        znorm_score(
            df_bidir, output_dir, gene, screen_names,
            pthrs=[0.05, 0.01, 0.001], score_type='LFC3D', gene_type=gene_type
        )

        df_lfc3d = pd.read_csv(f'{output_dir}/LFC3D/{gene_type}_{gene}_NonAggr_LFC3D.tsv', sep='\t')
        
        for screen_name in screen_names:
            average_split_bin_plots(
                df_lfc3d,
                workdir = output_dir,
                input_gene = gene,
                screen_name=screen_name, # BLANK FOR META #
                func='', # BLANK FOR NON AGGR #
                pthr=single_screen_pthr,
                score_type='LFC3D',
                aggregate_dir='LFC3D',
                save_type='svg'
                )
    
        for score_type in ['LFC', 'LFC3D']:
            df_pvals = pd.read_csv(f'{output_dir}/{score_type}/{gene_type}_{gene}_NonAggr_{score_type}.tsv', sep='\t')

            for screen_name in screen_names:
                df_hits_clust, distances, yvalues = clustering( 
                    df_struc, df_pvals,
                    output_dir, gene,
                    psig_columns=[f'{screen_name}_{score_type}_neg_05_psig',
                                f'{screen_name}_{score_type}_pos_05_psig',
                                f'{screen_name}_{score_type}_neg_01_psig',
                                f'{screen_name}_{score_type}_pos_01_psig',
                                f'{screen_name}_{score_type}_neg_001_psig',
                                f'{screen_name}_{score_type}_pos_001_psig'
                                ],
                    pthr_cutoffs=['p<0.05', 'p<0.05',
                                'p<0.01', 'p<0.01',
                                'p<0.001', 'p<0.001',
                                ],
                    screen_name=screen_name, score_type=score_type,
                    max_distances=20, merge_cols=['unipos', 'chain'],
                    atom_level= pdb_file if atom_level_naa == True else False
                )

                # PLOTTING #
                plot_clustering(
                    df_struc, df_pvals,
                    df_hits_clust, clustering_radius,
                    output_dir, gene,
                    distances, yvalues,
                    names=['Negative', 'Positive', 'Negative', 'Positive', 'Negative', 'Positive'],
                    psig_columns=[f'{screen_name}_{score_type}_neg_05_psig',
                                f'{screen_name}_{score_type}_pos_05_psig',
                                f'{screen_name}_{score_type}_neg_01_psig',
                                f'{screen_name}_{score_type}_pos_01_psig',
                                f'{screen_name}_{score_type}_neg_001_psig',
                                f'{screen_name}_{score_type}_pos_001_psig'
                                ],
                    pthr_cutoffs=['p<0.05', 'p<0.05',
                                'p<0.01', 'p<0.01',
                                'p<0.001', 'p<0.001',
                                ],                    
                    screen_name=screen_name, score_type=score_type,
                    merge_col=['unipos', 'chain'],
                    save_type='svg',
                    dendrogram_subplots_kwargs={'figsize':(15, 3.5)}
                )
        
        for screen_name in screen_names:
            # LFC vs LFC3D SCATTERPLOT #
            cutoff_str = str(single_screen_pthr).split('.')[1]
            df_lfc = pd.read_csv(f"{output_dir}/LFC/{gene_type}_{gene}_NonAggr_LFC.tsv", sep='\t')
            df_lfc_dis = pd.read_csv(f"{output_dir}/LFC/{gene_type}_{gene}_LFC_dis_wght.tsv", sep='\t')
            df_lfc3d = pd.read_csv(f"{output_dir}/LFC3D/{gene_type}_{gene}_NonAggr_LFC3D.tsv", sep='\t')
            df_lfc3d_dis = pd.read_csv(f"{output_dir}/LFC3D/{gene_type}_{gene}_LFC3D_dis_wght.tsv", sep='\t')

            df_dis_input = pd.DataFrame()
            df_dis_input['unipos'] = df_lfc_dis['unipos']

            df_dis_input = pd.concat([df_dis_input,df_lfc_dis.filter(regex=r'LFC$'),df_lfc3d_dis.filter(regex=r'LFC3D$'),df_lfc3d.filter(regex=fr'LFC3D_dis$|{cutoff_str}_psig$')])
            df_dis_input = df_dis_input.rename({f'{screen_name}_LFC3D_neg_{cutoff_str}_psig':f'{screen_name}_LFC3D_neg_psig',
                                                f'{screen_name}_LFC3D_pos_{cutoff_str}_psig':f'{screen_name}_LFC3D_pos_psig',
                                                },axis=1)
            lfc_lfc3d_scatter(
                df_input=df_dis_input,
                workdir=output_dir,
                input_gene=gene, screen_name=screen_name,
                pthr=single_screen_pthr,
            )
            os.rename(f'{output_dir}/characterization/plots/{gene}_LFC_LFC3D_scatter.png',
                    f'{output_dir}/characterization/plots/{gene}_LFC_LFC3D_scatter_{cutoff_str}_{screen_name}.png')
        
        # Load both LFC and lFC3D dataframes
        df_pvals_LFC3D = pd.read_csv(f'{output_dir}/LFC3D/{gene_type}_{gene}_NonAggr_LFC3D.tsv', sep='\t')
        df_pvals_LFC = pd.read_csv(f'{output_dir}/LFC/{gene_type}_{gene}_NonAggr_LFC.tsv', sep='\t')
        df_pvals = pd.concat([df_pvals_LFC3D, df_pvals_LFC.drop(columns=['unipos', 'unires', 'chain'])], axis=1)

        # Find union of LFC and LFC3D
        for screen_name in screen_names:
            for each_pthr in ['05','01','001']:
                df_pvals[f'{screen_name}_union_neg_{each_pthr}_psig'] = df_pvals[[f'{screen_name}_LFC_neg_{each_pthr}_psig', f'{screen_name}_LFC3D_neg_{each_pthr}_psig']].apply(lambda row: find_union(row, f'0.{each_pthr}'), axis=1)
                df_pvals[f'{screen_name}_union_pos_{each_pthr}_psig'] = df_pvals[[f'{screen_name}_LFC_pos_{each_pthr}_psig', f'{screen_name}_LFC3D_pos_{each_pthr}_psig']].apply(lambda row: find_union(row, f'0.{each_pthr}'), axis=1)

            df_hits_clust, distances, yvalues = clustering(
                df_struc, df_pvals,
                output_dir, gene,
                psig_columns=[f'{screen_name}_union_neg_05_psig',
                            f'{screen_name}_union_pos_05_psig',
                            f'{screen_name}_union_neg_01_psig',
                            f'{screen_name}_union_pos_01_psig',
                            f'{screen_name}_union_neg_001_psig',
                            f'{screen_name}_union_pos_001_psig'
                            ],
                pthr_cutoffs=['p<0.05', 'p<0.05',
                            'p<0.01', 'p<0.01',
                            'p<0.001', 'p<0.001',
                            ],                
                screen_name=screen_name, score_type='union',
                max_distances=20, merge_cols=['unipos', 'chain'],
                atom_level= pdb_file if atom_level_naa == True else False
            )

            # PLOTTING #
            plot_clustering(
                df_struc, df_pvals,
                df_hits_clust, clustering_radius,
                output_dir, gene,
                distances, yvalues,
                names=['Negative', 'Positive', 'Negative', 'Positive', 'Negative', 'Positive'],
                psig_columns=[f'{screen_name}_union_neg_05_psig',
                            f'{screen_name}_union_pos_05_psig',
                            f'{screen_name}_union_neg_01_psig',
                            f'{screen_name}_union_pos_01_psig',
                            f'{screen_name}_union_neg_001_psig',
                            f'{screen_name}_union_pos_001_psig'
                            ],
                pthr_cutoffs=['p<0.05', 'p<0.05',
                            'p<0.01', 'p<0.01',
                            'p<0.001', 'p<0.001',
                            ],                
                screen_name = screen_name, score_type='union',
                merge_col=['unipos', 'chain'],
                save_type='svg',
                dendrogram_subplots_kwargs={'figsize':(15, 3.5)}
            )
        return df_LFC_LFC3D
    
    df_LFC_LFC3D_original, df_LFC_LFC3D_alt, merged_df_LFC_LFC3D = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    
    if conservation_run:
        df_LFC_LFC3D_original = run_Clust3D_per_species(input_gene, [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is None], gene_type='Original', pdb_file=pdb_file) # For human
        df_LFC_LFC3D_alt = run_Clust3D_per_species(alt_gene_name, [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is not None], gene_type='Alternative', pdb_file=pdb_file) # For non-Human
        merged_df_LFC_LFC3D =  df_LFC_LFC3D_original.merge(df_LFC_LFC3D_alt, on='unipos', how='inner', suffixes=('', f'_{alt_gene_name}'))
        merged_df_LFC_LFC3D = merged_df_LFC_LFC3D.sort_index()
        gene = 'Merged'
        
    else:
        gene = input_gene
        merged_df_LFC_LFC3D =  run_Clust3D_per_species(input_gene, [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is None], gene_type='Original', pdb_file=pdb_file) # For human
    
    
    def merge_two_tsvs(gene_name, results_dir, file_pattern, original_gene, alternative_gene, priority_on_alternative=False, compression=False):
        if compression:
            original_tsv = os.path.join(results_dir,f'Original_{original_gene}_{file_pattern}.tsv.gz')
            alternative_tsv = os.path.join(results_dir,f'Alternative_{alternative_gene}_{file_pattern}.tsv.gz')
        else:
            original_tsv = os.path.join(results_dir,f'Original_{original_gene}_{file_pattern}.tsv')
            alternative_tsv = os.path.join(results_dir,f'Alternative_{alternative_gene}_{file_pattern}.tsv')

        if os.path.exists(original_tsv) and os.path.exists(alternative_tsv):
            df1 = pd.read_csv(original_tsv,sep='\t')
            df2 = pd.read_csv(alternative_tsv,sep='\t')
            if not priority_on_alternative:
                # Automatically detect shared columns
                common_cols = list(set(df1.columns) & set(df2.columns))
                merged = pd.merge(df1, df2, on=common_cols, how="outer")
            else:
                merged = pd.merge(df2, df1, how="outer")
            # Merge on all shared columns
            if compression:
                merged.to_csv(os.path.join(results_dir,f'{gene_name}_{file_pattern}.tsv.gz'), sep="\t", index=False, compression='gzip')
            else:
                merged.to_csv(os.path.join(results_dir,f'{gene_name}_{file_pattern}.tsv'), sep="\t", index=False)
        elif os.path.exists(original_tsv):
            df1 = pd.read_csv(original_tsv,sep='\t')
            if compression:
                df1.to_csv(os.path.join(results_dir,f'{gene_name}_{file_pattern}.tsv.gz'), sep="\t", index=False, compression='gzip')
            else:
                df1.to_csv(os.path.join(results_dir,f'{gene_name}_{file_pattern}.tsv'), sep="\t", index=False)
        else:
            df2 = pd.read_csv(alternative_tsv,sep='\t')
            if compression:
                df2.to_csv(os.path.join(results_dir,f'{gene_name}_{file_pattern}.tsv.gz'), sep="\t", index=False, compression='gzip')
            else:
                df2.to_csv(os.path.join(results_dir,f'{gene_name}_{file_pattern}.tsv'), sep="\t", index=False)

    merge_two_tsvs(input_gene,f'{output_dir}/LFC','LFC_bidirectional', input_gene, alt_gene_name, priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene,f'{output_dir}/LFC','LFC_dis_wght',input_gene, alt_gene_name,priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene,f'{output_dir}/LFC','NonAggr_LFC',input_gene, alt_gene_name,priority_on_alternative=priority_on_alternative)
    
    merge_two_tsvs(input_gene,f'{output_dir}/LFC3D','LFC3D_bidirectional',input_gene, alt_gene_name,priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene,f'{output_dir}/LFC3D','LFC3D_dis_wght',input_gene, alt_gene_name,priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene,f'{output_dir}/LFC3D','NonAggr_LFC3D',input_gene, alt_gene_name,priority_on_alternative=priority_on_alternative)    
    merge_two_tsvs(input_gene,f'{output_dir}/LFC3D','LFC_LFC3D_LFC3Dr',input_gene, alt_gene_name,priority_on_alternative=priority_on_alternative,compression=True)    
    
    if len(screen_names) > 1: # Meta-Aggregation automatically runs if the number of screens > 1
        # META-AGGREGATION ON LFC3D
        df_bidir_meta = average_split_meta(
            merged_df_LFC_LFC3D,
            output_dir, gene, screen_names,
            nRandom=nRandom,
            score_type='LFC3D',
            aggr_func_name=function_for_meta,
        )
        df_dis, _, _ = bin_meta(
            df_bidir_meta,
            output_dir, gene,
            score_type='LFC3D', aggr_func_name=function_for_meta,
        )
        
        znorm_meta(
            df_dis,
            output_dir, gene, screen_names,
            pthrs=[0.05, 0.01, 0.001], score_type='LFC3D', aggr_func_name=function_for_meta,
        )

        df_lfc3d = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_LFC3D.tsv', sep='\t')

        # for screen_name in screen_names:
        average_split_bin_plots(
            df_lfc3d,
            workdir = output_dir,
            input_gene = gene,
            screen_name='', # BLANK FOR META #
            func=function_for_meta, # BLANK FOR NON AGGR #
            pthr=multi_screen_pthr,
            score_type='LFC3D',
            aggregate_dir='meta-aggregate',
            save_type='svg'
            )

        # META-AGGREGATION ON LFC
        df_bidir_meta = average_split_meta(
            merged_df_LFC_LFC3D,
            output_dir, gene, screen_names,
            nRandom=nRandom,
            score_type='LFC',
            aggr_func_name=function_for_meta,
        )
        df_dis = bin_meta(
            df_bidir_meta,
            output_dir, gene,
            score_type='LFC', aggr_func_name=function_for_meta,
        )
        
        znorm_meta(
            df_bidir_meta,
            output_dir, gene, screen_names,
            pthrs=[0.05, 0.01, 0.001], score_type='LFC', aggr_func_name=function_for_meta,
        )

        df_lfc = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_LFC.tsv', sep='\t')

        for screen_name in screen_names:
            average_split_bin_plots(
                df_lfc,
                workdir = output_dir,
                input_gene = gene,
                screen_name='', # BLANK FOR META #
                func=function_for_meta, # BLANK FOR NON AGGR #
                pthr=multi_screen_pthr,
                score_type='LFC',
                aggregate_dir='meta-aggregate',
                save_type='svg'
                )
            
        ## CLUSTERING META-
        df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

        for score_type in ['LFC', 'LFC3D']:
            df_pvals = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_{score_type}.tsv', sep='\t')

            df_hits_clust, distances, yvalues = clustering(
                df_struc, df_pvals,
                output_dir, gene,
                # psig_columns=[f'{function_for_meta}_{score_type}_neg_{pthr_str_short}_psig',
                #             f'{function_for_meta}_{score_type}_pos_{pthr_str_short}_psig'],
                # pthr_cutoffs=[f'p<{pthr_str}', f'p<{pthr_str}'],
                psig_columns=[f'{function_for_meta}_{score_type}_neg_05_psig',
                            f'{function_for_meta}_{score_type}_pos_05_psig',
                            f'{function_for_meta}_{score_type}_neg_01_psig',
                            f'{function_for_meta}_{score_type}_pos_01_psig',
                            f'{function_for_meta}_{score_type}_neg_001_psig',
                            f'{function_for_meta}_{score_type}_pos_001_psig'
                            ],
                pthr_cutoffs=['p<0.05', 'p<0.05',
                            'p<0.01', 'p<0.01',
                            'p<0.001', 'p<0.001',
                            ],            
                screen_name='Meta', score_type=score_type,
                max_distances=20, merge_cols=['unipos', 'chain'],
                atom_level= pdb_file if atom_level_naa == True else False
            )

            # PLOTTING #
            plot_clustering(
                df_struc, df_pvals,
                df_hits_clust, clustering_radius,
                output_dir, gene,
                distances, yvalues,
                names=['Negative', 'Positive', 'Negative', 'Positive', 'Negative', 'Positive'],
                psig_columns=[f'{function_for_meta}_{score_type}_neg_05_psig',
                            f'{function_for_meta}_{score_type}_pos_05_psig',
                            f'{function_for_meta}_{score_type}_neg_01_psig',
                            f'{function_for_meta}_{score_type}_pos_01_psig',
                            f'{function_for_meta}_{score_type}_neg_001_psig',
                            f'{function_for_meta}_{score_type}_pos_001_psig'
                            ],
                pthr_cutoffs=['p<0.05', 'p<0.05',
                            'p<0.01', 'p<0.01',
                            'p<0.001', 'p<0.001',
                            ],            
                screen_name='Meta', score_type=score_type,
                merge_col=['unipos', 'chain'],
                save_type='svg',
                dendrogram_subplots_kwargs={'figsize':(15, 3.5)}
            )

        ## CLUSTERING ON UNION
        df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

        # Load both LFC and lFC3D dataframes
        df_pvals_LFC3D = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_LFC3D.tsv', sep='\t')
        df_pvals_LFC = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_LFC.tsv', sep='\t')
        df_pvals = pd.concat([df_pvals_LFC3D, df_pvals_LFC.drop(columns=['unipos', 'unires', 'chain'])], axis=1)

        # Find union of LFC and LFC3D
        for each_pthr in ['05','01','001']:        
            df_pvals[f'{function_for_meta}_union_neg_{each_pthr}_psig'] = df_pvals[[f'{function_for_meta}_LFC_neg_{each_pthr}_psig', f'{function_for_meta}_LFC3D_neg_{each_pthr}_psig']].apply(lambda row: find_union(row, f'0.{each_pthr}'), axis=1)
            df_pvals[f'{function_for_meta}_union_pos_{each_pthr}_psig'] = df_pvals[[f'{function_for_meta}_LFC_pos_{each_pthr}_psig', f'{function_for_meta}_LFC3D_pos_{each_pthr}_psig']].apply(lambda row: find_union(row, f'0.{each_pthr}'), axis=1)

        df_hits_clust, distances, yvalues = clustering(
            df_struc, df_pvals,
            output_dir, gene,
            psig_columns=[f'{function_for_meta}_union_neg_05_psig',
                        f'{function_for_meta}_union_pos_05_psig',
                        f'{function_for_meta}_union_neg_01_psig',
                        f'{function_for_meta}_union_pos_01_psig',
                        f'{function_for_meta}_union_neg_001_psig',
                        f'{function_for_meta}_union_pos_001_psig'
                        ],
            pthr_cutoffs=['p<0.05', 'p<0.05',
                        'p<0.01', 'p<0.01',
                        'p<0.001', 'p<0.001',
                        ],        
            screen_name='Meta', score_type='union',
            max_distances=20, merge_cols=['unipos', 'chain'],
            atom_level= pdb_file if atom_level_naa == True else False
        )

        # PLOTTING #
        plot_clustering(
            df_struc, df_pvals,
            df_hits_clust, clustering_radius,
            output_dir, gene,
            distances, yvalues,
            names=['Negative', 'Positive', 'Negative', 'Positive', 'Negative', 'Positive'],
            psig_columns=[f'{function_for_meta}_union_neg_05_psig',
                        f'{function_for_meta}_union_pos_05_psig',
                        f'{function_for_meta}_union_neg_01_psig',
                        f'{function_for_meta}_union_pos_01_psig',
                        f'{function_for_meta}_union_neg_001_psig',
                        f'{function_for_meta}_union_pos_001_psig'
                        ],
            pthr_cutoffs=['p<0.05', 'p<0.05',
                        'p<0.01', 'p<0.01',
                        'p<0.001', 'p<0.001',
                        ],        
            screen_name='Meta', score_type='union',
            merge_col=['unipos', 'chain'],
            save_type='svg',
            dendrogram_subplots_kwargs={'figsize':(15, 3.5)}
        )
   
        # CHARACTERIZATION
        df_domains = pd.read_csv(f'{output_dir}/sequence_structure/{input_gene}_{input_uniprot}_domains.tsv', sep='\t')
        df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

        for cutoff in [multi_screen_pthr]:
            for score_type in ['LFC', 'LFC3D']:
                # ENRICHMENT TEST #
                cutoff_str = str(cutoff).split('.')[1]
                if conservation_run:
                    df_meta = pd.read_csv(f'{output_dir}/meta-aggregate/Merged_MetaAggr_{score_type}.tsv', sep='\t')
                else:
                    df_meta = pd.read_csv(f'{output_dir}/meta-aggregate/{input_gene}_MetaAggr_{score_type}.tsv', sep='\t')
                

                input_df = pd.concat([df_meta, df_domains['Domain'], df_struc['pLDDT_dis']], axis=1)

                hit_columns = [f'{function_for_meta}_{score_type}_neg_{cutoff_str}_p', f'{function_for_meta}_{score_type}_pos_{cutoff_str}_p']
                input_df[hit_columns] = input_df[hit_columns].replace('-', np.nan).astype(float)

                results = enrichment_test(
                    input_df,
                    workdir=output_dir,
                    input_gene=input_gene,
                    hit_columns=hit_columns,
                    hit_threshold=cutoff,
                    feature_column='pLDDT_dis',
                    feature_values=['confident', 'low', 'very low'],
                    confidence_level=0.95,
                )

                plot_enrichment_test(
                    enrichment_results=results,
                    workdir=output_dir,
                    input_gene=input_gene,
                    hit_value=cutoff,
                    feature_values=['confident', 'low', 'very low'],
                )

                os.rename(f'{output_dir}/characterization/{input_gene}_enrichment_test.pickle',
                        f'{output_dir}/characterization/{input_gene}_enrichment_test_{score_type}_{cutoff_str}.pickle')
                os.rename(f'{output_dir}/characterization/plots/{input_gene}_enrichment_test.png',
                        f'{output_dir}/characterization/plots/{input_gene}_enrichment_test_{score_type}_{cutoff_str}.png')

                # BARPLOTS #
                colnames = [f'{function_for_meta}_{score_type}_neg_{cutoff_str}_psig', f'{function_for_meta}_{score_type}_pos_{cutoff_str}_psig']
                input_df = pd.concat([df_struc, df_meta[colnames], df_domains['Domain']], axis=1)

                hits_feature_barplot(
                    input_df,
                    workdir=output_dir,
                    input_gene=input_gene,
                    category_col='pLDDT_dis',
                    values_cols=colnames, values_vals=[f'p<0.{cutoff_str}', f'p<0.{cutoff_str}'], value_names=['NEG', 'POS'],
                    plot_type='Count', colors = ['darkred', 'darkblue'],
                )

                os.rename(f'{output_dir}/characterization/plots/{input_gene}_Count_pLDDT_dis_barplot.png',
                        f'{output_dir}/characterization/plots/{input_gene}_Count_pLDDT_dis_barplot_{score_type}_{cutoff_str}.png')

            # SCATTERPLOT #
            if conservation_run:
                df_lfc_dis = pd.read_csv(f"{output_dir}/meta-aggregate/Merged_LFC_dis_wght.tsv", sep='\t')
                df_lfc3d_dis = pd.read_csv(f"{output_dir}/meta-aggregate/Merged_LFC3D_dis_wght.tsv", sep='\t')
                df_lfc = pd.read_csv(f"{output_dir}/meta-aggregate/Merged_MetaAggr_LFC.tsv", sep='\t')
                df_lfc3d = pd.read_csv(f"{output_dir}/meta-aggregate/Merged_MetaAggr_LFC3D.tsv", sep='\t')
            else:
                df_lfc_dis = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_LFC_dis_wght.tsv", sep='\t')
                df_lfc3d_dis = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_LFC3D_dis_wght.tsv", sep='\t')
                df_lfc = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_MetaAggr_LFC.tsv", sep='\t')
                df_lfc3d = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_MetaAggr_LFC3D.tsv", sep='\t')
                
            df_dis_input = pd.DataFrame()
            df_dis_input['unipos'] = df_lfc_dis['unipos']
            for screen_name in screen_names:
                df_dis_input[f'Meta_LFC'] = df_lfc_dis[f'{function_for_meta}_LFC']
                df_dis_input[f'Meta_LFC3D'] = df_lfc3d_dis[f'{function_for_meta}_LFC3D']
                df_dis_input[f'Meta_LFC3D_dis'] = df_lfc3d_dis[f'{function_for_meta}_LFC3D_dis']

            df_dis_input[f'{function_for_meta}_LFC3D_neg_{cutoff_str}_psig'] = df_lfc3d[f'{function_for_meta}_LFC3D_neg_{cutoff_str}_psig']
            df_dis_input[f'{function_for_meta}_LFC3D_pos_{cutoff_str}_psig'] = df_lfc3d[f'{function_for_meta}_LFC3D_pos_{cutoff_str}_psig']

            df_dis_input = df_dis_input.rename(columns={
                f'{function_for_meta}_LFC3D_neg_{cutoff_str}_psig': f"Meta_LFC3D_neg_psig",
                f'{function_for_meta}_LFC3D_pos_{cutoff_str}_psig': f"Meta_LFC3D_pos_psig",
                })

            lfc_lfc3d_scatter(
                df_input=df_dis_input,
                workdir=output_dir,
                input_gene=input_gene, screen_name='Meta',
                pthr=cutoff,
            )
            os.rename(f'{output_dir}/characterization/plots/{input_gene}_LFC_LFC3D_scatter.png',
                    f'{output_dir}/characterization/plots/{input_gene}_LFC_LFC3D_scatter_{cutoff_str}_Meta.png')

        # SCATTERPLOT #
        if conservation_run:
            df_meta_wght = pd.read_csv(f"{output_dir}/meta-aggregate/Merged_{score_type}_dis_wght.tsv", sep='\t')
        else:
            df_meta_wght = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_{score_type}_dis_wght.tsv", sep='\t')
        df_input = pd.concat([df_struc, df_meta_wght], axis=1)

        df_input[f'{score_type}_wght'] = df_input[f'{function_for_meta}_{score_type}_wght'].abs() * 100
        df_input['direction'] = np.where(
            df_meta_wght[f'{function_for_meta}_{score_type}_wght'].astype(float) > 0, 'POS', np.where(df_meta_wght[f'{function_for_meta}_{score_type}_wght'].astype(float) < 0, 'NEG', 'ZERO')
        )
        df_input = df_input[df_input['direction'].isin(['NEG', 'POS'])]
        df_input = df_input[~df_input['bfactor_pLDDT'].isin(['-'])]
        df_input = df_input[~df_input['RSA'].isin(['-'])]
        df_input['bfactor_pLDDT'] = df_input['bfactor_pLDDT'].astype(float)
        df_input['RSA'] = df_input['RSA'].astype(float)

        pLDDT_RSA_scatter(
            df_input,
            workdir=output_dir,
            input_gene=input_gene,
            pLDDT_col='bfactor_pLDDT', RSA_col='RSA', size_col=f'{score_type}_wght', direction_col='direction',
            color_map = {'NEG': 'darkred', 'POS': 'darkblue'}
        )
        g2p_formatted_hit_cluster(output_dir, gene_list, screen_names, lfc_pthr=single_screen_pthr_str.split('.')[1], lfc3d_pthr=single_screen_pthr_str.split('.')[1], meta_pthr=multi_screen_pthr_str.split('.')[1], function_for_meta=function_for_meta, conservation=conservation_run, input_gene=input_gene)
    else:
        g2p_formatted_hit_cluster(output_dir, gene_list, screen_names, lfc_pthr=single_screen_pthr_str.split('.')[1], lfc3d_pthr=single_screen_pthr_str.split('.')[1], meta_pthr=multi_screen_pthr_str.split('.')[1], function_for_meta=False, conservation=conservation_run, input_gene=input_gene)

if __name__ == '__main__':
    config_yaml = sys.argv[1]
    config = load_config(config_yaml)
    beclust3d_path = config['beclust3d_path']
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), beclust3d_path)))

    from beclust3d.lfc3d.structure import sequence_structural_features
    from beclust3d.lfc3d.preprocess_data import parse_be_data
    from beclust3d.lfc3d.preprocess_data_plot import plot_rawdata
    from beclust3d.qc.hypothesis_tests import hypothesis_test
    from beclust3d.lfc3d.randomize_data import randomize_data
    from beclust3d.lfc3d.prioritize_sequence import prioritize_by_sequence
    from beclust3d.lfc3d.prioritize_sequence_plot import plot_screendata_sequence
    from beclust3d.lfc3d.randomize_sequence import randomize_sequence
    from beclust3d.lfc3d.clustering import clustering
    from beclust3d.lfc3d.clustering_plot import plot_clustering
    from beclust3d.lfc3d.calculate_lfc3d import calculate_lfc3d
    from beclust3d.aggregate.metaaggregate import average_split_meta, bin_meta, znorm_meta
    from beclust3d.aggregate.aggregate_plot import average_split_bin_plots
    from beclust3d.lfc3d.characterization import enrichment_test
    from beclust3d.lfc3d.characterization_plot import plot_enrichment_test, hits_feature_barplot, lfc_lfc3d_scatter, pLDDT_RSA_scatter
    from beclust3d.aggregate.nonaggregate import average_split_score, bin_score, znorm_score
    from beclust3d.lfc3d.conservation import conservation
    from beclust3d.helpers.visualization.g2p import g2p_formatted_hit_cluster
    
    ## REQUIRED
    input_gene = config['input_gene']
    input_uniprot = config['input_uniprot']
    input_chain = config['input_chain']
    screen_dir = config['screen_dir']
    screens = config['screens']
    output_dir= config['output_dir']
    
    conservation_run=config['conservation']['run']
    v_score_threshold = config['conservation']['v_score_threshold']
    alt_gene_name=config['conservation']['alt_gene_name']
    alt_uniprot_id=config['conservation']['alt_uniprot_id']
    alt_screen_start = config['conservation']['alt_screen_start']

    mut_col = config['database']['mut_col']
    val_col = config['database']['val_col']
    gene_col = config['database']['gene_col']
    edits_col = config['database']['edits_col']
    gRNA_col = config['database']['gRNA_col']
    mut_delimiter = config['database']['mut_delimiter']

    mut_categories = list()
    mut_categories.extend(config['mutation_category']['nonsense'])
    mut_categories.extend(config['mutation_category']['splice'])
    mut_categories.extend(config['mutation_category']['missense'])
    mut_categories.extend(config['mutation_category']['silent'])
    mut_categories.extend(config['mutation_category']['no_mutation'])
    mut_categories.extend(config['mutation_category']['intron'])
    
    # OPTIONAL
    user_fasta = config['user_fasta']
    user_pdb = config['user_pdb']
    user_dssp = config['user_dssp']
    function_for_lfc=config['function_for_lfc']
    function_for_lfc3d=config['function_for_lfc3d']    
    function_for_meta=config['function_for_meta']    
    nRandom = config['nRandom']
    single_screen_pthr = config['pthr']['single_screen']
    multi_screen_pthr = config['pthr']['multi_screen']
    structure_radius = config['structure_radius'] 
    clustering_radius = config['clustering_radius']
    qa_passed_only = config['qa']['qa_passed_only']
    qa_only = config['qa']['qa_only']
    qa_controls = config['qa']['controls']
    qa_cases = config['qa']['cases']
    priority_on_alternative = config['priority_on_alternative'],
    ppi_chain_gene_dict = config['ppi_chain_gene_dict']
    ppi_gene_edits_dict = config['ppi_gene_edits_dict']
    atom_level_naa = config['atom_level_naa']
    
    main(input_gene=input_gene, input_uniprot=input_uniprot, input_chain=input_chain, screen_dir=screen_dir,\
        screens=screens, mut_col=mut_col, val_col=val_col, gene_col=gene_col, edits_col=edits_col,\
        gRNA_col=gRNA_col, output_dir=output_dir, user_fasta=user_fasta, user_pdb=user_pdb, user_dssp=user_dssp, nRandom=nRandom,\
        single_screen_pthr=single_screen_pthr,multi_screen_pthr=multi_screen_pthr,
        structure_radius=structure_radius, clustering_radius=clustering_radius, function_for_lfc=function_for_lfc, function_for_lfc3d=function_for_lfc3d,\
        mut_categories=mut_categories,mut_delimiter=mut_delimiter, conservation_run=conservation_run, alt_gene_name=alt_gene_name,alt_uniprot_id=alt_uniprot_id,
        alt_screen_start=alt_screen_start,v_score_threshold=v_score_threshold,
        function_for_meta=function_for_meta, qa_passed_only=qa_passed_only, qa_only=qa_only, qa_controls=qa_controls, qa_cases=qa_cases,
        priority_on_alternative=priority_on_alternative, ppi_chain_gene_dict=ppi_chain_gene_dict, ppi_gene_edits_dict=ppi_gene_edits_dict,
        config_yaml=config_yaml,atom_level_naa=atom_level_naa
        )
