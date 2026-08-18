import os
import sys
import pandas as pd
import numpy as np
import yaml
import glob
import shutil
from datetime import datetime

def load_config(config_yaml):
	with open(config_yaml, "r") as file:
		config = yaml.safe_load(file)

	return config

def preprocess_ppi_partner(
	partner_dir, gene, uniprot, chain,
	screens, screen_dir, screen_names,
	mut_col, val_col, gene_col, edits_col, mut_categories, mut_delimiter,
	user_fasta, user_pdb, mutation_priority=None,
	conservation_run=False, alt_gene_name=None, alt_uniprot_id=None, alt_screen_start=None,
	v_score_threshold=3, muscle_path='muscle', priority_on_alternative=False,
):
	"""
	Lightweight preprocessing for a PPI partner chain: only produces
	{partner_dir}/screendata_sequence/{gene}_{screen_name}_protein_edits.tsv,
	the single file calculate_lfc3d reads as a cross-chain LFC lookup for
	that chain. Skips DSSP/radius/domains/burial, QA, and randomization
	entirely since none of them are used for that value.

	If conservation_run and alt_gene_name are set, screens whose name starts
	with alt_screen_start are treated as coming from a different species/ortholog
	(e.g. mouse Setdb1 screens for a human SETDB1 partner) -- the same
	gene-name switch + residue-map alignment main() applies to its own target
	gene, mirrored here so the partner's screendata is read under the right
	gene symbol and mouse edit positions get mapped onto the partner's own
	(human/PDB-numbered) residues before aggregating. The OUTPUT file is
	always written under the partner's own fixed `gene` name regardless,
	since that's the name ppi_gene_edits_dict/calculate_lfc3d look up by.
	"""
	structureid = f'PDB-{uniprot}' if user_pdb else f'AF-{uniprot}-F1-model_v6'
	df_struc_lite = sequence_structural_features_lite(
		partner_dir, gene, uniprot, structureid,
		target_chainid=chain,
		user_fasta=user_fasta, user_pdb=user_pdb,
	)

	input_dfs = [pd.read_csv(os.path.join(screen_dir, s), sep='\t') for s in screens]
	if mutation_priority:
		for df in input_dfs:
			df[mut_col] = df[mut_col].apply(lambda x: reduce_mutation_type(x, mut_delimiter, mutation_priority))

	# SPECIES-AWARE GENE NAME PER SCREEN, MIRRORING main()'s CONSERVATION SETUP #
	if conservation_run and alt_gene_name:
		_, df_residuemap = conservation(
			partner_dir, gene, alt_gene_name, uniprot, alt_uniprot_id,
			muscle_path=muscle_path,
		)
		conserv_dfs, gene_list = [], []
		for screen_name in screen_names:
			if priority_on_alternative or screen_name.startswith(alt_screen_start):
				conserv_dfs.append(df_residuemap)
				gene_list.append(alt_gene_name)
			else:
				conserv_dfs.append(None)
				gene_list.append(gene)
	else:
		conserv_dfs = [None] * len(screen_names)
		gene_list = [gene] * len(screen_names)

	parse_be_data(
		partner_dir, input_dfs, gene, screen_names,
		mut_col=mut_col, val_col=val_col, gene_col=gene_col, edits_col=edits_col,
		mut_categories=mut_categories, mut_delimiter=mut_delimiter,
		conserv_dfs=conserv_dfs, conserv_col='alternative_res_pos',
		gene_list=gene_list, v_score_threshold=v_score_threshold,
	)

	for screen_name, screen_gene, df_consrv in zip(screen_names, gene_list, conserv_dfs):
		df_control = pd.read_csv(f'{partner_dir}/screendata/{screen_gene}_{screen_name}_No_Mutation.tsv', sep='\t', index_col=0)
		df_dict = {}
		for mut in ['Missense', 'Silent', 'Nonsense']:
			filepath = f'{partner_dir}/screendata/{screen_gene}_{screen_name}_{mut}.tsv'
			if os.path.exists(filepath):
				df_dict[mut] = pd.read_csv(filepath, sep='\t', index_col=0)

		if df_consrv is not None:
			prioritize_by_sequence(
				df_dict, df_struc_lite, df_consrv, df_control,
				partner_dir, gene, screen_name, # ALWAYS the partner's own fixed name (not screen_gene) #
				target_res_pos='original_res_pos', alt_res_pos='alternative_res_pos', alt_res='alternative_res',
			)
		else:
			prioritize_by_sequence(
				df_dict, df_struc_lite, df_consrv, df_control,
				partner_dir, gene, screen_name,
			)

DEFAULT_FUNC_MAP = {'mean': np.mean, 'median': np.median, 'sum': np.sum, 'min': min, 'max': max}

def compute_blind_lfc3d(df_struc, target_chain, ppi_chain_gene_dict, ppi_edits_dict, lfc_colname, aggr_func):
	"""
	For each residue of the target chain, aggregate ONLY its cross-chain
	(listed partner) neighbor LFC values -- self and same-chain neighbors are
	always excluded, since there is no own-chain screen data to contribute in
	the first place. Returns (values, n_neighbors) lists, one entry per row
	of df_struc, in the same order.
	"""
	naa_pos_list = df_struc['Naa_pos'].fillna('-').tolist()
	naa_chain_list = df_struc['Naa_chain'].fillna('-').tolist()

	values, n_neighbors = [], []
	for naa_pos_str, naa_chain_str in zip(naa_pos_list, naa_chain_list):
		if naa_pos_str in ('-', '') or naa_chain_str in ('-', ''):
			values.append('-')
			n_neighbors.append(0)
			continue

		neighbor_vals = []
		for naa_chain, naa_pos in zip(naa_chain_str.split(';'), naa_pos_str.split(';')):
			if naa_chain == target_chain:
				continue # SAME-CHAIN NEIGHBORS EXCLUDED -- BLIND BY DEFINITION, NO SELF-DATA EXISTS #
			gene_identifier = ppi_chain_gene_dict.get(naa_chain)
			if gene_identifier is None:
				continue # CHAIN PRESENT IN THE PDB BUT NOT A LISTED PARTNER -- IGNORE #
			val = ppi_edits_dict.get(gene_identifier, {}).get(int(naa_pos) - 1)
			if val is not None and val != '-' and not (isinstance(val, float) and np.isnan(val)):
				neighbor_vals.append(float(val))

		if neighbor_vals:
			values.append(aggr_func(neighbor_vals))
			n_neighbors.append(len(neighbor_vals))
		else:
			values.append('-')
			n_neighbors.append(0)

	return values, n_neighbors

def split_neg_pos(values):
	"""
	Classify each per-residue combined value as neg (< 0) or pos (> 0), mutually
	exclusive, matching average_split_score's convention exactly (0.0 or missing
	falls into neither). Returns (neg_values, pos_values), same length as values.
	"""
	neg_values, pos_values = [], []
	for v in values:
		fv = float(v) if v != '-' else 0.0
		neg_values.append(fv if fv < 0 else '-')
		pos_values.append(fv if fv > 0 else '-')
	return neg_values, pos_values

def meta_split_neg_pos(per_screen_neg_cols, per_screen_pos_cols, df_out, meta_func):
	"""
	Meta-aggregate neg/pos separately across screens, matching average_split_meta's
	convention: per residue, gather every screen's neg value into one list and
	aggregate it, independently of gathering every screen's pos value into another
	list and aggregating that -- so a residue can end up with BOTH a Meta-neg and
	a Meta-pos value if different screens disagree in sign (e.g. screen A calls it
	negative, screen B calls it positive). Returns (meta_neg, meta_pos) lists.
	"""
	meta_neg, meta_pos = [], []
	for i in range(len(df_out)):
		neg_here = [df_out.at[i, c] for c in per_screen_neg_cols if df_out.at[i, c] != '-']
		pos_here = [df_out.at[i, c] for c in per_screen_pos_cols if df_out.at[i, c] != '-']
		meta_neg.append(meta_func(neg_here) if neg_here else '-')
		meta_pos.append(meta_func(pos_here) if pos_here else '-')
	return meta_neg, meta_pos

def write_bfactor_pdb(input_pdb_path, output_pdb_path, target_chain, value_by_unipos):
	"""
	Write a copy of the input PDB with the B-factor column of every ATOM record
	replaced: for the target chain, the residue's score (looked up by residue
	number == unipos, matching this codebase's PDB-numbering convention); 0.0
	for every other chain/atom and for target-chain residues with no score.
	"""
	from biopandas.pdb import PandasPdb

	ppdb = PandasPdb().read_pdb(input_pdb_path)
	atom_df = ppdb.df['ATOM'].copy()

	def _bfactor(row):
		if row['chain_id'] != target_chain:
			return 0.0
		val = value_by_unipos.get(int(row['residue_number']))
		if val is None or val == '-':
			return 0.0
		return float(val)

	atom_df['b_factor'] = atom_df.apply(_bfactor, axis=1)
	ppdb.df['ATOM'] = atom_df
	ppdb.to_pdb(path=output_pdb_path, records=['ATOM'], gz=False, append_newline=True)

def run_blind_target(
	target_gene, target_uniprot, target_chain, output_dir, partners,
	user_pdb=None, user_fasta=None, user_dssp=None,
	structure_radius=6.0, atom_level_naa=False,
	function_for_lfc='mean', function_for_lfc3d='mean', function_for_meta='mean',
	mut_delimiter_default=';', mutation_priority=None,
):
	"""
	Computes a purely partner-derived ("blind") LFC3D signal at the residues of
	target_chain -- for a chain with no base-editing screen data of its own (a
	structural scaffold chain), or one whose own data is deliberately ignored.

	Unlike main()'s pipeline, this never runs parse_be_data/hypothesis_test/
	prioritize_by_sequence/randomize_data/calculate_lfc3d for the target itself --
	only sequence_structural_features (for real radius-based neighbor topology)
	and preprocess_ppi_partner (for each partner's LFC lookup) are reused.

	IMPORTANT: this does NOT compute a significance-tested "hit". Only
	('cross', ...) neighbor sources ever apply here (no self/same-chain
	contribution is possible), and a partner's cross-chain LFC value is always
	its fixed real value, never re-randomized per permutation -- so a purely
	cross-chain-derived value would be identical across every randomization
	iteration, giving zero variance in the null distribution and making a
	z-score/p-value meaningless. This reports only the raw aggregated LFC3D
	value (split into neg/pos/overall) as a descriptive quantity.
	"""
	os.makedirs(output_dir, exist_ok=True)

	structureid = f'PDB-{target_uniprot}' if user_pdb else f'AF-{target_uniprot}-F1-model_v6'

	# FULL STRUCTURAL FEATURES FOR THE TARGET -- NEEDS THE REAL RADIUS-BASED
	# Naa_chain/Naa_pos NEIGHBOR TOPOLOGY, SO NOT THE _lite VARIANT #
	df_struc = sequence_structural_features(
		output_dir, target_gene, target_uniprot, structureid,
		target_chainid=target_chain,
		radius=structure_radius,
		user_fasta=user_fasta, user_pdb=user_pdb, user_dssp=user_dssp,
		atom_level_naa=atom_level_naa,
	)

	lfc_colname = f'{function_for_lfc}_Missense_LFC'
	ppi_chain_gene_dict = {}
	ppi_edits_dict_per_screen = {} # gene_identifier -> {screen_name: {row_idx: lfc_value}}
	all_screen_names = []

	for partner in partners:
		gene = partner['gene']
		uniprot = partner['uniprot']
		chain = partner['chain']
		screens = [s.strip() for s in partner['screens'].split(',')]
		screen_names = [s.split('.')[0] for s in screens]
		screen_dir = partner['screen_dir']
		mut_categories = partner['mut_categories']
		mut_delimiter = partner.get('mut_delimiter', mut_delimiter_default)

		gene_identifier = f'{gene}_chain_{chain}'
		partner_dir = os.path.join(output_dir, 'ppi_partners', gene_identifier)
		os.makedirs(partner_dir, exist_ok=True)

		preprocess_ppi_partner(
			partner_dir, gene, uniprot, chain,
			screens, screen_dir, screen_names,
			partner.get('mut_col', 'Mutation_type'), partner.get('val_col', 'sgRNA_score'),
			partner.get('gene_col', 'Gene'), partner.get('edits_col', 'Mutation_list'),
			mut_categories, mut_delimiter,
			partner.get('user_fasta'), partner.get('user_pdb', user_pdb),
			mutation_priority=partner.get('mutation_priority', mutation_priority),
			conservation_run=partner.get('conservation_run', False),
			alt_gene_name=partner.get('alt_gene_name'),
			alt_uniprot_id=partner.get('alt_uniprot_id'),
			alt_screen_start=partner.get('alt_screen_start'),
			v_score_threshold=partner.get('v_score_threshold', 3),
			muscle_path=partner.get('muscle_path', 'muscle'),
			priority_on_alternative=partner.get('priority_on_alternative', False),
		)

		ppi_chain_gene_dict[chain] = gene_identifier
		ppi_edits_dict_per_screen[gene_identifier] = {}
		for screen_name in screen_names:
			if screen_name not in all_screen_names:
				all_screen_names.append(screen_name)
			df_partner_edits = pd.read_csv(
				f'{partner_dir}/screendata_sequence/{gene}_{screen_name}_protein_edits.tsv', sep='\t')
			ppi_edits_dict_per_screen[gene_identifier][screen_name] = df_partner_edits[lfc_colname].to_dict()

	aggr_func = DEFAULT_FUNC_MAP[function_for_lfc3d]
	meta_func = DEFAULT_FUNC_MAP[function_for_meta]

	df_out = df_struc[['unipos', 'unires', 'chain']].copy()
	neg_cols, pos_cols = [], []
	for screen_name in all_screen_names:
		ppi_edits_dict = {gid: d.get(screen_name, {}) for gid, d in ppi_edits_dict_per_screen.items()}
		values, n_neighbors = compute_blind_lfc3d(
			df_struc, target_chain, ppi_chain_gene_dict, ppi_edits_dict, lfc_colname, aggr_func)
		neg_values, pos_values = split_neg_pos(values)

		col = f'{screen_name}_LFC3D_blind'
		df_out[col] = values
		df_out[f'{col}_nneighbors'] = n_neighbors
		df_out[f'{col}_neg'] = neg_values
		df_out[f'{col}_pos'] = pos_values
		neg_cols.append(f'{col}_neg')
		pos_cols.append(f'{col}_pos')

	# WHICH LEVEL TO REPORT/MAP ONTO THE PDB: THE SINGLE SCREEN'S VALUES IF THERE'S
	# ONLY ONE, OTHERWISE THE META-AGGREGATE ACROSS SCREENS (NOT MUTUALLY EXCLUSIVE --
	# A RESIDUE CAN HAVE BOTH A META-NEG AND META-POS IF SCREENS DISAGREE IN SIGN) #
	if len(all_screen_names) > 1:
		meta_neg, meta_pos = meta_split_neg_pos(neg_cols, pos_cols, df_out, meta_func)
		df_out['Meta_LFC3D_blind_neg'] = meta_neg
		df_out['Meta_LFC3D_blind_pos'] = meta_pos
		report_neg, report_pos = meta_neg, meta_pos
	else:
		report_neg, report_pos = df_out[neg_cols[0]].tolist(), df_out[pos_cols[0]].tolist()

	# OVERALL: abs(neg) + abs(pos) -- HIGHER (MORE "SIGNIFICANT") WHEN A RESIDUE HAS
	# EVIDENCE IN BOTH DIRECTIONS AT ONCE, NOT JUST ONE #
	def _overall(n, p):
		n_abs = abs(float(n)) if n != '-' else 0.0
		p_abs = abs(float(p)) if p != '-' else 0.0
		total = n_abs + p_abs
		return total if total > 0 else '-'
	report_overall = [_overall(n, p) for n, p in zip(report_neg, report_pos)]
	overall_colname = 'Meta_LFC3D_blind_overall' if len(all_screen_names) > 1 else f'{all_screen_names[0]}_LFC3D_blind_overall'
	df_out[overall_colname] = report_overall

	out_path = os.path.join(output_dir, f'{target_gene}_{target_chain}_blind_LFC3D.tsv')
	df_out.to_csv(out_path, sep='\t', index=False)
	print(f'Wrote {out_path}')

	# MAP NEG/POS/OVERALL ONTO THE INPUT PDB'S B-FACTOR COLUMN FOR THE TARGET CHAIN #
	if user_pdb:
		unipos_list = df_out['unipos'].tolist()
		for label, report_vals in [('neg', report_neg), ('pos', report_pos), ('overall', report_overall)]:
			value_by_unipos = {int(pos): val for pos, val in zip(unipos_list, report_vals)}
			pdb_out_path = os.path.join(output_dir, f'{target_gene}_{target_chain}_blind_LFC3D_{label}.pdb')
			write_bfactor_pdb(user_pdb, pdb_out_path, target_chain, value_by_unipos)
			print(f'Wrote {pdb_out_path}')

	return df_out

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
	muscle_path=kwargs['muscle_path']
	mutation_priority=kwargs['mutation_priority']

	if user_pdb:
		structureid = f'PDB-{input_uniprot}'
	else:
		structureid = f'AF-{input_uniprot}-F1-model_v6'
		
	single_screen_pthr_str = str(single_screen_pthr)
	multi_screen_pthr_str = str(multi_screen_pthr)
	pdb_file = os.path.join(output_dir, "sequence_structure", f"{structureid}_processed.pdb")

	def find_union(input, pthr_str): #TODO: this should be moved to helper
		if input[0] == f'p<{pthr_str}' or input[1] == f'p<{pthr_str}':
			return f'p<{pthr_str}'
		elif input[0] == f'p>={pthr_str}' or input[1] == f'p>={pthr_str}':
			return f'p>={pthr_str}'
		else:
			return '-'

	os.makedirs(output_dir, exist_ok=True)
	shutil.copy2(config_yaml, os.path.join(output_dir,os.path.basename(config_yaml)))
	print(f'All results will be saved in the following directory: {output_dir}')

	screens = [screen.strip() for screen in screens.split(',')]
	screen_names = [s.split('.')[0] for s in screens]
	input_dfs = [pd.read_csv(os.path.join(screen_dir,s), sep='\t') for s in screens]

	# COLLAPSE PER-EDIT MUTATION CATEGORY LISTS (e.g. 'Silent;Missense;') INTO A SINGLE CATEGORY PER GUIDE,
	# BEFORE hypothesis_test / parse_be_data DO EXACT-MATCH FILTERING AGAINST mut_categories #
	if mutation_priority:
		for df in input_dfs:
			df[mut_col] = df[mut_col].apply(lambda x: reduce_mutation_type(x, mut_delimiter, mutation_priority))

	df_residuemap = pd.DataFrame()
	conserv_dfs = list()
	gene_list = list()

	if conservation_run:
		_, df_residuemap = conservation(output_dir,
					input_gene, alt_gene_name,
					input_uniprot, alt_uniprot_id,
					muscle_path=muscle_path, 
		)

		if priority_on_alternative:
			for screen_name in screen_names:
				conserv_dfs.append(df_residuemap)
				gene_list.append(alt_gene_name)
		else:
			for screen_name in screen_names:
				if alt_gene_name and screen_name.startswith(alt_screen_start):
					conserv_dfs.append(df_residuemap)
					gene_list.append(alt_gene_name)
				else:
					conserv_dfs.append(None)
					gene_list.append(input_gene)
	else:
		for screen_name in screen_names:
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
		h2_ks_test_pd = pd.read_csv(f'{output_dir}/hypothesis_qa/KolmogorovSmirnov_hypothesis2.tsv',sep='\t')
		h2_ks_test_pd = h2_ks_test_pd.replace(-999,None)
		white_screen_list = h2_ks_test_pd[(h2_ks_test_pd[f"p_{'_'.join(qa_cases)}_vs_{'_'.join(qa_controls)}"]<0.05)&(h2_ks_test_pd['gene_name'].isin(gene_list))]['screenid'].to_list()
		print(f'original screen size: {len(screen_names)}, screen white list size: {len(white_screen_list)}, QA-passed screen size: {len(list(set(screen_names).intersection(white_screen_list)))}')
		screen_name_to_file = dict(zip([s.split('.')[0] for s in screens], screens))

		original_screen_names = screen_names
		passed_screen_names = list(set(original_screen_names).intersection(white_screen_list))
		failed_screen_names = [s for s in original_screen_names if s not in passed_screen_names]
		qa_status = 'PASSED' if passed_screen_names else 'FAILED'

		# TAG THE OUTPUT FOLDER WITH QA STATUS SO A FAILED GENE IS DISTINGUISHABLE WITHOUT RE-RUNNING #
		with open(os.path.join(output_dir, 'QA_STATUS.txt'), 'w') as f:
			f.write(f'qa_status: {qa_status}\n')
			f.write(f'n_screens_total: {len(original_screen_names)}\n')
			f.write(f'n_screens_passed: {len(passed_screen_names)}\n')
			f.write(f'screens_passed: {",".join(sorted(passed_screen_names))}\n')
			f.write(f'screens_failed: {",".join(sorted(failed_screen_names))}\n')

		screen_names = passed_screen_names
		input_dfs = [pd.read_csv(os.path.join(screen_dir,screen_name_to_file[s]), sep='\t') for s in screen_names]
		conserv_dfs = list()
		gene_list = list() # where we need only have genes passed QA, so OVERWRITE the previous raw gene list.
		for screen_name in screen_names:
			if alt_gene_name and screen_name.startswith(alt_screen_start):
				conserv_dfs.append(df_residuemap)
				gene_list.append(alt_gene_name)
			else:
				conserv_dfs.append(None)
				gene_list.append(input_gene)

		if qa_status == 'FAILED':
			print(f'QA FAILED: 0/{len(original_screen_names)} screens passed QA for gene {input_gene}. See {output_dir}/QA_STATUS.txt')
			sys.exit(1)

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
	
	df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')
	sanitary_check(df_struc, df_missense_list)
	
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
			ppi_gene_edits_dict=ppi_gene_edits_dict,
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
		
		# CHARACTERIZATION
		df_domains = pd.read_csv(f'{output_dir}/sequence_structure/{input_gene}_{input_uniprot}_domains.tsv', sep='\t')
		df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')
		df_struc['pLDDT_dis'] = df_struc['pLDDT_dis'].isin(['confident','high']).map(
			{True: 'True', False: 'False'}
		)

		if df_domains['Domain'].nunique() > 1:
			df_domains['Domain'] = df_domains['Domain'].notna().map(
				{True: 'True', False: 'False'}
		)
		cutoff_str = str(single_screen_pthr_str).split('.')[1]

		for screen_name in screen_names:
			for score_type, _df_ in zip(['LFC', 'LFC3D'],[df_lfc,df_lfc3d]):
				for each_feature in ['Domain','pLDDT_dis']:
					input_df = pd.concat([_df_, df_domains['Domain'], df_struc['pLDDT_dis']], axis=1)
					if len(input_df[each_feature].unique()) > 1:
						hit_columns = [f'{screen_name}_{score_type}_neg_{cutoff_str}_p', f'{screen_name}_{score_type}_pos_{cutoff_str}_p']
						input_df[hit_columns] = input_df[hit_columns].replace('-', np.nan).astype(float)

						enrichment_test_results = enrichment_test(
							input_df,
							workdir=output_dir,
							input_gene=input_gene,
							hit_columns=hit_columns,
							hit_threshold=single_screen_pthr,
							feature_column=each_feature,
							feature_values=['True'],
							confidence_level=0.95,
						)

						plot_enrichment_test(
							enrichment_results=enrichment_test_results,
							workdir=output_dir,
							input_gene=input_gene,
							hit_value=single_screen_pthr,
							log2=True,
							feature_values=['True','False']
						)

						os.rename(f'{output_dir}/characterization/{input_gene}_enrichment_test.pickle',
								f'{output_dir}/characterization/{input_gene}_enrichment_test_{screen_name}_{each_feature}_{score_type}_{cutoff_str}.pickle')
						os.rename(f'{output_dir}/characterization/plots/{input_gene}_enrichment_test.png',
								f'{output_dir}/characterization/plots/{input_gene}_enrichment_test_{screen_name}_{each_feature}_{score_type}_{cutoff_str}.png')

						# BARPLOTS #
						colnames = [f'{screen_name}_{score_type}_neg_{cutoff_str}_psig', f'{screen_name}_{score_type}_pos_{cutoff_str}_psig']
						input_df = pd.concat([df_struc, _df_[colnames], df_domains['Domain']], axis=1)

						hits_feature_barplot(
							input_df,
							workdir=output_dir,
							input_gene=input_gene,
							category_col=each_feature,
							score_type=score_type,
							values_cols=colnames, values_vals=[f'p<0.{cutoff_str}', f'p<0.{cutoff_str}'], value_names=['NEG', 'POS'],
							plot_type='Count',
						)

						os.rename(f'{output_dir}/characterization/plots/{input_gene}_Count_{each_feature}_barplot.png',
								f'{output_dir}/characterization/plots/{input_gene}_{screen_name}_Count_{each_feature}_barplot_{score_type}_{cutoff_str}.png')


		# SCATTERPLOT #
		# if conservation_run: # TODO: maybe later
		df_lfc_dis = pd.read_csv(f"{output_dir}/LFC/{gene_type}_{gene}_LFC_dis_wght.tsv", sep='\t')
		df_lfc3d_dis = pd.read_csv(f"{output_dir}/LFC3D/{gene_type}_{gene}_LFC3D_dis_wght.tsv", sep='\t')
		df_lfc = pd.read_csv(f"{output_dir}/LFC/{gene_type}_{gene}_NonAggr_LFC.tsv", sep='\t')
		df_lfc3d = pd.read_csv(f"{output_dir}/LFC3D/{gene_type}_{gene}_NonAggr_LFC3D.tsv", sep='\t')

		for screen_name in screen_names:
			df_dis_input = pd.DataFrame()
			df_dis_input['unipos'] = df_lfc_dis['unipos'].copy()
			df_dis_input[f'{screen_name}_LFC'] = df_lfc_dis[f'{screen_name}_LFC'].copy()
			df_dis_input[f'{screen_name}_LFC3D'] = df_lfc3d_dis[f'{screen_name}_LFC3D'].copy()
			df_dis_input[f'{screen_name}_LFC3D_dis'] = df_lfc3d_dis[f'{screen_name}_LFC3D_dis'].copy()
			df_dis_input[f'LFC3D_neg_psig'] = df_lfc3d[f'{screen_name}_LFC3D_neg_{cutoff_str}_psig'].copy()
			df_dis_input[f'LFC3D_pos_psig'] = df_lfc3d[f'{screen_name}_LFC3D_pos_{cutoff_str}_psig'].copy()
		
			lfc_lfc3d_scatter(
				df_input=df_dis_input,
				workdir=output_dir,
				input_gene=input_gene, screen_name=screen_name,
				pthr=single_screen_pthr,
			)
			os.rename(f'{output_dir}/characterization/plots/{input_gene}_LFC_LFC3D_scatter.png',
					f'{output_dir}/characterization/plots/{input_gene}_{screen_name}_LFC_LFC3D_scatter_{cutoff_str}.png')


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

		df_struc['pLDDT_dis'] = df_struc['pLDDT_dis'].isin(['confident','high']).map(
			{True: 'True', False: 'False'}
		)

		if df_domains['Domain'].nunique() > 1:
			df_domains['Domain'] = df_domains['Domain'].notna().map(
				{True: 'True', False: 'False'}
		)
			
		for cutoff in [multi_screen_pthr]:
			for score_type in ['LFC', 'LFC3D']:
				# ENRICHMENT TEST #
				cutoff_str = str(cutoff).split('.')[1]
				if conservation_run:
					df_meta = pd.read_csv(f'{output_dir}/meta-aggregate/Merged_MetaAggr_{score_type}.tsv', sep='\t')
				else:
					df_meta = pd.read_csv(f'{output_dir}/meta-aggregate/{input_gene}_MetaAggr_{score_type}.tsv', sep='\t')
				
				for each_feature in ['pLDDT_dis','Domain']:
					input_df = pd.concat([df_meta, df_domains['Domain'], df_struc['pLDDT_dis']], axis=1)

					if len(input_df[each_feature].unique()) > 1:
						hit_columns = [f'{function_for_meta}_{score_type}_neg_{cutoff_str}_p', f'{function_for_meta}_{score_type}_pos_{cutoff_str}_p']
						input_df[hit_columns] = input_df[hit_columns].replace('-', np.nan).astype(float)

						results = enrichment_test(
							input_df,
							workdir=output_dir,
							input_gene=input_gene,
							hit_columns=hit_columns,
							hit_threshold=cutoff,
							feature_column=each_feature,
							feature_values=['True'],
							confidence_level=0.95,
						)

						plot_enrichment_test(
							enrichment_results=results,
							workdir=output_dir,
							input_gene=input_gene,
							hit_value=cutoff,
							log2=True,
							feature_values=['True','False']

						)
						os.rename(f'{output_dir}/characterization/{input_gene}_enrichment_test.pickle',
								f'{output_dir}/characterization/{input_gene}_enrichment_test_{each_feature}_{score_type}_{cutoff_str}.pickle')
						os.rename(f'{output_dir}/characterization/plots/{input_gene}_enrichment_test.png',
								f'{output_dir}/characterization/plots/{input_gene}_enrichment_test_Meta_{each_feature}_{score_type}_{cutoff_str}.png')

						# BARPLOTS #
						colnames = [f'{function_for_meta}_{score_type}_neg_{cutoff_str}_psig', f'{function_for_meta}_{score_type}_pos_{cutoff_str}_psig']
						input_df = pd.concat([df_struc, df_meta[colnames], df_domains['Domain']], axis=1)

						hits_feature_barplot(
							input_df,
							workdir=output_dir,
							input_gene=input_gene,
							category_col=each_feature,
							score_type=score_type,
							values_cols=colnames, values_vals=[f'p<0.{cutoff_str}', f'p<0.{cutoff_str}'], value_names=['NEG', 'POS'],
							plot_type='Count',
						)

						os.rename(f'{output_dir}/characterization/plots/{input_gene}_Count_{each_feature}_barplot.png',
								f'{output_dir}/characterization/plots/{input_gene}_Count_Meta_{each_feature}_barplot_{score_type}_{cutoff_str}.png')

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
				

			for screen_name in ['Meta']:
				df_dis_input = pd.DataFrame()
				df_dis_input['unipos'] = df_lfc_dis['unipos']				
				df_dis_input[f'Meta_LFC'] = df_lfc_dis[f'{function_for_meta}_LFC']
				df_dis_input[f'Meta_LFC3D'] = df_lfc3d_dis[f'{function_for_meta}_LFC3D']
				df_dis_input[f'Meta_LFC3D_dis'] = df_lfc3d_dis[f'{function_for_meta}_LFC3D_dis']

				df_dis_input[f'Meta_LFC3D_neg_psig'] = df_lfc3d[f'{function_for_meta}_LFC3D_neg_{cutoff_str}_psig']
				df_dis_input[f'Meta_LFC3D_pos_psig'] = df_lfc3d[f'{function_for_meta}_LFC3D_pos_{cutoff_str}_psig']

				lfc_lfc3d_scatter(
					df_input=df_dis_input,
					workdir=output_dir,
					input_gene=input_gene, screen_name='Meta',
					pthr=cutoff,
				)
				os.rename(f'{output_dir}/characterization/plots/{input_gene}_LFC_LFC3D_scatter.png',
						f'{output_dir}/characterization/plots/{input_gene}_{screen_name}_LFC_LFC3D_scatter_{cutoff_str}.png')

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

	# REACHING HERE MEANS THE FULL PIPELINE RAN WITHOUT ERROR; QA-FAILED OR qa_only RUNS EXIT EARLIER AND NEVER WRITE THIS #
	with open(os.path.join(output_dir, 'RUN_COMPLETED.txt'), 'w') as f:
		f.write(f'status: SUCCESS\n')
		f.write(f'finished_at: {datetime.now().isoformat()}\n')
		f.write(f'input_gene: {input_gene}\n')
		f.write(f'screens: {",".join(screen_names)}\n')

if __name__ == '__main__':
	config_yaml = sys.argv[1]
	config = load_config(config_yaml)
	beclust3d_path = config['beclust3d_path']
	sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), beclust3d_path)))

	from beclust3d.lfc3d.structure import sequence_structural_features, sequence_structural_features_lite
	from beclust3d.lfc3d.preprocess_data import parse_be_data, sanitary_check
	from beclust3d.lfc3d.preprocess_data_helpers import reduce_mutation_type
	from beclust3d.lfc3d.preprocess_data_plot import plot_rawdata
	from beclust3d.qa.hypothesis_tests import hypothesis_test
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
	output_dir= config['output_dir']
	mode = config.get('mode', 'monomer') # monomer: single target gene (PPI neighbors, if any, are pre-built dirs in ppi_gene_edits_dict)
	                                      # complex: run the full pipeline once per gene in input_gene, lite-preprocessing every other chain as its PPI neighbor
	                                      # blind_target: input_gene has no screen data of its own (or it's ignored); LFC3D
	                                      # is computed purely from partners (config['partners']), never calls main() at all

	# THE FIELDS BELOW ARE ONLY EVER USED BY monomer/complex (main()'s full pipeline) -- READ DEFENSIVELY
	# VIA .get() SO A blind_target YAML (WHICH NEEDS NONE OF THIS) DOESN'T HAVE TO SUPPLY THEM #
	screen_dir = config.get('screen_dir')
	screens = config.get('screens')

	conservation_cfg = config.get('conservation') or {}
	conservation_run = conservation_cfg.get('run', False)
	v_score_threshold = conservation_cfg.get('v_score_threshold', 3)
	alt_gene_name = conservation_cfg.get('alt_gene_name')
	alt_uniprot_id = conservation_cfg.get('alt_uniprot_id')
	alt_screen_start = conservation_cfg.get('alt_screen_start')

	database_cfg = config.get('database') or {}
	mut_col = database_cfg.get('mut_col')
	val_col = database_cfg.get('val_col')
	gene_col = database_cfg.get('gene_col')
	edits_col = database_cfg.get('edits_col')
	gRNA_col = database_cfg.get('gRNA_col')
	mut_delimiter = database_cfg.get('mut_delimiter', ';')

	mutation_category_cfg = config.get('mutation_category') or {}
	mut_categories = list()
	mut_categories.extend(mutation_category_cfg.get('nonsense', []))
	mut_categories.extend(mutation_category_cfg.get('splice', []))
	mut_categories.extend(mutation_category_cfg.get('missense', []))
	mut_categories.extend(mutation_category_cfg.get('silent', []))
	mut_categories.extend(mutation_category_cfg.get('no_mutation', []))
	mut_categories.extend(mutation_category_cfg.get('intron', []))

	# OPTIONAL
	user_fasta = config.get('user_fasta')
	user_pdb = config.get('user_pdb')
	user_dssp = config.get('user_dssp')
	function_for_lfc = config.get('function_for_lfc', 'mean')
	function_for_lfc3d = config.get('function_for_lfc3d', 'mean')
	function_for_meta = config.get('function_for_meta', 'mean')
	nRandom = config.get('nRandom', 500)
	pthr_cfg = config.get('pthr') or {}
	single_screen_pthr = pthr_cfg.get('single_screen', 0.05)
	multi_screen_pthr = pthr_cfg.get('multi_screen', 0.05)
	structure_radius = config.get('structure_radius', 6.0)
	clustering_radius = config.get('clustering_radius', 6.0)
	qa_cfg = config.get('qa') or {}
	qa_passed_only = qa_cfg.get('qa_passed_only', False)
	qa_only = qa_cfg.get('qa_only', False)
	qa_controls = qa_cfg.get('controls', [])
	qa_cases = qa_cfg.get('cases', [])
	priority_on_alternative = config.get('priority_on_alternative', False)
	atom_level_naa = config.get('atom_level_naa', False)
	muscle_path = config.get('muscle_path', 'muscle')
	mutation_priority = config.get('mutation_priority') # optional; most-to-least-deleterious order for collapsing
	                                                      # delimiter-joined multi-category mut_col values (e.g. 'Silent;Missense;')

	# KWARGS SHARED ACROSS EVERY main() CALL, REGARDLESS OF MODE #
	common_kwargs = dict(
		screen_dir=screen_dir, screens=screens, mut_col=mut_col, val_col=val_col, gene_col=gene_col, edits_col=edits_col,
		gRNA_col=gRNA_col, user_fasta=user_fasta, user_pdb=user_pdb, user_dssp=user_dssp, nRandom=nRandom,
		single_screen_pthr=single_screen_pthr, multi_screen_pthr=multi_screen_pthr,
		structure_radius=structure_radius, clustering_radius=clustering_radius, function_for_lfc=function_for_lfc, function_for_lfc3d=function_for_lfc3d,
		mut_categories=mut_categories, mut_delimiter=mut_delimiter, conservation_run=conservation_run, alt_gene_name=alt_gene_name, alt_uniprot_id=alt_uniprot_id,
		alt_screen_start=alt_screen_start, v_score_threshold=v_score_threshold,
		function_for_meta=function_for_meta, qa_passed_only=qa_passed_only, qa_only=qa_only, qa_controls=qa_controls, qa_cases=qa_cases,
		priority_on_alternative=priority_on_alternative, config_yaml=config_yaml, atom_level_naa=atom_level_naa, muscle_path=muscle_path,
		mutation_priority=mutation_priority,
	)

	if mode == 'monomer':
		# PPI-RELATED CONFIG (ppi_chain_gene_dict, ppi_gene_edits_dict, partner_uniprot) IS IGNORED IN MONOMER MODE,
		# EVEN IF PRESENT IN THE YAML -- USE mode: complex (a single-gene list works fine) FOR ANY PPI-AWARE RUN #
		# NOTE: calculate_lfc3d distinguishes monomer vs PPI mode via `isinstance(ppi_chain_gene_dict, dict)`,
		# NOT truthiness -- an empty dict {} still counts as PPI mode and crashes on any cross-chain neighbor
		# found in a multi-chain PDB, so this must be None, not {} #
		ppi_chain_gene_dict, ppi_gene_edits_dict = None, {}

		main(input_gene=input_gene, input_uniprot=input_uniprot, input_chain=input_chain, output_dir=output_dir,
			ppi_chain_gene_dict=ppi_chain_gene_dict, ppi_gene_edits_dict=ppi_gene_edits_dict, **common_kwargs)

	elif mode == 'complex':
		gene_names = [g.strip() for g in input_gene.split(',')]
		uniprot_list = [u.strip() for u in input_uniprot.split(',')]
		chain_list = [c.strip() for c in input_chain.split(',')]
		assert len(gene_names) == len(uniprot_list) == len(chain_list), \
			'input_gene, input_uniprot, input_chain must list the same number of comma-separated entries in complex mode'
		gene_to_uniprot = dict(zip(gene_names, uniprot_list))
		# GENES THAT ARE PURE PPI NEIGHBORS (NEVER RUN AS A TARGET) SUPPLY THEIR UNIPROT HERE INSTEAD OF input_gene/input_uniprot #
		gene_to_uniprot.update(config.get('partner_uniprot') or {})

		ppi_chain_gene_dict_full = config['ppi_chain_gene_dict'] # chain -> gene, covers every chain in the complex
		screens_list = [s.strip() for s in screens.split(',')]
		screen_names = [s.split('.')[0] for s in screens_list]

		for target_gene, target_uniprot, target_chain in zip(gene_names, uniprot_list, chain_list):
			gene_output_dir = os.path.join(output_dir, target_gene)
			partner_items = [(ch, g) for ch, g in ppi_chain_gene_dict_full.items() if ch != target_chain]

			# LITE-PREPROCESS EVERY OTHER CHAIN AS A CROSS-CHAIN LFC LOOKUP FOR THIS TARGET GENE'S RUN #
			ppi_chain_gene_dict, ppi_gene_edits_dict = {}, {}
			for ch, g in partner_items:
				assert g in gene_to_uniprot, f'gene {g} (chain {ch}) has no matching entry in input_gene/input_uniprot or partner_uniprot'
				gene_identifier = f'{g}_chain_{ch}'
				partner_dir = os.path.join(gene_output_dir, 'ppi_partners', gene_identifier)
				os.makedirs(partner_dir, exist_ok=True)

				preprocess_ppi_partner(
					partner_dir, g, gene_to_uniprot[g], ch,
					screens_list, screen_dir, screen_names,
					mut_col, val_col, gene_col, edits_col, mut_categories, mut_delimiter,
					user_fasta, user_pdb, mutation_priority=mutation_priority,
					conservation_run=conservation_run, alt_gene_name=alt_gene_name, alt_uniprot_id=alt_uniprot_id,
					alt_screen_start=alt_screen_start, v_score_threshold=v_score_threshold,
					muscle_path=muscle_path, priority_on_alternative=priority_on_alternative,
				)
				ppi_chain_gene_dict[ch] = gene_identifier
				ppi_gene_edits_dict[gene_identifier] = partner_dir

			main(input_gene=target_gene, input_uniprot=target_uniprot, input_chain=target_chain, output_dir=gene_output_dir,
				ppi_chain_gene_dict=ppi_chain_gene_dict, ppi_gene_edits_dict=ppi_gene_edits_dict, **common_kwargs)

	elif mode == 'blind_target':
		# input_gene HAS NO SCREEN DATA OF ITS OWN (OR ITS OWN DATA IS DELIBERATELY IGNORED) --
		# NEVER CALLS main(): NO parse_be_data/hypothesis_test/prioritize_by_sequence/calculate_lfc3d FOR THE
		# TARGET ITSELF, ONLY sequence_structural_features (REAL, NOT _lite) + preprocess_ppi_partner PER PARTNER #
		run_blind_target(
			target_gene=input_gene, target_uniprot=input_uniprot, target_chain=input_chain,
			output_dir=output_dir, partners=config['partners'],
			user_pdb=user_pdb, user_fasta=user_fasta, user_dssp=user_dssp,
			structure_radius=structure_radius, atom_level_naa=atom_level_naa,
			function_for_lfc=function_for_lfc, function_for_lfc3d=function_for_lfc3d, function_for_meta=function_for_meta,
			mut_delimiter_default=mut_delimiter, mutation_priority=mutation_priority,
		)

	else:
		raise ValueError(f"Unknown mode '{mode}'; expected 'monomer', 'complex', or 'blind_target'")
