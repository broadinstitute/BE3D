# Functions Overview

This document summarizes the major functions used in the analysis pipeline. Each section provides a brief description and an example function call with all parameters shown clearly.

---

# BE-QA

## Quality Assessment Hypothesis Testing

### 1. `hypothesis_test`

**Description:** \
Runs Mann-Whitney U and Kolmogorov-Smirnov tests for Hypothesis 1 (case vs. control within a single screen) and Hypothesis 2 (case in one screen vs. controls pooled across all screens).

```python
hypothesis_test(
    workdir = 'PATH/TO/WORKING/DIRECTORY',   # output directory
    input_dfs = [pd.DataFrame()],            # one DataFrame per screen
    screen_names = ['screen_name_1'],        # screen identifier for each DataFrame in input_dfs
    cases = ['Nonsense', 'Splice Site'],     # mutation categories treated as the case group
    controls = ['Silent', 'No Mutation'],    # mutation categories treated as the control group
    # Optional
    comp_name = 'CaseVsControl',             # label used in output filenames and plot titles
    mut_col = 'Mutation category',           # mutation category column in input_dfs
    val_col = 'logFC',                       # numeric measurement column in input_dfs
    gene_col = 'Target Gene Symbol',         # gene identifier column in input_dfs
    save_type = 'png',                       # plot format ('png', 'pdf', 'svg', etc.)
)
```

Files are output to ```'[workdir]/hypothesis_qc'```

---

# BE-Clust3D

## Structure and Conservation

### 2. `sequence_structural_features`

**Description:** \
Queries UniProt, AlphaFold, and DSSP to generate a combined sequence-structure feature table.

```python
sequence_structural_features(
    workdir = 'PATH/TO/WORKING/DIRECTORY',   # output directory
    input_gene = 'GENE_NAME',                # gene name (e.g., 'DNMT3A', 'MEN1')
    input_uniprot = 'Q12345',                # UniProt accession ID for input_gene
    structureid = 'UNIQUE-ID',               # identifier used for naming output files
    target_chainid = 'A',                    # chain ID of input_gene in the PDB structure
    # Optional
    radius = 6.0,                            # neighbor count radius in Angstroms
    user_fasta = None,                       # path to user-supplied FASTA file; skips UniProt query
    user_pdb = None,                         # path to user-supplied PDB file; skips AlphaFold query
    user_dssp = None,                        # path to user-supplied DSSP file; skips DSSP locally
    domains_dict = None,                     # domain annotations e.g. {'ZnF': (1, 100), ...}
    atom_level_naa = False,                  # if True, counts neighbors at atom level rather than residue level
)
```

Files are output to ```'[workdir]/sequence_structure'```

---

### 3. `conservation`

**Description:** \
Aligns two protein sequences and generates per-residue conservation scores.

```python
conservation(
    workdir = 'PATH/TO/WORKING/DIRECTORY',   # output directory
    input_gene = 'GENE_NAME',                # gene name (e.g., 'DNMT3A', 'MEN1')
    alt_input_gene = 'ALT_GENE_NAME',        # alternate gene name (e.g., mouse ortholog or isoform)
    input_uniprot = 'Q12345',                # UniProt accession ID for input_gene
    alt_input_uniprot = 'P12345',            # UniProt accession ID for alt_input_gene
    # Optional
    alignment_filename = None,               # path to precomputed alignment file; skips MUSCLE entirely
    mode = 'run',                            # 'run' uses local MUSCLE; 'query' uses remote MUSCLE API
    title = None,                            # job title for remote API request; required if mode='query'
    email = None,                            # email for remote API request; required if mode='query'
    wait_time = 30,                          # seconds between API re-polls; only used if mode='query'
    muscle_path = 'muscle',                  # path to local MUSCLE executable; only used if mode='run'
    cons_dict = {'*': ('conserved', 3), ...} # alignment symbol to conservation label/score mapping
)
```

Files are output to ```'[workdir]/conservation'```

---

## Raw Data to LFC

### 4. `parse_be_data`

**Description:** \
Parses raw base editing screen data into per-mutation-type DataFrames for each screen.

```python
parse_be_data(
    workdir = 'PATH/TO/WORKING/DIRECTORY',              # output directory
    input_dfs = [pd.DataFrame()],                       # one DataFrame per screen
    input_gene = 'GENE_NAME',                           # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_names = ['screen_name_1'],                   # screen identifier for each DataFrame in input_dfs
    # Optional
    mut_col = 'Mutation category',                      # mutation category column in input_dfs
    val_col = 'logFC',                                  # numeric measurement column in input_dfs
    gene_col = 'Target Gene Symbol',                    # gene identifier column in input_dfs
    edits_col = 'Amino Acid Edits',                     # amino acid edits column in input_dfs (e.g., 'M1V,Q2Q')
    mut_categories = ["Nonsense", "Splice Site", ...],  # mutation categories to extract from mut_col
    mut_delimiter = ',',                                # delimiter used within edits_col
    conserv_dfs = [],                                   # conservation DataFrames from conservation()
    conserv_col = 'alt_res_pos',                        # residue position column in conserv_dfs to filter on
    v_score_threshold = 3,                              # minimum conservation score (-1, 1, 2, or 3)
    gene_list = False,                                  # if True, processes a list of genes
```

Files are output to ```'[workdir]/screendata'```

---

### 5. `plot_rawdata`

**Description:** \
Parses raw screen data and generates summary plots per mutation category for each screen.

```python
plot_rawdata(
    workdir = 'PATH/TO/WORKING/DIRECTORY',                              # output directory
    input_dfs = [pd.DataFrame()],                                       # one DataFrame per screen
    screen_names = ['screen_name_1'],                                   # screen identifier for each DataFrame in input_dfs
    # Optional
    mut_col = 'Mutation category',                                      # mutation category column in input_dfs
    val_col = 'logFC',                                                  # numeric measurement column in input_dfs
    gene_col = 'Target Gene Symbol',                                    # gene identifier column in input_dfs
    mut_categories = ["Nonsense", "Splice Site", "Missense", ...],      # mutation categories to plot from mut_col
    save_type = 'png',                                                  # plot format ('png', 'pdf', 'svg', etc.)
)
```

Files are output to ```'[workdir]/screendata/plots'```

---

### 6. `randomize_data`

**Description:** \
Randomizes mutation scores from a parsed screen DataFrame to create a baseline distribution.

```python
randomize_data(
    df_missense,                      # parsed mutation DataFrame from parse_be_data()
    workdir = 'PATH/TO/WORKING/DIRECTORY',   # output directory
    input_gene = 'GENE_NAME',         # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_name = 'screen_name_1',    # screen identifier for df_missense
    # Optional
    nRandom = 1000,                   # number of randomizations to perform
    val_colname = 'LFC',              # numeric measurement column in df_missense
    muttype = 'Missense',             # mutation category of df_missense
    seed = False,                     # if True, uses a fixed seed for reproducibility
)
```

Files are output to ```'[workdir]/screendata_rand'```

---

## LFC by Sequence to LFC3D

### 7. `prioritize_by_sequence`

**Description:** \
Takes in results across multiple edit types for a screen, and aggregates the edits for each residue with sequence and conservation information. 
    
```python
prioritize_by_sequence(
    df_dict,                               # dict of parsed DataFrames from parse_be_data() e.g. {'Missense': pd.DataFrame(), ...}
    df_struc,                              # structural feature DataFrame from sequence_structural_features()
    df_consrv,                             # conservation DataFrame from conservation(), or None
    df_control,                            # control/no-mutation LFC DataFrame, or None
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_name = 'screen_name_1',         # screen identifier for DataFrames in df_dict
    # Optional
    pthr = 0.05,                           # p-value threshold for significance labeling
    functions = [statistics.mean, min, max, sum],        # aggregation functions applied per residue
    function_names = ['mean', 'min', 'max', 'sum'],      # names corresponding to each function
    target_res_pos = 'human_res_pos',      # primary sequence residue position column in df_consrv
    alt_res_pos = 'mouse_res_pos',         # alternate sequence residue position column in df_consrv
    alt_res = 'mouse_res',                 # alternate sequence residue identity column in df_consrv
)
```

Files are output to ```'[workdir]/screendata_sequence'```

---

### 8. `randomize_sequence`

**Description:** \
Randomizes per-residue scores weighted by structural and conservation features to create a baseline distribution for significance testing.

```python
randomize_sequence(
    df_missense,                           # per-residue LFC DataFrame from prioritize_by_sequence()
    df_rand,                               # randomized per-guide LFC DataFrame from randomize_data()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_name = 'screen_name_1',         # screen identifier for df_missense
    # Optional
    nRandom = 1000,                        # number of randomizations to perform
    conservation = False,                  # if True, only aggregates conserved residues; match prioritize_by_sequence()
    muttype = 'Missense',                  # mutation category to randomize; match df_dict in prioritize_by_sequence()
    function_name = 'mean',                # aggregation function name; match function_names in prioritize_by_sequence()
    target_pos = 'unipos',                 # primary sequence residue position column in df_missense
    target_res = None,                     # primary sequence residue identity column in df_missense, or None
)
```

Files are output to ```'[workdir]/screendata_sequence_rand'```

---

### 9. `plot_screendata_sequence`

**Description:** \
Parse raw data and create plots for each input screen.

```python
plot_screendata_sequence(
    df_protein,                            # per-residue feature DataFrame from prioritize_by_sequence()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_name = 'screen_name_1',         # screen identifier for df_protein
    # Optional
    function_name = 'mean',                # aggregation function name to plot; match function_names in prioritize_by_sequence()
    muttype = 'Missense',                  # mutation category to plot; match df_dict in prioritize_by_sequence()
    save_type = 'png',                     # plot format ('png', 'pdf', 'svg', etc.)
)
```

Files are output to ```'[workdir]/screendata_sequence'```

---

### 10. `calculate_lfc3d`

**Description:** \
Calculates LFC3D scores by aggregating local structural neighborhood mutation effects across screens.

```python
calculate_lfc3d(
    df_struc,                              # structural feature DataFrame from sequence_structural_features()
    df_edits_list,                         # list of per-residue LFC DataFrames from prioritize_by_sequence()
    df_rand_list,                          # list of randomized per-residue LFC DataFrames from randomize_sequence()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_names = ['screen_name_1'],      # screen identifier for each DataFrame in df_edits_list and df_rand_list
    # Optional
    nRandom = 1000,                        # number of randomizations to perform
    muttype = 'Missense',                  # mutation category; match df_dict in prioritize_by_sequence()
    function_type_lfc = 'mean',            # per-residue LFC aggregation function name; match function_names in prioritize_by_sequence()
    function_type_lfc3d = 'mean',          # LFC3D aggregation function name; must be a key in func_map
    LFC_only = False,                      # if True, skips LFC3D calculation and outputs LFC scores only
    conserved_only = False,                # if True, aggregates conserved residues only; match randomize_sequence()
    gene_type = 'Human',                   # species or gene type label used in output naming
    target_gene_chain = 'A',              # chain ID of the target gene in the PDB structure
    ppi_chain_gene_dict = {},              # interacting gene to chain ID mapping (e.g., {'GENE1': 'B', ...})
    ppi_gene_edits_dict = {},              # interacting gene to edits dict mapping (e.g., {'GENE1': edits_dict, ...})
    func_map = {'mean': np.mean, ...},     # function name to callable mapping for LFC3D aggregation
)
```

Files are output to ```'[workdir]/LFC3D'```

---

## Non Aggregating for Single Screens

### 11. `average_split_score`

**Description:** \
Splits LFC or LFC3D scores into positive and negative components and aggregates randomized scores per screen.
    
```python
average_split_score(
    df_LFC_LFC3D,                          # per-residue mutation score DataFrame from calculate_lfc3d()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_names = ['screen_name_1'],      # screen identifiers for df_LFC_LFC3D
    # Optional
    score_type = 'LFC3D',                  # score type to split; 'LFC' or 'LFC3D'
    gene_type = 'Human',                   # species or gene type label used in output naming
)
```

Files are output to ```'[workdir]/[score_type]'```

---

### 12. `bin_score`

**Description:** \
Bins positive and negative LFC or LFC3D scores into percentile thresholds per screen.

```python
bin_score(
    df_bidir,                              # split positive/negative score DataFrame from average_split_score()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_names = ['screen_name_1'],      # screen identifiers for df_bidir
    # Optional
    score_type = 'LFC3D',                  # score type to bin; 'LFC' or 'LFC3D'
    gene_type = 'Human',                   # species or gene type label used in output naming
)
```

Files are output to ```'[workdir]/[score_type]'```

---

### 13. `znorm_score`

**Description:** \
Z-normalizes LFC or LFC3D scores against randomized control distributions and assigns significance labels at multiple p-value thresholds.

```python
znorm_score(
    df_bidir,                              # split positive/negative score DataFrame from average_split_score()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_names = ['screen_name_1'],      # screen identifiers for df_bidir
    # Optional
    score_type = 'LFC3D',                  # score type to normalize; 'LFC' or 'LFC3D'
    pthrs = [0.05, 0.01, 0.001],           # p-value thresholds for significance labeling
    gene_type = 'Human',                   # species or gene type label used in output naming
)
```

Files are output to ```'[workdir]/[score_type]'```

---

### 14. `average_split_bin_plots`

**Description:** \
Generates histograms and scatterplots for positive and negative scores with binning and significance thresholds.

```python
average_split_bin_plots(
    df_z,                                  # z-score DataFrame from znorm_score() or znorm_meta()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    # Optional
    pthr = 0.05,                           # p-value threshold for significance labeling
    screen_name = '',                      # '' for meta-aggregate, or screen identifier for per-screen output
    func = 'SUM',                          # '' for per-screen, or aggr_func_name from znorm_meta() for meta-aggregate
    score_type = 'LFC3D',                  # score type to plot; 'LFC' or 'LFC3D'
    aggregate_dir = 'meta-aggregate',      # subdirectory to save plots into
    save_type = 'png',                     # plot format ('png', 'pdf', 'svg', etc.)
)
```

Files are output to ```'[workdir]/[score_type]/plots'```

---

## Clustering

### 15. `clustering`

**Description:** \
Performs spatial agglomerative clustering of significant residues over a range of distance thresholds.

```python
clustering(
    df_struc,                              # structural feature DataFrame from sequence_structural_features()
    df_pvals,                              # significance label DataFrame from znorm_score(), znorm_meta(), or prioritize_by_sequence()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    # Optional
    max_distances = 25,                    # maximum clustering radius in Angstroms
    psig_columns = ['SUM_LFC3D_neg_05_psig', 'SUM_LFC3D_pos_05_psig'], # significance columns in df_pvals to cluster on
    pthr_cutoffs = ['p<0.05', 'p<0.05'],   # significance values in psig_columns to include in clustering
    screen_name = 'Meta',                  # screen identifier for output filenames
    score_type = 'LFC3D',                  # score type to cluster; 'LFC' or 'LFC3D'
    merge_cols = ['unipos', 'chain'],      # columns used to merge clustering results
    clustering_kwargs = {'n_clusters': None, 'metric': 'euclidean', 'linkage': 'single'} # None enables distance-threshold clustering
    atom_level = False,                    # if True, clusters at atom level rather than residue level
)
```

Files are output to ```'[workdir]/cluster_[score_type]'```

---

### 16. `plot_clustering`

**Description:** \
Generates line plots and dendrograms for clustering results at a specified distance threshold.

```python
plot_clustering(
    df_struc,                              # structural feature DataFrame from sequence_structural_features()
    df_pvals,                              # significance label DataFrame from znorm_score(), znorm_meta(), or prioritize_by_sequence()
    df_pvals_clust,                        # cluster label DataFrame from clustering()
    dist,                                  # clustering radius in Angstroms to plot results for
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    distances,                             # distances output from clustering()
    yvalues,                               # cluster counts output from clustering()
    # Optional
    psig_columns = ['SUM_LFC3D_neg_05_psig', 'SUM_LFC3D_pos_05_psig'], # significance columns in df_pvals; match clustering()
    names = ['Negative', 'Positive'],      # display names corresponding to each psig_column
    pthr_cutoffs = ['p<0.05', 'p<0.05'],   # significance values in psig_columns; match clustering()
    screen_name = 'Meta',                  # screen identifier for output filenames
    score_type = 'LFC3D',                  # score type to plot; 'LFC' or 'LFC3D'
    merge_col = ['unipos', 'chain'],       # columns used to merge clustering results
    clustering_kwargs = {'n_clusters': None, 'metric': 'euclidean', 'linkage': 'single'} # AgglomerativeClustering kwargs; match clustering()
    horizontal = False,                    # if True, renders plots with a horizontal layout
    line_subplots_kwargs = {'figsize': (10, 7)},       # kwargs for line plot figure
    dendrogram_subplots_kwargs = {'figsize': (15, 12)}, # kwargs for dendrogram figure
    save_type = 'png',                     # plot format ('png', 'pdf', 'svg', etc.)
)
```

Files are output to ```'[workdir]/cluster_[score_type]/plots'```

---

## Characterization

### 17. `enrichment_test`

**Description:** \
Performs enrichment tests (e.g., Fisher's exact test) for structural features.

```python
enrichment_test(
    df,
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    hit_columns,
    hit_threshold,
    feature_column,
    feature_values,

    # Optional
    confidence_level = 0.95,
)
```

Files are output to ```'[workdir]/characterization'```

---

### 18. `plot_enrichment_test`

**Description:** \
Plots enrichment test results as odds ratios with confidence intervals.

```python
plot_enrichment_test(
    enrichment_results,
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    hit_value,
    feature_values,

    # Optional
    padding = 0.5,
    save_type = 'png', # OUTPUT GRAPH SAVE TYPE (ie 'png', 'pdf', 'svg', etc)
)
```

Files are output to ```'[workdir]/characterization/plots''```

---

### 19. `lfc_lfc3d_scatter`

**Description:** \
Generates LFC vs LFC3D scatter plots colored by hit significance.

```python
lfc_lfc3d_scatter(
    df_input,
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_name = 'screen_name_1', # UNIQUE SCREEN IDENTIFIER FOR df_input

    # Optional
    pthr = 0.05, # P-VALUE CUTOFF TO Z-SCORE ON
    save_type = 'png', # OUTPUT GRAPH SAVE TYPE (ie 'png', 'pdf', 'svg', etc)
)
```

Files are output to ```'[workdir]/characterization/plots''```

---

### 20. `pLDDT_RSA_scatter`

**Description:** \
Generates scatter plot of RSA vs pLDDT scores, scaled by mutation weight.

```python
pLDDT_RSA_scatter(
    df_input,
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc

    # Optional
    pLDDT_col = 'bfactor_pLDDT',
    RSA_col = 'RSA',
    size_col = 'LFC3D_wght',
    direction_col = 'direction',
    color_map = {'NEG': 'darkred', 'POS': 'darkblue'},
    save_type = 'png', # OUTPUT GRAPH SAVE TYPE (ie 'png', 'pdf', 'svg', etc)
)
```

Files are output to ```'[workdir]/characterization/plots''```

---

### 21. `hits_feature_barplot`

**Description:** \
Generates bar plots of hit counts (or fractions) across different structural categories.

```python
hits_feature_barplot(
    df_input,
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    category_col,
    values_cols,
    values_vals,
    value_names,

    # Optional
    plot_type = 'Count',
    colors = ['darkred', 'darkblue'],
    save_type = 'png', # OUTPUT GRAPH SAVE TYPE (ie 'png', 'pdf', 'svg', etc)
)
```

Files are output to ```'[workdir]/characterization/plots'```

---

# BE-MetaClust3D

## Meta-Aggregation for Multiple Screens

Replaces Steps 10-13 under Non Aggregating for Single Screens

### 22. `average_split_meta`

**Description:** \
Aggregates scores across screens into a meta score, then splits into positive and negative components and averages randomized scores.

```python
average_split_meta(
    df_LFC_LFC3D,                          # per-residue mutation score DataFrame from calculate_lfc3d()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_names = ['screen_name_1'],      # screen identifiers for df_LFC_LFC3D
    # Optional
    score_type = 'LFC3D',                  # score type to aggregate and split; 'LFC' or 'LFC3D'
    nRandom = 500,                         # number of randomizations to perform
    aggr_func_name = 'SUM',                # aggregation function name; must be a key in func_map
    func_map = {'SUM': np.sum, ...},       # function name to callable mapping for aggregation
)
```

Files are output to ```'[workdir]/meta-aggregate'```

---

### 23. `bin_meta`

**Description:** \
Bins positive and negative meta-aggregated LFC or LFC3D scores into percentile thresholds.

```python
bin_meta(
    df_bidir_meta,                         # meta-aggregated split score DataFrame from average_split_meta()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    # Optional
    score_type = 'LFC3D',                  # score type to bin; 'LFC' or 'LFC3D'
    aggr_func_name = 'SUM',                # aggregation function name used in average_split_meta()
    quantiles = {'NEG_10p_v': 0.1, ...},  # percentile label to threshold value mapping
)
```

Files are output to ```'[workdir]/meta-aggregate'```

---

### 24. `znorm_meta`

**Description:** \
Z-normalizes meta-aggregated LFC or LFC3D scores against randomized control distributions and assigns significance labels at multiple p-value thresholds.
    
```python
znorm_meta(
    df_bidir_meta,                         # meta-aggregated split score DataFrame from average_split_meta()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_names = ['screen_name_1'],      # screen identifiers for df_bidir_meta
    # Optional
    score_type = 'LFC3D',                  # score type to normalize; 'LFC' or 'LFC3D'
    pthrs = [0.05, 0.01, 0.001],           # p-value thresholds for significance labeling
    aggr_func_name = 'SUM',                # aggregation function name used in average_split_meta()
)
```

Files are output to ```'[workdir]/meta-aggregate'```

---

### 14. `average_split_bin_plots`

**Description:** \
Generates histograms and scatterplots for positive and negative scores with binning and significance thresholds.

```python
average_split_bin_plots(
    df_z,                                  # z-score DataFrame from znorm_score() or znorm_meta()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    # Optional
    pthr = 0.05,                           # p-value threshold for significance labeling
    screen_name = '',                      # '' for meta-aggregate, or screen identifier for per-screen output
    func = 'SUM',                          # '' for per-screen, or aggr_func_name from znorm_meta() for meta-aggregate
    score_type = 'LFC3D',                  # score type to plot; 'LFC' or 'LFC3D'
    aggregate_dir = 'meta-aggregate',      # subdirectory to save plots into
    save_type = 'png',                     # plot format ('png', 'pdf', 'svg', etc.)
)
```

Files are output to ```'[workdir]/meta-aggregate/plots'```

---

# Notes

- All outputs are saved under specified working directories.
- Many functions allow customizations via optional parameters.
- Functions are modular and can be run screen-by-screen or in batch.
