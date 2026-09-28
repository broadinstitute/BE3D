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

Files are output to ```'[workdir]/hypothesis_qa'```

---

# BE-Clust3D

## Structure and Conservation

### 2. `sequence_structural_features`

**Description:** \
Queries AlphaFold (or reads a user PDB), UniProt, and DSSP to generate a combined sequence-structure feature table. By default the reference sequence (`unipos` 1..N) is taken from `target_chainid` in the structure itself, falling back to UniProt when that chain has gaps, does not start at residue 1, or differs in length from UniProt. The sequence used is saved to `sequence_structure/[structureid]_used_sequence.fasta`.

```python
sequence_structural_features(
    workdir = 'PATH/TO/WORKING/DIRECTORY',   # output directory
    input_gene = 'GENE_NAME',                # gene name (e.g., 'DNMT3A', 'MEN1')
    input_uniprot = 'Q12345',                # UniProt accession ID for input_gene
    structureid = 'UNIQUE-ID',               # identifier used for naming output files
    target_chainid = 'A',                    # chain ID of input_gene in the PDB structure
    # Optional
    radius = 6.0,                            # neighbor count radius in Angstroms
    user_fasta = None,                       # path to user-supplied FASTA file; overrides sequence_source
    user_pdb = None,                         # path to user-supplied PDB file; skips AlphaFold query
    user_dssp = None,                        # path to user-supplied DSSP file; skips DSSP locally
    domains_dict = None,                     # domain annotations e.g. {'ZnF': (1, 100), ...}
    atom_level_naa = False,                  # if True, counts neighbors at atom level rather than residue level
    sequence_source = 'structure',           # 'structure' (from target_chainid, UniProt fallback) or 'uniprot'
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
    alignment_filename = None,               # path to precomputed Clustal-format alignment (primary sequence first); skips UniProt query and MUSCLE entirely
    user_fasta = None,                       # path to local FASTA for input_gene; skips UniProt query; ignored if alignment_filename is given
                                             # (be3d_local.py passes the structure table's [structureid]_used_sequence.fasta here by default)
    alt_user_fasta = None,                   # path to local FASTA for alt_input_gene (e.g., RefSeq-only isoform); ignored if alignment_filename is given
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
    conserv_col = 'mouse_res_pos',                      # residue position column in conserv_dfs to filter on
    v_score_threshold = 3,                              # minimum conservation score (-1, 1, 2, or 3)
    gene_list = False,                                  # if True, processes a list of genes
    mutation_priority = None,                           # most-to-least deleterious category order used to collapse multi-category mut_col values (e.g., 'Silent;Missense;'); None uses mut_col as-is
)
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
    skip_no_coords = True,                 # if True, sets LFC/LFC3D (and randomized) scores to '-' at residues without xyz coordinates (x_coord == '-')
    gene_type = 'Human',                   # species or gene type label used in output naming
    target_gene_chain = 'A',               # chain ID of the target gene in the PDB structure
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
    quantiles = {'NEG 10th v': 0.1, ...},  # percentile label to threshold value mapping
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
    clustering_kwargs = {'n_clusters': None, 'metric': 'euclidean', 'linkage': 'single'}, # None enables distance-threshold clustering
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
    clustering_kwargs = {'n_clusters': None, 'metric': 'euclidean', 'linkage': 'single'}, # AgglomerativeClustering kwargs; match clustering()
    horizontal = False,                    # if True, renders plots with a horizontal layout
    line_subplots_kwargs = {'figsize': (6, 5)},        # kwargs for line plot figure
    dendrogram_subplots_kwargs = {'figsize': (12, 10)}, # kwargs for dendrogram figure
    save_type = 'png',                     # plot format ('png', 'pdf', 'svg', etc.)
    max_distance = None,                   # caps the dendrogram distance axis; cosmetic only, does not affect clustering
)
```

Files are output to ```'[workdir]/cluster_[score_type]/plots'```

---

## Characterization

### 17. `enrichment_test`

**Description:** \
Performs Fisher's exact test to assess enrichment of structural features among significant residues, and plots the results.

```python
enrichment_test(
    df,                                    # DataFrame containing per-residue significance scores and feature annotations
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    hit_columns,                           # significance score columns in df to test for enrichment
    hit_threshold,                         # threshold on hit_columns to define significant residues
    feature_column,                        # structural/functional feature column in df to test
    feature_values,                        # feature values within feature_column to test for enrichment
    # Optional
    confidence_level = 0.95,              # confidence level for odds ratio confidence intervals
)
```

Files are output to ```'[workdir]/characterization'```

---

### 18. `plot_enrichment_test`

**Description:** \
Plots enrichment test results as odds ratios with confidence intervals.

```python
plot_enrichment_test(
    enrichment_results,                    # enrichment results from enrichment_test()
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    hit_value,                             # significance threshold for highlighting odds ratios
    feature_values,                        # feature values tested; output from enrichment_test()
    # Optional
    padding = 0.5,                         # Y-axis padding above and below plotted points
    save_type = 'png',                     # plot format ('png', 'pdf', 'svg', etc.)
    log2 = False,                          # if True, plots odds ratios and confidence intervals on the log2 scale
)
```

Files are output to ```'[workdir]/characterization/plots''```

---

### 19. `lfc_lfc3d_scatter`

**Description:** \
Generates a scatter plot of LFC vs LFC3D scores, color-coded by significance.

```python
lfc_lfc3d_scatter(
    df_input,                              # DataFrame containing per-residue LFC, LFC3D, and significance columns
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    screen_name = 'screen_name_1',         # screen identifier for df_input
    # Optional
    pthr = 0.05,                           # p-value threshold for significance labeling
    save_type = 'png',                     # plot format ('png', 'pdf', 'svg', etc.)
    custom_palette = {'Not LFC3D Hit': 'grey', 'LFC3D Pos Hit': 'blue', ...}, # significance label to color mapping
)
```

Files are output to ```'[workdir]/characterization/plots''```

---

### 20. `pLDDT_RSA_scatter`

**Description:** \
Generates a scatter plot of RSA vs pLDDT scores, scaled by mutation effect weight.

```python
pLDDT_RSA_scatter(
    df_input,                              # DataFrame containing pLDDT, RSA, weight, and directionality columns
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    # Optional
    pLDDT_col = 'bfactor_pLDDT',          # pLDDT confidence score column in df_input
    RSA_col = 'RSA',                       # relative solvent accessibility column in df_input
    size_col = 'LFC3D_wght',              # point size column in df_input
    direction_col = 'direction',           # mutation effect direction column in df_input (e.g., 'NEG', 'POS')
    color_map = {'NEG': 'darkred', 'POS': 'darkblue'}, # direction to color mapping
    save_type = 'png',                     # plot format ('png', 'pdf', 'svg', etc.)
)
```

Files are output to ```'[workdir]/characterization/plots''```

---

### 21. `hits_feature_barplot`

**Description:** \
Generates bar plots of hit counts or fractions across structural feature categories.

```python
hits_feature_barplot(
    df_input,                              # DataFrame containing hit annotations and feature categories
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory
    input_gene = 'GENE_NAME',              # gene name (e.g., 'DNMT3A', 'MEN1')
    category_col,                          # feature category column in df_input to group by
    score_type,                            # score type label used in the plot title; 'LFC' or 'LFC3D'
    values_cols,                           # hit direction columns in df_input to plot
    values_vals,                           # values within values_cols that define a hit
    value_names,                           # display names for each hit category in the legend
    # Optional
    plot_type = 'Count',                   # 'Count' for raw counts or 'Fraction' for proportions
    color_map = {'NEG': 'darkred', 'POS': 'darkblue'}, # hit category to color mapping
    save_type = 'png',                     # plot format ('png', 'pdf', 'svg', etc.)
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
    quantiles = {'NEG 10th v': 0.1, ...},  # percentile label to threshold value mapping
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

# Helpers

## Structure, Preprocessing, and Visualization

### 25. `sequence_structural_features_lite`

**Description:** \
Lightweight variant of `sequence_structural_features()` for PPI partner chains used only as a cross-chain LFC lookup in `calculate_lfc3d()`. Maps the gene's UniProt sequence to its chain in the shared PDB structure and skips domains, DSSP, neighbor counting, and burial degree.

```python
sequence_structural_features_lite(
    workdir = 'PATH/TO/WORKING/DIRECTORY',   # output directory
    input_gene = 'GENE_NAME',                # gene name of the PPI partner
    input_uniprot = 'Q12345',                # UniProt accession ID for input_gene
    structureid = 'UNIQUE-ID',               # identifier used for naming output files
    target_chainid = 'B',                    # chain ID of input_gene in the shared PDB structure
    # Optional
    user_fasta = None,                       # path to user-supplied FASTA file; overrides sequence_source
    user_pdb = None,                         # path to user-supplied PDB/complex file; skips AlphaFold query
    sequence_source = 'structure',           # 'structure' or 'uniprot'; see sequence_structural_features()
)
```

Returns a DataFrame with columns ```['unipos', 'unires', 'chain']```. Files are output to ```'[workdir]/sequence_structure'```

---

### 26. `sanitary_check`

**Description:** \
Reports how many missense edits in each screen map onto the residues in the structure-sequence table, and warns when too many do not (usually a reference-sequence / screen-numbering mismatch).

```python
sanitary_check(
    df_struc,                              # structural feature DataFrame from sequence_structural_features()
    df_missense_list,                      # list of missense DataFrames from parse_be_data(), one per screen
    # Optional
    mute = True,                           # if False, prints mapped / unmapped missense edit counts and the unmapped edits
    screen_names = None,                   # screen identifiers for df_missense_list, used in the report and warnings
    warn_fraction = 0.05,                  # warn when more than this fraction of a screen's missense edits are unmapped
)
```

No files are written.

---

### 27. `check_screen_residues`

**Description:** \
Guardrail comparing every reference residue/position the screens edit (e.g. the `I257` of `I257V`) against the residue the PDB has at that position on `target_chainid`, and against the reference sequence. Any disagreement is printed as a `WARNING` per screen and written to a report; unresolved positions are listed but not counted as discrepancies. be3d_local.py runs it in every mode after `parse_be_data()`, controlled by the yaml key `on_residue_mismatch` (`warn` or `error`), and the notebooks display its reports with `show_residue_check()`.

```python
check_screen_residues(
    workdir = 'PATH/TO/WORKING/DIRECTORY', # output directory parse_be_data() wrote screendata/ into
    input_gene = 'GENE_NAME',              # gene whose own-species screens are checked
    screen_names = ['screen_name_1'],      # screen identifiers, as passed to parse_be_data()
    pdb_processed_file = 'PATH/TO/[structureid]_processed.pdb', # processed PDB from sequence_structural_features()
    target_chainid = 'A',                  # chain ID of input_gene in the PDB structure
    df_struc = pd.DataFrame(),             # residue table with 'unipos' and 'unires' (the reference sequence)
    # Optional
    gene_list = None,                      # per-screen gene symbol; screens not numbered on input_gene are skipped
    mut_categories = ('Missense', 'Silent', 'Nonsense'), # parse_be_data() tables with per-edit refAA/edit_pos
    on_mismatch = 'warn',                  # 'warn' reports and continues; 'error' raises ValueError after writing the report
)
```

Status per edited position: `pdb_mismatch` (PDB residue differs), `reference_mismatch` (unresolved in the PDB and the reference residue differs), `outside_reference` (beyond the reference sequence), `not_in_structure` (reference agrees, PDB has no residue; reported only). Files are output to ```'[workdir]/sequence_check'```

---

### 28. `reduce_mutation_type`

**Description:** \
Collapses a delimiter-joined multi-category mutation type (e.g., `'Silent;Missense;'`, one category per edit in the guide) into a single category, keeping whichever appears first in `priority_order`. Categories not in `priority_order` fall back to the first token; single-category values are returned unchanged. Typically applied per value to `mut_col` before `parse_be_data()`.

```python
df[mut_col] = df[mut_col].apply(lambda x: reduce_mutation_type(
    x,                                     # mutation type string from mut_col
    mut_delimiter = ';',                   # delimiter between per-edit categories
    priority_order = ['Nonsense', 'Splice Site', 'Missense', 'Silent', ...], # categories ordered most to least deleterious
))
```

No files are written.

---

### 29. `g2p_formatted_hit_cluster`

**Description:** \
Gathers LFC, LFC3D, and union hit and cluster labels (and scores) into TSVs formatted for hit cluster visualization on G2P.

```python
g2p_formatted_hit_cluster(
    results_dir = 'PATH/TO/WORKING/DIRECTORY', # output directory containing cluster_LFC, cluster_LFC3D, cluster_union, LFC, LFC3D
    gene_list = ['GENE_NAME'],             # gene name for each screen in screen_names
    screen_names = ['screen_name_1'],      # screen identifiers; paired with gene_list
    # Optional
    lfc_pthr = '05',                       # p-value threshold suffix for LFC hits ('05', '01', or '001')
    lfc3d_pthr = '05',                     # p-value threshold suffix for LFC3D hits ('05', '01', or '001')
    meta_pthr = '001',                     # p-value threshold suffix for meta-aggregate hits ('05', '01', or '001')
    dist = 6,                              # clustering radius in Angstroms whose cluster labels are exported
    function_for_meta = False,             # aggr_func_name from average_split_meta() (e.g., 'SUM', 'mean') to include meta-aggregate hits; False skips
    conservation = False,                  # if True, reads meta-aggregate files under gene name 'Merged'
    input_gene = None,                     # currently unused
)
```

Files are output to ```'[results_dir]/g2p_visualization'```

---

# Notes

- All outputs are saved under specified working directories.
- Many functions allow customizations via optional parameters.
- Functions are modular and can be run screen-by-screen or in batch.
