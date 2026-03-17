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
)

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
Aggregates mutation effects across edit types and sequence positions, combining
structural, conservation, and statistical features per residue.
    
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
Generates per-residue sequence-level plots for a single screen from prioritize_by_sequence() output.

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
Calculates LFC3D scores by aggregating local neighborhood mutation effects.

```python
calculate_lfc3d(
    df_struc, # OUTPUT DF FROM sequence_structural_features()
    df_edits_list, # LIST OF SCREEN DFs FROM prioritize_by_sequence()
    df_rand_list, # LIST OF SCREEN DFs FROM randomize_sequence()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS

    # Optional
    nRandom = 1000, # NUMBER OF RANDOMIZATIONS
    muttype = 'Missense', # ONE OF THE MUTATION TYPES IN df_dict IN prioritize_by_sequence()
    function_aggr = np.mean, # FUNCTION TO APPLY TO ALL LFC SCORES PER POSITION
    function_type = 'mean', # ASSOCIATED FUNCTION NAME
    LFC_only = False, # TRUE SKIPS THE LFC3D CALCULATION, FALSE KEEPS LFC3D CALCULATION
    conserved_only = False, # ONLY CONSIDER 'conserved' RESIDUES FOR SCORES TO AGGREGATE IF TRUE
)
```

Files are output to ```'[workdir]/LFC3D'```

---

## Non Aggregating for Single Screens

### 11. `average_split_score`

**Description:** \
Splits LFC/LFC3D scores into positive and negative components and aggregates randomized scores.

```python
average_split_score(
    df_LFC_LFC3D, # OUTPUT DF FROM calculate_lfc3d()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS

    # Optional
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
)
```

Files are output to ```'[workdir]/[score_type]'```

---

### 12. `bin_score`

**Description:** \
Bins positive and negative LFC3D scores into percentile thresholds.

```python
bin_score(
    df_bidir, # OUTPUT DF FROM average_split_score()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS

    # Optional
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
)
```

Files are output to ```'[workdir]/[score_type]'```

---

### 13. `znorm_score`

**Description:** \
Z-normalizes scores against randomized controls and labels significance.

```python
znorm_score(
    df_bidir, # OUTPUT DF FROM average_split_score()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS

    # Optional
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
    pthrs = [0.05, 0.01, 0.001], # LIST OF P-VALUE CUTOFF TO Z-SCORE ON
)
```

Files are output to ```'[workdir]/[score_type]'```

---

### 14. `average_split_bin_plots`

**Description:** \
Generates histograms, histplots, and scatterplots for positive and negative scores after binning.

```python
average_split_bin_plots(
    df_z, # OUTPUT DF FROM znorm_score()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc

    # Optional
    pthr = 0.05, # P-VALUE CUTOFF TO Z-SCORE ON
    screen_name = '', # '' IF META-AGGREGATE, OR A UNIQUE SCREEN IDENTIFIER FOR NON-AGGREGATE
    func = 'SUM', # '' IF NON-AGGREGATE, OR aggr_func_name FROM znorm_meta() IF META-AGGREGATE
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
    aggregate_dir = 'meta-aggregate', # DIRECTORY TO SAVE TO, SIMILAR TO score_type
    save_type = 'png', # OUTPUT GRAPH SAVE TYPE (ie 'png', 'pdf', 'svg', etc)
)
```

Files are output to ```'[workdir]/[score_type]/plots'```

---

## Clustering

### 15. `clustering`

**Description:** \
Performs spatial clustering of significant residues over a range of distance thresholds.

```python
clustering(
    df_struc, # OUTPUT DF FROM sequence_structural_features()
    df_pvals, # OUTPUT DF FROM znorm_score() OR znorm_meta() OR prioritize_by_sequence()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc

    # Optional
    max_distances = 25, # RANGE OF DISTANCES TO TRY CLUSTERING OVER
    psig_columns = ['SUM_LFC3D_neg_05_psig', 'SUM_LFC3D_pos_05_psig'], # COLUMNS IN df_pvals IDENTIFYING WHAT TO CLUSTER
    pthr_cutoffs = ['p<0.05', 'p<0.05'], # VALUES IN COLUMNS IN df_pvals IDENTIFYING WHAT TO CLUSTER
    screen_name = 'Meta', # IDENTIFIER FOR NAMING OUTPUT FILES
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
    merge_cols = ['unipos', 'chain'],
    clustering_kwargs = {
        'n_clusters': None,
        'metric': 'euclidean',
        'linkage': 'single'
    },
)
```

Files are output to ```'[workdir]/cluster_[score_type]'```

---

### 16. `plot_clustering`

**Description:** \
Plots clustering results including line plots and dendrograms.

```python
plot_clustering(
    df_struc, # OUTPUT DF FROM sequence_structural_features()
    df_pvals, # OUTPUT DF FROM znorm_score() OR znorm_meta() OR prioritize_by_sequence()
    df_pvals_clust, # OUTPUT DF FROM clustering()
    dist, # DISTANCE TO CLUSTER ON
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    distances, # OUTPUT FROM clustering()
    yvalues, # OUTPUT FROM clustering()

    # Optional
    psig_columns = ['SUM_LFC3D_neg_05_psig', 'SUM_LFC3D_pos_05_psig'], # COLUMNS IN df_pvals IDENTIFYING WHAT TO CLUSTER
    names = ['Negative', 'Positive'], # NAMES OF CONDITIONS IN df_pvals IDENTIFYING WHAT TO CLUSTER
    pthr_cutoffs = ['p<0.05', 'p<0.05'], # VALUES IN COLUMNS IN df_pvals IDENTIFYING WHAT TO CLUSTER
    screen_name = 'Meta', # IDENTIFIER FOR NAMING OUTPUT FILES
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
    merge_col = ['unipos', 'chain'],
    clustering_kwargs = {
        'n_clusters': None,
        'metric': 'euclidean',
        'linkage': 'single'
    },
    horizontal = False,
    line_subplots_kwargs = {'figsize': (10, 7)},
    dendogram_subplots_kwargs = {'figsize': (15, 12)},
    save_type = 'png', # OUTPUT GRAPH SAVE TYPE (ie 'png', 'pdf', 'svg', etc)
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
Aggregates scores across multiple screens into a meta score before splitting and averaging.

```python
average_split_meta(
    df_LFC_LFC3D, # OUTPUT DF FROM calculate_lfc3d()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS

    # Optional
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
    nRandom=1000, # NUMBER OF RANDOMIZATIONS
    aggr_func_name = 'SUM', # ASSOCIATED FUNCTION NAME
)
```

Files are output to ```'[workdir]/meta-aggregate'```

---

### 23. `bin_meta`

**Description:** \
Bins meta-aggregated LFC3D scores into percentile thresholds.

```python
bin_meta(
    df_bidir_meta, # OUTPUT DF FROM average_split_meta()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc

    # Optional
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
    aggr_func_name = 'SUM', # ASSOCIATED FUNCTION NAME FROM average_split_meta()
)
```

Files are output to ```'[workdir]/meta-aggregate'```

---

### 24. `znorm_meta`

**Description:** \
Z-normalizes meta-aggregated scores against randomized controls and labels significance.

```python
znorm_meta(
    df_bidir_meta, # OUTPUT DF FROM average_split_meta()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS

    # Optional
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
    pthrs = [0.05, 0.01, 0.001], # LIST OF P-VALUE CUTOFF TO Z-SCORE ON
    aggr_func_name = 'SUM', # ASSOCIATED FUNCTION NAME FROM average_split_meta()
)
```

Files are output to ```'[workdir]/meta-aggregate'```

---

### 14. `average_split_bin_plots`

**Description:** \
Generates histograms, histplots, and scatterplots for positive and negative scores after binning.

```python
average_split_bin_plots(
    df_z, # OUTPUT DF FROM znorm_meta()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc

    # Optional
    pthr=0.05,
    screen_name = '', # '' IF META-AGGREGATE, OR A UNIQUE SCREEN IDENTIFIER FOR NON-AGGREGATE
    func = 'SUM', # '' IF NON-AGGREGATE, OR aggr_func_name FROM znorm_meta() IF META-AGGREGATE
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
    aggregate_dir = 'meta-aggregate', # DIRECTORY TO SAVE TO, SIMILAR TO score_type
    save_type = 'png', # OUTPUT GRAPH SAVE TYPE (ie 'png', 'pdf', 'svg', etc)
)
```

Files are output to ```'[workdir]/meta-aggregate/plots'```

---

# Notes

- All outputs are saved under specified working directories.
- Many functions allow customizations via optional parameters.
- Functions are modular and can be run screen-by-screen or in batch.
