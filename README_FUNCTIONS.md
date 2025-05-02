# Functions Overview

This document summarizes the major functions used in the analysis pipeline. Each section provides a brief description and an example function call with all parameters shown clearly.

---

# Structure and Conservation

## 1. `sequence_structural_features`

**Description:**
Queries UniProt, AlphaFold, DSSP, and domain features.
Generates a combined sequence-structure feature table.

```python
sequence_structural_features(
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    input_uniprot = 'Q12345',
    structureid = 'UNIQUE-ID',

    # Optional
    chains = ['A'], # CHAIN OF input_gene
    radius = 6.0, # RADIUS OF LFC3D CALCULATION
    user_uniprot = None, # RELATIVE/PATH/FROM/WORKING/DIR/TO/USER/PROVIDED/.FASTA/FILE
    user_pdb = None, # RELATIVE/PATH/FROM/WORKING/DIR/TO/USER/PROVIDED/.PDB/FILE
    user_dssp = None, # RELATIVE/PATH/FROM/WORKING/DIR/TO/USER/PROVIDED/.DSSP/FILE
    domains_dict = None, # DICT OF DOMAINS ie {'ZnF':(1,100), ...}
)
```

Files are output to ```'[workdir]/sequence_structure'```

---

## 2. `conservation`

**Description:**
Generates dataframes of sequence conservation by aligning sequences across species or isoforms.

```python
conservation(
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    alt_input_gene = 'ALT_GENE_NAME',
    input_uniprot = 'Q12345',
    alt_input_uniprot = 'P12345',

    # Optional
    alignment_filename = None, # RELATIVE/PATH/FROM/WORKING/DIR/TO/USER/PROVIDED/.ALIGN/FILE
    mode = 'run', # 'run' USES LOCAL PACKAGES WHILE 'query' USES THE MUSCLE API
                # WHICH METHOD WORKS MAY VARY DEPENDING ON THE MACHINE
    title = None, # ONLY REQUIRED IF MODE='QUERY', TITLE OF JOB
    email = None, # ONLY REQUIRED IF MODE='QUERY'
    wait_time = 30, # ONLY REQUIRED IF MODE='QUERY'
)
```

Files are output to ```'[workdir]/conservation'```

---

# Raw Data to LFC

## 3. `hypothesis_test`

**Description:**
Conducts hypothesis 1 (within screen) and hypothesis 2 (across screens) statistical tests.

```python
hypothesis_test(
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_dfs = [pd.DataFrame()], # LIST OF DFs, ONE FOR EACH SCREEN
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS
    cases = ['Nonsense', 'Splice Site'], # DELETERIOUS TYPES OF MUTATIONS UNDER mut_col
    controls = ['Silent', 'No Mutation'], # DELETERIOUS TYPES OF MUTATIONS UNDER mut_col

    # Optional
    comp_name = 'CaseVsControl', # FOR NAMING HYPOTHESIS DF COLUMNS
    mut_col = 'Mutation category', # MUTATION CATEGORY COLUMN IN input_dfs
    val_col = 'logFC', # SCORE COLUMN IN input_dfs
    gene_col = 'Target Gene Symbol', # GENE NAME COLUMN IN input_dfs
    save_type = 'png', # OUTPUT GRAPH SAVE TYPE (ie 'png', 'pdf', 'svg', etc)
)
```

Files are output to ```'[workdir]/hypothesis_qc'```

---

## 4. `parse_be_data`

**Description:**
Parses raw base editing screen data into DataFrames for each mutation type.

```python
parse_be_data(
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_dfs = [pd.DataFrame()], # LIST OF DFs, ONE FOR EACH SCREEN
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS

    # Optional
    mut_col = 'Mutation category', # MUTATION CATEGORY COLUMN IN input_dfs (ie 'Missense', 'Silent', etc)
    val_col = 'logFC', # SCORE COLUMN IN input_dfs
    gene_col = 'Target Gene Symbol', # GENE NAME COLUMN IN input_dfs
    edits_col = 'Amino Acid Edits', # EDITS CATEGORY COLUMN IN input_dfs (ie 'M1V,Q2Q', etc)
    mut_categories = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"], # CATEGORIES IN mut_col
    mut_delimiter = ',', # DELIMITER IN edits_col
    conserv_dfs = [], # CONSERVATION DFs FROM conservation()
    conserv_col = 'alt_res_pos', # ALTERNATE RES POS ALIGNED WITH MAIN SEQUENCE
)
```

Files are output to ```'[workdir]/screendata'```

---

## 5. `plot_rawdata`

**Description:**
Parses raw screen data and generates summary plots per mutation category.

```python
plot_rawdata(
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_dfs = [pd.DataFrame()], # LIST OF DFs, ONE FOR EACH SCREEN
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS

    # Optional
    mut_col = 'Mutation category', # MUTATION CATEGORY COLUMN IN input_dfs (ie 'Missense', 'Silent', etc)
    val_col = 'logFC', # SCORE COLUMN IN input_dfs
    gene_col = 'Target Gene Symbol', # GENE NAME COLUMN IN input_dfs
    mut_categories = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"], # CATEGORIES IN mut_col
    save_type = 'png', # OUTPUT GRAPH SAVE TYPE (ie 'png', 'pdf', 'svg', etc)
)
```

Files are output to ```'[workdir]/screendata/plots'```

---

## 6. `randomize_data`

**Description:**
Randomizes missense mutation scores to create a baseline distribution.

```python
randomize_data(
    df_missense, # OUTPUT DF FROM parse_be_data()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_name = 'screen_name_1', # UNIQUE SCREEN IDENTIFIER FOR df_missense

    # Optional
    nRandom = 1000, # NUMBER OF RANDOMIZATIONS
    val_colname = 'LFC', # SCORE COLUMN IN df_missense 
    muttype = 'Missense', # MUTATION TYPE OF df_missense
    seed = False, # SEED CREATES REPRODUCIBLE RANDOMIZATIONS
)
```

Files are output to ```'[workdir]/screendata_rand'```

---

# LFC by Sequence to LFC3D

## 7. `prioritize_by_sequence`

**Description:**
Aggregates mutation effects across edit types, sequence positions, and conservation features.

```python
prioritize_by_sequence(
    df_dict, # DICTIONARY OF OUTPUT DF FROM prioritize_by_sequence() {'Missense':pd.DataFrame(), ...}
    df_struc, # OUTPUT DF FROM sequence_structural_features()
    df_consrv, # OUTPUT DF FROM conservation() OR 'NONE'
    df_control, # OUTPUT DF FROM prioritize_by_sequence() OR SEPARATE CONTROL POPULATION DF
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_name = 'screen_name_1', # UNIQUE SCREEN IDENTIFIER FOR ITEMS IN df_dict

    # Optional
    pthr = 0.05, # P-VALUE CUTOFF TO Z-SCORE ON
    functions = [statistics.mean, min, max], # LIST OF FUNCTIONS TO APPLY TO ALL SCORES PER POSITION
    function_names = ['mean', 'min', 'max'], # LIST OF ASSOCIATED FUNCTION NAMES
    target_res_pos = 'original_res_pos', # ORIGINAL RES POS OF THE MAIN SEQUENCE
    target_res = 'original_res', # ORIGINAL RES OF THE MAIN SEQUENCE
)
```

Files are output to ```'[workdir]/screendata_sequence'```

---

## 8. `randomize_sequence`

**Description:**
Randomizes scores based on structural sequence and conservation information.

```python
randomize_sequence(
    df_missense, # OUTPUT DF FROM prioritize_by_sequence()
    df_rand, # OUTPUT DF FROM randomize_data()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_name = 'screen_name_1', # UNIQUE SCREEN IDENTIFIER FOR df_missense

    # Optional
    nRandom = 1000, # NUMBER OF RANDOMIZATIONS
    conservation = False, # TRUE IF df_consrv IS NONE IN prioritize_by_sequence()
    muttype = 'Missense', # ONE OF THE MUTATION TYPES IN df_dict IN prioritize_by_sequence()
    function_name = 'mean', # ONE OF THE FUNCTIONS APPLIED IN prioritize_by_sequence()
    target_pos = 'unipos', # ORIGINAL RES POS OF THE MAIN SEQUENCE
    target_res = None, # ORIGINAL RES OF THE MAIN SEQUENCE
)
```

Files are output to ```'[workdir]/screendata_sequence_rand'```

---

## 9. `calculate_lfc3d`

**Description:**
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

# Non Aggregating for Single Screens

## 10. `average_split_score`

**Description:**
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

## 11. `bin_score`

**Description:**
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

## 12. `znorm_score`

**Description:**
Z-normalizes scores against randomized controls and labels significance.

```python
znorm_score(
    df_bidir, # OUTPUT DF FROM average_split_score()
    neg_stats_list, # OUTPUT NEG STATS DF FROM bin_score()
    pos_stats_list, # OUTPUT POS STATS DF FROM bin_score()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc
    screen_names = ['screen_name_1'], # LIST OF UNIQUE SCREEN IDENTIFIERS

    # Optional
    pthrs = [0.05, 0.01, 0.001], # LIST OF P-VALUE CUTOFF TO Z-SCORE ON
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
)
```

Files are output to ```'[workdir]/[score_type]'```

---

## 16. `average_split_bin_plots`

**Description:**
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

# Meta Aggregating for Multiple Screens

## 13. `average_split_meta`

**Description:**
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
    aggr_func = np.sum, # FUNCTIONS TO APPLY TO ALL SCORES PER POSITION TO AGGREGATE
    aggr_func_name = 'SUM', # ASSOCIATED FUNCTION NAME
)
```

Files are output to ```'[workdir]/meta-aggregate'```

---

## 14. `bin_meta`

**Description:**
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

## 15. `znorm_meta`

**Description:**
Z-normalizes meta-aggregated scores against randomized controls and labels significance.

```python
znorm_meta(
    df_bidir_meta, # OUTPUT DF FROM average_split_meta()
    neg_stats, # OUTPUT NEG STATS DF FROM bin_meta()
    pos_stats, # OUTPUT POS STATS DF FROM bin_meta()
    workdir = 'PATH/TO/WORKING/DIRECTORY',
    input_gene = 'GENE_NAME', # DNMT3A, MEN1, etc

    # Optional
    score_type = 'LFC3D', # 'LFC' OR 'LFC3D'
    pthrs = [0.05, 0.01, 0.001], # LIST OF P-VALUE CUTOFF TO Z-SCORE ON
    aggr_func_name = 'SUM', # ASSOCIATED FUNCTION NAME FROM average_split_meta()
)
```

Files are output to ```'[workdir]/meta-aggregate'```

---

## 16. `average_split_bin_plots`

**Description:**
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

# Clustering

## 17. `clustering`

**Description:**
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

## 18. `plot_clustering`

**Description:**
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

# Characterization

## 19. `enrichment_test`

**Description:**
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

## 20. `plot_enrichment_test`

**Description:**
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

## 21. `lfc_lfc3d_scatter`

**Description:**
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

## 22. `pLDDT_RSA_scatter`

**Description:**
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

## 23. `hits_feature_barplot`

**Description:**
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

# Notes

- All outputs are saved under specified working directories.
- Many functions allow customizations via optional parameters.
- Functions are modular and can be run screen-by-screen or in batch.
