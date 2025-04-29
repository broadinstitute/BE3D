# Functions Overview

This document summarizes the major functions used in the analysis pipeline. Each section provides a brief description and an example function call.

---

# Structure and Conservation

## 1. `sequence_structural_features`

**Description:**\
Queries UniProt, AlphaFold, DSSP, and domain features. Generates a combined sequence-structure feature table.

```python
sequence_structural_features(workdir, input_gene, input_uniprot, structureid)
```

---

## 2. `conservation`

**Description:**\
Generates dataframes of sequence conservation by aligning sequences across species or isoforms.

```python
conservation(workdir, input_gene, alt_input_gene, input_uniprot, alt_input_uniprot)
```

---

## 3. `hypothesis_test`

**Description:**\
Conducts hypothesis 1 (within screen) and hypothesis 2 (across screens) statistical tests.

```python
hypothesis_test(workdir, input_dfs, screen_names, cases, controls, comp_name)
```

---

# Raw Data to LFC

## 4. `parse_be_data`

**Description:**\
Parses raw base editing screen data into DataFrames for each mutation type.

```python
parse_be_data(workdir, input_dfs, input_gene, screen_names)
```

---

## 5. `plot_rawdata`

**Description:**\
Parses raw screen data and generates summary plots per mutation category.

```python
plot_rawdata(workdir, input_dfs, screen_names)
```

---

## 6. `randomize_data`

**Description:**\
Randomizes missense mutation scores to create a baseline distribution.

```python
randomize_data(df_missense, workdir, input_gene, screen_name)
```

---

# LFC by Sequence to LFC3D

## 7. `prioritize_by_sequence`

**Description:**\
Aggregates mutation effects across edit types, sequence positions, and conservation features.

```python
prioritize_by_sequence(df_struc, df_consrv, df_control, workdir, input_gene, screen_name, df_dict)
```

---

## 8. `randomize_sequence`

**Description:**\
Randomizes scores based on structural sequence and conservation information.

```python
randomize_sequence(df_missense, df_rand, workdir, input_gene, screen_name)
```

---

## 9. `calculate_lfc3d`

**Description:**\
Calculates LFC3D scores by aggregating local neighborhood mutation effects.

```python
calculate_lfc3d(df_str_cons, df_edits_list, df_rand_list, workdir, input_gene, screen_names)
```

---

# Non Aggregating for Single Screens

## 10. `average_split_score`

**Description:**\
Splits LFC/LFC3D scores into positive and negative components and aggregates randomized scores.

```python
average_split_score(df_LFC_LFC3D, workdir, input_gene, screen_names)
```

---

## 11. `bin_score`

**Description:**\
Bins positive and negative LFC3D scores into percentile thresholds.

```python
bin_score(df_bidir, workdir, input_gene, screen_names)
```

---

## 12. `znorm_score`

**Description:**\
Z-normalizes scores against randomized controls and labels significance.

```python
znorm_score(df_bidir, neg_stats_list, pos_stats_list, workdir, input_gene, screen_names)
```

---

## 16. `average_split_bin_plots`

**Description:**\
Generates histograms, histplots, and scatterplots for positive and negative scores after binning.

```python
average_split_bin_plots(df_Z, workdir, input_gene)
```

# Meta Aggregating for Multiple Screens

## 13. `average_split_meta`

**Description:**\
Aggregates scores across multiple screens into a meta score before splitting and averaging.

```python
average_split_meta(df_LFC_LFC3D, workdir, input_gene, screen_names)
```

---

## 14. `bin_meta`

**Description:**\
Bins meta-aggregated LFC3D scores into percentile thresholds.

```python
bin_meta(df_bidir_meta, workdir, input_gene)
```

---

## 15. `znorm_meta`

**Description:**\
Z-normalizes meta-aggregated scores against randomized controls and labels significance.

```python
znorm_meta(df_bidir_meta, neg_stats, pos_stats, workdir, input_gene)
```

---

## 16. `average_split_bin_plots`

**Description:**\
Generates histograms, histplots, and scatterplots for positive and negative scores after binning.

```python
average_split_bin_plots(df_Z, workdir, input_gene)
```

---

# Clustering

## 17. `clustering`

**Description:**\
Performs spatial clustering of significant residues over a range of distance thresholds.

```python
clustering(df_struc, df_pvals, workdir, input_gene)
```

---

## 18. `plot_clustering`

**Description:**\
Plots clustering results including line plots and dendrograms.

```python
plot_clustering(df_struc, df_pvals, df_pvals_clust, dist, workdir, input_gene, distances, yvalues)
```

---

# Characterization

## 19. `enrichment_test`

**Description:**\
Performs enrichment tests (e.g., Fisher's exact test) for structural features.

```python
enrichment_test(df, workdir, input_gene, hit_columns, hit_threshold, feature_column, feature_values)
```

---

## 20. `plot_enrichment_test`

**Description:**\
Plots enrichment test results as odds ratios with confidence intervals.

```python
plot_enrichment_test(enrichment_results, workdir, input_gene, hit_value, feature_values)
```

---

## 21. `lfc_lfc3d_scatter`

**Description:**\
Generates LFC vs LFC3D scatter plots colored by hit significance.

```python
lfc_lfc3d_scatter(df_input, workdir, input_gene, screen_name)
```

---

## 22. `pLDDT_RSA_scatter`

**Description:**\
Generates scatter plot of RSA vs pLDDT scores, scaled by mutation weight.

```python
pLDDT_RSA_scatter(df_input, workdir, input_gene)
```

---

## 23. `hits_feature_barplot`

**Description:**\
Generates bar plots of hit counts (or fractions) across different structural categories.

```python
hits_feature_barplot(df_input, workdir, input_gene, category_col, values_cols, values_vals, value_names)
```

---

# Notes

- All outputs are saved under specified working directories.
- Many functions allow customizations via optional parameters.
- Functions are modular and can be run screen-by-screen or in batch.
